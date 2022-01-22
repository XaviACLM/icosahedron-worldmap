import os
from functools import cached_property
import pickle

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from osgeo import gdal
from PIL import Image

from util import PHI, EPSILON, ICOSAHEDRON_VERTICES, ICOSAHEDRON_FACES
from util import complete_base, x_rotation_matrix, force_extension, simplex2_indices


Point3 = np.ndarray
Point2 = np.ndarray

class GeoData:
    def __init__(self, data):
        self.data = data
        self.height,self.width = data.shape
    def at(self, p: Point3):
        a = np.arctan2(p[0],p[2])-EPSILON
        b = np.arcsin(p[1])
        a = int((a/(2*np.pi)+0.5)*self.width)
        b = int((b/np.pi+0.5)*self.height)
        return self.data[b,a]
    def at_vect(self, coords):
        a = np.arctan2(coords[0],coords[2])-EPSILON
        b = np.arcsin(coords[1])
        a = ((a/(2*np.pi)+0.5)*self.width).astype(int)
        b = ((b/np.pi+0.5)*self.height).astype(int)
        return self.data[b,a]
    

class SphereTriangle:
    def __init__(self, p1: Point3, p2: Point3, p3: Point3):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.d2 = p2-p1
        self.d3 = p3-p1
    def at(self, x, y):
        p = self.p1 + self.d2*x + self.d3*y
        return p/np.linalg.norm(p)
    def at_vect(self, coords):
        v = np.expand_dims(self.p1,axis=1)+np.outer(self.d2,coords[0])+np.outer(self.d3,coords[1])
        v /= np.linalg.norm(v,axis=0,keepdims=True)
        return v
    def show_data(self, data, resolution, offset):
        tarray = np.zeros((resolution,resolution))
        for sx in range(resolution):
            for sy in range(sx):
                x,y = sx/resolution, sy/resolution
                p = self.at(x,y)
                tarray[sx-sy//2,sy] = data.at(p)+offset
        plt.imshow(tarray)
        plt.show()


class PlaneTriangle:
    def __init__(self, p1: Point2, p2: Point2, p3: Point2, reduce:float=0):
        c = (p1+p2+p3)/3
        p1,p2,p3 = (c*reduce + p*(1-reduce) for p in (p1,p2,p3))
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.d2 = p2-p1
        self.d3 = p3-p1
    def at_vect(self, coords):
        #coords.shape = (2,k)
        return np.expand_dims(self.p1,axis=1)+np.outer(self.d2,coords[0])+np.outer(self.d3,coords[1])


def project_triangle(array, data, sphereTriangle, screenTriangle, resolution):
    RES_FACTOR = 1.1
    sp0 = screenTriangle[0]
    sd1,sd2 = screenTriangle[1]-sp0,screenTriangle[2]-sp0
    #TODO vectorize
    for pidx,px in enumerate(np.linspace(0,1,num=int(resolution*RES_FACTOR))):
        for py in np.linspace(0,1-px,num=int(resolution*RES_FACTOR-pidx)):
            height_val = data.at(sphereTriangle.at(px,py))
            sp = sp0+px*sd1+py*sd2
            if height_val>0:
                array[int(resolution*sp[0]),int(resolution*sp[1])]=height_val+2000
                #array[int(resolution*sp[0]),int(resolution*sp[1])]=1000*(height_val%10000<5000)+2000
            else:
                array[int(resolution*sp[0]),int(resolution*sp[1])]=1000

    
class IcosahedralTiling:
    def __init__(self, transform):
        self.transform = transform
        self.vertices = list(map(lambda v:transform@v, ICOSAHEDRON_VERTICES))
        self.faces = ICOSAHEDRON_FACES
        
    @classmethod
    def parametrized(cls, direction, rotation):
        rmat = x_rotation_matrix(rotation)
        smat = np.stack(complete_base(direction),axis=1)
        mat = smat@rmat
        if np.linalg.det(mat)<0: mat = -mat
        return cls(mat)
    
    @force_extension("tiling")
    def save(self, filename):
        with open(filename,"wb") as f:
            pickle.dump(self.transform,f)
            
    @classmethod
    @force_extension("tiling")
    def load(cls, filename):
        with open(filename,"rb") as f:
            transform = pickle.load(f)
        return cls(transform)
        

class IcosahedralProjection:
    def __init__(self, tiling, data, copy=None):
        self.tiling = tiling
        self.data = data
        self.copy = copy
    @cached_property
    def cost_graph(self):
        #dual graph of icosahedron. Note self.faces is in the specific order for this to work
        graph = nx.dodecahedral_graph()
        SAMPLES = 100
        for f1,f2 in graph.edges:
            cost = 0
            v1,v2 = map(self.tiling.vertices.__getitem__,set(self.tiling.faces[f1]).intersection(set(self.tiling.faces[f2])))

            vs = np.linspace(v1,v2,num=SAMPLES,endpoint=False).T
            vs /= np.linalg.norm(vs,axis=0,keepdims=True)
            if self.copy is None:
                profile = self.data.at_vect(vs)
            else:
                profile = self.copy.data.at_vect(vs)
            cost = np.mean(profile>(430/(9055+430)+0.01))

            is_weeaboo = (vs[0]>0.4)*(vs[1]<-0.545)*(vs[2]<-0.5)
            cost += 10000*np.sum(is_weeaboo)
            
            graph[f1][f2]["weight"] = cost
            del vs
            del profile
        return graph
    @cached_property
    def tree(self):
        return nx.algorithms.tree.maximum_spanning_tree(self.cost_graph)
    @cached_property
    def cost(self):
        cost = 0
        for v1,v2 in self.cost_graph.edges:
            if (v1,v2) not in self.tree.edges:
                cost += self.cost_graph[v1][v2]["weight"]
        return cost
    @cached_property
    def vertex_positions(self):
        remaining_faces = list(range(20))
        face_index = 0
        face = self.tiling.faces[face_index]
        positions = dict()
        positions[face_index] = dict()
        positions[face_index][face[0]] = np.asarray((0,0))
        positions[face_index][face[1]] = np.asarray((1,0))
        positions[face_index][face[2]] = np.asarray((0.5,np.sqrt(3)/2))

        explored = [0]
        to_explore = list(self.tree[0])
        parent = dict()
        for fcdx in to_explore: parent[fcdx] = face_index
        
        while to_explore:
            face_index = to_explore.pop(0)
            if face_index in explored: continue
            explored.append(face_index)
            face = self.tiling.faces[face_index]
            anterior_face_index = parent[face_index]
            anterior_face = self.tiling.faces[anterior_face_index]

            ip = [i for i in anterior_face if i in face]
            ap = [i for i in anterior_face if i not in face][0]
            pp = [i for i in face if i not in anterior_face][0]
            positions[face_index]=dict()
            positions[face_index][ip[0]] = positions[anterior_face_index][ip[0]].copy()
            positions[face_index][ip[1]] = positions[anterior_face_index][ip[1]].copy()
            nnp = positions[anterior_face_index][ip[0]]\
                 +positions[anterior_face_index][ip[1]]\
                 -positions[anterior_face_index][ap]
            positions[face_index][pp] = nnp.copy()

            
            for neighbor in self.tree[face_index]:
                if neighbor not in explored:
                    to_explore.append(neighbor)
                    parent[neighbor] = face_index
        mnx=float("inf"); mxx=-float("inf")
        mny=float("inf"); mxy=-float("inf")
        for face in positions.values():
            for vertex in face.values():
                mnx = min(mnx,vertex[0]); mxx = max(mxx,vertex[0])
                mny = min(mny,vertex[1]); mxy = max(mxy,vertex[1])
        for j,face in positions.items():
            for i,vertex in face.items():
                positions[j][i] = np.asarray((positions[j][i][0]-mnx,positions[j][i][1]-mny))
        self.max_x = mxx-mnx
        self.max_y = mxy-mny
        return positions
    def draw(self,resolution,shrink_coeff):
        self.vertex_positions
        width = int(resolution*self.max_x+1)
        height = int(resolution*self.max_y+1)
        array = np.zeros((width,height))

        RES_FACTOR = 1.001
        pos_indices = simplex2_indices(int(resolution*RES_FACTOR))/(resolution*RES_FACTOR)
        for face_index in range(20):
            vertex_indices = self.tiling.faces[face_index]
            globe_vertices = map(self.tiling.vertices.__getitem__,vertex_indices)
            plane_vertices = map(self.vertex_positions[face_index].__getitem__,vertex_indices)
            st = SphereTriangle(*globe_vertices)
            pt = PlaneTriangle(*plane_vertices, reduce=shrink_coeff)

            #japan schwapan finder (max:0.5773502691896257)
            #sphere_indices = st.at_vect(pos_indices)
            #height_data = self.data.at_vect(sphere_indices)
            #height_data += 1000*((sphere_indices[0]>0.4)*(sphere_indices[1]<-0.545)*(sphere_indices[2]<-0.5))

            height_data = self.data.at_vect(st.at_vect(pos_indices))
            height_data = height_data*0.98+0.02 #vis. land, icosahedron unfolding
            #height_data += (height_data>0)+0.1 #vis. land, icosahedron unfolding
            #height_data += 2*((height_data>0)+1) #vis. land, icosahedron unfolding
            #height_data += 1000*((height_data>0)+1) #vis. land, icosahedron unfolding
            screen_indices = (pt.at_vect(pos_indices)*resolution).astype(int)
            array[screen_indices[0],screen_indices[1]]=height_data
        return array


heightmap_upper_pad_shape = None
heightmap_shape = None
def load_height_array():
    if "height_array.bin" in os.listdir("./data_arrays"):
        with open("./data_arrays/height_array.bin","rb") as f:
            return pickle.load(f)
    else:
        height_data_path = "./mn30_grd/mn30_grd/"
        gd = gdal.Open(height_data_path)
        array = gd.ReadAsArray()
        pad = np.zeros((int(array.shape[0]*(0.5-84/180)),array.shape[1]))
        global heightmap_upper_pad_shape
        heightmap_upper_pad_shape = pad.shape
        height_array = np.concatenate((pad,array))
        global heightmap_shape
        heightmap_shape = height_array.shape

        with open("./data_arrays/height_array.bin","wb") as f:
            pickle.dump(height_array, f)
        return height_array

def load_pop_array():
    if "pop_array.bin" in os.listdir("./data_arrays"):
        with open("./data_arrays/pop_array.bin","rb") as f:
            return pickle.load(f)
    else:
        pop_data_path = "./population/ppp_2020_1km_Aggregated.tif"
        pop_gd = gdal.Open(pop_data_path)
        pop_array = pop_gd.ReadAsArray()
        pop_array *= pop_array>=0 # tif shenanigans
        pop_array = 1.1**np.log(1+pop_array)-1#1.1**np.log(1+pop_array)
        upper_pad = np.zeros(heightmap_upper_pad_shape)
        lower_pad = np.zeros((heightmap_shape[0]-heightmap_upper_pad_shape[0]-pop_array.shape[0],heightmap_shape[1]))
        pop_array = np.concatenate((upper_pad,pop_array,lower_pad))

        with open("./data_arrays/pop_array.bin","wb") as f:
            pickle.dump(pop_array, f)
        return pop_array

if __name__=="__main__":pass

height_array = load_height_array()
pop_array = load_pop_array()
height_array-=np.min(height_array)
height_array/=np.max(height_array)
pop_array-=np.min(pop_array)
pop_array/=np.max(pop_array)

#TODO: missing pixels btw triangles
#TODO: pretty up and look into printing
#TODO: better unfold copying
height_data = GeoData(height_array)
pop_data = GeoData(pop_array)


"""
QT_TRIES = 3000
lowest_cost = float("inf")
best_projection = None
for _ in range(QT_TRIES):
    print(_,end=", ")
    direction, rotation = np.random.normal(0,1,3),np.random.uniform(0,1)
    tiling = IcosahedralTiling.parametrized(direction,rotation)
    projection = IcosahedralProjection(tiling,data)
    if projection.cost < lowest_cost:
        lowest_cost = projection.cost
        best_projection = projection
        print("\nProjection cost: ",projection.cost)
best_projection.tiling.save("best")
"""

tiling = IcosahedralTiling.load("best")
height_projection = IcosahedralProjection(tiling, height_data)
pop_projection = IcosahedralProjection(tiling, pop_data, copy=height_projection)

height_drawing = height_projection.draw(1000,0)
pop_drawing = pop_projection.draw(1000,0)
plt.imshow(height_drawing)
plt.show()
plt.imshow(pop_drawing)
plt.show()
plt.imsave("cozo.png",np.stack((height_drawing,pop_drawing,np.zeros(height_drawing.shape)),axis=2))
#plt.imsave("higresfig.png",drawing)
exit()

res=15000 #seems to be about max res
face_index = 4#8#5#1
shrink_coeff = 0.01
RES_FACTOR = 1.001
pos_indices = simplex2_indices(int(res*RES_FACTOR))/(res*RES_FACTOR)
vertex_indices = tiling.faces[face_index]
globe_vertices = list(map(tiling.vertices.__getitem__,vertex_indices))
plane_vertices = map(np.asarray,((0,0),(1,0),(0.5,np.sqrt(3)/2)))
st = SphereTriangle(*globe_vertices)
pt = PlaneTriangle(*plane_vertices, reduce=shrink_coeff)
array= np.zeros((res+1,res+1))
t_height_data = height_data.at_vect(st.at_vect(pos_indices))
t_height_data += 1000*((t_height_data>0)+1) #vis. land, icosahedron unfolding
screen_indices = (pt.at_vect(pos_indices)*res).astype(int)
array[screen_indices[0],screen_indices[1]]=t_height_data
plt.imshow(array)
plt.show()
plt.imsave("hrfig.png",array)

