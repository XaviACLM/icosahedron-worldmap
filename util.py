import numpy as np

EPSILON = 1e-10

PHI = (1+np.sqrt(5))/2

ICOSAHEDRON_VERTICES = list(map(lambda l: np.asarray(l)/np.linalg.norm(np.asarray(l))
                                ,[(0,+1,+PHI),(0,+1,-PHI),(0,-1,+PHI),(0,-1,-PHI),
                                  (+1,+PHI,0),(+1,-PHI,0),(-1,+PHI,0),(-1,-PHI,0),
                                  (+PHI,0,+1),(+PHI,0,-1),(-PHI,0,+1),(-PHI,0,-1),
                                  ]))
"""
#not gonna work again unless de-normalized
faces = []
for i1,p1 in enumerate(ICOSAHEDRON_VERTICES):
    for i2,p2 in enumerate(ICOSAHEDRON_VERTICES[:i1]):
        for i3,p3 in enumerate(ICOSAHEDRON_VERTICES[:i2]):
            l1,l2,l3 = np.linalg.norm(p1-p2),np.linalg.norm(p2-p3),np.linalg.norm(p3-p1)
            l = (l1+l2+l3)/2
            a = np.sqrt(l*(l-l1)*(l-l2)*(l-l3))
            if abs(a-1.732)<0.01:
                faces.append((i1,i2,i3))
print(faces)
"""

ICOSAHEDRON_FACES = [(6,4,0),
(10,6,0),
(11,10,6),
(11,6,1),
(11,3,1),
(11,7,3),
(11,10,7),
(10,7,2),
(10,2,0),
(8,2,0),
(8,4,0),
(9,8,4),
(9,8,5),
(8,5,2),
(7,5,2),
(7,5,3),
(9,5,3),
(9,3,1),
(9,4,1),
(6,4,1)]
def complete_base(p):
    p /= np.linalg.norm(p)
    n = p.size
    base = [p]
    for i in range(1,n):
        pp = np.random.uniform(-1,1,3)
        for op in base:
            pp -= op*np.dot(op,pp)
        base.append(pp/np.linalg.norm(pp))
    return base

def x_rotation_matrix(angle):
    return np.asarray(((1,0,0),(0,np.cos(angle),-np.sin(angle)),(0,np.sin(angle),np.cos(angle))))

def simplex2_indices(n):
    n+=1
    qt = n*(n+1)//2
    array = np.zeros((2,qt))
    indexer = 0
    header = 0
    current_len = n
    while current_len!=0:
        array[0,header:header+current_len] = indexer
        array[1,header:header+current_len] = np.arange(current_len)
        header += current_len
        current_len -= 1
        indexer += 1
    return array

def force_extension(extension):
    if not extension[0]==".": extension = "."+extension
    def decorator(f):
        def decorated_f(self, filename):
            if filename.endswith(extension):
                return f(self,filename)
            else:
                return f(self,filename+extension)
        decorated_f.__name__ = f.__name__+"_forcedExtension"
        return decorated_f
    return decorator
