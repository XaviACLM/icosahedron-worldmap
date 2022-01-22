# icosahedron-worldmap
I've always wanted a huge printed map, but putting up a projection with significant global distortion would feel weird. This is my attempt to plot something nice on a dymaxion-like projection.

Of note:

- Regarding the necessary data: the program needs data on height and population density to run. By default, it will look for these in ./mn30_grd/mn30_grd/ and ./population/ppp_2020_1km_Aggregated.tif, respectively. These names should suggest where to find the datasets themselves (they are not in the repo, of course). 

- Regarding the reading of this data: gdal is used, through the python package osgeo. The installation process might be a bit finnicky. In the script, the data is immediately cast to array format and then a backup of this array is made, to be read later if the program is ever re-executed (as the reading of the data is by far the slowest part). Note that this creates two files of about 7Gb each, in ./data_arrays. These can be safely deleted and will be re-generated whenever the program is executed again.

- No triangle slicing. This is both to limit the scope of the project and because I think it looks ugly. It seems that with this constraint, the best (least land separated) way to unfold an icosahedral projection always involves slicing Japan in half. To placate the weeaboo I cohabitate with, this was specifically avoided, at the cost of slicing off the lower tip of the american continent and a little bit of Sakhalin island.

![](https://i.imgur.com/maqM53u.jpg)

- Rather fast vectorized computation. This 15x15mpx image took about 6.5s to generate:

![](https://i.imgur.com/F0goIlw.jpg)

- Also regarding the image above, maximum resolution (at least for terrain) seems to be at about 1km:3px.

Both maps above are only with height data. Population data is also available, and both can be shown in tandem:

![](https://i.imgur.com/8VYUlou.png)
