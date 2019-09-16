# vtugeo

vtugeo is a small python tool which convert  2D (x-y) vtu unstructured grids to georeferenced geodataframe structures and vice versa

## Features
* Converts VTU grids consisting of triangles to georeferenc geodataframes having element and node ids as attributes
* Optional Output as shapefiles (if converting from vtu to geo)
* Optionally, problem of limit digit length in dbase column header (shapefiles) is intercepted by comparing stripped string with additional list of full strings
* Problem of counter-clockwise order of nodes in vtu and clockwise in shp is recognized and intercepted

## Dependencies

* Python 3
* numpy
* geopandas
* meshio
* shapely

## Limitations
* io stream problems of meshio. I suggest to switch around between "vtu-ascii" and "vtu-binary"
* Currently, only scalar fields are supported



## Code Example
```python
import vtugeo
#convert from vtu to geodataframe
gdf_test=vtugeo.vtu_to_geo(fctn_inpt='.\Examples\Songhua.vtu',Output=None)
#convert the geodataframe structure back to vtu
vtugeo.geo_to_vtu(fctn_inpt=gdf_test,Output='.\Examples\outputmsh.vtu')
```

## Authors

* **Erik Nixdorf**

## Acknowledgments

* Thx to all who helped to improve the code


