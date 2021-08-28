# delaunay_voronoi
Delaunay triangulation and Voronoi tessellation in Python.

2020, Xavier R. Hoffmann <xrhoffmann@gmail.com>

## Basics

Given a set of points in an Euclidean 2D space, we compute the [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) (DT)
using the [Bowyer-Watson algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm). 
We compute the [Voronoi tessellation](https://en.wikipedia.org/wiki/Voronoi_diagram) (VT) as 
the dual graph of the DT.

For clarity, we define the DT graph as a set of **vertices** connected by **edges**. On the 
other hand, the VT network is defined as a set of **nodes** connected by **links**, with 
additional **arrows** pointing outward.

We only use core Python. The examples use `matplotlib` for plotting.

Coordinates and distances are assumed Euclidean. Algorithm and implementation details can be found in XXX.

## Usage
* Initialize an instance of the `DelVor` class with the `x` and `y` coordinates.
```
from delaunayvoronoi import delvor

x = <my_x_coordinates>
y = <my_y_coordinates>
my_instance = delvor.DelVor(x=x, y=y)
```
* Compute the Delaunay triangluation. 
```
vertices, edges = my_instance.compute_delaunay()
```
This returns the original points (`vertices`)
and the connections between them (`edges`). 
* Compute the Voronoi triangulation.
```
nodes, links, arrows = my_instance.compute_voronoi()
```
This returns all the circumcenters (`nodes`),
the connections between them (`links`) and the lines that 
extends towards (`arrows`).
* Define limits of the bounding box and prepare elements for plotting.
```
lower_left = (x_min, y_min)
upper_right = (x_max, y_max)
bounding_box = (lower_left, upper_right)
p_del, e_del, p_vor, e_vor = a.prepare_plot(bbox=bounding_box)
```


### Representation


## Example
