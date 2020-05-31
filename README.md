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

## Usage