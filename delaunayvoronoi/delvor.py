"""Perform Delaunay triangulation and/or Voronoi tessellation.

2020, Xavier R. Hoffmann <xrhoffmann@gmail.com>
"""

import math
from collections import Counter, defaultdict
from typing import Dict, Sequence, Tuple


class DelVor:
    # TODO document
    """Compute triangulation and/or tessellation."""

    def __init__(self, *, x: Sequence[float], y: Sequence[float]) -> None:
        """Constructor for class instance.

        Args:
            x: x-coordinates of points, 1d-array.
            y: y-coordinates of points, 1d-array.

        Raises:
            ValueError: If x and y have different length.
            ValueError: If less than 3 points.
        """
        if len(x) != len(y):
            err = f"Coordinates x ({len(x)}) and y ({len(y)}) must have same length."
            raise ValueError(err)
        elif len(x) < 3:
            err = f"Unsufficient number of points ({len(x)}). Must be >= 3."
            raise ValueError(err)
        else:
            self.coord = tuple(zip(x, y))
            self._triangulation = False
            self._tessellation = False

    def __repr__(self) -> str:
        """Representation."""
        return f"DelVor({len(self.coord)} points)"

    def _compute_bbox(self) -> None:
        """Construct bounding box and scale."""
        x, y = zip(*self.coord)
        self._xmin = min(x) - 1.0
        self._ymin = min(y) - 1.0
        self._xmax = max(x) + 1.0
        self._ymax = max(y) + 1.0
        self._scale = 0.1 * math.sqrt(
            (self._xmax - self._xmin) ** 2 + (self._ymax - self._ymin) ** 2
        )

    def _make_supertriangle(self) -> None:
        """Create vertices of supertriangle."""
        self._vertices = {
            -3: (self._xmin - (self._ymax - self._ymin), self._ymin),
            -2: (self._xmax + (self._ymax - self._ymin), self._ymin),
            -1: (
                0.5 * (self._xmin + self._xmax),
                self._ymax + 0.5 * (self._xmax - self._xmin),
            ),
        }

    def _bisector(self, *, v1: int, v2: int) -> Tuple[float, float, float, float]:
        """Compute bisecting line between two vertices.

        Args:
            v1:
            v2:

        Returns:
            px:
            py:
            vx:
            vy:
        """
        px = 0.5 * (self._vertices[v2][0] + self._vertices[v1][0])
        py = 0.5 * (self._vertices[v2][1] + self._vertices[v1][1])
        vx = self._vertices[v1][0] - self._vertices[v2][0]
        vy = self._vertices[v1][1] - self._vertices[v2][1]
        mod = math.sqrt(vx ** 2 + vy ** 2)
        vx, vy = vy / mod, -vx / mod
        return px, py, vx, vy

    def _circumference(
        self, *, triangle: Tuple[int, int, int]
    ) -> Tuple[Tuple[float, float], float]:
        """Compute circumcircle of a triangle.

        Args:
            triangle: Vertices of triangle in form (v1, v2, v3).

        Returns:
            center: Center of circumcircle.
            radius: Radius of circumcircle

        """
        line1 = self._bisector(v1=triangle[0], v2=triangle[1])
        line2 = self._bisector(v1=triangle[0], v2=triangle[2])
        if line1[2] == 0:
            # line1 horizontal
            x0 = line1[0]
            y0 = line2[1] + (x0 - line2[0]) * line2[3] / line2[2]
        elif line2[2] == 0:
            # line2 horizontal
            x0 = line2[0]
            y0 = line1[1] + (x0 - line1[0]) * line1[3] / line1[2]
        else:
            a1 = line1[3] / line1[2]
            b1 = line1[1] - line1[0] * a1
            a2 = line2[3] / line2[2]
            b2 = line2[1] - line2[0] * a2
            x0 = (b2 - b1) / (a1 - a2)
            y0 = a1 * x0 + b1
        radius = math.sqrt(
            (x0 - self._vertices[triangle[0]][0]) ** 2
            + (y0 - self._vertices[triangle[0]][1]) ** 2
        )
        center = (x0, y0)
        return center, radius

    @staticmethod
    def euclidean_distance(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        """Euclidean distance between two points.

        Args:
            p1: Coordinates in form (x, y).
            p2: Coordinates in form (x, y).

        Returns:
            Euclidean distance between input points.

        """
        d = (p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2
        return math.sqrt(d)

    def _arrow_vector(
        self, *, edge: Tuple[int, int], triangle: Tuple[int, int, int]
    ) -> Tuple[float, float]:
        """Arrow extending to infinity."""
        # TODO document
        v1, v2 = edge
        v3 = [x for x in triangle if x not in edge][0]

        # compute arrow
        x0 = 0.5 * (self._vertices[v1][0] + self._vertices[v2][0])
        y0 = 0.5 * (self._vertices[v1][1] + self._vertices[v2][1])
        xc, yc = self._triangles[triangle][0]
        vcx = xc - x0
        vcy = yc - y0
        mod = math.sqrt(vcx ** 2 + vcy ** 2)
        vcx *= self._scale / mod
        vcy *= self._scale / mod

        # invert direction?
        v3x = self._vertices[v3][0] - x0
        v3y = self._vertices[v3][1] - y0
        cosine = vcx * v3x + vcy * v3y
        if cosine > 0:
            vcx *= -1
            vcy *= -1
        return vcx, vcy

    def _boundary_arrow(
        self, *, node: Tuple[float, float], vector: Tuple[float, float]
    ) -> Tuple[float, float]:
        """Extend arrows to bounding box."""
        if vector[0] < 0:
            # left boundary
            bound_x = self._xmin
            bound_y = node[1] + (bound_x - node[0]) * vector[1] / vector[0]
            if bound_y < self._ymin:
                bound_y = self._ymin
                bound_x = node[0] + (bound_y - node[1]) * vector[0] / vector[1]
            elif bound_y > self._ymax:
                bound_y = self._ymax
                bound_x = node[0] + (bound_y - node[1]) * vector[0] / vector[1]
        elif vector[0] > 0:
            # right boundary
            bound_x = self._xmax
            bound_y = node[1] + (bound_x - node[0]) * vector[1] / vector[0]
            if bound_y < self._ymin:
                bound_y = self._ymin
                bound_x = node[0] + (bound_y - node[1]) * vector[0] / vector[1]
            elif bound_y > self._ymax:
                bound_y = self._ymax
                bound_x = node[0] + (bound_y - node[1]) * vector[0] / vector[1]
        else:
            # vertical line
            bound_x = node[0]
            if vector[1] > 0:
                bound_y = self._ymax
            else:
                bound_y = self._ymin
        return bound_x, bound_y

    def compute_delaunay(
        self
    ) -> Tuple[Dict[int, Tuple[float, float]], Tuple[Tuple[int, int], ...]]:
        # TODO document
        """Delaunay triangulation."""
        if not self._triangulation:
            # construct bbox
            self._compute_bbox()

            # make supertriangle
            self._make_supertriangle()
            super_vertices = tuple(self._vertices.keys())

            # add first vertex
            self._vertices[0] = self.coord[0]
            self._triangles = {}
            for triangle in [(-3, -2, 0), (-3, -1, 0), (-2, -1, 0)]:
                self._triangles[triangle] = self._circumference(triangle=triangle)

            # iterate over all vertices
            for i, vertex_coord in enumerate(self.coord[1:]):
                vertex_id = i + 1
                self._vertices[vertex_id] = vertex_coord
                # find bad triangles
                bad_triangles = []
                edges = []
                for triangle, circumcircle in self._triangles.items():
                    dist = self.euclidean_distance(
                        self._vertices[vertex_id], circumcircle[0]
                    )
                    if dist < circumcircle[1]:
                        bad_triangles.append(triangle)
                        edges += [
                            (triangle[i], triangle[j])
                            for i, j in [(0, 1), (0, 2), (1, 2)]
                        ]
                # remove bad triangles
                for triangle in bad_triangles:
                    del self._triangles[triangle]
                # add new triangles
                new_triangles = [
                    (edge[0], edge[1], vertex_id)
                    for edge, count in Counter(edges).items()
                    if count == 1
                ]
                for triangle in new_triangles:
                    self._triangles[triangle] = self._circumference(triangle=triangle)

            # remove super vertices
            for vertex in super_vertices:
                del self._vertices[vertex]
            # remove bad triangles (including super vertices)
            bad_triangles = [
                triangle
                for triangle in self._triangles
                if any(vertex in super_vertices for vertex in triangle)
            ]
            for triangle in bad_triangles:
                del self._triangles[triangle]

            # compute edges
            self._edges = defaultdict(list)
            for triangle in self._triangles.keys():
                for i, j in [(0, 1), (0, 2), (1, 2)]:
                    self._edges[(triangle[i], triangle[j])].append(triangle)

            # assign triangulation
            self._triangulation = True

        return self._vertices, tuple(self._edges.keys())

    def compute_voronoi(
        self
    ) -> Tuple[
        Dict[Tuple[int, int, int], Tuple[float, float]],
        Tuple[Tuple[Tuple[int, int, int], Tuple[int, int, int]], ...],
        Dict[Tuple[int, int, int], Tuple[float, float]],
    ]:
        # TODO document
        """Voronoi tessellation."""
        if not self._triangulation:
            # compute delaunay if necessary
            _, _ = self.compute_delaunay()
        if not self._tessellation:
            # assign nodes
            self._nodes = {idx: circum[0] for idx, circum in self._triangles.items()}
            # assign links
            self._links = []
            self._arrows = {}
            for edge, triangles in self._edges.items():
                if len(triangles) == 2:
                    # interior edges
                    self._links.append((triangles[0], triangles[1]))
                elif len(triangles) == 1:
                    # outer edges
                    vcx, vcy = self._arrow_vector(edge=edge, triangle=triangles[0])
                    self._arrows[triangles[0]] = (vcx, vcy)
                else:
                    err = [
                        f"Edge {edge} pertains to {len(triangles)} triangles.",
                        "Should be 1 or 2.",
                    ]
                    raise ValueError(" ".join(err))

            # assign tessellation
            self._tessellation = True

        return self._nodes, (*self._links,), self._arrows

    def prepare_plot(
        self, bbox: Tuple[Tuple[float, float], Tuple[float, float]] = None
    ) -> Tuple[
        Tuple[Tuple[float, float], ...],
        Tuple[Tuple[Tuple[float, float], Tuple[float, float]], ...],
        Tuple[Tuple[float, float], ...],
        Tuple[Tuple[Tuple[float, float], Tuple[float, float]], ...],
    ]:
        """Compute elements for plotting."""
        if not self._tessellation:
            _, _, _ = self.compute_voronoi()
        if bbox is not None:
            self._xmin = bbox[0][0]
            self._ymin = bbox[0][1]
            self._xmax = bbox[1][0]
            self._ymax = bbox[1][1]

        # triangulation
        points_delaunay = tuple(self._vertices.values())
        edges_delaunay = tuple(
            (self._vertices[edge[0]], self._vertices[edge[1]]) for edge in self._edges
        )

        # tessellation
        # separate interior and exterior nodes
        int_nodes = {
            triangle: center
            for triangle, center in self._nodes.items()
            if self._xmin <= center[0] <= self._xmax
            and self._ymin <= center[1] <= self._ymax
        }
        ext_nodes = {
            triangle: center
            for triangle, center in self._nodes.items()
            if triangle in set(self._nodes.keys()).difference(set(int_nodes.keys()))
        }
        points_voronoi = tuple(int_nodes.values())

        # add interior edges and list pending arrows
        edges_voronoi = []
        arrows = [
            (self._nodes[node], vector)
            for node, vector in self._arrows.items()
            if node in int_nodes
        ]
        for link in self._links:
            if link[0] in ext_nodes:
                vx = self._nodes[link[0]][0] - self._nodes[link[1]][0]
                vy = self._nodes[link[0]][1] - self._nodes[link[1]][1]
                arrows.append((self._nodes[link[1]], (vx, vy)))
            elif link[1] in ext_nodes:
                vx = self._nodes[link[1]][0] - self._nodes[link[0]][0]
                vy = self._nodes[link[1]][1] - self._nodes[link[0]][1]
                arrows.append((self._nodes[link[0]], (vx, vy)))
            else:
                # interior node
                edges_voronoi.append((self._nodes[link[0]], self._nodes[link[1]]))

        # extend arrows
        for arrow in arrows:
            bounds = self._boundary_arrow(node=arrow[0], vector=arrow[1])
            edges_voronoi.append((arrow[0], bounds))

        return points_delaunay, edges_delaunay, points_voronoi, tuple(edges_voronoi)
