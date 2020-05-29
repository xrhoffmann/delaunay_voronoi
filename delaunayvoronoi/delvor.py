"""Perform Delaunay triangulation and/or Voronoi tessellation.

2020, Xavier R. Hoffmann <xrhoffmann@gmail.com>
"""

from collections import Counter, defaultdict

import numpy as np


class DelVor:
    # TODO document
    """Compute triangulation and/or tessellation."""

    def __init__(self, *, x, y, buffer=1):
        # TODO document
        """Constructor."""
        if len(x) != len(y):
            err = "Coordinates x and y should have the same length."
            raise ValueError(err)
        elif len(x) < 3:
            err = "At least 3 points are necessary."
            raise ValueError(err)
        elif buffer < 0:
            err = f"Buffer ({buffer}) should be >= 0."
            raise ValueError(err)
        else:
            self.coord = tuple(zip(x, y))
            self.buffer = buffer
            self.triangulation = None
            self.tessellation = None
            self._xmin = None
            self._ymin = None
            self._xmax = None
            self._ymax = None
            self._vertices = None
            self._edges = None
            self._triangles = None
            self._nodes = None
            self._links = None

    def __repr__(self):
        """Representation."""
        return f"DelVor({len(self.points)} points)"

    def _compute_bbox(self):
        """Construct bounding box."""
        x, y = zip(*self.coord)
        self._xmin = min(x) - self.buffer
        self._ymin = min(y) - self.buffer
        self._xmax = max(x) + self.buffer
        self._ymax = max(y) + self.buffer

    def _make_supertriangle(self):
        """Create vertices of supertriangle."""
        self._vertices = {
            -3: (self._xmin - (self._ymax - self._ymin), self._ymin),
            -2: (self._xmax + (self._ymax - self._ymin), self._ymin),
            -1: (
                0.5 * (self._xmin + self._xmax),
                self._ymax + 0.5 * (self._xmax - self._xmin),
            ),
        }

    def _bisector(self, *, v1, v2):
        """Compute bisecting line between two vertices."""
        px = 0.5 * (self._vertices[v2][0] + self._vertices[v1][0])
        py = 0.5 * (self._vertices[v2][1] + self._vertices[v1][1])
        vx = self._vertices[v1][0] - self._vertices[v2][0]
        vy = self._vertices[v1][1] - self._vertices[v2][1]
        mod = np.sqrt(vx ** 2 + vy ** 2)
        vx /= mod
        vy /= mod
        dx = vy
        dy = -vx
        return px, py, dx, dy

    def _circumference(self, triangle):
        """Compute circumcircle of a triangle."""
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
        radius = np.sqrt(
            (x0 - self._vertices[triangle[0]][0]) ** 2
            + (y0 - self._vertices[triangle[0]][1]) ** 2
        )
        return (x0, y0), radius

    @staticmethod
    def euclidean_distance(p1, p2):
        """Euclidean distance between two points."""
        d = (p1[0] - p2[0]) ** 2 + (p2[1] - p1[1]) ** 2
        return np.sqrt(d)

    def compute_delaunay(self):
        # TODO document
        """Delaunay triangulation."""
        # result: add _bbox and triangulation = {dict}, (tuple)
        # {dict} = node
        if self.triangulation is None:
            # construct bbox
            self._compute_bbox()

            # make supertriangle
            self._make_supertriangle()
            super_vertices = tuple(self._vertices.keys())

            # add first vertex
            self._vertices[0] = self.coord[0]
            self._triangles = {}
            for tr in [(-3, -2, 0), (-3, -1, 0), (-2, -1, 0)]:
                self._triangles[tr] = self._circumference(tr)

            # iterate over all vertices
            for i, vertex_coord in enumerate(self.coord[1:]):
                vertex_id = i + 1
                self._vertices[vertex_id] = vertex_coord
                # find bad triangles
                bad_triangles = []
                edges = []
                # TODO optimize with numpy vectoritzation?
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
                    self._triangles[triangle] = self._circumference(triangle)

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
            self.triangulation = self._vertices, set(self._edges.keys())

        return self.triangulation

    def compute_voronoi(self):
        # TODO document
        """Voronoi tessellation."""
        if self.triangulation is None:
            _, _ = self.compute_delaunay()
        if self.tessellation is None:
            # assign nodes
            self._nodes = {idx: circum[0] for idx, circum in self._triangles.items()}
            # assign links

            # assign tessellation
            self.tessellation = self._nodes, self._links

        return self.tessellation
