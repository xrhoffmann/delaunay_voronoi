"""Perform Delaunay triangulation and/or Voronoi tessellation.

2020, Xavier R. Hoffmann <xrhoffmann@gmail.com>
"""

from collections import Counter, defaultdict
from typing import Sequence, Tuple

import numpy as np


class DelVor:
    # TODO document
    """Compute triangulation and/or tessellation."""

    def __init__(
        self, *, x: Sequence[float], y: Sequence[float], buffer: float = 1
    ) -> None:
        # TODO document
        """Constructor."""
        if len(x) != len(y):
            err = f"Coordinates x ({len(x)}) and y ({len(y)}) must have same length."
            raise ValueError(err)
        elif len(x) < 3:
            err = f"Unsufficient number of points ({len(x)}). Must be >= 3."
            raise ValueError(err)
        elif buffer < 0:
            err = f"Buffer ({buffer}) should be >= 0."
            raise ValueError(err)
        else:
            self.coord = tuple(zip(x, y))
            self.buffer = buffer
            self._triangulation = None
            self._tessellation = None

    def __repr__(self) -> str:
        """Representation."""
        return f"DelVor({len(self.coord)} points)"

    def _compute_bbox(self) -> None:
        """Construct bounding box."""
        x, y = zip(*self.coord)
        self._xmin = min(x) - self.buffer
        self._ymin = min(y) - self.buffer
        self._xmax = max(x) + self.buffer
        self._ymax = max(y) + self.buffer

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
        """Compute bisecting line between two vertices."""
        px = 0.5 * (self._vertices[v2][0] + self._vertices[v1][0])
        py = 0.5 * (self._vertices[v2][1] + self._vertices[v1][1])
        vx = self._vertices[v1][0] - self._vertices[v2][0]
        vy = self._vertices[v1][1] - self._vertices[v2][1]
        mod = np.sqrt(vx ** 2 + vy ** 2)
        return px, py, vy / mod, -vx / mod

    def _circumference(
        self, *, triangle: Tuple[int, int, int]
    ) -> Tuple[Tuple[float, float], float]:
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
    def euclidean_distance(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        """Euclidean distance between two points."""
        d = (p1[0] - p2[0]) ** 2 + (p2[1] - p1[1]) ** 2
        return np.sqrt(d)

    def compute_delaunay(self):
        # TODO document
        """Delaunay triangulation."""
        # result: add _bbox and triangulation = {dict}, (tuple)
        # {dict} = node
        if self._triangulation is None:
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
            self._triangulation = 1

        return self._vertices, set(self._edges.keys())

    def compute_voronoi(self):
        # TODO document
        """Voronoi tessellation."""
        if self._triangulation is None:
            # compute delaunay if necessary
            _, _ = self.compute_delaunay()
        if self._tessellation is None:
            # assign nodes
            self._nodes = {idx: circum[0] for idx, circum in self._triangles.items()}
            # assign links
            scale = 0.1 * np.sqrt(
                (self._xmax - self._xmin) ** 2 + (self._ymax - self._ymin) ** 2
            )
            self._links = []
            self._arrows = {}
            for edge, triangles in self._edges.items():
                if len(triangles) == 2:
                    # interior edges
                    self._links.append(tuple(triangles))
                elif len(triangles) == 1:
                    # outer edges
                    v1, v2 = edge
                    v3 = [x for x in triangles[0] if x not in edge][0]

                    # compute arrow
                    x0 = 0.5 * (self._vertices[v1][0] + self._vertices[v2][0])
                    y0 = 0.5 * (self._vertices[v1][1] + self._vertices[v2][1])
                    xc, yc = self._triangles[triangles[0]][0]
                    vcx = xc - x0
                    vcy = yc - y0
                    mod = np.sqrt(vcx ** 2 + vcy ** 2)
                    vcx *= scale / mod
                    vcy *= scale / mod

                    # invert direction?
                    v3x = self._vertices[v3][0] - x0
                    v3y = self._vertices[v3][1] - y0
                    mod3 = np.sqrt(v3x ** 2 + v3y ** 2)
                    cosine = (vcx * v3x + vcy * v3y) / (scale * mod3)
                    if cosine > 0:
                        vcx *= -1
                        vcy *= -1
                    self._arrows[triangles[0]] = (vcx, vcy)
                else:
                    err = [
                        f"Edge {edge} pertains to {len(triangles)} triangles.",
                        "Should be 1 or 2.",
                    ]
                    raise ValueError(" ".join(err))

            # assign tessellation
            self._tessellation = 1

        return self._nodes, set(self._links), self._arrows
