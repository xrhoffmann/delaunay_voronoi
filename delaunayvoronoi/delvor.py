"""Perform Delaunay triangulation and/or Voronoi tessellation.

2020, Xavier R. Hoffmann <xrhoffmann@gmail.com>
"""


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
        self._vertices = {
            -3: (self._xmin - (self._ymax - self._ymin), self._ymin),
            -2: (self._xmax + (self._ymax - self._ymin), self._ymin),
            -1: (
                0.5 * (self._xmin + self._xmax),
                self._ymax + 0.5 * (self._xmax - self._xmin),
            ),
        }

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
            _super_vertices = self._vertices.keys()

            # add first vertex

            # iterate over all vertices

            # remove super vertices
            for vertex in _super_vertices:
                pass

        return self.triangulation

    def compute_tessellation(self):
        # TODO document
        """Voronoi tessellation."""
        if self.triangulation is None:
            pass
            # self.compute_delaunay()
        pass
