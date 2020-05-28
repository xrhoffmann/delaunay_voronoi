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

    def _construct_bbox(self):
        """Construct bounding box."""
        x, y = zip(*self.coord)
        self._bbox = (
            (min(x) - self.buffer, min(y) - self.buffer),
            (max(x) + self.buffer, max(y) + self.buffer),
        )

    def compute_delaunay(self):
        # TODO document
        """Delaunay triangulation."""
        # result: add _bbox and triangulation = {dict}, (tuple)
        # {dict} = node
        if self.triangulation is None:
            # construct bbox
            self._construct_bbox()
        return self.triangulation

    def compute_tessellation(self):
        # TODO document
        """Voronoi tessellation."""
        if self.triangulation is None:
            pass
            # self.compute_delaunay()
        pass
