"""Perform Delaunay triangulation and/or Voronoi tessellation.

2020, Xavier R. Hoffmann <xrhoffmann@gmail.com>
"""


class DelVor:
    """Compute triangulation and/or tessellation."""

    def __init__(self, *, x, y):
        """Constructor."""
        if len(x) != len(y):
            err = "Coordinates x and y should have the same length."
            raise ValueError(err)
        elif len(x) < 3:
            err = "At least 3 points are necessary."
            raise ValueError(err)
        else:
            pass

    def __repr__(self):
        """Representation."""
        pass
