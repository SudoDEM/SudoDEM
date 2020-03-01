# pygts - python package for the manipulation of triangulated surfaces
#
#   Copyright (C) 2009 Thomas J. Duck
#   All rights reserved.
#
#   Thomas J. Duck <tom.duck@dal.ca>
#   Department of Physics and Atmospheric Science,
#   Dalhousie University, Halifax, Nova Scotia, Canada, B3H 3J5
#
# NOTICE
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Library General Public
#   License as published by the Free Software Foundation; either
#   version 2 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Library General Public License for more details.
#
#   You should have received a copy of the GNU Library General Public
#   License along with this library; if not, write to the
#   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#   Boston, MA 02111-1307, USA.

"""A package for constructing and manipulating triangulated surfaces.

PyGTS is a python binding for the GNU Triangulated Surface (GTS) 
Library, which may be used to build, manipulate, and perform
computations on triangulated surfaces.

The following geometric primitives are provided:

  * Point - a point in 3D space
  * Vertex - a Point in 3D space that may be used to define a Segment
  * Segment - a line defined by two Vertex end-points
  * Edge - a Segment that may be used to define the edge of a Triangle
  * Triangle - a triangle defined by three Edges
  * Face - a Triangle that may be used to define a face on a Surface
  * Surface - a surface composed of Faces

A tetrahedron is assembled from these primitives as follows.  First,
create Vertices for each of the tetrahedron's points:

  .. code-block:: python

    import gts
    
    v1 = gts.Vertex(1,1,1)
    v2 = gts.Vertex(-1,-1,1)
    v3 = gts.Vertex(-1,1,-1)
    v4 = gts.Vertex(1,-1,-1)

Next, connect the four vertices to create six unique Edges:

  .. code-block:: python

    e1 = gts.Edge(v1,v2)
    e2 = gts.Edge(v2,v3)
    e3 = gts.Edge(v3,v1)
    e4 = gts.Edge(v1,v4)
    e5 = gts.Edge(v4,v2)
    e6 = gts.Edge(v4,v3)

The four triangular faces are composed using three edges each:

  .. code-block:: python

    f1 = gts.Face(e1,e2,e3)
    f2 = gts.Face(e1,e4,e5)
    f3 = gts.Face(e2,e5,e6)
    f4 = gts.Face(e3,e4,e6)

Finally, the surface is assembled from the faces:

  .. code-block:: python

    s = gts.Surface()
    for face in [f1,f2,f3,f4]:
        s.add(face)

Some care must be taken in the orientation of the faces.  In the above
example, the surface normals are pointing inward, and so the surface
technically defines a void, rather than a solid.  To create a 
tetrahedron with surface normals pointing outward, use the following
instead:

  .. code-block:: python

    f1.revert()
    s = Surface()
    for face in [f1,f2,f3,f4]:
        if not face.is_compatible(s):
            face.revert()
        s.add(face)

Once the Surface is constructed, there are many different operations that
can be performed.  For example, the volume can be calculated using:

  .. code-block:: python

    s.volume()

The difference between two Surfaces s1 and s2 is given by:

  .. code-block:: python

    s3 = s2.difference(s1)

Etc.

It is also possible to read in GTS data files and plot surfaces to
the screen.  See the example programs packaged with PyGTS for
more information.
"""

from pygts import *
