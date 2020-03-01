/* pygts - python package for the manipulation of triangulated surfaces
 *
 *   Copyright (C) 2009 Thomas J. Duck
 *   All rights reserved.
 *
 *   Thomas J. Duck <tom.duck@dal.ca>
 *   Department of Physics and Atmospheric Science,
 *   Dalhousie University, Halifax, Nova Scotia, Canada, B3H 3J5
 *
 * NOTICE
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.
 *
 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, write to the
 *   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 *   Boston, MA 02111-1307, USA.
 */

#ifndef __PYGTS_EDGE_H__
#define __PYGTS_EDGE_H__

typedef struct _PygtsObject PygtsEdge;

#define PYGTS_EDGE(obj) ((PygtsEdge*)obj)

#define PYGTS_EDGE_AS_GTS_EDGE(o) (GTS_EDGE(PYGTS_OBJECT(o)->gtsobj))

extern PyTypeObject PygtsEdgeType;

gboolean pygts_edge_check(PyObject* o);
gboolean pygts_edge_is_ok(PygtsEdge *e);

PygtsEdge* pygts_edge_new(GtsEdge *e);


/*-------------------------------------------------------------------------*/
/* Parent GTS triangle for GTS edges */

/* Define a GtsTriangle subclass that can be readily identified as the parent 
 * of an encapsulated GtsEdge.  The pygts_parent_triangle_class() function 
 * is defined at the bottom, and is what ultimately allows the distinction 
 * to be made.  This capability is used for edge replacement operations.
 */
typedef struct _GtsTriangle PygtsParentTriangle;

#define PYGTS_PARENT_TRIANGLE(obj) GTS_OBJECT_CAST(obj,\
					      GtsTriangle,\
					      pygts_parent_triangle_class())

#define PYGTS_IS_PARENT_TRIANGLE(obj)(gts_object_is_from_class(obj,\
                                             pygts_parent_triangle_class()))

GtsTriangleClass* pygts_parent_triangle_class(void);


/* GTS edges in parent triangles */

typedef struct _GtsEdge PygtsParentEdge;

#define PYGTS_PARENT_EDGE(obj) GTS_OBJECT_CAST(obj,\
						 GtsEdge,\
						 pygts_parent_edge_class())

#define PYGTS_IS_PARENT_EDGE(obj)(gts_object_is_from_class(obj,\
                                             pygts_parent_edge_class()))

GtsEdgeClass* pygts_parent_edge_class(void);

#endif /* __PYGTS_EDGE_H__ */
