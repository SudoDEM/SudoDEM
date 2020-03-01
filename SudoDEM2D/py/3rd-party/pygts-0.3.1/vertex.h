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

#ifndef __PYGTS_VERTEX_H__
#define __PYGTS_VERTEX_H__

typedef struct _PygtsObject PygtsVertex;


#define PYGTS_VERTEX(o)							\
  ( PyObject_TypeCheck((PyObject*)o, &PygtsVertexType) ?		\
    (PygtsVertex*)o :							\
    pygts_vertex_from_sequence((PyObject*)o) )

#define PYGTS_VERTEX_AS_GTS_VERTEX(o)					\
  ( PyObject_TypeCheck((PyObject*)o, &PygtsVertexType) ?		\
    GTS_VERTEX(PYGTS_OBJECT(o)->gtsobj) :				\
    GTS_VERTEX(PYGTS_OBJECT(PYGTS_VERTEX(o))->gtsobj) )

extern PyTypeObject PygtsVertexType;

gboolean pygts_vertex_check(PyObject* o);
gboolean pygts_vertex_is_ok(PygtsVertex *v);

PygtsVertex* pygts_vertex_new(GtsVertex *f);
PygtsVertex* pygts_vertex_from_sequence(PyObject *tuple);


/*-------------------------------------------------------------------------*/
/* Parent GTS segment for GTS vertices */

/* Define a GtsSegment subclass that can be readily identified as the parent 
 * of an encapsulated GtsVertex.  The pygts_parent_segment_class() function 
 * is defined at the bottom, and is what ultimately allows the distinction 
 * to be made.  This capability is used for vertex replacement operations.
 */
typedef struct _GtsSegment PygtsParentSegment;

#define PYGTS_PARENT_SEGMENT(obj) GTS_OBJECT_CAST(obj,\
					     GtsSegment,\
					     pygts_parent_segment_class())

#define PYGTS_IS_PARENT_SEGMENT(obj)(gts_object_is_from_class(obj,\
                                             pygts_parent_segment_class()))

GtsSegmentClass* pygts_parent_segment_class(void);


/* GTS vertices in parent segments */

typedef struct _GtsVertex PygtsParentVertex;

#define PYGTS_PARENT_VERTEX(obj) GTS_OBJECT_CAST(obj,\
						 GtsVertex,\
						 pygts_parent_vertex_class())

#define PYGTS_IS_PARENT_VERTEX(obj)(gts_object_is_from_class(obj,\
                                             pygts_parent_vertex_class()))

GtsVertexClass *pygts_parent_vertex_class(void);

#endif /* __PYGTS_VERTEX_H__ */
