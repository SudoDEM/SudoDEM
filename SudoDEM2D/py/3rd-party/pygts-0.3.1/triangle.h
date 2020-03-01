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

#ifndef __PYGTS_TRIANGLE_H__
#define __PYGTS_TRIANGLE_H__

typedef struct _PygtsObject PygtsTriangle;

#define PYGTS_TRIANGLE(obj) ((PygtsTriangle*)obj)

#define PYGTS_TRIANGLE_AS_GTS_TRIANGLE(o) \
  (GTS_TRIANGLE(PYGTS_OBJECT(o)->gtsobj))

extern PyTypeObject PygtsTriangleType;

gboolean pygts_triangle_check(PyObject* o);
gboolean pygts_triangle_is_ok(PygtsTriangle *t);

PygtsTriangle* pygts_triangle_new(GtsTriangle *t);

int pygts_triangle_compare(GtsTriangle* t1,GtsTriangle* t2);



/* Replacement for gts_triangle_is_ok().  The problem is that sometimes the 
 * "reserved" variable is set in a face by gts, and so this function fails.  
 * e.g., The error occurs when gts_triangle_is_ok() is called during 
 * iteration over faces in a surface.  This function ignores that check so
 * that there is no failure when PYGTS_DEBUG is set.  A bug report should be 
 * submitted.
 */
gboolean pygts_gts_triangle_is_ok(GtsTriangle *t);

#endif /* __PYGTS_TRIANGLE_H__ */
