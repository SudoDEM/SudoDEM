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

#ifndef __PYGTS_POINT_H__
#define __PYGTS_POINT_H__

typedef struct _PygtsObject PygtsPoint;

#define PYGTS_POINT(o) ( PyObject_TypeCheck((PyObject*)o, &PygtsPointType) ? \
			 (PygtsPoint*)o :				\
			 pygts_point_from_sequence((PyObject*)o) )

#define PYGTS_POINT_AS_GTS_POINT(o) (GTS_POINT(PYGTS_OBJECT(o)->gtsobj))

extern PyTypeObject PygtsPointType;

gboolean pygts_point_check(PyObject* o);
gboolean pygts_point_is_ok(PygtsPoint *o);

PygtsPoint* pygts_point_from_sequence(PyObject *tuple);
int pygts_point_compare(GtsPoint* p1,GtsPoint* p2);

gint pygts_point_rotate(GtsPoint* p,gdouble dx,gdouble dy,gdouble dz,gdouble a);
gint pygts_point_scale(GtsPoint* p, gdouble dx, gdouble dy, gdouble dz);
gint pygts_point_translate(GtsPoint* p, gdouble dx, gdouble dy, gdouble dz);

#endif /* __PYGTS_POINT_H__ */
