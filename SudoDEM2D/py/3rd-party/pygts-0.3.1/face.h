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

#ifndef __PYGTS_FACE_H__
#define __PYGTS_FACE_H__

#ifndef gts_face_is_unattached
#define gts_face_is_unattached(f) ((f)->surfaces == NULL ? TRUE : FALSE)
#endif

typedef struct _PygtsObject PygtsFace;

#define PYGTS_FACE(obj) ((PygtsFace*)obj)

#define PYGTS_FACE_AS_GTS_FACE(o) (GTS_FACE(PYGTS_OBJECT(o)->gtsobj))

extern PyTypeObject PygtsFaceType;

gboolean pygts_face_check(PyObject* o);
gboolean pygts_face_is_ok(PygtsFace *f);

PygtsFace* pygts_face_new(GtsFace *f);

#endif /* __PYGTS_FACE_H__ */
