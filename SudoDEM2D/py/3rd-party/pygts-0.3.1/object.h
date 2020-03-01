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

#ifndef __PYGTS_OBJECT_H__
#define __PYGTS_OBJECT_H__


typedef struct _PygtsObject PygtsObject;
typedef struct _PygtsMethods PygtsMethods;

#define PYGTS_OBJECT(obj) ((PygtsObject*)obj)

struct _PygtsObject {
  PyObject_HEAD
  GtsObject *gtsobj;         /* Encapsulated GtsObject */
  GtsObject *gtsobj_parent;  /* A parent object to ensure persistence */
};

extern PyTypeObject PygtsObjectType;
extern PygtsMethods PygtsObjectMethods;

gboolean pygts_object_check(PyObject* o);
gboolean pygts_object_is_ok(PygtsObject *o);

extern GHashTable *obj_table; /* GtsObject key, associated PyObject value */
void pygts_object_register(PygtsObject *o);
void pygts_object_deregister(PygtsObject *o);

#endif /* __PYGTS_OBJECT_H__ */
