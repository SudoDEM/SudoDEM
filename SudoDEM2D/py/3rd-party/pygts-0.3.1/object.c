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

#include "pygts.h"


/*-------------------------------------------------------------------------*/
/* Methods exported to python */

static PyObject*
is_unattached(PygtsObject *self, PyObject *args, PyObject *kwds)
{
  /* Objects are unattached by default */
  Py_INCREF(Py_False);
  return Py_False;
}


/* Methods table */
static PyMethodDef methods[] = {
  {"is_unattached", (PyCFunction)is_unattached,
   METH_NOARGS,
   "True if this Object o is not attached to another Object.\n"
   "Otherwise False.\n"
   "\n"
   "Trace: o.is_unattached().\n"
  }, 

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Attributes exported to python */

static PyObject *
id(PygtsObject *self, void *closure)
{
  if( self->gtsobj == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "GTS object does not exist!");
    return NULL;
  }
  /* Use the pointer of the gtsobj */
  return Py_BuildValue("i",(long)(self->gtsobj));
}


/* Methods table */
static PyGetSetDef getset[] = {
    {"id", (getter)id, NULL, "GTS object id", NULL},
    {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Python type methods */

static void
dealloc(PygtsObject* self)
{
  /* De-register entry from the object table */
  pygts_object_deregister(self);

  if(self->gtsobj_parent!=NULL) {
    /* Free the parent; GTS will free the child unless it is attached
     * to something else.
     */
    gts_object_destroy(self->gtsobj_parent);
    self->gtsobj_parent=NULL;
  }
  else {
    /* We have the only reference, and so it is safe to destroy the gtsobj 
     * (unless it was never created in the first place).
     */
    if(self->gtsobj!=NULL) {
      gts_object_destroy(self->gtsobj);
      self->gtsobj=NULL;
    }
  }
  self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PygtsObject *self;

  /* Chain up object allocation */
  self = PYGTS_OBJECT(type->tp_alloc(type, 0));
  if( self == NULL ) return NULL;

  /* Object initialization */
  self->gtsobj = NULL;
  self->gtsobj_parent = NULL;

  return (PyObject *)self;
}


static int
init(PygtsObject *self, PyObject *args, PyObject *kwds)
{
  if( self->gtsobj == NULL ) {
    PyErr_SetString(PyExc_RuntimeError, "Cannot create abstract Object");
    return -1;
  }
  
  return 0;
}


static int 
compare(PygtsObject *o1, PygtsObject *o2)
{
  if(o1->gtsobj==o2->gtsobj) {
    return 0;
  }
  else {
    return -1;
  }
}


/* Methods table */
PyTypeObject PygtsObjectType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /* ob_size */
  "gts.Object"  ,            /* tp_name */
  sizeof(PygtsObject),       /* tp_basicsize */
  0,                         /* tp_itemsize */
  (destructor)dealloc,       /* tp_dealloc */
  0,                         /* tp_print */
  0,                         /* tp_getattr */
  0,                         /* tp_setattr */
  (cmpfunc)compare,          /* tp_compare */
  0,                         /* tp_repr */
  0,                         /* tp_as_number */
  0,                         /* tp_as_sequence */
  0,                         /* tp_as_mapping */
  0,                         /* tp_hash */
  0,                         /* tp_call */
  0,                         /* tp_str */
  0,                         /* tp_getattro */
  0,                         /* tp_setattro */
  0,                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,     /* tp_flags */
  "Base object",             /* tp_doc */
  0,		             /* tp_traverse */
  0,		             /* tp_clear */
  0,		             /* tp_richcompare */
  0,		             /* tp_weaklistoffset */
  0,		             /* tp_iter */
  0,		             /* tp_iternext */
  methods,                   /* tp_methods */
  0,                         /* tp_members */
  getset,                    /* tp_getset */
  0,                         /* tp_base: attached in pygts.c */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)init,            /* tp_init */
  0,                         /* tp_alloc */
  (newfunc)new               /* tp_new */
};


/*-------------------------------------------------------------------------*/
/* Pygts functions */

gboolean 
pygts_object_check(PyObject* o)
{
  if(! PyObject_TypeCheck(o, &PygtsObjectType)) {
    return FALSE;
  }
  else {
#if PYGTS_DEBUG
    return pygts_object_is_ok(PYGTS_OBJECT(o));
#else
    return TRUE;
#endif
  }
}


gboolean 
pygts_object_is_ok(PygtsObject *o)
{
  g_return_val_if_fail(o->gtsobj!=NULL,FALSE);
  g_return_val_if_fail(g_hash_table_lookup(obj_table,o->gtsobj)!=NULL,FALSE);
  return TRUE;
}


/*-------------------------------------------------------------------------*/
/* Object table functions */

GHashTable *obj_table; /* GtsObject key, associated PyObject value */

void
pygts_object_register(PygtsObject *o)
{
  if( g_hash_table_lookup(obj_table,o->gtsobj) == NULL ) {
    g_hash_table_insert(obj_table,o->gtsobj,o);
  }
}


void
pygts_object_deregister(PygtsObject *o)
{
  if(o->gtsobj!=NULL) {
    if(g_hash_table_lookup(obj_table,o->gtsobj)==o) {
      g_hash_table_remove(obj_table,o->gtsobj);
    }
  }
}

