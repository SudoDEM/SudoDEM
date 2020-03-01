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


#if PYGTS_DEBUG
  #define SELF_CHECK if(!pygts_edge_check((PyObject*)self)) {         \
                       PyErr_SetString(PyExc_RuntimeError,            \
                       "problem with self object (internal error)");  \
		       return NULL;                                   \
                     }
#else
  #define SELF_CHECK
#endif


/*-------------------------------------------------------------------------*/
/* Methods exported to python */

static PyObject*
is_ok(PygtsEdge *self, PyObject *args)
{
  if(pygts_edge_is_ok(self)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
is_unattached(PygtsEdge *self, PyObject *args)
{
  guint n;

  SELF_CHECK

  /* Check for attachments other than to the gtsobj_parent */
  n = g_slist_length(PYGTS_EDGE_AS_GTS_EDGE(self)->triangles);
  if( n > 1 ) {
    Py_INCREF(Py_False);
    return Py_False;
  }
  else if( n == 1 ){
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    PyErr_SetString(PyExc_RuntimeError,"Edge lost parent (internal error)");
    return NULL;
  }
}


/* replace() works, but can break Triangles and so is disabled */

/* static PyObject* */
/* replace(PygtsEdge *self, PyObject *args) */
/* { */
/*   PyObject *e2_; */
/*   PygtsEdge *e2; */
/*   GSList *parents=NULL, *i, *cur; */

/* #if PYGTS_DEBUG */
/*   if(!pygts_edge_check((PyObject*)self)) { */
/*     PyErr_SetString(PyExc_TypeError, */
/* 		    "problem with self object (internal error)"); */
/*     return NULL; */
/*   } */
/* #endif */

/*   /\* Parse the args *\/   */
/*   if(! PyArg_ParseTuple(args, "O", &e2_) ) { */
/*     return NULL; */
/*   } */

/*   /\* Convert to PygtsObjects *\/ */
/*   if(!pygts_edge_check(e2_)) { */
/*     PyErr_SetString(PyExc_TypeError,"expected an Edge"); */
/*     return NULL; */
/*   } */
/*   e2 = PYGTS_EDGE(e2_); */

/*   if(PYGTS_OBJECT(self)->gtsobj!=PYGTS_OBJECT(e2)->gtsobj) { */
/*     /\* (Ignore self-replacement) *\/ */

/*     /\* Detach and save any parent triangles *\/ */
/*     i = GTS_EDGE(PYGTS_OBJECT(self)->gtsobj)->triangles; */
/*     while(i!=NULL) { */
/*       cur = i; */
/*       i = i->next; */
/*       if(PYGTS_IS_PARENT_TRIANGLE(cur->data)) { */
/* 	GTS_EDGE(PYGTS_OBJECT(self)->gtsobj)->triangles =  */
/* 	  g_slist_remove_link(GTS_EDGE(PYGTS_OBJECT(self)->gtsobj)->triangles, */
/* 			      cur); */
/* 	parents = g_slist_prepend(parents,cur->data); */
/* 	g_slist_free_1(cur); */
/*       } */
/*     } */

/*     /\* Perform the replace operation *\/ */
/*     gts_edge_replace(GTS_EDGE(PYGTS_OBJECT(self)->gtsobj), */
/* 		     GTS_EDGE(PYGTS_OBJECT(e2)->gtsobj)); */

/*     /\* Reattach the parent segments *\/ */
/*     i = parents; */
/*     while(i!=NULL) { */
/*       GTS_EDGE(PYGTS_OBJECT(self)->gtsobj)->triangles =  */
/* 	g_slist_prepend(GTS_EDGE(PYGTS_OBJECT(self)->gtsobj)->triangles, */
/* 			i->data); */
/*       i = i->next; */
/*     } */
/*     g_slist_free(parents); */
/*   } */

/* #if PYGTS_DEBUG */
/*   if(!pygts_edge_check((PyObject*)self)) { */
/*     PyErr_SetString(PyExc_TypeError, */
/* 		    "problem with self object (internal error)"); */
/*     return NULL; */
/*   } */
/* #endif */

/*   Py_INCREF(Py_None); */
/*   return Py_None; */
/* } */


static PyObject*
face_number(PygtsEdge *self, PyObject *args)
{
  PyObject *s_;
  GtsSurface *s;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &s_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_surface_check(s_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Surface");
    return NULL;
  }
  s = PYGTS_SURFACE_AS_GTS_SURFACE(s_);

  return Py_BuildValue("i",
		       gts_edge_face_number(PYGTS_EDGE_AS_GTS_EDGE(self),s));
}


static PyObject*
belongs_to_tetrahedron(PygtsEdge *self, PyObject *args)
{
  SELF_CHECK

  if(gts_edge_belongs_to_tetrahedron(PYGTS_EDGE_AS_GTS_EDGE(self))) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
is_boundary(PygtsEdge *self, PyObject *args)
{
  PyObject *s_;
  GtsSurface *s;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &s_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_surface_check(s_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Surface");
    return NULL;
  }
  s = PYGTS_SURFACE_AS_GTS_SURFACE(s_);

  /* Make the call and return */
  if(gts_edge_is_boundary(PYGTS_EDGE_AS_GTS_EDGE(self),s)!=NULL) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
contacts(PygtsEdge *self, PyObject *args)
{
  SELF_CHECK

  return Py_BuildValue("i",
		       gts_edge_is_contact(PYGTS_EDGE_AS_GTS_EDGE(self)));
}


/* Methods table */
static PyMethodDef methods[] = {
  {"is_ok", (PyCFunction)is_ok,
   METH_NOARGS,
   "True if this Edge e is not degenerate or duplicate.\n"
   "False otherwise.  Degeneracy implies e.v1.id == e.v2.id.\n"
   "\n"
   "Signature: e.is_ok()\n"
  },  

  {"is_unattached", (PyCFunction)is_unattached,
   METH_NOARGS,
   "True if this Edge e is not part of any Triangle.\n"
   "\n"
   "Signature: e.is_unattached()\n"
  },

/* Edge replace() method works but results in Triangles that are "not ok";
 * i.e., they have edges that don't connect.  We don't want that problem
 * and so the method has been disabled.
 */
/*
  {"replace", (PyCFunction)replace,
   METH_VARARGS,
   "Replaces this Edge e1 with Edge e2 in all Triangles that have e1.\n"
   "Edge e1 itself is left unchanged.\n"
   "\n"
   "Signature: e1.replace(e2).\n"
  },
*/

  {"face_number", (PyCFunction)face_number,
   METH_VARARGS,
   "Returns number of faces using this Edge e on Surface s.\n"
   "\n"
   "Signature: e.face_number(s)\n"
  },

  {"is_boundary", (PyCFunction)is_boundary,
   METH_VARARGS,
   "Returns True if this Edge e is a boundary on Surface s.\n"
   "Otherwise False.\n"
   "\n"
   "Signature: e.is_boundary(s)\n"
  },

  {"belongs_to_tetrahedron", (PyCFunction)belongs_to_tetrahedron,
   METH_NOARGS,
   "Returns True if this Edge e belongs to a tetrahedron.\n"
   "Otherwise False.\n"
   "\n"
   "Signature: e.belongs_to_tetrahedron()\n"
  },

  {"contacts", (PyCFunction)contacts,
   METH_NOARGS,
   "Returns number of sets of connected triangles share this Edge e\n"
   "as a contact Edge.\n"
   "\n"
   "Signature: e.contacts()\n"
  },

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Python type methods */

static GtsObject* parent(GtsEdge *e1);

static PyObject *
new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyObject *o;
  PygtsObject *obj;
  GtsEdge *tmp;
  GtsObject *edge=NULL;
  PyObject *v1_,*v2_;
  PygtsVertex *v1,*v2;
  guint alloc_gtsobj = TRUE;
  guint N;

  /* Parse the args */
  if(kwds) {
    o = PyDict_GetItemString(kwds,"alloc_gtsobj");
    if(o==Py_False) {
      alloc_gtsobj = FALSE;
    }
    if(o!=NULL) {
      PyDict_DelItemString(kwds, "alloc_gtsobj");
    }
  }
  if(kwds) {
    Py_INCREF(Py_False);
    PyDict_SetItemString(kwds,"alloc_gtsobj", Py_False);
  }

  /* Allocate the gtsobj (if needed) */
  if( alloc_gtsobj ) {

    /* Parse the args */
    if( (N = PyTuple_Size(args)) < 2 ) {
      PyErr_SetString(PyExc_TypeError,"expected two Vertices");
      return NULL;
    }
    v1_ = PyTuple_GET_ITEM(args,0);
    v2_ = PyTuple_GET_ITEM(args,1);

    /* Convert to PygtsObjects */
    if(!pygts_vertex_check(v1_)) {
      PyErr_SetString(PyExc_TypeError,"expected two Vertices");
      return NULL;
    }
    if(!pygts_vertex_check(v2_)) {
      PyErr_SetString(PyExc_TypeError,"expected two Vertices");
      return NULL;
    }
    v1 = PYGTS_VERTEX(v1_);
    v2 = PYGTS_VERTEX(v2_);

    /* Error check */
    if(PYGTS_OBJECT(v1)->gtsobj == PYGTS_OBJECT(v2)->gtsobj) {
      PyErr_SetString(PyExc_ValueError,"Vertices given are the same");
      return NULL;
    }

    /* Create the GtsEdge */
    edge = GTS_OBJECT(gts_edge_new(gts_edge_class(),
				   GTS_VERTEX(v1->gtsobj),
				   GTS_VERTEX(v2->gtsobj)));
    if( edge == NULL )  {
      PyErr_SetString(PyExc_MemoryError, "could not create Edge");
      return NULL;
    }

    /* Check for duplicate */
    tmp = gts_edge_is_duplicate(GTS_EDGE(edge));
    if( tmp != NULL ) {
      gts_object_destroy(edge);
      edge = GTS_OBJECT(tmp);
    }

    /* If corresponding PyObject found in object table, we are done */
    if( (obj=g_hash_table_lookup(obj_table,edge)) != NULL ) {
      Py_INCREF(obj);
      return (PyObject*)obj;
    }
  }

  /* Chain up */
  obj = PYGTS_OBJECT(PygtsSegmentType.tp_new(type,args,kwds));

  if( alloc_gtsobj ) {
    obj->gtsobj = edge;

    /* Create a parent GtsTriangle */
    if( (obj->gtsobj_parent = parent(GTS_EDGE(obj->gtsobj))) == NULL ) {
      gts_object_destroy(obj->gtsobj);
      obj->gtsobj = NULL;
      return NULL;
    }

    pygts_object_register(PYGTS_OBJECT(obj));
  }

  return (PyObject*)obj;
}


static int
init(PygtsEdge *self, PyObject *args, PyObject *kwds)
{
  gint ret;

  /* Chain up */
  if( (ret=PygtsSegmentType.tp_init((PyObject*)self,args,kwds)) != 0 ){
    return ret;
  }

  return 0;
}


/* Methods table */
PyTypeObject PygtsEdgeType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "gts.Edge",              /* tp_name */
    sizeof(PygtsEdge),       /* tp_basicsize */
    0,                       /* tp_itemsize */
    0,                       /* tp_dealloc */
    0,                       /* tp_print */
    0,                       /* tp_getattr */
    0,                       /* tp_setattr */
    0,                       /* tp_compare */
    0,                       /* tp_repr */
    0,                       /* tp_as_number */
    0,                       /* tp_as_sequence */
    0,                       /* tp_as_mapping */
    0,                       /* tp_hash */
    0,                       /* tp_call */
    0,                       /* tp_str */
    0,                       /* tp_getattro */
    0,                       /* tp_setattro */
    0,                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
      Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Edge object",           /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    methods,                 /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    (initproc)init,          /* tp_init */
    0,                       /* tp_alloc */
    (newfunc)new             /* tp_new */
};


/*-------------------------------------------------------------------------*/
/* Pygts functions */

gboolean 
pygts_edge_check(PyObject* o)
{
  if(! PyObject_TypeCheck(o, &PygtsEdgeType)) {
    return FALSE;
  }
  else {
#if PYGTS_DEBUG
    return pygts_edge_is_ok(PYGTS_EDGE(o));
#else
    return TRUE;
#endif
  }
}

gboolean
pygts_edge_is_ok(PygtsEdge *e)
{
  GSList *parent;
  PygtsObject *obj;

  obj = PYGTS_OBJECT(e);

  if(!pygts_segment_is_ok(PYGTS_SEGMENT(e))) return FALSE;

  /* Check for a valid parent */
  g_return_val_if_fail(obj->gtsobj_parent!=NULL,FALSE);
  g_return_val_if_fail(PYGTS_IS_PARENT_TRIANGLE(obj->gtsobj_parent),FALSE);
  parent = g_slist_find(GTS_EDGE(obj->gtsobj)->triangles,
			obj->gtsobj_parent);
  g_return_val_if_fail(parent!=NULL,FALSE);

  return TRUE;
}


static GtsObject* 
parent(GtsEdge *e1) {
  GtsVertex *v1,*v2,*v3;
  GtsPoint *p1,*p2;
  GtsEdge *e2, *e3;
  GtsTriangle *p;

  /* Create a third vertex for the triangle */
  v1 = GTS_SEGMENT(e1)->v1;
  v2 = GTS_SEGMENT(e1)->v2;
  p1 = GTS_POINT(v1);
  p2 = GTS_POINT(v2);
  v3 = gts_vertex_new(pygts_parent_vertex_class(),
		      p1->x+p2->x,p1->y+p2->y,p1->z+p2->z);
  if( v3 == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create Vertex");
    return NULL;
  }

  /* Create another two edges */
  if( (e2 = gts_edge_new(pygts_parent_edge_class(),v2,v3)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create Edge");
    return NULL;
  }
  if( (e3 = gts_edge_new(pygts_parent_edge_class(),v3,v1)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create Edge");
    gts_object_destroy(GTS_OBJECT(e2));
    return NULL;
  }

  /* Create and return the parent */
  if( (p = gts_triangle_new(pygts_parent_triangle_class(),e1,e2,e3)) 
      == NULL ) {
    gts_object_destroy(GTS_OBJECT(e2));
    gts_object_destroy(GTS_OBJECT(e3));
    PyErr_SetString(PyExc_MemoryError, "could not create Triangle");
    return NULL;
  }

  return GTS_OBJECT(p);
}


PygtsEdge *
pygts_edge_new(GtsEdge *e)
{
  PyObject *args, *kwds;
  PygtsObject *edge;

  /* Check for Edge in the object table */
  if( (edge = PYGTS_OBJECT(g_hash_table_lookup(obj_table,GTS_OBJECT(e)))) 
      !=NULL ) {
    Py_INCREF(edge);
    return PYGTS_EDGE(edge);
  }

  /* Build a new Edge */
  args = Py_BuildValue("OO",Py_None,Py_None);
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_False);
  edge = PYGTS_EDGE(PygtsEdgeType.tp_new(&PygtsEdgeType, args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( edge == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Edge");
    return NULL;
  }
  edge->gtsobj = GTS_OBJECT(e);

  /* Attach the parent */
  if( (edge->gtsobj_parent = parent(e)) == NULL ) {
    Py_DECREF(edge);
    return NULL;
  }

  /* Register and return */
  pygts_object_register(edge);
  return PYGTS_EDGE(edge);
}


GtsTriangleClass*
pygts_parent_triangle_class(void)
{
  static GtsTriangleClass *klass = NULL;
  GtsObjectClass *super = NULL;

  if (klass == NULL) {

    super = GTS_OBJECT_CLASS(gts_triangle_class());

    GtsObjectClassInfo pygts_parent_triangle_info = {
      "PygtsParentTriangle",
      sizeof(PygtsParentTriangle),
      sizeof(GtsTriangleClass),
      (GtsObjectClassInitFunc)(super->info.class_init_func),
      (GtsObjectInitFunc)(super->info.object_init_func),
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new(gts_object_class(),
				 &pygts_parent_triangle_info);
  }

  return klass;
}


GtsEdgeClass*
pygts_parent_edge_class(void)
{
  static GtsEdgeClass *klass = NULL;
  GtsObjectClass *super = NULL;

  if (klass == NULL) {

    super = GTS_OBJECT_CLASS(pygts_parent_segment_class());

    GtsObjectClassInfo pygts_parent_edge_info = {
      "PygtsParentEdge",
      sizeof(PygtsParentEdge),
      sizeof(GtsEdgeClass),
      (GtsObjectClassInitFunc)(super->info.class_init_func),
      (GtsObjectInitFunc)(super->info.object_init_func),
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new(gts_object_class(),
				 &pygts_parent_edge_info);
  }

  return klass;
}
