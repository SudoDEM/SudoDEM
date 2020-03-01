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
  #define SELF_CHECK if(!pygts_face_check((PyObject*)self)) {         \
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
is_ok(PygtsFace *self, PyObject *args)
{
  if(pygts_face_is_ok(self)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
is_unattached(PygtsFace *self, PyObject *args)
{
  guint n;

  /* Check for attachments other than to the gtsobj_parent */
  n = g_slist_length(PYGTS_FACE_AS_GTS_FACE(self)->surfaces);
  if( n > 1 ) {
    Py_INCREF(Py_False);
    return Py_False;
  }
  else if( n == 1 ){
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    PyErr_SetString(PyExc_RuntimeError, "Face lost parent (internal error)");
    return NULL;
  }
}


static PyObject*
neighbor_number(PygtsFace *self, PyObject *args)
{
  PyObject *s_=NULL;
  PygtsSurface *s=NULL;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &s_) )
    return NULL;

  /* Convert to PygtsObjects */
  if( pygts_surface_check(s_) ) {
    s = PYGTS_SURFACE(s_);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "expected a Surface");
      return NULL;
  }

  return Py_BuildValue("i",
      gts_face_neighbor_number(PYGTS_FACE_AS_GTS_FACE(self),
			       PYGTS_SURFACE_AS_GTS_SURFACE(s)));
}


static PyObject*
neighbors(PygtsFace *self, PyObject *args)
{
  PyObject *s_=NULL;
  PygtsSurface *s=NULL;
  guint i,N;
  PyObject *tuple;
  GSList *faces,*f;
  PygtsFace *face;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &s_) )
    return NULL;

  /* Convert to PygtsObjects */
  if( pygts_surface_check(s_) ) {
    s = PYGTS_SURFACE(s_);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "expected a Surface");
      return NULL;
  }

  N = gts_face_neighbor_number(PYGTS_FACE_AS_GTS_FACE(self),
			       PYGTS_SURFACE_AS_GTS_SURFACE(s));

  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Could not create tuple");
    return NULL;
  }

  /* Get the neighbors */
  faces = gts_face_neighbors(PYGTS_FACE_AS_GTS_FACE(self),
			     PYGTS_SURFACE_AS_GTS_SURFACE(s));
  f = faces;
			     
  for(i=0;i<N;i++) {
    if( (face = pygts_face_new(GTS_FACE(f->data))) == NULL ) {
      Py_DECREF(tuple);
      return NULL;
    }
    PyTuple_SET_ITEM(tuple, i, (PyObject*)face);
    f = g_slist_next(f);
  }  

  return (PyObject*)tuple;
}


static PyObject*
is_compatible(PygtsFace *self, PyObject *args)
{
  PyObject *o1_=NULL;
  GtsEdge *e=NULL;
  PygtsTriangle *t=NULL;
  PygtsSurface *s=NULL;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &o1_) )
    return NULL;

  /* Convert to PygtsObjects */
  if( pygts_triangle_check(o1_) ) {
    t = PYGTS_TRIANGLE(o1_);
  }
  else {
    if( pygts_surface_check(o1_) ) {
      s = PYGTS_SURFACE(o1_);
    }
    else {
      PyErr_SetString(PyExc_TypeError, "expected a Triangle or Surface");
      return NULL;
    }
  }

  if(t!=NULL) {
    if( (e = gts_triangles_common_edge(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				       PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t))) 
	== NULL ) {
      PyErr_SetString(PyExc_RuntimeError, "Faces do not share common edge");
      return NULL;
    }
    if(gts_triangles_are_compatible(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				    PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t),
				    e)==TRUE) {

      Py_INCREF(Py_True);
      return Py_True;
    }
  }
  else {
    if(gts_face_is_compatible(PYGTS_FACE_AS_GTS_FACE(self),
			      PYGTS_SURFACE_AS_GTS_SURFACE(s))==TRUE) {
      Py_INCREF(Py_True);
      return Py_True;
    }
  }
  Py_INCREF(Py_False);
  return Py_False;
}


static PyObject*
is_on(PygtsFace *self, PyObject *args)
{
  PyObject *s_=NULL;
  PygtsSurface *s=NULL;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &s_) )
    return NULL;

  /* Convert to PygtsObjects */
  if( pygts_surface_check(s_) ) {
    s = PYGTS_SURFACE(s_);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "expected a Surface");
      return NULL;
  }

  if( gts_face_has_parent_surface(PYGTS_FACE_AS_GTS_FACE(self),
				  PYGTS_SURFACE_AS_GTS_SURFACE(s)) ) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


/* Methods table */
static PyMethodDef methods[] = {
  {"is_ok", (PyCFunction)is_ok,
   METH_NOARGS,
   "True if this Face f is non-degenerate and non-duplicate.\n"
   "False otherwise.\n"
   "\n"
   "Signature: f.is_ok()\n"
  },  

  {"is_unattached", (PyCFunction)is_unattached,
   METH_NOARGS,
   "True if this Face f is not part of any Surface.\n"
   "\n"
   "Signature: f.is_unattached().\n"
  },

  {"neighbor_number", (PyCFunction)neighbor_number,
   METH_VARARGS,
   "Returns the number of neighbors of Face f belonging to Surface s.\n"
   "\n"
   "Signature: f.neighbor_number(s).\n"
  },

  {"neighbors", (PyCFunction)neighbors,
   METH_VARARGS,
   "Returns a tuple of neighbors of this Face f belonging to Surface s.\n"
   "\n"
   "Signature: f.neighbors(s).\n"
  },

  {"is_compatible", (PyCFunction)is_compatible,
   METH_VARARGS,
   "True if Face f is compatible with all neighbors in Surface s.\n"
   "False otherwise.\n"
   "\n"
   "Signature: f.is_compatible(s).\n"
  },  

  {"is_on", (PyCFunction)is_on,
   METH_VARARGS,
   "True if this Face f is on Surface s.  False otherwise.\n"
   "\n"
   "Signature: f.is_on(s).\n"
  },  

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Python type methods */


static GtsObject * parent(GtsFace *face);

static PyObject *
new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyObject *o;
  PygtsObject *obj;
  guint alloc_gtsobj = TRUE;
  PyObject *o1_,*o2_,*o3_;
  GtsVertex *v1=NULL, *v2=NULL, *v3=NULL;
  GtsEdge *e1=NULL,*e2=NULL,*e3=NULL,*e;
  GtsSegment *s1,*s2,*s3;
  gboolean flag=FALSE;  /* Flag when the args are gts.Point objects */
  GtsFace *f;
  GtsTriangle *t;
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
    if( (N = PyTuple_Size(args)) < 3 ) {
      PyErr_SetString(PyExc_TypeError,"expected three Edges or three Vertices");
      return NULL;
    }
    o1_ = PyTuple_GET_ITEM(args,0);
    o2_ = PyTuple_GET_ITEM(args,1);
    o3_ = PyTuple_GET_ITEM(args,2);

    /* Convert to PygtsObjects */
    if( pygts_edge_check(o1_) ) {
      e1 = PYGTS_EDGE_AS_GTS_EDGE(o1_);
    }
    else {
      if( pygts_vertex_check(o1_) ) {
	v1 = PYGTS_VERTEX_AS_GTS_VERTEX(o1_);
	flag = TRUE;
      }
    }

    if( pygts_edge_check(o2_) ) {
      e2 = PYGTS_EDGE_AS_GTS_EDGE(o2_);
    }
    else {
      if( pygts_vertex_check(o2_) ) {
	v2 = PYGTS_VERTEX_AS_GTS_VERTEX(o2_);
	flag = TRUE;
      }
    }

    if( pygts_edge_check(o3_) ) {
      e3 = PYGTS_EDGE_AS_GTS_EDGE(o3_);
    }
    else {
      if(pygts_vertex_check(o3_)) {
	v3 = PYGTS_VERTEX_AS_GTS_VERTEX(o3_);
	flag = TRUE;
      }
    }
    
    /* Check for three edges or three vertices */
    if( !((e1!=NULL && e2!=NULL && e3!=NULL) ||
	  (v1!=NULL && v2!=NULL && v3!=NULL)) ) {
      PyErr_SetString(PyExc_TypeError,
		      "three Edge or three Vertex objects expected");
      return NULL;
    }

    if(flag) {

      /* Create gts edges */
      if( (e1 = gts_edge_new(gts_edge_class(),v1,v2)) == NULL ) {
	PyErr_SetString(PyExc_MemoryError, "could not create Edge");
	return NULL;
      }
      if( (e2 = gts_edge_new(gts_edge_class(),v2,v3)) == NULL ) {
	PyErr_SetString(PyExc_MemoryError, "could not create Edge");
	gts_object_destroy(GTS_OBJECT(e1));
	return NULL;
      }
      if( (e3 = gts_edge_new(gts_edge_class(),v3,v1)) == NULL ) {
	PyErr_SetString(PyExc_MemoryError, "could not create Edge");
	gts_object_destroy(GTS_OBJECT(e1));
	gts_object_destroy(GTS_OBJECT(e2));
	return NULL;
      }
      
      /* Check for duplicates */
      if( (e = gts_edge_is_duplicate(e1)) != NULL ) {
	gts_object_destroy(GTS_OBJECT(e1));
	e1 = e;
      }
      if( (e = gts_edge_is_duplicate(e2)) != NULL ) {
	gts_object_destroy(GTS_OBJECT(e2));
	e2 = e;
      }
      if( (e = gts_edge_is_duplicate(e3)) != NULL ) {
	gts_object_destroy(GTS_OBJECT(e3));
	e3 = e;
      }
    }
  
    /* Check that edges connect */
    s1 = GTS_SEGMENT(e1);
    s2 = GTS_SEGMENT(e2);
    s3 = GTS_SEGMENT(e3);
    if( !((s1->v1==s3->v2 && s1->v2==s2->v1 && s2->v2==s3->v1) ||
	  (s1->v1==s3->v2 && s1->v2==s2->v2 && s2->v1==s3->v1) ||
	  (s1->v1==s3->v1 && s1->v2==s2->v1 && s2->v2==s3->v2) ||
	  (s1->v2==s3->v2 && s1->v1==s2->v1 && s2->v2==s3->v1) ||
	  (s1->v1==s3->v1 && s1->v2==s2->v2 && s2->v1==s3->v2) ||
	  (s1->v2==s3->v2 && s1->v1==s2->v2 && s2->v1==s3->v1) ||
	  (s1->v2==s3->v1 && s1->v1==s2->v1 && s2->v2==s3->v2) ||
	  (s1->v2==s3->v1 && s1->v1==s2->v2 && s2->v1==s3->v2)) ) {
      PyErr_SetString(PyExc_RuntimeError, "Edges in face must connect");
      if(!g_hash_table_lookup(obj_table,GTS_OBJECT(e1))) {
	gts_object_destroy(GTS_OBJECT(e1));
      }
      if(!g_hash_table_lookup(obj_table,GTS_OBJECT(e1))) {
	gts_object_destroy(GTS_OBJECT(e2));
      }
      if(!g_hash_table_lookup(obj_table,GTS_OBJECT(e1))) {
	gts_object_destroy(GTS_OBJECT(e3));
      }
      return NULL;
    }

    /* Create the GtsFace */
    if( (f = gts_face_new(gts_face_class(),e1,e2,e3)) == NULL )  {
      PyErr_SetString(PyExc_MemoryError, "could not create Face");
      if(!g_hash_table_lookup(obj_table,GTS_OBJECT(e1))) {
	gts_object_destroy(GTS_OBJECT(e1));
      }
      if(!g_hash_table_lookup(obj_table,GTS_OBJECT(e1))) {
	gts_object_destroy(GTS_OBJECT(e2));
      }
      if(!g_hash_table_lookup(obj_table,GTS_OBJECT(e1))) {
	gts_object_destroy(GTS_OBJECT(e3));
      }
      return NULL;
    }

    /* Check for duplicate */
    t = gts_triangle_is_duplicate(GTS_TRIANGLE(f));
    if( t != NULL ) {
      gts_object_destroy(GTS_OBJECT(f));
      if(!GTS_IS_FACE(t)) {
	PyErr_SetString(PyExc_TypeError, "expected a Face (internal error)");
      }
      f = GTS_FACE(t);
    }

    /* If corresponding PyObject found in object table, we are done */
    if( (obj=g_hash_table_lookup(obj_table,GTS_OBJECT(f))) != NULL ) {
      Py_INCREF(obj);
      return (PyObject*)obj;
    }
  }
  
  /* Chain up */
  obj = PYGTS_OBJECT(PygtsTriangleType.tp_new(type,args,kwds));

  if( alloc_gtsobj ) {

    obj->gtsobj = GTS_OBJECT(f);

    /* Create the parent GtsSurface */
    if( (obj->gtsobj_parent = parent(GTS_FACE(obj->gtsobj))) == NULL ) {
      gts_object_destroy(obj->gtsobj);
      obj->gtsobj = NULL;
      return NULL;
    }

    pygts_object_register(PYGTS_OBJECT(obj));
  }

  return (PyObject*)obj;
}


static int
init(PygtsFace *self, PyObject *args, PyObject *kwds)
{
  gint ret;

  /* Chain up */
  if( (ret=PygtsTriangleType.tp_init((PyObject*)self,args,kwds)) != 0 ){
    return ret;
  }

  return 0;
}


/* Methods table */
PyTypeObject PygtsFaceType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "gts.Face",              /* tp_name */
    sizeof(PygtsFace),       /* tp_basicsize */
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
    "Face object",           /* tp_doc */
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
pygts_face_check(PyObject* o)
{
  if(! PyObject_TypeCheck(o, &PygtsFaceType)) {
    return FALSE;
  }
  else {
#if PYGTS_DEBUG
    return pygts_face_is_ok(PYGTS_FACE(o));
#else
    return TRUE;
#endif
  }
}


gboolean 
pygts_face_is_ok(PygtsFace *f)
{
  GSList *parent;
  PygtsObject *obj;

  obj = PYGTS_OBJECT(f);

  if(!pygts_triangle_is_ok(PYGTS_TRIANGLE(f))) return FALSE;

  /* Check for a valid parent */
  g_return_val_if_fail(obj->gtsobj_parent!=NULL,FALSE);
  g_return_val_if_fail(GTS_IS_SURFACE(obj->gtsobj_parent),FALSE);
  parent = g_slist_find(GTS_FACE(obj->gtsobj)->surfaces,
			obj->gtsobj_parent);
  g_return_val_if_fail(parent!=NULL,FALSE);

  return TRUE;
}


static GtsObject *
parent(GtsFace *face) {
  GtsSurface *p;

  p = gts_surface_new(gts_surface_class(), gts_face_class(),
		      gts_edge_class(), gts_vertex_class());

  if( p == NULL )  {
    PyErr_SetString(PyExc_MemoryError, "could not create parent");
    return NULL;
  }
  gts_surface_add_face(p,face);

  return GTS_OBJECT(p);
}


PygtsFace *
pygts_face_new(GtsFace *f)
{
  PyObject *args, *kwds;
  PygtsObject *face;

  /* Check for Face in the object table */
  if( (face=PYGTS_OBJECT(g_hash_table_lookup(obj_table,GTS_OBJECT(f))))
      != NULL ) {
    Py_INCREF(face);
    return PYGTS_FACE(face);
  }

  /* Build a new Face */
  args = Py_BuildValue("OOO",Py_None,Py_None,Py_None);
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_False);
  face = PYGTS_OBJECT(PygtsFaceType.tp_new(&PygtsFaceType, args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( face == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create Face");
    return NULL;
  }
  face->gtsobj = GTS_OBJECT(f);

  /* Attach the parent */
  if( (face->gtsobj_parent = parent(f)) == NULL ) {
    Py_DECREF(face);
    return NULL;
  }
  
  /* Register and return */
  pygts_object_register(face);
  return PYGTS_FACE(face);
}
