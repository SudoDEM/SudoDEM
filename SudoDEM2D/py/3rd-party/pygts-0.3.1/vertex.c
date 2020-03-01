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
  #define SELF_CHECK if(!pygts_vertex_check((PyObject*)self)) {      \
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
is_ok(PygtsVertex *self, PyObject *args)
{
  if(pygts_vertex_is_ok(self)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
is_unattached(PygtsVertex *self, PyObject *args)
{
  guint n;

  SELF_CHECK

  /* Check for attachments other than to the gtsobj_parent */
  n = g_slist_length(PYGTS_VERTEX_AS_GTS_VERTEX(self)->segments);
  if( n > 1 ) {
    Py_INCREF(Py_False);
    return Py_False;
  }
  else {
    Py_INCREF(Py_True);
    return Py_True;
  }
}


static PyObject*
is_boundary(PygtsVertex* self, PyObject *args)
{
  PyObject *s_;
  PygtsObject *s;

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
  s = PYGTS_OBJECT(s_);

  if( gts_vertex_is_boundary(PYGTS_VERTEX_AS_GTS_VERTEX(self),
			     PYGTS_SURFACE_AS_GTS_SURFACE(s)) ) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
contacts(PygtsVertex* self, PyObject *args)
{
  PyObject *sever_=NULL;
  gboolean sever=FALSE;
  guint n;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "|O", &sever_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if( sever_ != NULL ) {
    if(!PyBool_Check(sever_)) {
      PyErr_SetString(PyExc_TypeError,"expected a Boolean");
      return NULL;
    }
    if( sever_ == Py_True ) {
      sever = TRUE;
    }
  }

  n = gts_vertex_is_contact(PYGTS_VERTEX_AS_GTS_VERTEX(self),sever);
  return Py_BuildValue("i",n);
}


static PyObject*
is_connected(PygtsVertex *self, PyObject *args)
{
  PyObject *v_;
  PygtsVertex *v;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &v_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_vertex_check(v_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Vertex");
    return NULL;
  }
  v = PYGTS_VERTEX(v_);

  if( gts_vertices_are_connected(PYGTS_VERTEX_AS_GTS_VERTEX(self),
				 PYGTS_VERTEX_AS_GTS_VERTEX(v)) != NULL ) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
replace(PygtsVertex *self, PyObject *args)
{
  PyObject *p2_;
  PygtsVertex *p2;
  GSList *parents=NULL, *i, *cur;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &p2_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_vertex_check(p2_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Vertex");
    return NULL;
  }
  p2 = PYGTS_VERTEX(p2_);

  if( self != p2 ) {
    /* (Ignore self-replacement) */

    /* Detach and save any parent segments */
    i = PYGTS_VERTEX_AS_GTS_VERTEX(self)->segments;
    while(i!=NULL) {
      cur = i;
      i = g_slist_next(i);
      if(PYGTS_IS_PARENT_SEGMENT(cur->data)) {
	PYGTS_VERTEX_AS_GTS_VERTEX(self)->segments = 
	  g_slist_remove_link(PYGTS_VERTEX_AS_GTS_VERTEX(self)->segments,
			      cur);
	parents = g_slist_prepend(parents,cur->data);
	g_slist_free_1(cur);
      }
    }

    /* Perform the replace operation */
    gts_vertex_replace(PYGTS_VERTEX_AS_GTS_VERTEX(self),
		       PYGTS_VERTEX_AS_GTS_VERTEX(p2));

    /* Reattach the parent segments */
    i = parents;
    while(i!=NULL) {
      PYGTS_VERTEX_AS_GTS_VERTEX(self)->segments = 
	g_slist_prepend(PYGTS_VERTEX_AS_GTS_VERTEX(self)->segments,i->data);
      i = g_slist_next(i);
    }
    g_slist_free(parents);
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
neighbors(PygtsVertex* self, PyObject *args)
{
  PyObject *s_=NULL;
  GtsSurface *s=NULL;
  GSList *vertices,*v;
  PygtsVertex *vertex;
  PyObject *tuple;
  guint n,N;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "|O", &s_) ) {
    return NULL;
  }

  /* Convert */
  if( s_ != NULL ) {
    if(!pygts_surface_check(s_)) {
      PyErr_SetString(PyExc_TypeError,"expected a Surface");
      return NULL;
    }
    s = PYGTS_SURFACE_AS_GTS_SURFACE(s_);
  }

  /* Get the neighbors */
  vertices = gts_vertex_neighbors(PYGTS_VERTEX_AS_GTS_VERTEX(self),
				  NULL,s);
  N = g_slist_length(vertices);

  /* Create the tuple */
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create tuple");
    return NULL;
  }

  /* Put PygtsVertex objects into the tuple */
  v = vertices;
  for(n=0;n<N;n++) {

    /* Skip this vertex if it is a parent */
    while( v!=NULL && PYGTS_IS_PARENT_VERTEX(GTS_VERTEX(v->data)) ) {
      v = g_slist_next(v);      
    }
    if( v==NULL ) break;

    if( (vertex = pygts_vertex_new(GTS_VERTEX(v->data))) == NULL ) {
      Py_DECREF((PyObject*)tuple);
      return NULL;
    }

    PyTuple_SET_ITEM(tuple, n, (PyObject*)vertex);
    
    v = g_slist_next(v);
  }

  if(_PyTuple_Resize(&tuple,n)!=0) {
    Py_DECREF(tuple);
    return NULL;
  }

  return tuple;
}


static PyObject*
faces(PygtsVertex* self, PyObject *args)
{
  PyObject *s_=NULL;
  GtsSurface *s=NULL;
  GSList *faces,*f;
  PygtsFace *face;
  PyObject *tuple;
  guint n,N;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "|O", &s_) ) {
    return NULL;
  }

  /* Convert */
  if( s_ != NULL ) {
    if(!pygts_surface_check(s_)) {
      PyErr_SetString(PyExc_TypeError,"expected a Surface");
      return NULL;
    }
    s = PYGTS_SURFACE_AS_GTS_SURFACE(s_);
  }

  /* Get the faces */
  faces = gts_vertex_faces(PYGTS_VERTEX_AS_GTS_VERTEX(self),s,NULL);
  N = g_slist_length(faces);

  /* Create the tuple */
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"expected a tuple");
    return NULL;
  }

  /* Put PygtsVertex objects into the tuple */
  f = faces;
  for(n=0;n<N;n++) {

    if( (face = pygts_face_new(GTS_FACE(f->data))) == NULL ) {
      Py_DECREF(tuple);
      return NULL;
    }

    PyTuple_SET_ITEM(tuple, n, (PyObject*)face);
    
    f = g_slist_next(f);
  }

  return tuple;
}


static PyObject*
encroaches(PygtsVertex *self, PyObject *args)
{
  PyObject *e_;
  PygtsEdge *e;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &e_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_edge_check(e_)) {
    PyErr_SetString(PyExc_TypeError,"expected an Edge");
    return NULL;
  }
  e = PYGTS_EDGE(e_);

  if(gts_vertex_encroaches_edge(PYGTS_VERTEX_AS_GTS_VERTEX(self),
				PYGTS_EDGE_AS_GTS_EDGE(e))) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
triangles(PygtsVertex *self, PyObject *args)
{
  GSList *triangles, *t;
  PygtsTriangle *triangle;
  guint i,N;
  PyObject *tuple;

  SELF_CHECK

  triangles = gts_vertex_triangles(PYGTS_VERTEX_AS_GTS_VERTEX(self),NULL);
  N = g_slist_length(triangles);

  /* Create the tuple */
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create tuple");
    return NULL;
  }

  /* Put PygtsVertex objects into the tuple */
  t = triangles;
  for(i=0;i<N;i++) {

    if( (triangle = pygts_triangle_new(GTS_TRIANGLE(t->data))) == NULL ) {
      Py_DECREF(tuple);
      return NULL;
    }

    PyTuple_SET_ITEM(tuple, i, (PyObject*)triangle);
    
    t = g_slist_next(t);
  }

  return tuple;
}


/* Methods table */
static PyMethodDef methods[] = {

  {"is_ok", (PyCFunction)is_ok,
   METH_NOARGS,
   "True if this Vertex v is OK.  False otherwise.\n"
   "This method is useful for unit testing and debugging.\n"
   "\n"
   "Signature: v.is_ok().\n"
  },  

  {"is_unattached", (PyCFunction)is_unattached,
   METH_NOARGS,
   "True if this Vertex v is not the endpoint of any Segment.\n"
   "\n"
   "Signature: v.is_unattached().\n"
  },

  {"is_boundary", (PyCFunction)is_boundary,
   METH_VARARGS,
   "True if this Vertex v is used by a boundary Edge of Surface s.\n"
   "\n"
   "Signature: v.is_boundary().\n"
  },

  {"contacts", (PyCFunction)contacts,
   METH_VARARGS,
   "Returns the number of sets of connected Triangles sharing this\n"
   "Vertex v.\n"
   "\n"
   "Signature: v.contacts().\n"
   "\n"
   "If sever is True (default: False) and v is a contact vertex then\n"
   "the vertex is replaced in each Triangle with clones.\n"
  },

  {"is_connected", (PyCFunction)is_connected,
   METH_VARARGS,
   "Return True if this Vertex v1 is connected to Vertex v2\n"
   "by a Segment.\n"
   "\n"
   "Signature: v1.is_connected().\n"
  },

  {"replace", (PyCFunction)replace,
   METH_VARARGS,
   "Replaces this Vertex v1 with Vertex v2 in all Segments that have v1.\n"
   "Vertex v1 itself is left unchanged.\n"
   "\n"
   "Signature: v1.replace(v2).\n"
  },

  {"neighbors", (PyCFunction)neighbors,
   METH_VARARGS,
   "Returns a tuple of Vertices attached to this Vertex v\n"
   "by a Segment.\n"
   "\n"
   "If a Surface s is given, only Vertices on s are considered.\n"
   "\n"
   "Signature: v.neighbors() or v.neighbors(s).\n"
  },

  {"faces", (PyCFunction)faces,
   METH_VARARGS,
   "Returns a tuple of Faces that have this Vertex v.\n"
   "\n"
   "If a Surface s is given, only Vertices on s are considered.\n"
   "\n"
   "Signature: v.faces() or v.faces(s).\n"
  },

  {"encroaches", (PyCFunction)encroaches,
   METH_VARARGS,
   "Returns True if this Vertex v is strictly contained in the\n"
   "diametral circle of Edge e.  False otherwise.\n"
   "\n"
   "Only the projection onto the x-y plane is considered.\n"
   "\n"
   "Signature: v.encroaches(e)\n"
  },

  {"triangles", (PyCFunction)triangles,
   METH_NOARGS,
   "Returns a list of Triangles that have this Vertex v.\n"
   "\n"
   "Signature: v.triangles()\n"
  },

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Python type methods */

static GtsObject * parent(GtsVertex *v1);

static PyObject *
new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyObject *o;
  PygtsObject *obj;
  guint alloc_gtsobj = TRUE;

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

  /* Chain up */
  obj = PYGTS_OBJECT(PygtsPointType.tp_new(type,args,kwds));

  /* Allocate the gtsobj (if needed) */
  if( alloc_gtsobj ) {
    obj->gtsobj = GTS_OBJECT(gts_vertex_new(gts_vertex_class(),0,0,0));
    if( obj->gtsobj == NULL )  {
      PyErr_SetString(PyExc_MemoryError, "could not create Vertex");
      return NULL;
    }

    /* Create the parent GtsSegment */
    if( (obj->gtsobj_parent=parent(GTS_VERTEX(obj->gtsobj))) == NULL ) {
      gts_object_destroy(obj->gtsobj);
      obj->gtsobj = NULL;
      return NULL;
    }

    pygts_object_register(obj);
  }

  return (PyObject*)obj;
}


static int
init(PygtsVertex *self, PyObject *args, PyObject *kwds)
{
  gint ret;

  /* Chain up */
  if( (ret=PygtsPointType.tp_init((PyObject*)self,args,kwds)) != 0 ) {
    return ret;
  }

  return 0;
}


/* Methods table */
PyTypeObject PygtsVertexType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "gts.Vertex",            /* tp_name */
    sizeof(PygtsVertex),     /* tp_basicsize */
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
    "Vertex object",         /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    methods,                 /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base: attached in pygts.c */
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
pygts_vertex_check(PyObject* o)
{
  gboolean check = FALSE;
  guint i,N;
  PyObject *obj;

  /* Check for a Vertex */
  if( PyObject_TypeCheck(o, &PygtsVertexType) ) {
    check = TRUE;
  }

  /* Convert list into tuple */
  if(PyList_Check(o)) {
    o = PyList_AsTuple(o);
  }
  else {
    Py_INCREF(o);
  }

  /* Check for a tuple of floats */
  if( PyTuple_Check(o) ) {
    if( (N = PyTuple_Size(o)) <= 3 ) {
      check = TRUE;
      for(i=0;i<N;i++) {
	obj = PyTuple_GET_ITEM(o,i);
	if(!PyFloat_Check(obj) && !PyInt_Check(obj)) {
	  check = FALSE;
	}
      }
    }
  }
  Py_DECREF(o);

  if( !check ) {
    return FALSE;
  }
  else {
#if PYGTS_DEBUG
    if( PyObject_TypeCheck(o, &PygtsVertexType) ) {
      return pygts_vertex_is_ok(PYGTS_VERTEX(o));
    }
#endif
    return TRUE;
  }
}


gboolean 
pygts_vertex_is_ok(PygtsVertex *v)
{
  GSList *parent;
  PygtsObject *obj;

  obj = PYGTS_OBJECT(v);

  if(!pygts_point_is_ok(PYGTS_POINT(v))) return FALSE;

  /* Check for a valid parent */
  g_return_val_if_fail(obj->gtsobj_parent!=NULL,FALSE);
  g_return_val_if_fail(PYGTS_IS_PARENT_SEGMENT(obj->gtsobj_parent),FALSE);
  parent = g_slist_find(GTS_VERTEX(obj->gtsobj)->segments,
			obj->gtsobj_parent);
  g_return_val_if_fail(parent!=NULL,FALSE);

  return TRUE;
}


static GtsObject *
parent(GtsVertex *v1) {
  GtsPoint *p1;
  GtsVertex *v2;
  GtsSegment *p;

  /* Create another Vertex */
  p1 = GTS_POINT(v1);
  if( (v2 = gts_vertex_new(pygts_parent_vertex_class(),p1->x,p1->y,p1->z+1)) 
    == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create parent");
    return NULL;
  }

  /* Create and return the parent */
  if( (p = gts_segment_new(pygts_parent_segment_class(),v1,v2))
      == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create parent");
    gts_object_destroy(GTS_OBJECT(v2));
    return NULL;
  }

  return GTS_OBJECT(p);
}


PygtsVertex *
pygts_vertex_new(GtsVertex *v)
{
  PyObject *args, *kwds;
  PygtsObject *vertex;

  /* Check for Vertex in the object table */
  if( (vertex = PYGTS_OBJECT(g_hash_table_lookup(obj_table,GTS_OBJECT(v)))) 
      !=NULL ) {
    Py_INCREF(vertex);
    return PYGTS_VERTEX(vertex);
  }

  /* Build a new Vertex */
  args = Py_BuildValue("ddd",0,0,0);
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_False);
  vertex = PYGTS_VERTEX(PygtsVertexType.tp_new(&PygtsVertexType, args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( vertex == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Vertex");
    return NULL;
  }
  vertex->gtsobj = GTS_OBJECT(v);

  /* Attach the parent */
  if( (vertex->gtsobj_parent=parent(v)) == NULL ) {
	Py_DECREF(vertex);
	return NULL;
  }

  /* Register and return */
  pygts_object_register(vertex);
  return PYGTS_VERTEX(vertex);
}


PygtsVertex *
pygts_vertex_from_sequence(PyObject *tuple) {
  guint i,N;
  gdouble x=0,y=0,z=0;
  PyObject *obj;
  GtsVertex *v;
  PygtsVertex *vertex;

  /* Convert list into tuple */
  if(PyList_Check(tuple)) {
    tuple = PyList_AsTuple(tuple);
  }
  else {
    Py_INCREF(tuple);
  }
  if(!PyTuple_Check(tuple)) {
    Py_DECREF(tuple);
    PyErr_SetString(PyExc_TypeError,"expected a list or tuple of vertices");
    return NULL;
  }

  /* Get the tuple size */
  if( (N = PyTuple_Size(tuple)) > 3 ) {
    PyErr_SetString(PyExc_RuntimeError,
		    "expected a list or tuple of up to three floats");
    Py_DECREF(tuple);
    return NULL;
  }

  /* Get the coordinates */
  for(i=0;i<N;i++) {
    obj = PyTuple_GET_ITEM(tuple,i);

    if(!PyFloat_Check(obj) && !PyInt_Check(obj)) {
      PyErr_SetString(PyExc_TypeError,"expected a list or tuple of floats");
      Py_DECREF(tuple);
      return NULL;
    }
    if(i==0) {
      if(PyFloat_Check(obj)) x = PyFloat_AsDouble(obj);
      else  x = (double)PyInt_AsLong(obj);
    }
    if(i==1) {
      if(PyFloat_Check(obj)) y = PyFloat_AsDouble(obj);
      else  y = (double)PyInt_AsLong(obj);
    }
    if(i==2) {
      if(PyFloat_Check(obj)) z = PyFloat_AsDouble(obj);
      else  z = (double)PyInt_AsLong(obj);
    }
  }
  Py_DECREF(tuple);

  /* Create the vertex */  
  if( (v = gts_vertex_new(gts_vertex_class(),x,y,z)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Vertex");
  }
  if( (vertex = pygts_vertex_new(v)) == NULL ) {
    gts_object_destroy(GTS_OBJECT(v));
    return NULL;
  }

  return vertex;
}


GtsSegmentClass*
pygts_parent_segment_class(void)
{
  static GtsSegmentClass *klass = NULL;
  GtsObjectClass *super = NULL;

  if (klass == NULL) {

    super = GTS_OBJECT_CLASS(gts_segment_class());

    GtsObjectClassInfo pygts_parent_segment_info = {
      "PygtsParentSegment",
      sizeof(PygtsParentSegment),
      sizeof(GtsSegmentClass),
      (GtsObjectClassInitFunc)(super->info.class_init_func),
      (GtsObjectInitFunc)(super->info.object_init_func),
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new(gts_object_class(),
				 &pygts_parent_segment_info);
  }

  return klass;
}


GtsVertexClass*
pygts_parent_vertex_class(void)
{
  static GtsVertexClass *klass = NULL;
  GtsObjectClass *super = NULL;

  if (klass == NULL) {

    super = GTS_OBJECT_CLASS(gts_vertex_class());

    GtsObjectClassInfo pygts_parent_vertex_info = {
      "PygtsParentVertex",
      sizeof(PygtsParentVertex),
      sizeof(GtsVertexClass),
      (GtsObjectClassInitFunc)(super->info.class_init_func),
      (GtsObjectInitFunc)(super->info.object_init_func),
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new(gts_object_class(),
				 &pygts_parent_vertex_info);
  }

  return klass;
}
