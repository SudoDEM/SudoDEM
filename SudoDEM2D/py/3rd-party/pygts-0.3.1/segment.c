/* pygts - python segment for the manipulation of triangulated surfaces
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
  #define SELF_CHECK if(!pygts_segment_check((PyObject*)self)) {      \
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
is_ok(PygtsSegment *self, PyObject *args)
{
  if(pygts_segment_is_ok(self)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
intersects(PygtsSegment *self, PyObject *args)
{
  PyObject *s_;
  GtsSegment *s;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &s_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_segment_check(s_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Segment");
    return NULL;
  }
  s = PYGTS_SEGMENT_AS_GTS_SEGMENT(s_);

  return Py_BuildValue("i",
      gts_segments_are_intersecting(PYGTS_SEGMENT_AS_GTS_SEGMENT(self),s));
}


static PyObject*
connects(PygtsSegment *self, PyObject *args)
{
  PyObject *v1_,*v2_;
  GtsVertex *v1,*v2;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "OO", &v1_, &v2_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_vertex_check(v1_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Vertex");
    return NULL;
  }
  v1 = PYGTS_VERTEX_AS_GTS_VERTEX(v1_);

  if(!pygts_vertex_check(v2_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Vertex");
    return NULL;
  }
  v2 = PYGTS_VERTEX_AS_GTS_VERTEX(v2_);

  if(gts_segment_connect(PYGTS_SEGMENT_AS_GTS_SEGMENT(self),v1,v2)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
touches(PygtsSegment *self, PyObject *args)
{
  PyObject *s_;
  GtsSegment *s;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &s_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_segment_check(s_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Segment");
    return NULL;
  }
  s = PYGTS_SEGMENT_AS_GTS_SEGMENT(s_);

  if(gts_segments_touch(PYGTS_SEGMENT_AS_GTS_SEGMENT(self),s)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
midvertex(PygtsSegment *self, PyObject *args)
{
  PygtsVertex *vertex;
  GtsVertex *v;

  SELF_CHECK

  v = gts_segment_midvertex(PYGTS_SEGMENT_AS_GTS_SEGMENT(self),
			    gts_vertex_class());

  if( (vertex = pygts_vertex_new(v)) == NULL ) {
    return NULL;
  }

  return (PyObject*)vertex;
}


static PyObject*
intersection(PygtsSegment *self, PyObject *args)
{
  PyObject *t_,*boundary_=NULL;
  PygtsTriangle *t;
  gboolean boundary=TRUE;
  GtsVertex *v;
  PygtsObject *vertex;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O|O", &t_, &boundary_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_triangle_check(t_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Triangle and boolean");
    return NULL;
  }
  t = PYGTS_TRIANGLE(t_);

  if( boundary_ != NULL ) {
    if(PyBool_Check(boundary_)==FALSE) {
      PyErr_SetString(PyExc_TypeError,"expected a Triangle and boolean");
      return NULL;
    }
    if( boundary_ == Py_False ){  /* Default TRUE */
      boundary = FALSE;
    }
  }

  v = GTS_VERTEX( gts_segment_triangle_intersection(
		      PYGTS_SEGMENT_AS_GTS_SEGMENT(self),
		      PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t),
		      boundary,
		      GTS_POINT_CLASS(gts_vertex_class())) );

  if( v == NULL ) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  if( (vertex = pygts_vertex_new(v)) == NULL ) {
    return NULL;
  }

  return (PyObject *)vertex;
}


/* Methods table */
static PyMethodDef methods[] = {
  {"is_ok", (PyCFunction)is_ok,
   METH_NOARGS,
   "True if this Segment s is not degenerate or duplicate.\n"
   "False otherwise.  Degeneracy implies s.v1.id == s.v2.id.\n"
   "\n"
   "Signature: s.is_ok().\n"
  },  

  {"intersects", (PyCFunction)intersects,
   METH_VARARGS,
   "Checks if this Segment s1 intersects with Segment s2.\n"
   "Returns 1 if they intersect, 0 if an endpoint of one Segment lies\n"
   "on the other Segment, -1 otherwise\n"
   "\n"
   "Signature: s1.intersects(s2).\n"
  },  

  {"connects", (PyCFunction)connects,
   METH_VARARGS,
   "Returns True if this Segment s1 connects Vertices v1 and v2.\n"
   "False otherwise.\n"
   "\n"
   "Signature: s1.connects(v1,v2).\n"
  },  

  {"touches", (PyCFunction)touches,
   METH_VARARGS,
   "Returns True if this Segment s1 touches Segment s2\n"
   "(i.e., they share a common Vertex).  False otherwise.\n"
   "\n"
   "Signature: s1.touches(s2).\n"
  },  

  {"midvertex", (PyCFunction)midvertex,
   METH_NOARGS,
   "Returns a new Vertex at the mid-point of this Segment s.\n"
   "\n"
   "Signature: s.midvertex().\n"
  },  

  {"intersection", (PyCFunction)intersection,
   METH_VARARGS,
   "Returns the intersection of Segment s with Triangle t\n"
   "\n"
   "This function is geometrically robust in the sense that it will\n"
   "return None if s and t do not intersect and will return a\n"
   "Vertex if they do. However, the point coordinates are subject\n"
   "to round-off errors.  None will be returned if s is contained\n"
   "in the plane defined by t.\n"
   "\n"
   "Signature: s.intersection(t) or s.intersection(t,boundary).\n"
   "\n"
   "If boundary is True (default), the boundary of s is taken into\n"
   "account.\n"
   "\n"
   "Returns a summit of t (if boundary is True), one of the endpoints\n"
   "of s, a new Vertex at the intersection of s with t, or None if\n"
   "s and t don't intersect.\n"
  },  

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Attributes exported to python */

static PyObject *
get_v1(PygtsSegment *self, void *closure)
{
  PygtsObject *v1;

  SELF_CHECK

  if( (v1=pygts_vertex_new(PYGTS_SEGMENT_AS_GTS_SEGMENT(self)->v1)) == NULL ) {
    return NULL;
  }

  return (PyObject *)v1;
}


static PyObject *
get_v2(PygtsSegment *self, void *closure)
{
  PygtsObject *v2;

  SELF_CHECK

  if( (v2=pygts_vertex_new(PYGTS_SEGMENT_AS_GTS_SEGMENT(self)->v2)) == NULL ) {
    return NULL;
  }

  return (PyObject *)v2;
}


/* Methods table */
static PyGetSetDef getset[] = {
    {"v1", (getter)get_v1, NULL, "Vertex 1", NULL},
    {"v2", (getter)get_v2, NULL, "Vertex 2", NULL},
    {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Python type methods */

static PyObject *
new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyObject *o;
  PygtsObject *obj;
  GtsSegment *tmp;
  GtsObject *segment=NULL;
  PyObject *v1_=NULL,*v2_=NULL;
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
      PyErr_SetString(PyExc_ValueError,"Vertices are identical");
      return NULL;
    }

    /* Create the GtsSegment */
    segment = GTS_OBJECT(gts_segment_new(gts_segment_class(),
					 GTS_VERTEX(v1->gtsobj),
					 GTS_VERTEX(v2->gtsobj)));
    if( segment == NULL )  {
      PyErr_SetString(PyExc_MemoryError, "could not create Segment");
      return NULL;
    }

    /* Check for duplicate */
    tmp = gts_segment_is_duplicate(GTS_SEGMENT(segment));
    if( tmp != NULL ) {
      gts_object_destroy(segment);
      segment = GTS_OBJECT(tmp);
    }

    /* If corresponding PyObject found in object table, we are done */
    if( (obj=g_hash_table_lookup(obj_table,segment)) != NULL ) {
      Py_INCREF(obj);
      return (PyObject*)obj;
    }
  }  

  /* Chain up */
  obj = PYGTS_OBJECT(PygtsObjectType.tp_new(type,args,kwds));

  if( alloc_gtsobj ) {
    obj->gtsobj = segment;
    pygts_object_register(PYGTS_OBJECT(obj));
  }

  return (PyObject*)obj;
}


static int
init(PygtsSegment *self, PyObject *args, PyObject *kwds)
{
  gint ret;

  /* Chain up */
  if( (ret=PygtsObjectType.tp_init((PyObject*)self,args,kwds)) != 0 ){
    return ret;
  }

  return 0;
}


static int 
compare(PygtsSegment *s1_, PygtsSegment *s2_)
{
  GtsSegment *s1, *s2;

#if PYGTS_DEBUG
  pygts_segment_check((PyObject*)s1_);
  pygts_segment_check((PyObject*)s2_);
#endif

  s1 = PYGTS_SEGMENT_AS_GTS_SEGMENT(s1_);
  s2 = PYGTS_SEGMENT_AS_GTS_SEGMENT(s2_);
  
  return pygts_segment_compare(s1,s2);
  
}


/* Methods table */
PyTypeObject PygtsSegmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "gts.Segment",           /* tp_name */
    sizeof(PygtsSegment),    /* tp_basicsize */
    0,                       /* tp_itemsize */
    0,                       /* tp_dealloc */
    0,                       /* tp_print */
    0,                       /* tp_getattr */
    0,                       /* tp_setattr */
    (cmpfunc)compare,        /* tp_compare */
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
    "Segment object",        /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    methods,                 /* tp_methods */
    0,                       /* tp_members */
    getset,                  /* tp_getset */
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
pygts_segment_check(PyObject* o)
{
  if(! PyObject_TypeCheck(o, &PygtsSegmentType)) {
    return FALSE;
  }
  else {
#if PYGTS_DEBUG
    return pygts_segment_is_ok(PYGTS_SEGMENT(o));
#else
    return TRUE;
#endif
  }
}


gboolean 
pygts_segment_is_ok(PygtsSegment *s)
{
  if(!pygts_object_is_ok(PYGTS_OBJECT(s))) return FALSE;
  return gts_segment_is_ok(PYGTS_SEGMENT_AS_GTS_SEGMENT(s));
}


PygtsSegment * 
pygts_segment_new(GtsSegment *s)
{
  PyObject *args, *kwds;
  PygtsObject *segment;

  /* Check for Segment in the object table */
  if( (segment=PYGTS_OBJECT(g_hash_table_lookup(obj_table,GTS_OBJECT(s))))
      != NULL ) {
    Py_INCREF(segment);
    return PYGTS_FACE(segment);
  }

  /* Build a new Segment */
  args = Py_BuildValue("OO",Py_None,Py_None);
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_False);
  segment = PYGTS_SEGMENT(PygtsSegmentType.tp_new(&PygtsSegmentType, 
						  args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( segment == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Segment");
    return NULL;
  }
  segment->gtsobj = GTS_OBJECT(s);

  /* Register and return */
  pygts_object_register(segment);
  return PYGTS_SEGMENT(segment);
}


int 
pygts_segment_compare(GtsSegment* s1,GtsSegment* s2)
{
  if( (pygts_point_compare(GTS_POINT(s1->v1),GTS_POINT(s2->v1))==0 &&
       pygts_point_compare(GTS_POINT(s1->v2),GTS_POINT(s2->v2))==0) || 
      (pygts_point_compare(GTS_POINT(s1->v1),GTS_POINT(s2->v2))==0 &&
       pygts_point_compare(GTS_POINT(s1->v2),GTS_POINT(s2->v1))==0) ) {
    return 0;
  }
  return -1;
}
