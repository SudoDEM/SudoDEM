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
  #define SELF_CHECK if(!pygts_triangle_check((PyObject*)self)) {         \
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
is_ok(PygtsTriangle *self, PyObject *args)
{
  if(pygts_triangle_is_ok(self)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
area(PygtsTriangle *self, PyObject *args)
{
  SELF_CHECK

  return Py_BuildValue("d",
      gts_triangle_area(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)));
}


static PyObject*
perimeter(PygtsTriangle *self, PyObject *args)
{
  SELF_CHECK

  return Py_BuildValue("d",
      gts_triangle_perimeter(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)));
}


static PyObject*
quality(PygtsTriangle *self, PyObject *args)
{
  SELF_CHECK

  return Py_BuildValue("d",
      gts_triangle_quality(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)));
}


static PyObject*
normal(PygtsTriangle *self, PyObject *args)
{
  gdouble x,y,z;

  SELF_CHECK

  gts_triangle_normal(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),&x,&y,&z);
  return Py_BuildValue("ddd",x,y,z);
}


static PyObject*
revert(PygtsTriangle *self, PyObject *args)
{
  SELF_CHECK

  gts_triangle_revert(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
orientation(PygtsTriangle *self, PyObject *args)
{
  SELF_CHECK

  return Py_BuildValue("d",
      gts_triangle_orientation(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)));
}


static PyObject*
angle(PygtsTriangle* self, PyObject *args)
{
  PyObject *t_;
  PygtsTriangle *t;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &t_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_triangle_check(t_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Triangle");
    return NULL;
  }
  t = PYGTS_TRIANGLE(t_);

  return Py_BuildValue("d",
      gts_triangles_angle(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
			  PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t)));
}


static PyObject*
is_compatible(PygtsTriangle *self, PyObject *args)
{
  PyObject *t2_;
  PygtsTriangle *t2;
  GtsEdge *e;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &t2_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_triangle_check(t2_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Triangle");
    return NULL;
  }
  t2 = PYGTS_TRIANGLE(t2_);

  /* Get the common edge */
  if( (e = gts_triangles_common_edge(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				     PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t2)))
      == NULL ) {
    PyErr_SetString(PyExc_RuntimeError,"Triangles do not share common edge");
    return NULL;
  }

  if( gts_triangles_are_compatible(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				   PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t2),e) ) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
common_edge(PygtsTriangle *self, PyObject *args)
{
  PyObject *t2_;
  PygtsTriangle *t2;
  GtsEdge *e;
  PygtsEdge *edge;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &t2_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_triangle_check(t2_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Triangle");
    return NULL;
  }
  t2 = PYGTS_TRIANGLE(t2_);

  /* Get the common edge */
  if( (e = gts_triangles_common_edge(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				     PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t2))) 
      == NULL ) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  if( (edge = pygts_edge_new(GTS_EDGE(e))) == NULL ) {
    return NULL;
  }

  return (PyObject*)edge;
}


static PyObject*
opposite(PygtsTriangle *self, PyObject *args)
{
  PyObject *o_;
  PygtsEdge *e=NULL;
  PygtsVertex *v=NULL;
  GtsVertex *vertex=NULL,*v1,*v2,*v3;
  GtsEdge *edge=NULL;
  GtsTriangle *triangle;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &o_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(pygts_edge_check(o_)) {
    e = PYGTS_TRIANGLE(o_);
  }
  else {
    if(pygts_vertex_check(o_)) {
      v = PYGTS_TRIANGLE(o_);
    }
    else {
      PyErr_SetString(PyExc_TypeError,"expected an Edge or a Vertex");
      return NULL;
    }
  }

  /* Error check */
  triangle = PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self);
  if( e!=NULL ) {
    edge = PYGTS_EDGE_AS_GTS_EDGE(e);
    if(! ((triangle->e1==edge)||(triangle->e2==edge)||(triangle->e3==edge)) ) {
      PyErr_SetString(PyExc_RuntimeError,"Edge not in Triangle");
      return NULL;
    }
  }
  else {
    vertex = PYGTS_VERTEX_AS_GTS_VERTEX(v);
    gts_triangle_vertices(triangle,&v1,&v2,&v3);
    if(! ((vertex==v1)||(vertex==v2)||(vertex==v3)) ) {
      PyErr_SetString(PyExc_RuntimeError,"Vertex not in Triangle");
      return NULL;
    }
  }

  /* Get the opposite and return */
  if( e!=NULL) {
    vertex = gts_triangle_vertex_opposite(triangle, edge);
    if( (v = pygts_vertex_new(vertex)) == NULL ) {
      return NULL;
    }
    return (PyObject*)v;
  }
  else{
    edge = gts_triangle_edge_opposite(triangle, vertex);
    if( (e = pygts_edge_new(edge)) == NULL ) {
      return NULL;
    }
    return (PyObject*)e;
  }
}


static PyObject *
vertices(PygtsTriangle *self,PyObject *args)
{
  GtsVertex *v1_,*v2_,*v3_;
  PygtsObject *v1,*v2,*v3;

  SELF_CHECK

  /* Get the vertices */
  gts_triangle_vertices(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
			&v1_, &v2_, &v3_);

  if( (v1 = pygts_vertex_new(v1_)) == NULL ) {
    return NULL;
  }

  if( (v2 = pygts_vertex_new(v2_)) == NULL ) {
    Py_DECREF(v1);
    return NULL;
  }

  if( (v3 = pygts_vertex_new(v3_)) == NULL ) {
    Py_DECREF(v1);
    Py_DECREF(v2);
    return NULL;
  }

  return Py_BuildValue("OOO",v1,v2,v3);
}


static PyObject *
vertex(PygtsTriangle *self,PyObject *args)
{
  GtsVertex *v1_;
  PygtsObject *v1;

  SELF_CHECK

  /* Get the vertices */
  v1_ = gts_triangle_vertex(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self));

  if( (v1 = pygts_vertex_new(v1_)) == NULL ) {
    return NULL;
  }

  return (PyObject*)v1;
}


static PyObject *
circumcenter(PygtsTriangle *self,PyObject *args)
{
  PygtsVertex *v;
  GtsVertex *vertex;

  SELF_CHECK

  /* Get the Vertex */
  vertex = GTS_VERTEX(
	       gts_triangle_circumcircle_center(
	           PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
		   GTS_POINT_CLASS(gts_vertex_class())));

  if( vertex == NULL ) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  if( (v = pygts_vertex_new(vertex)) == NULL ) {
    return NULL;
  }

  return (PyObject*)v;
}


static PyObject *
is_stabbed(PygtsTriangle *self,PyObject *args)
{
  PyObject *p_;
  PygtsVertex *p;
  GtsObject *obj;
  PygtsVertex *vertex;
  PygtsEdge *edge;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &p_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_point_check(p_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Point");
    return NULL;
  }
  p = PYGTS_POINT(p_);

  obj = gts_triangle_is_stabbed(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				GTS_POINT(PYGTS_OBJECT(p)->gtsobj),
				NULL);

  if( obj == NULL ) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  if(GTS_IS_VERTEX(obj)) {
    if( (vertex = pygts_vertex_new(GTS_VERTEX(obj))) == NULL ) {
      return NULL;
    }
    return (PyObject*)vertex;
  }

  if(GTS_IS_EDGE(obj)) {
    if( (edge = pygts_edge_new(GTS_EDGE(obj))) == NULL ) {
      return NULL;
    }
    return (PyObject*)edge;
  }

  Py_INCREF(self);
  return (PyObject*)self;
}


static PyObject *
interpolate_height(PygtsTriangle *self,PyObject *args)
{
  PyObject *p_;
  PygtsPoint *p;
  GtsPoint point;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &p_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_point_check(p_)) {
    PyErr_SetString(PyExc_TypeError,"expected a Point");
    return NULL;
  }
  p = PYGTS_POINT(p_);

  point.x = PYGTS_POINT_AS_GTS_POINT(p)->x;
  point.y = PYGTS_POINT_AS_GTS_POINT(p)->y;

  gts_triangle_interpolate_height(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self),
				  &point);

  return Py_BuildValue("d",point.z);
}


/* Methods table */
static PyMethodDef methods[] = {
  {"is_ok", (PyCFunction)is_ok,
   METH_NOARGS,
   "True if this Triangle t is non-degenerate and non-duplicate.\n"
   "False otherwise.\n"
   "\n"
   "Signature: t.is_ok()\n"
  },  

  {"area", (PyCFunction)area,
   METH_NOARGS,
   "Returns the area of Triangle t.\n"
   "\n"
   "Signature: t.area()\n"
  },  

  {"perimeter", (PyCFunction)perimeter,
   METH_NOARGS,
   "Returns the perimeter of Triangle t.\n"
   "\n"
   "Signature: t.perimeter()\n"
  },  

  {"quality", (PyCFunction)quality,
   METH_NOARGS,
   "Returns the quality of Triangle t.\n"
   "\n"
   "The quality of a triangle is defined as the ratio of the square\n"
   "root of its surface area to its perimeter relative to this same\n"
   "ratio for an equilateral triangle with the same area.  The quality\n"
   "is then one for an equilateral triangle and tends to zero for a\n"
   "very stretched triangle."
   "\n"
   "Signature: t.quality()\n"
  },  

  {"normal", (PyCFunction)normal,
   METH_NOARGS,
   "Returns a tuple of coordinates of the oriented normal of Triangle t\n"
   "as the cross-product of two edges, using the left-hand rule.  The\n"
   "normal is not normalized.  If this triangle is part of a closed and\n" 
   "oriented surface, the normal points to the outside of the surface.\n"
   "\n"
   "Signature: t.normal()\n"
  },  

  {"revert", (PyCFunction)revert,
   METH_NOARGS,
   "Changes the orientation of triangle t, turning it inside out.\n"
   "\n"
   "Signature: t.revert()\n"
  },  

  {"orientation", (PyCFunction)orientation,
   METH_NOARGS,
   "Determines orientation of the plane (x,y) projection of Triangle t\n"
   "\n"
   "Signature: t.orientation()\n"
   "\n"
   "Returns a positive value if Points p1, p2 and p3 in Triangle t\n"
   "appear in counterclockwise order, a negative value if they appear\n"
   "in clockwise order and zero if they are colinear.\n"
  },  

  {"angle", (PyCFunction)angle,
   METH_VARARGS,
   "Returns the angle (radians) between Triangles t1 and t2\n"
   "\n"
   "Signature: t1.angle(t2)\n"
  },  

  {"is_compatible", (PyCFunction)is_compatible,
   METH_VARARGS,
   "True if this triangle t1 and other t2 are compatible;\n"
   "otherwise False.\n"
   "\n"
   "Checks if this triangle t1 and other t2, which share a common\n"
   "Edge, can be part of the same surface without conflict in the\n"
   "surface normal orientation.\n"
   "\n"
   "Signature: t1.is_compatible(t2)\n"
  },  

  {"common_edge", (PyCFunction)common_edge,
   METH_VARARGS,
   "Returns Edge common to both this Triangle t1 and other t2.\n"
   "Returns None if the triangles do not share an Edge.\n"
   "\n"
   "Signature: t1.common_edge(t2)\n"
  },  

  {"opposite", (PyCFunction)opposite,
   METH_VARARGS,
   "Returns Vertex opposite to Edge e or Edge opposite to Vertex v\n"
   "for this Triangle t.\n"
   "\n"
   "Signature: t.opposite(e) or t.opposite(v)\n"
  },  

  {"vertices", (PyCFunction)vertices,
   METH_NOARGS,
   "Returns the three oriented set of vertices in Triangle t.\n"
   "\n"
   "Signature: t.vertices()\n"
  },  

  {"vertex", (PyCFunction)vertex,
   METH_NOARGS,
   "Returns the Vertex of this Triangle t not in t.e1.\n"
   "\n"
   "Signature: t.vertex()\n"
  },  

  {"circumcenter", (PyCFunction)circumcenter,
   METH_NOARGS,
   "Returns a Vertex at the center of the circumscribing circle of\n"
   "this Triangle t, or None if the circumscribing circle is not\n"
   "defined.\n"
   "\n"
   "Signature: t.circumcircle_center()\n"
  },  

  {"is_stabbed", (PyCFunction)is_stabbed,
   METH_VARARGS,
   "Returns the component of this Triangle t that is stabbed by a\n"
   "ray projecting from Point p to z=infinity.  The result\n"
   "can be this Triangle t, one of its Edges or Vertices, or None.\n"
   "If the ray is contained in the plan of this Triangle then None is\n"
   "also returned.\n"
   "\n"
   "Signature: t.is_stabbed(p)\n"
  },  

  {"interpolate_height", (PyCFunction)interpolate_height,
   METH_VARARGS,
   "Returns the height of the plane defined by Triangle t at Point p.\n"
   "Only the x- and y-coordinates of p are considered.\n"
   "\n"
   "Signature: t.interpolate_height(p)\n"
  },  

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Attributes exported to python */

static PyObject *
get_e1(PygtsTriangle *self, void *closure)
{
  PygtsEdge *e1;

  SELF_CHECK

  if( (e1=pygts_edge_new(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)->e1)) == NULL ) {
    return NULL;
  }

  return (PyObject *)e1;
}


static PyObject *
get_e2(PygtsTriangle *self, void *closure)
{
  PygtsEdge *e2;

  SELF_CHECK

  if( (e2=pygts_edge_new(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)->e2)) == NULL ) {
    return NULL;
  }

  return (PyObject *)e2;
}


static PyObject *
get_e3(PygtsTriangle *self, void *closure)
{
  PygtsEdge *e3;

  SELF_CHECK

  if( (e3=pygts_edge_new(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(self)->e3)) == NULL ) {
    return NULL;
  }

  return (PyObject *)e3;
}


/* Methods table */
static PyGetSetDef getset[] = {
    {"e1", (getter)get_e1, NULL, "Edge 1", NULL},

    {"e2", (getter)get_e2, NULL, "Edge 2", NULL},

    {"e3", (getter)get_e3, NULL, "Edge 3", NULL},

    {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Python type methods */

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
  GtsTriangle *t,*t_;
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
      PyErr_SetString(PyExc_TypeError,"expected three Edges or three Vertices");
      return NULL;
    }
    if( (v1==v2 || v2==v3 || v1==v3) && v1!=NULL ) {
      PyErr_SetString(PyExc_ValueError,"three Vertices must be different");
      return NULL;
    }

    /* Get gts edges */
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

    /* Check that edges connect with common vertices */
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
      PyErr_SetString(PyExc_RuntimeError,
		      "Edges in triangle must connect");
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

    /* Create the GtsTriangle */
    if( (t = gts_triangle_new(gts_triangle_class(),e1,e2,e3)) == NULL )  {
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
    t_ = gts_triangle_is_duplicate(GTS_TRIANGLE(t));
    if( t_ != NULL ) {
      gts_object_destroy(GTS_OBJECT(t));
      t = t_;
    }

    /* If corresponding PyObject found in object table, we are done */
    if( (obj=g_hash_table_lookup(obj_table,GTS_OBJECT(t))) != NULL ) {
      Py_INCREF(obj);
      return (PyObject*)obj;
    }
  }
  
  /* Chain up */
  obj = PYGTS_OBJECT(PygtsObjectType.tp_new(type,args,kwds));

  if( alloc_gtsobj ) {
    obj->gtsobj = GTS_OBJECT(t);
    pygts_object_register(PYGTS_OBJECT(obj));
  }

  return (PyObject*)obj;
}


static int
init(PygtsTriangle *self, PyObject *args, PyObject *kwds)
{
  gint ret;

  /* Chain up */
  if( (ret=PygtsObjectType.tp_init((PyObject*)self,args,kwds)) != 0 ){
    return ret;
  }

#if PYGTS_DEBUG
  if(!pygts_triangle_check((PyObject*)self)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "problem with self object (internal error)");
    return -1;
  }
#endif

  return 0;
}


static int
compare(PyObject *o1, PyObject *o2)
{
  GtsTriangle *t1, *t2;

  if( !(pygts_triangle_check(o1) && pygts_triangle_check(o2)) ) {
    return -1;
  }
  t1 = PYGTS_TRIANGLE_AS_GTS_TRIANGLE(o1);
  t2 = PYGTS_TRIANGLE_AS_GTS_TRIANGLE(o2);
  
  return pygts_triangle_compare(t1,t2);  
}


/* Methods table */
PyTypeObject PygtsTriangleType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "gts.Triangle",          /* tp_name */
    sizeof(PygtsTriangle),   /* tp_basicsize */
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
    "Triangle object",       /* tp_doc */
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
pygts_triangle_check(PyObject* o)
{
  if(! PyObject_TypeCheck(o, &PygtsTriangleType)) {
    return FALSE;
  }
  else {
#if PYGTS_DEBUG
    return pygts_triangle_is_ok(PYGTS_TRIANGLE(o));
#else
    return TRUE;
#endif
  }
}


gboolean 
pygts_triangle_is_ok(PygtsTriangle *t)
{
  if(!pygts_object_is_ok(PYGTS_OBJECT(t))) return FALSE;
  return pygts_gts_triangle_is_ok(PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t));
}


PygtsTriangle *
pygts_triangle_new(GtsTriangle *t)
{
  PyObject *args, *kwds;
  PygtsObject *triangle;

  /* Check for Triangle in the object table */
  if( (triangle = PYGTS_OBJECT(g_hash_table_lookup(obj_table,GTS_OBJECT(t)))) 
      !=NULL ) {
    Py_INCREF(triangle);
    return PYGTS_TRIANGLE(triangle);
  }

  /* Build a new Triangle */
  args = Py_BuildValue("OOO",Py_None,Py_None,Py_None);
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_False);
  triangle = PYGTS_OBJECT(PygtsTriangleType.tp_new(&PygtsTriangleType,
						   args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( triangle == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Triangle");
    return NULL;
  }
  triangle->gtsobj = GTS_OBJECT(t);

  /* Register and return */
  pygts_object_register(triangle);
  return PYGTS_TRIANGLE(triangle);
}


int 
pygts_triangle_compare(GtsTriangle* t1,GtsTriangle* t2)
{
  if( (pygts_segment_compare(GTS_SEGMENT(t1->e1),GTS_SEGMENT(t2->e1))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e2),GTS_SEGMENT(t2->e2))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e3),GTS_SEGMENT(t2->e3))==0) ||
      (pygts_segment_compare(GTS_SEGMENT(t1->e1),GTS_SEGMENT(t2->e3))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e2),GTS_SEGMENT(t2->e1))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e3),GTS_SEGMENT(t2->e2))==0) ||
      (pygts_segment_compare(GTS_SEGMENT(t1->e1),GTS_SEGMENT(t2->e2))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e2),GTS_SEGMENT(t2->e3))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e3),GTS_SEGMENT(t2->e1))==0) ||

      (pygts_segment_compare(GTS_SEGMENT(t1->e1),GTS_SEGMENT(t2->e3))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e2),GTS_SEGMENT(t2->e2))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e3),GTS_SEGMENT(t2->e1))==0) ||
      (pygts_segment_compare(GTS_SEGMENT(t1->e1),GTS_SEGMENT(t2->e2))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e2),GTS_SEGMENT(t2->e1))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e3),GTS_SEGMENT(t2->e3))==0) ||
      (pygts_segment_compare(GTS_SEGMENT(t1->e1),GTS_SEGMENT(t2->e1))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e2),GTS_SEGMENT(t2->e3))==0 &&
       pygts_segment_compare(GTS_SEGMENT(t1->e3),GTS_SEGMENT(t2->e2))==0) ) {
    return 0;
  }
  return -1;
}


/**
 * gts_triangle_is_ok:
 * @t: a #GtsTriangle.
 *
 * Returns: %TRUE if @t is a non-degenerate, non-duplicate triangle,
 * %FALSE otherwise.
 */
gboolean pygts_gts_triangle_is_ok (GtsTriangle * t)
{
  g_return_val_if_fail (t != NULL, FALSE);
  g_return_val_if_fail (t->e1 != NULL, FALSE);
  g_return_val_if_fail (t->e2 != NULL, FALSE);
  g_return_val_if_fail (t->e3 != NULL, FALSE);
  g_return_val_if_fail (t->e1 != t->e2 && t->e1 != t->e3 && t->e2 != t->e3, 
                        FALSE);
  g_return_val_if_fail (gts_segments_touch (GTS_SEGMENT (t->e1), 
                                            GTS_SEGMENT (t->e2)),
                        FALSE);
  g_return_val_if_fail (gts_segments_touch (GTS_SEGMENT (t->e1), 
                                            GTS_SEGMENT (t->e3)), 
                        FALSE);
  g_return_val_if_fail (gts_segments_touch (GTS_SEGMENT (t->e2), 
                                            GTS_SEGMENT (t->e3)), 
                        FALSE);
  g_return_val_if_fail (GTS_SEGMENT (t->e1)->v1 != GTS_SEGMENT (t->e1)->v2, 
                        FALSE);
  g_return_val_if_fail (GTS_SEGMENT (t->e2)->v1 != GTS_SEGMENT (t->e2)->v2, 
                        FALSE);
  g_return_val_if_fail (GTS_SEGMENT (t->e3)->v1 != GTS_SEGMENT (t->e3)->v2, 
                        FALSE);
  /*  g_return_val_if_fail (GTS_OBJECT (t)->reserved == NULL, FALSE); */
  g_return_val_if_fail (!gts_triangle_is_duplicate (t), FALSE);
  return TRUE;
}
