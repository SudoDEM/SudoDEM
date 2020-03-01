/* pygts - python point for the manipulation of triangulated surfaces
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
  #define SELF_CHECK if(!pygts_point_check((PyObject*)self)) {      \
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
is_ok(PygtsPoint *self, PyObject *args)
{
  if(pygts_point_is_ok(self)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
set(PygtsPoint *self, PyObject *args)
{
  gdouble x=0,y=0,z=0;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "|ddd",  &x,&y,&z)) {
    return NULL;
  }

  gts_point_set(PYGTS_POINT_AS_GTS_POINT(self),x,y,z);

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
coords(PygtsPoint *self, PyObject *args)
{
  SELF_CHECK

  return Py_BuildValue("ddd",PYGTS_POINT_AS_GTS_POINT(self)->x,
		       PYGTS_POINT_AS_GTS_POINT(self)->y,
		       PYGTS_POINT_AS_GTS_POINT(self)->z);
}


static PyObject*
is_in_rectangle(PygtsPoint* self, PyObject *args)
{
  PyObject *o1_,*o2_;
  PygtsPoint *p1, *p2;
  gboolean flag = FALSE;
  gdouble x,y;

  SELF_CHECK

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "OO", &o1_, &o2_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!(pygts_point_check(o1_) && pygts_point_check(o2_))) {
    PyErr_SetString(PyExc_TypeError,"expected two Points");
    return NULL;
  }

  p1 = PYGTS_POINT(o1_);
  p2 = PYGTS_POINT(o2_);

  /* Test if point *may* be on rectangle perimeter */
  x = PYGTS_POINT_AS_GTS_POINT(self)->x;
  y = PYGTS_POINT_AS_GTS_POINT(self)->y;
  if( PYGTS_POINT_AS_GTS_POINT(p1)->x == x ||
      PYGTS_POINT_AS_GTS_POINT(p1)->y == y ||
      PYGTS_POINT_AS_GTS_POINT(p2)->x == x ||
      PYGTS_POINT_AS_GTS_POINT(p2)->y == y ) {
    flag = TRUE;
  }

  if( gts_point_is_in_rectangle(PYGTS_POINT_AS_GTS_POINT(self), 
				PYGTS_POINT_AS_GTS_POINT(p1), 
				PYGTS_POINT_AS_GTS_POINT(p2)) ) {
    if(flag) {
      return Py_BuildValue("i",0);
    }
    else {
      return Py_BuildValue("i",1);
    }
  }
  else {
    if( flag && 
	gts_point_is_in_rectangle(PYGTS_POINT_AS_GTS_POINT(self), 
				  PYGTS_POINT_AS_GTS_POINT(p2), 
				  PYGTS_POINT_AS_GTS_POINT(p1))) {
      return Py_BuildValue("i",0);
    }
    else {
      return Py_BuildValue("i",-1);
    }
  }
}


static PyObject*
distance(PygtsPoint* self, PyObject *args)
{
  PyObject *o_;
  PygtsPoint *p=NULL;
  PygtsSegment *s=NULL;
  PygtsTriangle *t=NULL;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &o_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(pygts_point_check(o_)) {
    p = PYGTS_POINT(o_);
  }
  else {
    if(pygts_segment_check(o_)) {
      s = PYGTS_SEGMENT(o_);
    }
    else {
      if(pygts_triangle_check(o_)) {
	t = PYGTS_TRIANGLE(o_);
      }
      else {
	PyErr_SetString(PyExc_TypeError,
			"expected a Point, Segment or Triangle");
	return NULL;
      }
    }
  }

  if(p!=NULL) {
    return Py_BuildValue("d",
        gts_point_distance(PYGTS_POINT_AS_GTS_POINT(self),
			   PYGTS_POINT_AS_GTS_POINT(p)));
  }
  else {
    if(s!=NULL) {
      return Py_BuildValue("d",
          gts_point_segment_distance(PYGTS_POINT_AS_GTS_POINT(self),
				     PYGTS_SEGMENT_AS_GTS_SEGMENT(s) )
			   );
    }
    else {
      return Py_BuildValue("d",
          gts_point_triangle_distance(PYGTS_POINT_AS_GTS_POINT(self),
				      PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t) )
			   );
    }
  }
}


static PyObject*
distance2(PygtsPoint* self, PyObject *args)
{
  PyObject *o_;
  PygtsPoint *p=NULL;
  PygtsSegment *s=NULL;
  PygtsTriangle *t=NULL;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O", &o_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(pygts_point_check(o_)) {
    p = PYGTS_POINT(o_);
  }
  else {
    if(pygts_segment_check(o_)) {
      s = PYGTS_SEGMENT(o_);
    }
    else {
      if(pygts_triangle_check(o_)) {
	t = PYGTS_TRIANGLE(o_);
      }
      else {
	PyErr_SetString(PyExc_TypeError,
			"expected a Point, Segment or Triangle");
	return NULL;
      }
    }
  }

  if(p!=NULL) {
    return Py_BuildValue("d",
        gts_point_distance2(PYGTS_POINT_AS_GTS_POINT(self),
			    PYGTS_POINT_AS_GTS_POINT(p)));
  }
  else {
    if(s!=NULL) {
      return Py_BuildValue("d",
        gts_point_segment_distance2(PYGTS_POINT_AS_GTS_POINT(self),
				    PYGTS_SEGMENT_AS_GTS_SEGMENT(s) )
			   );
    }
    else {
      return Py_BuildValue("d",
          gts_point_triangle_distance2(PYGTS_POINT_AS_GTS_POINT(self),
				       PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t) )
			   );
    }
  }
}


static PyObject*
orientation_3d(PygtsPoint* self, PyObject *args)
{
  PyObject *p1_,*p2_,*p3_;
  PygtsPoint *p1,*p2,*p3;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "OOO", &p1_, &p2_, &p3_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_point_check(p1_)) {
    PyErr_SetString(PyExc_TypeError,"expected three Points");
    return NULL;
  }
  if(!pygts_point_check(p2_)) {
    PyErr_SetString(PyExc_TypeError,"expected three Points");
    return NULL;
  }
  if(!pygts_point_check(p3_)) {
    PyErr_SetString(PyExc_TypeError,"expected three Points");
    return NULL;
  }
  p1 = PYGTS_POINT(p1_);
  p2 = PYGTS_POINT(p2_);
  p3 = PYGTS_POINT(p3_);

  return Py_BuildValue("d",
             gts_point_orientation_3d(PYGTS_POINT_AS_GTS_POINT(p1),
				      PYGTS_POINT_AS_GTS_POINT(p2),
				      PYGTS_POINT_AS_GTS_POINT(p3),
				      PYGTS_POINT_AS_GTS_POINT(self)));
}


static PyObject*
orientation_3d_sos(PygtsPoint* self, PyObject *args)
{
  PyObject *p1_,*p2_,*p3_;
  PygtsPoint *p1,*p2,*p3;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "OOO", &p1_, &p2_, &p3_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(!pygts_point_check(p1_)) {
    PyErr_SetString(PyExc_TypeError,"expected three Points");
    return NULL;
  }
  if(!pygts_point_check(p2_)) {
    PyErr_SetString(PyExc_TypeError,"expected three Points");
    return NULL;
  }
  if(!pygts_point_check(p3_)) {
    PyErr_SetString(PyExc_TypeError,"expected three Points");
    return NULL;
  }
  p1 = PYGTS_POINT(p1_);
  p2 = PYGTS_POINT(p2_);
  p3 = PYGTS_POINT(p3_);

  return Py_BuildValue("i",gts_point_orientation_3d_sos(
                             PYGTS_POINT_AS_GTS_POINT(p1),
			     PYGTS_POINT_AS_GTS_POINT(p2),
			     PYGTS_POINT_AS_GTS_POINT(p3),
			     PYGTS_POINT_AS_GTS_POINT(self)));
}


static PyObject*
is_in_circle(PygtsPoint* self, PyObject *args)
{
  PyObject *o1_=NULL,*o2_=NULL,*o3_=NULL;
  PygtsPoint *p1=NULL, *p2=NULL, *p3=NULL;
  PygtsTriangle *t=NULL;
  gdouble result;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "O|OO", &o1_, &o2_, &o3_) ) {
    return NULL;
  }
  if( (o2_==NULL && o3_!=NULL) || (o2_!=NULL && o3_==NULL) ) {
    PyErr_SetString(PyExc_TypeError,"expected three Points or one Triangle");
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(o2_==NULL && o3_==NULL) {
    if(!pygts_triangle_check(o1_)) {
      PyErr_SetString(PyExc_TypeError,"expected three Points or one Triangle");
      return NULL;
    }
    t = PYGTS_TRIANGLE(o1_);
  } 
  else {
    if(!pygts_point_check(o1_)) {
      PyErr_SetString(PyExc_TypeError,"expected three Points or one Triangle");
      return NULL;
    }
    if(!pygts_point_check(o2_)) {
      PyErr_SetString(PyExc_TypeError,"expected three Points or one Triangle");
      return NULL;
    }
    if(!pygts_point_check(o3_)) {
      PyErr_SetString(PyExc_TypeError,"expected three Points or one Triangle");
      return NULL;
    }
    p1 = PYGTS_POINT(o1_);
    p2 = PYGTS_POINT(o2_);
    p3 = PYGTS_POINT(o3_);
  }


  if(t!=NULL){
    result=gts_point_in_triangle_circle(PYGTS_POINT_AS_GTS_POINT(self),
					PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t));
  }
  else {
    result = gts_point_in_circle(PYGTS_POINT_AS_GTS_POINT(self),
				 PYGTS_POINT_AS_GTS_POINT(p1), 
				 PYGTS_POINT_AS_GTS_POINT(p2), 
				 PYGTS_POINT_AS_GTS_POINT(p3));
  }
  if(result>0) return Py_BuildValue("i",1);
  if(result==0) return Py_BuildValue("i",0);
  return Py_BuildValue("i",-1);
}


static PyObject*
is_in(PygtsPoint* self, PyObject *args)
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

  return Py_BuildValue("i",
      gts_point_is_in_triangle(PYGTS_POINT_AS_GTS_POINT(self),
			       PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t)));
}


static PyObject*
is_inside(PygtsPoint* self, PyObject *args)
{
  PyObject *s_;
  PygtsSurface *s;
  GNode *tree;
  gboolean is_open=FALSE, ret;

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
  s = PYGTS_SURFACE(s_);

  /* Error check */
  if(!gts_surface_is_closed(PYGTS_SURFACE_AS_GTS_SURFACE(s))) {
    PyErr_SetString(PyExc_RuntimeError,"Surface is not closed");
    return NULL;
  }

  /* Determing is_open parameter; note the meaning is different from the 
   * error check above.
   */
  if( gts_surface_volume(PYGTS_SURFACE_AS_GTS_SURFACE(s))<0. ) {
    is_open = TRUE;
  }

  /* Construct the tree */
  if((tree=gts_bb_tree_surface(PYGTS_SURFACE_AS_GTS_SURFACE(s))) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create GTree");
    return NULL;
  }
  
  /* Make the call */
  ret = gts_point_is_inside_surface(PYGTS_POINT_AS_GTS_POINT(self), tree,
				    is_open);

  g_node_destroy(tree);

  if(ret) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}


static PyObject*
closest(PygtsPoint* self, PyObject *args)
{
  PyObject *o1_,*o2_;
  PygtsPoint *p=NULL;
  PygtsSegment *s=NULL;
  PygtsTriangle *t=NULL;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTuple(args, "OO", &o1_, &o2_) ) {
    return NULL;
  }

  /* Convert to PygtsObjects */
  if(pygts_segment_check(o1_)) {
    s = PYGTS_SEGMENT(o1_);
  }
  else {
    if(pygts_triangle_check(o1_)) {
      t = PYGTS_TRIANGLE(o1_);
    }
    else {
      PyErr_SetString(PyExc_TypeError,
		      "expected a Segment or Triangle, and a Point");
      return NULL;
    }
  }
  if(pygts_point_check(o2_)) {
    p = PYGTS_POINT(o2_);
  }
  else {
	PyErr_SetString(PyExc_TypeError,
			"expected a Segment or Triangle, and a Point");
	return NULL;
  }

  if(s!=NULL) {
    gts_point_segment_closest(PYGTS_POINT_AS_GTS_POINT(p),
			      PYGTS_SEGMENT_AS_GTS_SEGMENT(s),
			      PYGTS_POINT_AS_GTS_POINT(self));
  }
  else {
    gts_point_triangle_closest(PYGTS_POINT_AS_GTS_POINT(p),
			       PYGTS_TRIANGLE_AS_GTS_TRIANGLE(t),
			       PYGTS_POINT_AS_GTS_POINT(self));
  }

  Py_INCREF(self);
  return (PyObject*)self;
}


/* Helper for rotate() */
gint
pygts_point_rotate(GtsPoint* p, gdouble dx, gdouble dy, gdouble dz, gdouble a)
{
  GtsMatrix *m;
  GtsVector v;

  v[0] = dx; v[1] = dy; v[2] = dz;
  if( (m = gts_matrix_rotate(NULL,v,a)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create matrix");
    return -1;
  }
  gts_point_transform(p,m);
  gts_matrix_destroy(m);

  return 0;
}


static PyObject*
rotate(PygtsPoint* self, PyObject *args, PyObject *keywds)
{
  static char *kwlist[] = {"dx", "dy", "dz", "a", NULL};
  gdouble dx=0,dy=0,dz=0,a=0;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTupleAndKeywords(args, keywds,"|dddd", kwlist,
				   &dx, &dy, &dz, &a) ) {
    return NULL;
  }

  if(pygts_point_rotate(PYGTS_POINT_AS_GTS_POINT(self),dx,dy,dz,a)==-1)
    return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}


/* Helper for scale() */
gint
pygts_point_scale(GtsPoint* p, gdouble dx, gdouble dy, gdouble dz)
{
  GtsMatrix *m;
  GtsVector v;

  v[0] = dx; v[1] = dy; v[2] = dz;
  if( (m = gts_matrix_scale(NULL,v)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create matrix");
    return -1;
  }
  gts_point_transform(p,m);
  gts_matrix_destroy(m);

  return 0;
}


static PyObject*
scale(PygtsPoint* self, PyObject *args, PyObject *keywds)
{
  static char *kwlist[] = {"dx", "dy", "dz", NULL};
  gdouble dx=1,dy=1,dz=1;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTupleAndKeywords(args, keywds,"|ddd", kwlist,
				   &dx, &dy, &dz) ) {
    return NULL;
  }

  if(pygts_point_scale(PYGTS_POINT_AS_GTS_POINT(self),dx,dy,dz)==-1)
    return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}


/* Helper for translate() */
gint
pygts_point_translate(GtsPoint* p, gdouble dx, gdouble dy, gdouble dz)
{
  GtsMatrix *m;
  GtsVector v;

  v[0] = dx; v[1] = dy; v[2] = dz;
  if( (m = gts_matrix_translate(NULL,v)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create matrix");
    return -1;
  }  
  gts_point_transform(p,m);
  gts_matrix_destroy(m);

  return 0;
}


static PyObject*
translate(PygtsPoint* self, PyObject *args, PyObject *keywds)
{
  static char *kwlist[] = {"dx", "dy", "dz", NULL};
  gdouble dx=0,dy=0,dz=0;

  SELF_CHECK

  /* Parse the args */
  if(! PyArg_ParseTupleAndKeywords(args, keywds,"|ddd", kwlist,
				   &dx, &dy, &dz) ) {
    return NULL;
  }

  if(pygts_point_translate(PYGTS_POINT_AS_GTS_POINT(self),dx,dy,dz)==-1)
    return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}


/* Methods table */
static PyMethodDef methods[] = {

  {"is_ok", (PyCFunction)is_ok,
   METH_NOARGS,
   "True if this Point p is OK.  False otherwise.\n"
   "This method is useful for unit testing and debugging.\n"
   "\n"
   "Signature: p.is_ok().\n"
  },  

  {"set", (PyCFunction)set,
   METH_VARARGS,
   "Sets x, y, and z coordinates of this Point p.\n"
   "\n"
   "Signature: p.set(x,y,z)\n"
  },  

  {"coords", (PyCFunction)coords,
   METH_VARARGS,
   "Returns a tuple of the x, y, and z coordinates for this Point p.\n"
   "\n"
   "Signature: p.coords(x,y,z)\n"
  },  

  {"is_in_rectangle", (PyCFunction)is_in_rectangle,
   METH_VARARGS,
   "True if this Point p is in box with bottom-left and upper-right\n"
   "Points p1 and p2.\n"
   "\n"
   "Signature: p.is_in_rectange(p1,p2)\n"
  },

  {"distance", (PyCFunction)distance,
   METH_VARARGS,
   "Returns Euclidean distance between this Point p and other Point p2,\n"
   "Segment s, or Triangle t."
   "\n"
   "Signature: p.distance(p2), p.distance(s) or p.distance(t)\n"
  },  

  {"distance2", (PyCFunction)distance2,
   METH_VARARGS,
   "Returns squared Euclidean distance between Point p and Point p2,\n"
   "Segment s, or Triangle t.\n"
   "\n"
   "Signature: p.distance2(p2), p.distance2(s), or p.distance2(t)\n"
  },  

  {"orientation_3d", (PyCFunction)orientation_3d,
   METH_VARARGS,
   "Determines if this Point p is above, below or on plane of 3 Points\n"
   "p1, p2 and p3.\n"
   "\n"
   "Signature: p.orientation_3d(p1,p2,p3)\n"
   "\n"
   "Below is defined so that p1, p2 and p3 appear in counterclockwise\n"
   "order when viewed from above the plane.\n"
   "\n"
   "The return value is positive if p4 lies below the plane, negative\n"
   "if p4 lies above the plane, and zero if the four points are\n"
   "coplanar.  The value is an approximation of six times the signed\n"
   "volume of the tetrahedron defined by the four points.\n"
  },  

  {"orientation_3d_sos", (PyCFunction)orientation_3d_sos,
   METH_VARARGS,
   "Determines if this Point p is above, below or on plane of 3 Points\n"
   "p1, p2 and p3.\n"
   "\n"
   "Signature: p.orientation_3d_sos(p1,p2,p3)\n"
   "\n"
   "Below is defined so that p1, p2 and p3 appear in counterclockwise\n"
   "order when viewed from above the plane.\n"
   "\n"
   "The return value is +1 if p4 lies below the plane, and -1 if p4\n"
   "lies above the plane.  Simulation of Simplicity (SoS) is used to\n"
   "break ties when the orientation is degenerate (i.e. the point lies\n"
   "on the plane definedby p1, p2 and p3)."
  },  

  {"is_in_circle", (PyCFunction)is_in_circle,
   METH_VARARGS,
   "Tests if this Point p is inside or outside circumcircle.\n"
   "The planar projection (x,y) of Point p is tested against the\n"
   "circumcircle defined by the planar projection of p1, p2 and p3,\n"
   "or alternatively the Triangle t\n"
   "\n"
   "Signature: p.in_circle(p1,p2,p3) or p.in_circle(t) \n"
   "\n"
   "Returns +1 if p lies inside, -1 if p lies outside, and 0 if p lies\n"
   "on the circle.  The Points p1, p2, and p3 must be in\n"
   "counterclockwise order, or the sign of the result will be reversed.\n"
  },

  {"is_in", (PyCFunction)is_in,
   METH_VARARGS,
   "Tests if this Point p is inside or outside Triangle t.\n"
   "The planar projection (x,y) of Point p is tested against the\n"
   "planar projection of Triangle t.\n"
   "\n"
   "Signature: p.in_circle(p1,p2,p3) or p.in_circle(t) \n"
   "\n"
   "Returns a +1 if p lies inside, -1 if p lies outside, and 0\n"
   "if p lies on the triangle.\n"
  },

  {"is_inside", (PyCFunction)is_inside,
   METH_VARARGS,
   "True if this Point p is inside or outside Surface s.\n"
   "False otherwise.\n"
   "\n"
   "Signature: p.in_inside(s)\n"
  },

  {"closest", (PyCFunction)closest,
   METH_VARARGS,
   "Set the coordinates of Point p to the Point on Segment s\n"
   "or Triangle t closest to the Point p2\n"
   "\n"
   "Signature: p.closest(s,p2) or p.closest(t,p2)\n"
   "\n"
   "Returns the (modified) Point p.\n"
  },

  {"rotate", (PyCFunction)rotate,
   METH_VARARGS | METH_KEYWORDS,
   "Rotates Point p around vector dx,dy,dz by angle a.\n"
   "The sense of the rotation is given by the right-hand-rule.\n"
   "\n"
   "Signature: p.rotate(dx=0,dy=0,dz=0,a=0)\n"
  },

  {"scale", (PyCFunction)scale,
   METH_VARARGS | METH_KEYWORDS,
   "Scales Point p by vector dx,dy,dz.\n"
   "\n"
   "Signature: p.scale(dx=1,dy=1,dz=1)\n"
  },

  {"translate", (PyCFunction)translate,
   METH_VARARGS | METH_KEYWORDS,
   "Translates Point p by vector dx,dy,dz.\n"
   "\n"
   "Signature: p.translate(dx=0,dy=0,dz=0)\n"
  },

  {NULL}  /* Sentinel */
};


/*-------------------------------------------------------------------------*/
/* Attributes exported to python */

static PyObject *
getx(PygtsPoint *self, void *closure)
{
  SELF_CHECK

  return Py_BuildValue("d",PYGTS_POINT_AS_GTS_POINT(self)->x);
}


static int
setx(PygtsPoint *self, PyObject *value, void *closure)
{
  if(PyFloat_Check(value)) {
    PYGTS_POINT_AS_GTS_POINT(self)->x = PyFloat_AsDouble(value);
  }
  else if(PyInt_Check(value)) {
    PYGTS_POINT_AS_GTS_POINT(self)->x = (gdouble)PyInt_AsLong(value);
  }
  else {
    PyErr_SetString(PyExc_TypeError,"expected a float");
    return -1;
  }
  return 0;
}


static PyObject *
gety(PygtsPoint *self, void *closure)
{
  SELF_CHECK

  return Py_BuildValue("d",PYGTS_POINT_AS_GTS_POINT(self)->y);
}


static int
sety(PygtsPoint *self, PyObject *value, void *closure)
{
  if(PyFloat_Check(value)) {
    PYGTS_POINT_AS_GTS_POINT(self)->y = PyFloat_AsDouble(value);
  }
  else if(PyInt_Check(value)) {
    PYGTS_POINT_AS_GTS_POINT(self)->y = (gdouble)PyInt_AsLong(value);
  }
  else {
    PyErr_SetString(PyExc_TypeError,"expected a float");
    return -1;
  }
  return 0;
}


static PyObject *
getz(PygtsPoint *self, void *closure)
{
  SELF_CHECK

  return Py_BuildValue("d",PYGTS_POINT_AS_GTS_POINT(self)->z);
}


static int
setz(PygtsPoint *self, PyObject *value, void *closure)
{
  if(PyFloat_Check(value)) {
    PYGTS_POINT_AS_GTS_POINT(self)->z = PyFloat_AsDouble(value);
  }
  else if(PyInt_Check(value)) {
    PYGTS_POINT_AS_GTS_POINT(self)->z = (gdouble)PyInt_AsLong(value);
  }
  else {
    PyErr_SetString(PyExc_TypeError,"expected a float");
    return -1;
  }
  return 0;
}


/* Methods table */
static PyGetSetDef getset[] = {
    {"x", (getter)getx, (setter)setx, "x value", NULL},
    {"y", (getter)gety, (setter)sety, "y value", NULL},
    {"z", (getter)getz, (setter)setz, "z value", NULL},
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
  obj = PYGTS_OBJECT(PygtsObjectType.tp_new(type,args,kwds));

  /* Allocate the gtsobj (if needed) */
  if( alloc_gtsobj ) {
    obj->gtsobj = GTS_OBJECT(gts_point_new(gts_point_class(),0,0,0));
    if( obj->gtsobj == NULL )  {
      PyErr_SetString(PyExc_MemoryError, "could not create Point");
      return NULL;
    }

    pygts_object_register(obj);
  }

  return (PyObject*)obj;
}


static int
init(PygtsPoint *self, PyObject *args, PyObject *kwds)
{
  PygtsObject *obj;
  gdouble x=0,y=0,z=0;
  guint a;
  gint ret;
  static char *kwlist[] = {"x", "y", "z", "alloc_gtsobj", NULL};

  obj = PYGTS_OBJECT(self);

  /* Parse the args */
  if(! PyArg_ParseTupleAndKeywords(args, kwds, "|dddi", kwlist, &x,&y,&z,&a)) {
    return -1;
  }

  /* Initialize */
  gts_point_set(GTS_POINT(obj->gtsobj),x,y,z);

  /* Chain up */
  if( (ret=PygtsObjectType.tp_init((PyObject*)self,args,kwds)) != 0 ) {
    return ret;
  }

  return 0;
}


static int 
compare(PygtsPoint *p1_, PygtsPoint *p2_)
{
  GtsPoint *p1, *p2;

#if PYGTS_DEBUG
  pygts_point_check((PyObject*)p1_);
  pygts_point_check((PyObject*)p2_);
#endif

  p1 = PYGTS_POINT_AS_GTS_POINT(p1_);
  p2 = PYGTS_POINT_AS_GTS_POINT(p2_);

  return pygts_point_compare(p1,p2);
}


/* Methods table */
PyTypeObject PygtsPointType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /* ob_size */
  "gts.Point"  ,             /* tp_name */
  sizeof(PygtsPoint),        /* tp_basicsize */
  0,                         /* tp_itemsize */
  0,                         /* tp_dealloc */
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
  "Point object",            /* tp_doc */
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
pygts_point_check(PyObject* o)
{
  gboolean check = FALSE;
  guint i,N;
  PyObject *obj;

  /* Check for a Point */
  if( PyObject_TypeCheck(o, &PygtsPointType) ) {
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
    if( PyObject_TypeCheck(o, &PygtsPointType) ) {
      return pygts_point_is_ok(PYGTS_POINT(o));
    }
#endif
    return TRUE;
  }
}


gboolean 
pygts_point_is_ok(PygtsPoint *p)
{
  return pygts_object_is_ok(PYGTS_OBJECT(p));
}


PygtsPoint *
pygts_point_new(GtsPoint *p)
{
  PyObject *args, *kwds;
  PygtsObject *point;

  /* Check for Point in the object table */
  if( (point = PYGTS_OBJECT(g_hash_table_lookup(obj_table,GTS_OBJECT(p)))) 
      !=NULL ) {
    Py_INCREF(point);
    return PYGTS_POINT(point);
  }

  /* Build a new Point */
  args = Py_BuildValue("ddd",0,0,0);
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_False);
  point = PYGTS_POINT(PygtsPointType.tp_new(&PygtsPointType, args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( point == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Point");
    return NULL;
  }
  point->gtsobj = GTS_OBJECT(p);

  /* Register and return */
  pygts_object_register(point);
  return PYGTS_POINT(point);
}


PygtsPoint *
pygts_point_from_sequence(PyObject *tuple) {
  guint i,N;
  gdouble x=0,y=0,z=0;
  PyObject *obj;
  GtsPoint *p;
  PygtsPoint *point;

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
  if( (p = gts_point_new(gts_point_class(),x,y,z)) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Point");
  }
  if( (point = pygts_point_new(p)) == NULL ) {
    gts_object_destroy(GTS_OBJECT(p));
    return NULL;
  }

  return point;
}


int 
pygts_point_compare(GtsPoint* p1,GtsPoint* p2)
{
  double r1,r2;

  if( (p1->x==p2->x) && (p1->y==p2->y) && (p1->z==p2->z) ) {
    return 0;
  }

  /* Compare distances from origin */
  r1 = sqrt(pow(p1->x,2) + pow(p1->y,2) + pow(p1->z,2));
  r2 = sqrt(pow(p2->x,2) + pow(p2->y,2) + pow(p2->z,2));
  if(r1<r2) return -1;
  if(r1>r2) return 1;

  /* Compare horizontal distances from origin */
  r1 = sqrt(pow(p1->x,2) + pow(p1->y,2));
  r2 = sqrt(pow(p2->x,2) + pow(p2->y,2));
  if(r1<r2) return -1;
  if(r1>r2) return 1;

  /* Compare x */
  r1 = p1->x;
  r2 = p2->x;
  if(r1<r2) return -1;
  if(r1>r2) return 1;

  /* Compare y */
  r1 = p1->y;
  r2 = p2->y;
  if(r1<r2) return -1;
  if(r1>r2) return 1;

  /* Compare z */
  r1 = p1->z;
  r2 = p2->z;
  if(r1<r2) return -1;
  return 1;
}
