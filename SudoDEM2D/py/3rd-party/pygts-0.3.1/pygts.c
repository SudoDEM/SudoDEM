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

#if PYGTS_HAS_NUMPY
  #include "numpy/arrayobject.h"
#endif


static PyObject*
merge(PyObject *self, PyObject *args)
{
  PyObject *tuple, *obj;
  guint i,N;
  GList *vertices=NULL,*v;
  gdouble epsilon;
  PygtsVertex *vertex;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "Od", &tuple, &epsilon) ) {
    return NULL;
  }
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

  /* Assemble the GList */
  N = PyTuple_Size(tuple);
  for(i=N-1;i>0;i--) {
    obj = PyTuple_GET_ITEM(tuple,i);
    if(!pygts_vertex_check(obj)) {
      Py_DECREF(tuple);
      g_list_free(vertices);
      PyErr_SetString(PyExc_TypeError,"expected a list or tuple of vertices");
      return NULL;
    }
    vertices = g_list_prepend(vertices,PYGTS_VERTEX_AS_GTS_VERTEX(obj));
  }
  Py_DECREF(tuple);

  /* Make the call */
  vertices = pygts_vertices_merge(vertices,epsilon,NULL);

  /* Assemble the return tuple */
  N = g_list_length(vertices);
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create tuple");
    return NULL;
  }
  v = vertices;
  for(i=0;i<N;i++) {
    if( (vertex = PYGTS_VERTEX(g_hash_table_lookup(obj_table,
						   GTS_OBJECT(v->data))
			       )) ==NULL ) {
      PyErr_SetString(PyExc_RuntimeError,
		      "could not get object from table (internal error)");
      g_list_free(vertices);
      return NULL;
    }
    Py_INCREF(vertex);
    PyTuple_SET_ITEM(tuple,i,(PyObject*)vertex);
    v = g_list_next(v);
  }

  g_list_free(vertices);

  return tuple;
}


static PyObject*
vertices(PyObject *self, PyObject *args)
{
  PyObject *tuple, *obj;
  guint i,N;
  GSList *segments=NULL,*vertices=NULL,*v;
  PygtsVertex *vertex;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &tuple) ) {
    return NULL;
  }
  if(PyList_Check(tuple)) {
    tuple = PyList_AsTuple(tuple);
  }
  else {
    Py_INCREF(tuple);
  }
  if(!PyTuple_Check(tuple)) {
    Py_DECREF(tuple);
    PyErr_SetString(PyExc_TypeError,"expected a list or tuple of Segments");
    return NULL;
  }

  /* Assemble the GSList */
  N = PyTuple_Size(tuple);
  for(i=0;i<N;i++) {
    obj = PyTuple_GET_ITEM(tuple,N-1-i);
    if(!pygts_segment_check(obj)) {
      Py_DECREF(tuple);
      g_slist_free(segments);
      PyErr_SetString(PyExc_TypeError,"expected a list or tuple of Segments");
      return NULL;
    }
    segments = g_slist_prepend(segments,PYGTS_SEGMENT_AS_GTS_SEGMENT(obj));
  }
  Py_DECREF(tuple);

  /* Make the call */
  vertices = gts_vertices_from_segments(segments);
  g_slist_free(segments);

  /* Assemble the return tuple */
  N = g_slist_length(vertices);
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create tuple");
    return NULL;
  }
  v = vertices;
  for(i=0;i<N;i++) {
    if( (vertex = pygts_vertex_new(GTS_VERTEX(v->data))) == NULL ) {
	Py_DECREF(tuple);
	g_slist_free(vertices);
	return NULL;
    }
    PyTuple_SET_ITEM(tuple,i,(PyObject*)vertex);
    v = g_slist_next(v);
  }

  g_slist_free(vertices);

  return tuple;
}


static PyObject*
segments(PyObject *self, PyObject *args)
{
  PyObject *tuple, *obj;
  guint i,n,N;
  GSList *segments=NULL,*vertices=NULL,*s;
  PygtsSegment *segment;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &tuple) ) {
    return NULL;
  }
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

  /* Assemble the GSList */
  N = PyTuple_Size(tuple);
  for(i=0;i<N;i++) {
    obj = PyTuple_GET_ITEM(tuple,i);
    if(!pygts_vertex_check(obj)) {
      Py_DECREF(tuple);
      g_slist_free(vertices);
      PyErr_SetString(PyExc_TypeError,"expected a list or tuple of vertices");
      return NULL;
    }
    vertices = g_slist_prepend(vertices,PYGTS_VERTEX_AS_GTS_VERTEX(obj));
  }
  Py_DECREF(tuple);

  /* Make the call */
  if( (segments = gts_segments_from_vertices(vertices)) == NULL ) {
    PyErr_SetString(PyExc_RuntimeError,"could not retrieve segments");
    return NULL;
  }
  g_slist_free(vertices);

  /* Assemble the return tuple */
  N = g_slist_length(segments);
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create tuple");
    return NULL;
  }

  s = segments;
  n=0;
  while(s!=NULL) {

    /* Skip any parent segments */
    if(PYGTS_IS_PARENT_SEGMENT(s->data) || PYGTS_IS_PARENT_EDGE(s->data)) {
      s = g_slist_next(s);
      segment = NULL;
      continue;
    }

    /* Fill in the tuple */
    if(GTS_IS_EDGE(s->data)) {
      segment = PYGTS_SEGMENT(pygts_edge_new(GTS_EDGE(s->data)));
    }
    else {
      segment = pygts_segment_new(GTS_SEGMENT(s->data));
    }
    if( segment == NULL ) {
      Py_DECREF(tuple);
      g_slist_free(segments);
      return NULL;
    }
    PyTuple_SET_ITEM(tuple,n,(PyObject*)segment);
    s = g_slist_next(s);
    n += 1;
  }

  g_slist_free(segments);

  if(_PyTuple_Resize(&tuple,n)!=0) {
    Py_DECREF(tuple);
    return NULL;
  }

  return tuple;
}


static PyObject*
triangles(PyObject *self, PyObject *args)
{
  PyObject *tuple, *obj;
  guint i,n,N;
  GSList *edges=NULL,*triangles=NULL,*t;
  PygtsTriangle *triangle;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &tuple) ) {
    return NULL;
  }
  if(PyList_Check(tuple)) {
    tuple = PyList_AsTuple(tuple);
  }
  else {
    Py_INCREF(tuple);
  }
  if(!PyTuple_Check(tuple)) {
    Py_DECREF(tuple);
    PyErr_SetString(PyExc_TypeError,"expected a list or tuple of edges");
    return NULL;
  }

  /* Assemble the GSList */
  N = PyTuple_Size(tuple);
  for(i=0;i<N;i++) {
    obj = PyTuple_GET_ITEM(tuple,i);
    if(!pygts_edge_check(obj)) {
      Py_DECREF(tuple);
      g_slist_free(edges);
      PyErr_SetString(PyExc_TypeError,"expected a list or tuple of edges");
      return NULL;
    }
    edges = g_slist_prepend(edges,PYGTS_EDGE_AS_GTS_EDGE(obj));
  }
  Py_DECREF(tuple);

  /* Make the call */
  if( (triangles = gts_triangles_from_edges(edges)) == NULL ) {
    PyErr_SetString(PyExc_RuntimeError,"could not retrieve triangles");
    return NULL;
  }
  g_slist_free(edges);

  /* Assemble the return tuple */
  N = g_slist_length(triangles);
  if( (tuple=PyTuple_New(N)) == NULL) {
    PyErr_SetString(PyExc_MemoryError,"could not create tuple");
    return NULL;
  }

  t = triangles;
  n=0;
  while(t!=NULL) {

    /* Skip any parent triangles */
    if(PYGTS_IS_PARENT_TRIANGLE(t->data)) {
      t = g_slist_next(t);
      triangle = NULL;
      continue;
    }

    /* Fill in the tuple */
    if(GTS_IS_FACE(t->data)) {
      triangle = PYGTS_TRIANGLE(pygts_face_new(GTS_FACE(t->data)));
    }
    else {
      triangle = pygts_triangle_new(GTS_TRIANGLE(t->data));
    }
    if( triangle == NULL ) {
      Py_DECREF(tuple);
      g_slist_free(triangles);
      return NULL;
    }
    PyTuple_SET_ITEM(tuple,n,(PyObject*)triangle);
    t = g_slist_next(t);
    n += 1;
  }

  g_slist_free(triangles);

  if(_PyTuple_Resize(&tuple,n)!=0) {
    Py_DECREF(tuple);
    return NULL;
  }

  return tuple;
}


static PyObject*
triangle_enclosing(PyObject *self, PyObject *args)
{
  PyObject *tuple, *obj;
  guint i,N;
  GSList *points=NULL;
  GtsTriangle *t;
  PygtsTriangle *triangle;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &tuple) ) {
    return NULL;
  }
  if(PyList_Check(tuple)) {
    tuple = PyList_AsTuple(tuple);
  }
  else {
    Py_INCREF(tuple);
  }
  if(!PyTuple_Check(tuple)) {
    Py_DECREF(tuple);
    PyErr_SetString(PyExc_TypeError,"expected a list or tuple of points");
    return NULL;
  }

  /* Assemble the GSList */
  N = PyTuple_Size(tuple);
  for(i=0;i<N;i++) {
    obj = PyTuple_GET_ITEM(tuple,i);
    if(!pygts_point_check(obj)) {
      Py_DECREF(tuple);
      g_slist_free(points);
      PyErr_SetString(PyExc_TypeError,"expected a list or tuple of points");
      return NULL;
    }
    points = g_slist_prepend(points,PYGTS_POINT_AS_GTS_POINT(obj));
  }
  Py_DECREF(tuple);

  /* Make the call */
  t = gts_triangle_enclosing(gts_triangle_class(),points,1.0);
  g_slist_free(points);

  if(t==NULL) {
    PyErr_SetString(PyExc_RuntimeError,"could not compute triangle");
    return NULL;
  }

  if( (triangle = pygts_triangle_new(t)) == NULL ) {
    return NULL;
  }

  return (PyObject*)triangle;
}


static PyObject*
pygts_read(PygtsSurface *self, PyObject *args)
{
  PyObject *f_;
  FILE *f;
  GtsFile *fp;
  guint lineno;
  GtsSurface *s;
  PygtsSurface *surface;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "O", &f_) )
    return NULL;

  /* Convert to PygtsObjects */
  if(!PyFile_Check(f_)) {
    PyErr_SetString(PyExc_TypeError,"expected a File");
    return NULL;
  }
  f = PyFile_AsFile(f_);

  if(feof(f)) {
    PyErr_SetString(PyExc_EOFError,"End of File");
    return NULL;
  }

  /* Create a temporary surface to read into */
  if( (s = gts_surface_new(gts_surface_class(), gts_face_class(),
			   gts_edge_class(), gts_vertex_class())) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Surface");
    return NULL;
  }

  /* Read from the file */
  fp = gts_file_new(f);
  if( (lineno = gts_surface_read(s,fp)) != 0 ) {
    PyErr_SetString(PyExc_RuntimeError,fp->error);
    gts_file_destroy(fp);
    return NULL;
  }
  gts_file_destroy(fp);
  if( (surface = pygts_surface_new(s)) == NULL )  {
    gts_object_destroy(GTS_OBJECT(s));
    return NULL;
  }

  /* Clean up the surface */
  pygts_edge_cleanup(PYGTS_SURFACE_AS_GTS_SURFACE(surface));
  pygts_face_cleanup(PYGTS_SURFACE_AS_GTS_SURFACE(surface));

  return (PyObject*)surface;
}


static PyObject*
disk(PyObject *self, PyObject *args)
{
  PyObject *kwds;
  PygtsSurface *surface;
  guint geodesation_order;

  /* Parse the args */  
  if(! PyArg_ParseTuple(args, "i", &geodesation_order) )
    return NULL;

  /* Chain up object allocation */
  args = Py_BuildValue("()");
  kwds = Py_BuildValue("{s:O}","alloc_gtsobj",Py_True);
  surface = PYGTS_SURFACE(PygtsSurfaceType.tp_new(&PygtsSurfaceType, 
						  args, kwds));
  Py_DECREF(args);
  Py_DECREF(kwds);
  if( surface == NULL ) {
    PyErr_SetString(PyExc_MemoryError, "could not create Surface");
    return NULL;
  }

  gts_surface_generate_disk(PYGTS_SURFACE_AS_GTS_SURFACE(surface),
			      geodesation_order);

  pygts_object_register(PYGTS_OBJECT(surface));
  return (PyObject*)surface;
}


#if PYGTS_HAS_NUMPY

/* Helper for pygts_iso to fill f with a layer of data from scalar */
static void isofunc(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
  PyArrayObject *scalars = (PyArrayObject *)data;
  int i, j;

  for (i = 0; i < scalars->dimensions[0]; i++) {
    for (j = 0; j < scalars->dimensions[1]; j++) {
      f[i][j] = *(gdouble *)(scalars->data + i*scalars->strides[0] + \
			     j*scalars->strides[1] + k*scalars->strides[2]);
    }
  }
}

#define ISO_CLEANUP \
  if (scalars) { Py_DECREF(scalars); } \
  if (extents) { Py_DECREF(extents); }

static PyObject*
isosurface(PyObject *self, PyObject *args, PyObject *kwds)
{
  double isoval[1];
  PyObject *Oscalars = NULL, *Oextents = NULL;
  PyArrayObject *scalars = NULL, *extents = NULL;
  GtsCartesianGrid g;
  GtsSurface *s;
  PygtsSurface *surface;
  char *method = "cubes";
  
  static char *kwlist[] = {"scalars", "isoval", "method", "extents", NULL};

  if(!PyArg_ParseTupleAndKeywords(args, kwds, "Od|sO", kwlist, 
				  &Oscalars, isoval, &method, &Oextents)) {
    return NULL;
  }
  
  if(!(scalars = (PyArrayObject *) 
       PyArray_ContiguousFromObject(Oscalars, PyArray_DOUBLE, 3, 3))) {
    ISO_CLEANUP;
    return NULL;
  }

  if(Oextents && 
     (!(extents =  (PyArrayObject *)
	PyArray_ContiguousFromObject(Oextents, PyArray_DOUBLE, 1, 1)))) {
    ISO_CLEANUP;
    return NULL;
  }

  if(extents && extents->dimensions[0] < 6) {
    PyErr_SetString(PyExc_ValueError, "extents must have at least 6 elements");
    ISO_CLEANUP;
    return NULL;
  }
  
  if(extents) {
    int s = extents->strides[0];
    g.x = *(gdouble*)(extents->data + 0*s);
    g.nx = scalars->dimensions[0];
    g.dx = (*(gdouble*)(extents->data + 1*s) - \
	    *(gdouble*)(extents->data + 0* s))/(g.nx-1);

    g.y = *(gdouble*)(extents->data + 2*s);
    g.ny = scalars->dimensions[1];
    g.dy = (*(gdouble*)(extents->data + 3*s) - \
	    *(gdouble*)(extents->data + 2*s))/(g.ny-1);

    g.z = *(gdouble*)(extents->data + 4*s);
    g.nz = scalars->dimensions[2];
    g.dz = (*(gdouble*)(extents->data + 5*s) - \
	    *(gdouble*)(extents->data + 4*s))/(g.nz-1);
  }
  else {
    g.x = -1.0;
    g.nx = scalars->dimensions[0];
    g.dx = 2.0/(scalars->dimensions[0]-1);
    g.y = -1.0;
    g.ny = scalars->dimensions[1];
    g.dy = 2.0/(scalars->dimensions[1]-1);
    g.z = -1.0;
    g.nz = scalars->dimensions[2];
    g.dz = 2.0/(scalars->dimensions[2]-1);
  }

  /* Create the surface */
  if((s = gts_surface_new(gts_surface_class(), gts_face_class(),
			  gts_edge_class(), gts_vertex_class())) == NULL ) {
    PyErr_SetString(PyExc_MemoryError,"could not create Surface");
    return NULL;
  }

  /* Make the call */
  switch(method[0]) {
  case 'c': /* cubes */
    gts_isosurface_cartesian(s, g, isofunc, scalars, isoval[0]);
    break;
  case 't': /* tetra */
    gts_isosurface_tetra(s, g, isofunc, scalars, isoval[0]);
    /* *** ATTENTION ***
     * Isosurface produced is "inside-out", and so we must revert it.
     * This is a bug in GTS.
     */
    gts_surface_foreach_face(s, (GtsFunc)gts_triangle_revert, NULL);
    /* *** ATTENTION *** */
    break;
  case 'b': /* tetra bounded */
    gts_isosurface_tetra_bounded(s, g, isofunc, scalars, isoval[0]);
    /* *** ATTENTION ***
     * Isosurface produced is "inside-out", and so we must revert it.
     * This is a bug in GTS.
     */
    gts_surface_foreach_face(s, (GtsFunc)gts_triangle_revert, NULL);
    /* *** ATTENTION *** */
    break;
  case 'd': /* tetra bcl*/
    gts_isosurface_tetra_bcl(s, g, isofunc, scalars, isoval[0]);
    /* *** ATTENTION ***
     * Isosurface produced is "inside-out", and so we must revert it.
     * This is a bug in GTS.
     */
    gts_surface_foreach_face(s, (GtsFunc)gts_triangle_revert, NULL);
    /* *** ATTENTION *** */
    break;
  default:
    PyErr_SetString(PyExc_ValueError, "unknown method");
    ISO_CLEANUP;
    return NULL;
  }    

  ISO_CLEANUP;

  if( (surface = pygts_surface_new(s)) == NULL )  {
    gts_object_destroy(GTS_OBJECT(s));
    return NULL;
  }

  return (PyObject*)surface;
}

#endif /* PYGTS_HAS_NUMPY */


static PyMethodDef gts_methods[] = {

  {"read", (PyCFunction)pygts_read,
   METH_VARARGS,
   "Returns the data read from File f as a Surface.\n"
   "The File data must be in GTS format (e.g., as written using\n"
   "Surface.write())\n"
   "\n"
   "Signature: read(f)\n"
  },

  { "disk", disk, METH_VARARGS,
    "Returns a unit disk generated by recursive subdivision.\n"
    "First approximation is an isocahedron; each level of refinement\n"
    "(geodesation_order) increases the number of triangles by a factor\n"
    "of 4.\n"
    "\n"
    "Signature: disk(geodesation_order)\n"
  },

#if PYGTS_HAS_NUMPY
  {"isosurface",  (PyCFunction)isosurface, 
   METH_VARARGS|METH_KEYWORDS,
   "Adds to surface new faces defining the isosurface data[x,y,z] = c\n"
   "\n"
   "Signature: isosurface(data, c, ...)\n"
   "\n"
   "data is a 3D numpy array.\n"
   "c    is the isovalue defining the surface\n"
   "\n"
   "Keyword arguments:\n"
   "extents= [xmin, xmax, ymin, ymax, zmin, zmax]\n"
   "         A numpy array defining the extent of the data cube.\n" 
   "         Default is the cube with corners at (-1,-1,-1) and (1,1,1)\n"
   "         Data is assumed to be regularly sampled in the cube.\n"
   "method=  ['cube'|'tetra'|'dual'|'bounded']\n"
   "         String (only the first character counts) specifying the\n"
   "         method.\n"
   "         cube    -- marching cubes (default)\n"
   "         tetra   -- marching tetrahedra\n"
   "         dual    -- maching tetrahedra producing dual 'body-centred'\n"
   "                    faces relative to 'tetra'\n"
   "         bounded -- marching tetrahedra ensuring the surface is\n"
   "                    bounded by adding a border of large negative\n"
   "                    values around the domain.\n"
   "\n"
   "By convention, the normals to the surface are pointing towards\n"
   "positive values of data[x,y,z] - c.\n"
  },
#endif /* PYGTS_HAS_NUMPY */

  { "merge", merge, METH_VARARGS,
    "Merges list of Vertices that are within a box of side-length\n"
    "epsilon of each other.\n"
    "\n"
    "Signature: merge(list,epsilon)\n"
  },

  { "vertices", vertices, METH_VARARGS,
    "Returns tuple of Vertices from a list or tuple of Segments.\n"
    "\n"
    "Signature: vertices(list)\n"
  },

  { "segments", segments, METH_VARARGS,
    "Returns tuple of Segments from a list or tuple of Vertices.\n"
    "\n"
    "Signature: segments(list)\n"
  },

  { "triangles", triangles, METH_VARARGS,
    "Returns tuple of Triangles from a list or tuple of Edges.\n"
    "\n"
    "Signature: triangles(list)\n"
  },

  { "triangle_enclosing", triangle_enclosing, METH_VARARGS,
    "Returns a Triangle that encloses the plane projection of a list\n"
    "or tuple of Points.  The Triangle is equilateral and encloses a\n"
    "rectangle defined by the maximum and minimum x and y coordinates\n"
    "of the points.\n"
    "\n"
    "Signature: triangles(list)\n"
  },


  {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_gts(void) 
{
  PyObject* m;

  /* Allocate the object table */
  if( (obj_table=g_hash_table_new(NULL,NULL)) == NULL ) return;

  /* Set class base types and make ready (i.e., inherit methods) */
  if (PyType_Ready(&PygtsObjectType) < 0) return;

  PygtsPointType.tp_base = &PygtsObjectType;
  if (PyType_Ready(&PygtsPointType) < 0) return;

  PygtsVertexType.tp_base = &PygtsPointType;
  if (PyType_Ready(&PygtsVertexType) < 0) return;

  PygtsSegmentType.tp_base = &PygtsObjectType;
  if (PyType_Ready(&PygtsSegmentType) < 0) return;

  PygtsEdgeType.tp_base = &PygtsSegmentType;
  if (PyType_Ready(&PygtsEdgeType) < 0) return;

  PygtsTriangleType.tp_base = &PygtsObjectType;
  if (PyType_Ready(&PygtsTriangleType) < 0) return;

  PygtsFaceType.tp_base = &PygtsTriangleType;
  if (PyType_Ready(&PygtsFaceType) < 0) return;

  PygtsSurfaceType.tp_base = &PygtsObjectType;
  if (PyType_Ready(&PygtsSurfaceType) < 0) return;


  /* Initialize the module */
  m = Py_InitModule3("_gts", gts_methods,"Gnu Triangulated Surface Library");
  if (m == NULL) return;

#if PYGTS_HAS_NUMPY
  /* Make sure Surface.iso can work with numpy arrays */
  import_array()
#endif

  /* Add new types to python */
  Py_INCREF(&PygtsObjectType);
  PyModule_AddObject(m, "Object", (PyObject *)&PygtsObjectType);

  Py_INCREF(&PygtsPointType);
  PyModule_AddObject(m, "Point", (PyObject *)&PygtsPointType);

  Py_INCREF(&PygtsVertexType);
  PyModule_AddObject(m, "Vertex", (PyObject *)&PygtsVertexType);

  Py_INCREF(&PygtsSegmentType);
  PyModule_AddObject(m, "Segment", (PyObject *)&PygtsSegmentType);

  Py_INCREF(&PygtsEdgeType);
  PyModule_AddObject(m, "Edge", (PyObject *)&PygtsEdgeType);

  Py_INCREF(&PygtsTriangleType);
  PyModule_AddObject(m, "Triangle", (PyObject *)&PygtsTriangleType);

  Py_INCREF(&PygtsFaceType);
  PyModule_AddObject(m, "Face", (PyObject *)&PygtsFaceType);

  Py_INCREF(&PygtsSurfaceType);
  PyModule_AddObject(m, "Surface", (PyObject *)&PygtsSurfaceType);
}
