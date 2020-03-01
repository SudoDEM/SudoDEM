/* pygts - python package for the manipulation of triangulated surfaces
 *
 *   Copyright (C) 1999 Stéphane Popinet
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


/*
 * Below are functions for cleaning up duplicated edges and faces on
 * a surface.  This file contains modified functions from the GTS 
 * distribution.
 */

#include "pygts.h"


/**
 * Original documentation from GTS's vertex.c:
 *
 * gts_vertices_merge:
 * @vertices: a list of #GtsVertex.
 * @epsilon: half the size of the bounding box to consider for each vertex.
 * @check: function called for each pair of vertices about to be merged
 * or %NULL.
 *
 * For each vertex v in @vertices look if there are any vertex of
 * @vertices contained in a box centered on v of size 2*@epsilon. If
 * there are and if @check is not %NULL and returns %TRUE, replace
 * them with v (using gts_vertex_replace()), destroy them and remove
 * them from list.  This is done efficiently using Kd-Trees.
 *
 * Returns: the updated list of vertices.  
 */
/* This function is modified from the original in GTS in order to avoid 
 * deallocating any objects referenced by the live-objects table.  The 
 * approach is similar to what is used for replace() in vertex.c.
 */
GList*
pygts_vertices_merge(GList* vertices, gdouble epsilon,
		     gboolean (* check) (GtsVertex *, GtsVertex *))
{
  GPtrArray *array;
  GList *i, *next;
  GNode *kdtree;
  GtsVertex *v;
  GtsBBox *bbox;
  GSList *selected, *j;
  GtsVertex *sv;
  PygtsObject *obj;
  PygtsVertex *vertex=NULL;
  GSList *parents=NULL, *ii,*cur;

  g_return_val_if_fail(vertices != NULL, 0);

  array = g_ptr_array_new();
  i = vertices;
  while (i) {
    g_ptr_array_add(array, i->data);
    i = g_list_next(i);
  }
  kdtree = gts_kdtree_new(array, NULL);
  g_ptr_array_free(array, TRUE);

  i = vertices;
  while(i) {
    v = i->data;
    if (!GTS_OBJECT(v)->reserved) { /* Do something only if v is active */

      /* build bounding box */
      bbox = gts_bbox_new(gts_bbox_class(), v, 
			  GTS_POINT(v)->x - epsilon,
			  GTS_POINT(v)->y - epsilon,
			  GTS_POINT(v)->z - epsilon,
			  GTS_POINT(v)->x + epsilon,
			  GTS_POINT(v)->y + epsilon,
			  GTS_POINT(v)->z + epsilon);

      /* select vertices which are inside bbox using kdtree */
      j = selected = gts_kdtree_range(kdtree, bbox, NULL);
      while(j) {
        sv = j->data;
        if( sv!=v && !GTS_OBJECT(sv)->reserved && (!check||(*check)(sv, v)) ) {
          /* sv is not v and is active */
	  if( (obj = g_hash_table_lookup(obj_table,GTS_OBJECT(sv))) !=NULL ) {
	    vertex = PYGTS_VERTEX(obj);
	    /* Detach and save any parent segments */
	    ii = sv->segments;
	    while(ii!=NULL) {
	      cur = ii;
	      ii = g_slist_next(ii);
	      if(PYGTS_IS_PARENT_SEGMENT(cur->data)) {
		sv->segments = g_slist_remove_link(sv->segments, cur);
		parents = g_slist_prepend(parents,cur->data);
		g_slist_free_1(cur);
	      }
	    } 
	  }

          gts_vertex_replace(sv, v);
          GTS_OBJECT(sv)->reserved = sv; /* mark sv as inactive */

	  /* Reattach the parent segments */
	  if( vertex != NULL ) {
	    ii = parents;
	    while(ii!=NULL) {
	      sv->segments = g_slist_prepend(sv->segments, ii->data);
	      ii = g_slist_next(ii);
	    }
	    g_slist_free(parents);
	    parents = NULL;
	  }
	  vertex = NULL;
        }
        j = g_slist_next(j);
      }
      g_slist_free(selected);
      gts_object_destroy(GTS_OBJECT(bbox));
    }
    i = g_list_next(i);
  }
  gts_kdtree_destroy(kdtree);


  /* destroy inactive vertices and removes them from list */

  /* we want to control vertex destruction */
  gts_allow_floating_vertices = TRUE;

  i = vertices;
  while (i) {
    v = i->data;
    next = g_list_next(i);
    if(GTS_OBJECT(v)->reserved) { /* v is inactive */
      if( g_hash_table_lookup(obj_table,GTS_OBJECT(v))==NULL ) {
	gts_object_destroy(GTS_OBJECT(v));
      }
      else {
	GTS_OBJECT(v)->reserved = 0;
      }
      vertices = g_list_remove_link(vertices, i);
      g_list_free_1(i);
    }
    i = next;
  }
  gts_allow_floating_vertices = FALSE; 

  return vertices;
}


static void 
build_list(gpointer data, GSList ** list)
{
  *list = g_slist_prepend(*list, data);
}


static void
build_list1(gpointer data, GList ** list)
{
  *list = g_list_prepend(*list, data);
}


void 
pygts_vertex_cleanup(GtsSurface *s, gdouble threshold)
{
  GList * vertices = NULL;

  /* merge vertices which are close enough */
  /* build list of vertices */
  gts_surface_foreach_vertex(s, (GtsFunc) build_list1, &vertices);

  /* merge vertices: we MUST update the variable vertices because this function
     modifies the list (i.e. removes the merged vertices). */
  vertices = pygts_vertices_merge(vertices, threshold, NULL);

  /* free the list */
  g_list_free(vertices);
}


void 
pygts_edge_cleanup(GtsSurface *s)
{
  GSList *edges = NULL;
  GSList *i, *ii, *cur, *parents=NULL;
  PygtsEdge *edge;
  GtsEdge *e, *duplicate;

  g_return_if_fail(s != NULL);

  /* build list of edges */
  gts_surface_foreach_edge(s, (GtsFunc)build_list, &edges);

  /* remove degenerate and duplicate edges.
     Note: we could use gts_edges_merge() to remove the duplicates and then
     remove the degenerate edges but it is more efficient to do everything 
     at once (and it's more pedagogical too ...) */

  /* We want to control manually the destruction of edges */
  gts_allow_floating_edges = TRUE;

  i = edges;
  while(i) {
    e = i->data;
    if(GTS_SEGMENT(e)->v1 == GTS_SEGMENT(e)->v2) {
      /* edge is degenerate */
      if( !g_hash_table_lookup(obj_table,GTS_OBJECT(e)) ) {
	/* destroy e */
	gts_object_destroy(GTS_OBJECT(e));
      }
    }
    else {
      if((duplicate = gts_edge_is_duplicate(e))) {

	/* Detach and save any parent triangles */
	if( (edge = PYGTS_EDGE(g_hash_table_lookup(obj_table,GTS_OBJECT(e))))
	    !=NULL ) {
	  ii = e->triangles;
	  while(ii!=NULL) {
	    cur = ii;
	    ii = g_slist_next(ii);
	    if(PYGTS_IS_PARENT_TRIANGLE(cur->data)) {
	      e->triangles = g_slist_remove_link(e->triangles, cur);
	      parents = g_slist_prepend(parents,cur->data);
	      g_slist_free_1(cur);
	    }
	  } 
	}

	/* replace e with its duplicate */
	gts_edge_replace(e, duplicate);

	/* Reattach the parent segments */
	if( edge != NULL ) {
	  ii = parents;
	  while(ii!=NULL) {
	    e->triangles = g_slist_prepend(e->triangles, ii->data);
	    ii = g_slist_next(ii);
	  }
	  g_slist_free(parents);
	  parents = NULL;
	}

	if( !g_hash_table_lookup(obj_table,GTS_OBJECT(e)) ) {
	  /* destroy e */
	  gts_object_destroy(GTS_OBJECT (e));
	}
      }
    }
    i = g_slist_next(i);
  }
  
  /* don't forget to reset to default */
  gts_allow_floating_edges = FALSE;

  /* free list of edges */
  g_slist_free (edges);
}


void 
pygts_face_cleanup(GtsSurface * s)
{
  GSList *triangles = NULL;
  GSList * i;

  g_return_if_fail(s != NULL);

  /* build list of triangles */
  gts_surface_foreach_face(s, (GtsFunc) build_list, &triangles);
  
  /* remove duplicate and degenerate triangles */
  i = triangles;
  while(i) {
    GtsTriangle * t = i->data;
    if (!gts_triangle_is_ok(t)) {
      /* destroy t, its edges (if not used by any other triangle)
	 and its corners (if not used by any other edge) */
      if( g_hash_table_lookup(obj_table,GTS_OBJECT(t))==NULL ) {
	gts_object_destroy(GTS_OBJECT(t));
      }
      else {
	gts_surface_remove_face(PYGTS_SURFACE_AS_GTS_SURFACE(s),GTS_FACE(t));
      }
    }
    i = g_slist_next(i);
  }
  
  /* free list of triangles */
  g_slist_free(triangles);
}


/* old main program (below) - retained as an example of how to use the
 * functions above.
 */

/* cleanup - using a given threshold merge vertices which are too close.
   Eliminate degenerate and duplicate edges.
   Eliminate duplicate triangles . */
/* static int  */
/* main (int argc, char * argv[]) */
/* { */
/*   GtsSurface * s, * m; */
/*   GList * vertices = NULL; */
/*   gboolean verbose = FALSE, sever = FALSE, boundary = FALSE; */
/*   gdouble threshold; */
/*   int c = 0; */
/*   GtsFile * fp; */
/*   gboolean (* check) (GtsVertex *, GtsVertex *) = NULL; */

/*   if (!setlocale (LC_ALL, "POSIX")) */
/*     g_warning ("cannot set locale to POSIX"); */

/*   s = gts_surface_new (gts_surface_class (), */
/* 		       gts_face_class (), */
/* 		       gts_edge_class (), */
/* 		       gts_vertex_class ()); */

/*   /\* parse options using getopt *\/ */
/*   while (c != EOF) { */
/* #ifdef HAVE_GETOPT_LONG */
/*     static struct option long_options[] = { */
/*       {"2D", no_argument, NULL, 'c'}, */
/*       {"boundary", no_argument, NULL, 'b'}, */
/*       {"merge", required_argument, NULL, 'm'}, */
/*       {"sever", no_argument, NULL, 's'}, */
/*       {"help", no_argument, NULL, 'h'}, */
/*       {"verbose", no_argument, NULL, 'v'}, */
/*       { NULL } */
/*     }; */
/*     int option_index = 0; */
/*     switch ((c = getopt_long (argc, argv, "hvsm:bc", */
/* 			      long_options, &option_index))) { */
/* #else /\* not HAVE_GETOPT_LONG *\/ */
/*     switch ((c = getopt (argc, argv, "hvsm:bc"))) { */
/* #endif /\* not HAVE_GETOPT_LONG *\/ */
/*     case 'c': /\* 2D *\/ */
/*       check = check_boundaries; */
/*       break; */
/*     case 'b': /\* boundary *\/ */
/*       boundary = TRUE; */
/*       break; */
/*     case 's': /\* sever *\/ */
/*       sever = TRUE; */
/*       break; */
/*     case 'm': { /\* merge *\/ */
/*       FILE * fptr = fopen (optarg, "rt"); */
/*       GtsFile * fp; */

/*       if (fptr == NULL) { */
/* 	fprintf (stderr, "cleanup: cannot open file `%s' for merging\n", */
/* 		 optarg); */
/* 	return 1; /\* failure *\/ */
/*       } */
/*       m = gts_surface_new (gts_surface_class (), */
/* 			   gts_face_class (), */
/* 			   gts_edge_class (), */
/* 			   s->vertex_class); */
/*       fp = gts_file_new (fptr); */
/*       if (gts_surface_read (m, fp)) { */
/* 	fprintf (stderr, "cleanup: file `%s' is not a valid GTS file\n", */
/* 		 optarg); */
/* 	fprintf (stderr, "%s:%d:%d: %s\n", */
/* 		 optarg, fp->line, fp->pos, fp->error); */
/* 	return 1; /\* failure *\/ */
/*       } */
/*       gts_file_destroy (fp); */
/*       fclose (fptr); */
/*       gts_surface_merge (s, m); */
/*       gts_object_destroy (GTS_OBJECT (m)); */
/*       break; */
/*     } */
/*     case 'v': /\* verbose *\/ */
/*       verbose = TRUE; */
/*       break; */
/*     case 'h': /\* help *\/ */
/*       fprintf (stderr, */
/* 	       "Usage: cleanup [OPTION] THRESHOLD < FILE\n" */
/* 	       "Merge vertices of the GTS surface FILE if they are closer than THRESHOLD,\n" */
/* 	       "eliminate degenerate, duplicate edges and duplicate triangles.\n" */
/* 	       "\n" */
/* 	       "  -c      --2D        2D boundary merging\n" */
/* 	       "  -b      --boundary  only consider boundary vertices for merging\n" */
/* 	       "  -s      --sever     sever \"contact\" vertices\n" */
/* 	       "  -m FILE --merge     merge surface FILE\n" */
/* 	       "  -v      --verbose   print statistics about the surface\n" */
/* 	       "  -h      --help      display this help and exit\n" */
/* 	       "\n" */
/* 	       "Report bugs to %s\n", */
/* 	       GTS_MAINTAINER); */
/*       return 0; /\* success *\/ */
/*       break; */
/*     case '?': /\* wrong options *\/ */
/*       fprintf (stderr, "Try `cleanup --help' for more information.\n"); */
/*       return 1; /\* failure *\/ */
/*     } */
/*   } */

/*   if (optind >= argc) { /\* missing threshold *\/ */
/*     fprintf (stderr, */
/* 	     "cleanup: missing THRESHOLD\n" */
/* 	     "Try `cleanup --help' for more information.\n"); */
/*     return 1; /\* failure *\/ */
/*   } */

/*   threshold = atof (argv[optind]); */

/*   if (threshold < 0.0) { /\* threshold must be positive *\/ */
/*      fprintf (stderr, */
/* 	     "cleanup: THRESHOLD must be >= 0.0\n" */
/* 	     "Try `cleanup --help' for more information.\n"); */
/*     return 1; /\* failure *\/ */
/*   } */

/*   /\* read surface in *\/ */
/*   m = gts_surface_new (gts_surface_class (), */
/* 		       gts_face_class (), */
/* 		       gts_edge_class (), */
/* 		       s->vertex_class); */
/*   fp = gts_file_new (stdin); */
/*   if (gts_surface_read (m, fp)) { */
/*     fputs ("cleanup: file on standard input is not a valid GTS file\n", */
/* 	   stderr); */
/*     fprintf (stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error); */
/*     return 1; /\* failure *\/ */
/*   } */
/*   gts_surface_merge (s, m); */
/*   gts_object_destroy (GTS_OBJECT (m)); */

/*   /\* if verbose on print stats *\/ */
/*   if (verbose) */
/*     gts_surface_print_stats (s, stderr); */
 
/*   /\* merge vertices which are close enough *\/ */
/*   /\* build list of vertices *\/ */
/*   gts_surface_foreach_vertex (s, boundary ? (GtsFunc) build_list2 : (GtsFunc) build_list1, */
/* 			      &vertices); */
/*   /\* merge vertices: we MUST update the variable vertices because this function */
/*      modifies the list (i.e. removes the merged vertices). *\/ */
/*   vertices = gts_vertices_merge (vertices, threshold, check); */

/*   /\* free the list *\/ */
/*   g_list_free (vertices); */

/*   /\* eliminate degenerate and duplicate edges *\/ */
/*   edge_cleanup (s); */
/*   /\* eliminate duplicate triangles *\/ */
/*   triangle_cleanup (s); */

/*   if (sever) */
/*     gts_surface_foreach_vertex (s, (GtsFunc) vertex_cleanup, NULL); */

/*   /\* if verbose on print stats *\/ */
/*   if (verbose) { */
/*     GtsBBox * bb = gts_bbox_surface (gts_bbox_class (), s); */
/*     gts_surface_print_stats (s, stderr); */
/*     fprintf (stderr, "# Bounding box: [%g,%g,%g] [%g,%g,%g]\n", */
/* 	     bb->x1, bb->y1, bb->z1, */
/* 	     bb->x2, bb->y2, bb->z2); */
/*   } */

/*   /\* write surface *\/ */
/*   gts_surface_write (s, stdout); */

/*   return 0; /\* success *\/ */
/* } */
