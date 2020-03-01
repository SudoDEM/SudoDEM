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

#ifndef __PYGTS_SEGMENT_H__
#define __PYGTS_SEGMENT_H__

typedef struct _PygtsObject PygtsSegment;

#define PYGTS_SEGMENT(obj) ((PygtsSegment*)obj)

#define PYGTS_SEGMENT_AS_GTS_SEGMENT(o) (GTS_SEGMENT(PYGTS_OBJECT(o)->gtsobj))

extern PyTypeObject PygtsSegmentType;

gboolean pygts_segment_check(PyObject* o);
gboolean pygts_segment_is_ok(PygtsSegment *t);

PygtsSegment* pygts_segment_new(GtsSegment *f);

int pygts_segment_compare(GtsSegment* s1,GtsSegment* s2);

#endif /* __PYGTS_SEGMENT_H__ */
