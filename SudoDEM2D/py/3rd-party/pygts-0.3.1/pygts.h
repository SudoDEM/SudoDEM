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

#ifndef __PYGTS_H__
#define __PYGTS_H__

#ifndef PYGTS_DEBUG
#define PYGTS_DEBUG 1
#endif /* PYGTS_DEBUG */

#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <structmember.h>

/* Defined for arrayobject.h which is only included where needed */
#define PY_ARRAY_UNIQUE_SYMBOL PYGTS

#include <glib.h>
#include <gts.h>

#include "object.h"
#include "point.h"
#include "vertex.h"
#include "segment.h"
#include "edge.h"
#include "triangle.h"
#include "face.h"
#include "surface.h"

#include "cleanup.h"

#endif /* __PYGTS_H__ */
