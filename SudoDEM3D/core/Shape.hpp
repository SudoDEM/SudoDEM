/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/lib/multimethods/Indexable.hpp>
#include<sudodem/core/Dispatcher.hpp>

#define BV_FUNCTOR_CACHE

class BoundFunctor;

class Shape: public Serializable, public Indexable {
	public:
		~Shape() {}; // vtable
		#ifdef BV_FUNCTOR_CACHE
			shared_ptr<BoundFunctor> boundFunctor;
		#endif

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Shape,Serializable,"Geometry of a body",
		((Vector3r,color,Vector3r(1,1,1),,"Color for rendering (normalized RGB)."))
		((bool,wire,false,,"Whether this Shape is rendered using color surfaces, or only wireframe (can still be overridden by global config of the renderer)."))
		((bool,highlight,false,,"Whether this Shape will be highlighted when rendered.")),
		/*ctor*/,
		/*py*/ SUDODEM_PY_TOPINDEXABLE(Shape)
	);
	REGISTER_INDEX_COUNTER(Shape);
};
REGISTER_SERIALIZABLE(Shape);

