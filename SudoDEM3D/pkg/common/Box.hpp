/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#ifndef BOX_H
#define BOX_H
#include <sudodem/pkg/common/Dispatching.hpp>
#include <sudodem/pkg/common/Aabb.hpp>
#include<sudodem/core/Shape.hpp>
#include<sudodem/pkg/common/GLDrawFunctors.hpp>

class Box: public Shape{
	public:
		Box(const Vector3r& _extents): extents(_extents){}
		virtual ~Box () {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(Box,Shape,"Box (cuboid) particle geometry. (Avoid using in new code, prefer :yref:`Facet` instead.",
		((Vector3r,extents,,,"Half-size of the cuboid")),
		/* ctor */ createIndex();
	);
	REGISTER_CLASS_INDEX(Box,Shape);
};
REGISTER_SERIALIZABLE(Box);

class Bo1_Box_Aabb : public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);
	FUNCTOR1D(Box);
	SUDODEM_CLASS_BASE_DOC(Bo1_Box_Aabb,BoundFunctor,"Create/update an :yref:`Aabb` of a :yref:`Box`.");
};

REGISTER_SERIALIZABLE(Bo1_Box_Aabb);

class Gl1_Box : public GlShapeFunctor{
	public :
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
	RENDERS(Box);
	SUDODEM_CLASS_BASE_DOC(Gl1_Box,GlShapeFunctor,"Renders :yref:`Box` object");
};
#endif
