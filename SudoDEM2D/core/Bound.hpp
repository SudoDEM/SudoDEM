/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include<sudodem/lib/base/Math.hpp>
#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/lib/multimethods/Indexable.hpp>
#include<sudodem/core/Dispatcher.hpp>

/*! Interface for approximate body locations in space

	Note: the min and max members refer to shear coordinates, in periodic
	and sheared space, not cartesian coordinates in the physical space.

*/

class Bound: public Serializable, public Indexable{
	public:
	SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(Bound,Serializable,"Object bounding part of space taken by associated body; might be larger, used to optimalize collision detection",
		((int,lastUpdateIter,0,Attr::readonly,"record iteration of last reference position update |yupdate|"))
		((Vector2r,refPos,Vector2r(NaN,NaN),Attr::readonly,"Reference position, updated at current body position each time the bound dispatcher update bounds |yupdate|"))
// 		((bool,isBounding,false,Attr::readonly,"A flag used to tell when the body moves out of bounds - only used if oriVerlet striding is active :yref:`BoundDispatcher::updatingDispFactor`>0 |yupdate|"))
		((Real,sweepLength,0, Attr::readonly,"The length used to increase the bounding boxe size, can be adjusted on the basis of previous displacement if :yref:`BoundDispatcher::targetInterv`>0. |yupdate|"))
		((Vector3r,color,Vector3r(1,1,1),,"Color for rendering this object"))
		((Vector2r,min,Vector2r(NaN,NaN),(Attr::noSave | Attr::readonly),"Lower corner of box containing this bound (and the :yref:`Body` as well)"))
		((Vector2r,max,Vector2r(NaN,NaN),(Attr::noSave | Attr::readonly),"Upper corner of box containing this bound (and the :yref:`Body` as well)"))
		,
		/*deprec*/,
		/* init */,
		/* ctor*/,
		/*py*/
		SUDODEM_PY_TOPINDEXABLE(Bound)
	);
	REGISTER_INDEX_COUNTER(Bound);
};
REGISTER_SERIALIZABLE(Bound);
