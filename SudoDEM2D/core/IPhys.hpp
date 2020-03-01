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

class IPhys : public Serializable, public Indexable
{
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(IPhys,Serializable,"Physical (material) properties of :yref:`interaction<Interaction>`.",
		/*attrs*/,
		/*ctor*/,
		/*py*/SUDODEM_PY_TOPINDEXABLE(IPhys)
	);
	REGISTER_INDEX_COUNTER(IPhys);
};
REGISTER_SERIALIZABLE(IPhys);


