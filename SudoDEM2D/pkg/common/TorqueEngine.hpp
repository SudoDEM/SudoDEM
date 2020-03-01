/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <sudodem/core/PartialEngine.hpp>
#include <sudodem/core/Scene.hpp>

class TorqueEngine: public PartialEngine{
	public:
		virtual void action() {
			FOREACH(const Body::id_t id, ids){
			// check that body really exists?
			scene->forces.addTorque(id,moment);
		}
	}
	SUDODEM_CLASS_BASE_DOC_ATTRS(TorqueEngine,PartialEngine,"Apply given torque (momentum) value at every subscribed particle, at every step.",
		((Real,moment,0.0,,"Torque value to be applied."))
	);
};
REGISTER_SERIALIZABLE(TorqueEngine);
