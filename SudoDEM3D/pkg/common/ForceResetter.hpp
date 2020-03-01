#pragma once

#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/core/Scene.hpp>

class Scene;
class ForceResetter: public GlobalEngine{
	public:
		virtual void action() {
			scene->forces.reset(scene->iter);
			if(scene->trackEnergy) scene->energy->resetResettables();
		}
	SUDODEM_CLASS_BASE_DOC(ForceResetter,GlobalEngine,"Reset all forces stored in Scene::forces (``O.forces`` in python). Typically, this is the first engine to be run at every step. In addition, reset those energies that should be reset, if energy tracing is enabled.");
};
REGISTER_SERIALIZABLE(ForceResetter);


