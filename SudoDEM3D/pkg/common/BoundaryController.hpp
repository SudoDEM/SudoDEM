#pragma once
#include<sudodem/core/GlobalEngine.hpp>
class BoundaryController: public GlobalEngine{
	virtual void action() {
		{ throw std::runtime_error("BoundaryController must not be used in simulations directly (BoundaryController::action called)."); }
	}
	SUDODEM_CLASS_BASE_DOC(BoundaryController,GlobalEngine,"Base for engines controlling boundary conditions of simulations. Not to be used directly.");
};
REGISTER_SERIALIZABLE(BoundaryController);
