// 2014 Burak ER <burak.er@btu.edu.tr>

#pragma once
#include<sudodem/pkg/common/InsertionSortCollider.hpp>


/*!
	Adaptive Integration Sort Collider:

	Changing the Integrator dependence from Newton Integrator to Arbitrary Integrators. Arbitrary integrators should use the Integrator interface.

*/



#ifdef ISC_TIMING
	#define ISC_CHECKPOINT(cpt) timingDeltas->checkpoint(cpt)
#else
	#define ISC_CHECKPOINT(cpt)
#endif

class Integrator;

class GeneralIntegratorInsertionSortCollider: public InsertionSortCollider{

	// we need this to find out about current maxVelocitySq
	shared_ptr<Integrator> integrator;
	// if False, no type of striding is used

	public:

	virtual bool isActivated(); //override this function to change NewtonIntegrator dependency.

	virtual void action(); //override this function to change behaviour with the NewtonIntegrator dependency.

	SUDODEM_CLASS_BASE_DOC(GeneralIntegratorInsertionSortCollider,InsertionSortCollider," This class is the adaptive version of the InsertionSortCollider and changes the NewtonIntegrator dependency of the collider algorithms to the Integrator interface which is more general.");
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(GeneralIntegratorInsertionSortCollider);
