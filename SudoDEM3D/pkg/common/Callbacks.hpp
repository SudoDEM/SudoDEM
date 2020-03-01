// 2010 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include<sudodem/lib/serialization/Serializable.hpp>

class Interaction;
class Body;
class Scene;

class IntrCallback: public Serializable{
	public:
	virtual ~IntrCallback() {}; // vtable
	typedef void(*FuncPtr)(IntrCallback*,Interaction*);
	// should be set at every step by InteractionLoop
	Scene* scene;
	/*
	At the beginning of each timestep, perform initialization and
	return pointer to the static member function that does the actual work.
	Returned value might be NULL, in which case the callback will be deactivated during that timestep.
	*/
	virtual FuncPtr stepInit(){ throw std::runtime_error("Called IntrCallback::stepInit() of the base class?"); }
	SUDODEM_CLASS_BASE_DOC(IntrCallback,Serializable,"Abstract callback object which will be called for every (real) :yref:`Interaction` after the interaction has been processed by :yref:`InteractionLoop`.\n\nAt the beginning of the interaction loop, ``stepInit`` is called, initializing the object; it returns either ``NULL`` (to deactivate the callback during this time step) or pointer to function, which will then be passed (1) pointer to the callback object itself and (2) pointer to :yref:`Interaction`.\n\n.. note::\n\t(NOT YET DONE) This functionality is accessible from python by passing 4th argument to :yref:`InteractionLoop` constructor, or by appending the callback object to :yref:`InteractionLoop::callbacks`.\n");
};
REGISTER_SERIALIZABLE(IntrCallback);

#ifdef SUDODEM_BODY_CALLBACKS
	class BodyCallback: public Serializable{
		public:
		virtual ~BodyCallback() {}; // vtable
		typedef void(*FuncPtr)(BodyCallback*,Body*);
		// set at every step, before stepInit() is called
		Scene* scene;
		virtual FuncPtr stepInit(){ throw std::runtime_error("Called BodyCallback::stepInit() of the base class?"); }
		SUDODEM_CLASS_BASE_DOC(BodyCallback,Serializable,"Abstract callback object which will be called for every :yref:`Body` after being processed by :yref:`NewtonIntegrator`. See :yref:`IntrCallback` for details.");
	};
	REGISTER_SERIALIZABLE(BodyCallback);
#endif
