
#pragma once
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/core/PartialEngine.hpp>
//#include<sudodem/lib/base/Math.hpp>

struct KinematicEngine;

struct CombinedKinematicEngine: public PartialEngine{
	virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(CombinedKinematicEngine,PartialEngine,"Engine for applying combined displacements on pre-defined bodies. Constructed using ``+`` operator on regular :yref:`KinematicEngines<KinematicEngine>`. The ``ids`` operated on are those of the first engine in the combination (assigned automatically).",
		((vector<shared_ptr<KinematicEngine> >,comb,,,"Kinematic engines that will be combined by this one, run in the order given."))
		, /* ctor */
		, /* py */  .def("__add__",&CombinedKinematicEngine::appendOne)
	);
	// exposed as operator + in python
	static const shared_ptr<CombinedKinematicEngine> appendOne(const shared_ptr<CombinedKinematicEngine>& self, const shared_ptr<KinematicEngine>& other){ self->comb.push_back(other); return self; }
	static const shared_ptr<CombinedKinematicEngine> fromTwo(const shared_ptr<KinematicEngine>& first, const shared_ptr<KinematicEngine>& second);
};
REGISTER_SERIALIZABLE(CombinedKinematicEngine);

struct KinematicEngine: public PartialEngine{
	virtual void action();
	virtual void apply(const vector<Body::id_t>& ids){ LOG_ERROR("KinematicEngine::apply called, derived class ("<<getClassName()<<") did not override that method?"); }
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(KinematicEngine,PartialEngine,"Abstract engine for applying prescribed displacement.\n\n.. note:: Derived classes should override the ``apply`` with given list of ``ids`` (not ``action`` with :yref:`PartialEngine.ids`), so that they work when combined together; :yref:`velocity<State.vel>` and :yref:`angular velocity<State.angVel>` of all subscribed bodies is reset before the ``apply`` method is called, it should therefore only increment those quantities.",
		/* attrs*/, /* ctor */, /* py */ .def("__add__",&CombinedKinematicEngine::fromTwo)
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(KinematicEngine);


struct TranslationEngine: public KinematicEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	void postLoad(TranslationEngine&){ translationAxis.normalize(); }
	SUDODEM_CLASS_BASE_DOC_ATTRS(TranslationEngine,KinematicEngine,"This engine is the base class for different engines, which require any kind of motion.",
		((Real,velocity,,,"Velocity [m/s]"))
		((Vector3r,translationAxis,,Attr::triggerPostLoad,"Direction [Vector3]"))
	);
};
REGISTER_SERIALIZABLE(TranslationEngine);

struct HarmonicMotionEngine: public KinematicEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	SUDODEM_CLASS_BASE_DOC_ATTRS(HarmonicMotionEngine,KinematicEngine,"This engine implements the harmonic oscillation of bodies. http://en.wikipedia.org/wiki/Simple_harmonic_motion#Dynamics_of_simple_harmonic_motion",
		((Vector3r,A,Vector3r::Zero(),,"Amplitude [m]"))
		((Vector3r,f,Vector3r::Zero(),,"Frequency [hertz]"))
		((Vector3r,fi,Vector3r(Mathr::PI/2.0, Mathr::PI/2.0, Mathr::PI/2.0),,"Initial phase [radians]. By default, the body oscillates around initial position."))
	);
};
REGISTER_SERIALIZABLE(HarmonicMotionEngine);

struct RotationEngine: public KinematicEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	void postLoad(RotationEngine&){ rotationAxis.normalize(); }
	SUDODEM_CLASS_BASE_DOC_ATTRS(RotationEngine,KinematicEngine,"Engine applying rotation (by setting angular velocity) to subscribed bodies. If :yref:`rotateAroundZero<RotationEngine.rotateAroundZero>` is set, then each body is also displaced around :yref:`zeroPoint<RotationEngine.zeroPoint>`.",
		((Real,angularVelocity,0,,"Angular velocity. [rad/s]"))
		((Vector3r,rotationAxis,Vector3r::UnitX(),Attr::triggerPostLoad,"Axis of rotation (direction); will be normalized automatically."))
		((bool,rotateAroundZero,false,,"If True, bodies will not rotate around their centroids, but rather around ``zeroPoint``."))
		((Vector3r,zeroPoint,Vector3r::Zero(),,"Point around which bodies will rotate if ``rotateAroundZero`` is True"))
	);
};
REGISTER_SERIALIZABLE(RotationEngine);

struct HelixEngine:public RotationEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	SUDODEM_CLASS_BASE_DOC_ATTRS(HelixEngine,RotationEngine,"Engine applying both rotation and translation, along the same axis, whence the name HelixEngine",
		((Real,linearVelocity,0,,"Linear velocity [m/s]"))
		((Real,angleTurned,0,,"How much have we turned so far. |yupdate| [rad]"))
	);
};
REGISTER_SERIALIZABLE(HelixEngine);

struct InterpolatingHelixEngine: public HelixEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	SUDODEM_CLASS_BASE_DOC_ATTRS(InterpolatingHelixEngine,HelixEngine,"Engine applying spiral motion, finding current angular velocity by linearly interpolating in times and velocities and translation by using slope parameter. \n\n The interpolation assumes the margin value before the first time point and last value after the last time point. If wrap is specified, time will wrap around the last times value to the first one (note that no interpolation between last and first values is done).",
		((vector<Real>,times,,,"List of time points at which velocities are given; must be increasing [s]"))
		((vector<Real>,angularVelocities,,,"List of angular velocities; manadatorily of same length as times. [rad/s]"))
		((bool,wrap,false,,"Wrap t if t>times_n, i.e. t_wrapped=t-N*(times_n-times_0)"))
		((Real,slope,0,,"Axial translation per radian turn (can be negative) [m/rad]"))
		((size_t,_pos,0,(Attr::hidden),"holder of interpolation state, should not be touched by user"))
	);
};
REGISTER_SERIALIZABLE(InterpolatingHelixEngine);

struct HarmonicRotationEngine: public RotationEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	SUDODEM_CLASS_BASE_DOC_ATTRS(HarmonicRotationEngine,RotationEngine,"This engine implements the harmonic-rotation oscillation of bodies. http://en.wikipedia.org/wiki/Simple_harmonic_motion#Dynamics_of_simple_harmonic_motion ; please, set dynamic=False for bodies, droven by this engine, otherwise amplitude will be 2x more, than awaited.",
		((Real,A,0,,"Amplitude [rad]"))
		((Real,f,0,,"Frequency [hertz]"))
		((Real,fi,Mathr::PI/2.0,,"Initial phase [radians]. By default, the body oscillates around initial position."))
	);
};
REGISTER_SERIALIZABLE(HarmonicRotationEngine);

struct ServoPIDController: public TranslationEngine{
  virtual void apply(const vector<Body::id_t>& ids);
  SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(ServoPIDController,TranslationEngine,"PIDController servo-engine for applying prescribed force on bodies. http://en.wikipedia.org/wiki/PID_controller",
    ((Real,maxVelocity,0.0,,"Velocity [m/s]"))
    ((Vector3r,axis,Vector3r::Zero(),,"Unit vector along which apply the velocity [-]"))
    ((Real,target,0.0,,"Target value for the controller [N]"))
    ((Vector3r,current,Vector3r::Zero(),,"Current value for the controller [N]"))
    ((Real,kP,0.0,,"Proportional gain/coefficient for the PID-controller [-]"))
    ((Real,kI,0.0,,"Integral gain/coefficient for the PID-controller [-]"))
    ((Real,kD,0.0,,"Derivative gain/coefficient for the PID-controller [-]"))
    ((Real,iTerm,0.0,,"Integral term [N]"))
    ((Real,curVel,0.0,,"Current applied velocity [m/s]"))
    ((Real,errorCur,0.0,,"Current error [N]"))
    ((Real,errorPrev,0.0,,"Previous error [N]"))
    ((long,iterPeriod,100.0,,"Periodicity criterion of velocity correlation [-]"))
    ((long,iterPrevStart,-1.0,,"Previous iteration of velocity correlation [-]"))
    /* attrs*/,
    /* ctor */,
    /* py */
  )
  DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ServoPIDController);

struct BicyclePedalEngine: public KinematicEngine{
	virtual void apply(const vector<Body::id_t>& ids);
	void postLoad(BicyclePedalEngine&){ rotationAxis.normalize(); }
	SUDODEM_CLASS_BASE_DOC_ATTRS(BicyclePedalEngine,KinematicEngine,"Engine applying the linear motion of ``bicycle pedal`` e.g. moving points around the axis without rotation",
		((Real,angularVelocity,0,,"Angular velocity. [rad/s]"))
		((Vector3r,rotationAxis,Vector3r::UnitX(),Attr::triggerPostLoad,"Axis of rotation (direction); will be normalized automatically."))
		((Real,radius,-1.0,,"Rotation radius. [m]"))
		((Real,fi,Mathr::PI/2.0,,"Initial phase [radians]"))
	);
};
REGISTER_SERIALIZABLE(BicyclePedalEngine);
