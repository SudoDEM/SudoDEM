// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/lib/multimethods/Indexable.hpp>
#include<sudodem/core/Dispatcher.hpp>

// delete later and remove relevant code, to not support old State.blockedDOFs=['x','y','rz'] syntax anymore
//#define SUDODEM_DEPREC_DOF_LIST

class State: public Serializable, public Indexable{
	public:
		/// linear motion (references to inside se3)
		Vector2r& pos;
		/// rotational motion (reference to inside se3)
		Rotationr& ori;

		//! mutex for updating the parameters from within the interaction loop (only used rarely)
		boost::mutex updateMutex;

		// bits for blockedDOFs
		enum {DOF_NONE=0,DOF_X=1,DOF_Y=2,DOF_RZ=4};
		//! shorthand for all DOFs blocked
		static const unsigned DOF_ALL=DOF_X|DOF_Y|DOF_RZ;
		//! shorthand for all displacements blocked
		static const unsigned DOF_XY=DOF_X|DOF_Y;
		//! shorthand for all rotations blocked
		//static const unsigned DOF_RZ=DOF_RX|DOF_RY|DOF_RZ;

		//! Return DOF_* constant for given axis∈{0,1} and rotationalDOF∈{false(default),true}; e.g. axisDOF(0,true)==DOF_RX
		static unsigned axisDOF(int axis, bool rotationalDOF=false){return 1<<(axis+(rotationalDOF?3:0));}
		//! set DOFs according to two Vector3r arguments (blocked is when disp[i]==1.0 or rot[i]==1.0)
		void setDOFfromVector3r(Vector3r disp_rot);//disp_rot(X,Y,RZ)
		//! Getter of blockedDOFs for list of strings (e.g. DOF_X | DOR_RX | DOF_RZ → 'xXZ')
		std::string blockedDOFs_vec_get() const;
		//! Setter of blockedDOFs from string ('xXZ' → DOF_X | DOR_RX | DOF_RZ)
		#ifdef SUDODEM_DEPREC_DOF_LIST
			void blockedDOFs_vec_set(const python::object&);
		#else
			void blockedDOFs_vec_set(const std::string& dofs);
		#endif

		//! Return displacement (current-reference position)
		Vector2r displ() const {return pos-refPos;}
		//! Return rotation (current-reference orientation, as Vector3r)//debug for 2D
		Real rot() const { return ori.angle(); }

		// python access functions: pos and ori are references to inside Se3r and cannot be pointed to directly
		Vector2r pos_get() const {return pos;}
		void pos_set(const Vector2r p) {pos=p;}
		Rotationr ori_get() const {return ori; }
		void ori_set(const Rotationr o){ori=o;}

	SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(State,Serializable,"State of a body (spatial configuration, internal variables).",
		((Se2r,se2,Se2r(Vector2r::Zero(),Rotationr::Identity()),,"Position and orientation as one object."))
		((Vector2r,vel,Vector2r::Zero(),,"Current linear velocity."))
		((Real,mass,0,,"Mass of this body"))
		((Real,angVel,0.0,,"Current angular velocity"))
		((Real,angMom,0.0,,"Current angular momentum"))
		((Real,inertia,0.0,,"Inertia of associated body, in local coordinate system."))
		((Vector2r,refPos,Vector2r::Zero(),,"Reference position"))
		((Rotationr,refOri,Rotationr::Identity(),,"Reference orientation"))
		((unsigned,blockedDOFs,,,"[Will be overridden]"))
		((bool,isDamped,true,,"Damping in :yref:`Newtonintegrator` can be deactivated for individual particles by setting this variable to FALSE. E.g. damping is inappropriate for particles in free flight under gravity but it might still be applicable to other particles in the same simulation."))
		((Real,densityScaling,1,,"|yupdate| see :yref:`GlobalStiffnessTimeStepper::targetDt`.")),
		/* additional initializers */
			((pos,se2.position))
			((ori,se2.rotation)),
		/* ctor */,
		/*py*/
		SUDODEM_PY_TOPINDEXABLE(State)
		.add_property("blockedDOFs",&State::blockedDOFs_vec_get,&State::blockedDOFs_vec_set,"Degress of freedom where linear/angular velocity will be always constant (equal to zero, or to an user-defined value), regardless of applied force/torque. String that may contain 'xyzXYZ' (translations and rotations).")
		// references must be set using wrapper funcs
		.add_property("pos",&State::pos_get,&State::pos_set,"Current position.")
		.add_property("ori",&State::ori_get,&State::ori_set,"Current orientation.")
		.def("displ",&State::displ,"Displacement from :yref:`reference position<State.refPos>` (:yref:`pos<State.pos>` - :yref:`refPos<State.refPos>`)")
		.def("rot",&State::rot,"Rotation from :yref:`reference orientation<State.refOri>` (as rotation vector)")
	);
	REGISTER_INDEX_COUNTER(State);
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(State);
