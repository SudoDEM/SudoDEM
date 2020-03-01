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
		Vector3r& pos;
		/// rotational motion (reference to inside se3)
		Quaternionr& ori;

		//! mutex for updating the parameters from within the interaction loop (only used rarely)
		boost::mutex updateMutex;

		// bits for blockedDOFs
		enum {DOF_NONE=0,DOF_X=1,DOF_Y=2,DOF_Z=4,DOF_RX=8,DOF_RY=16,DOF_RZ=32};
		//! shorthand for all DOFs blocked
		static const unsigned DOF_ALL=DOF_X|DOF_Y|DOF_Z|DOF_RX|DOF_RY|DOF_RZ;
		//! shorthand for all displacements blocked
		static const unsigned DOF_XYZ=DOF_X|DOF_Y|DOF_Z;
		//! shorthand for all rotations blocked
		static const unsigned DOF_RXRYRZ=DOF_RX|DOF_RY|DOF_RZ;

		//! Return DOF_* constant for given axis∈{0,1,2} and rotationalDOF∈{false(default),true}; e.g. axisDOF(0,true)==DOF_RX
		static unsigned axisDOF(int axis, bool rotationalDOF=false){return 1<<(axis+(rotationalDOF?3:0));}
		//! set DOFs according to two Vector3r arguments (blocked is when disp[i]==1.0 or rot[i]==1.0)
		void setDOFfromVector3r(Vector3r disp,Vector3r rot=Vector3r::Zero());
		//! Getter of blockedDOFs for list of strings (e.g. DOF_X | DOR_RX | DOF_RZ → 'xXZ')
		std::string blockedDOFs_vec_get() const;
		//! Setter of blockedDOFs from string ('xXZ' → DOF_X | DOR_RX | DOF_RZ)
		#ifdef SUDODEM_DEPREC_DOF_LIST
			void blockedDOFs_vec_set(const python::object&);
		#else
			void blockedDOFs_vec_set(const std::string& dofs);
		#endif

		//! Return displacement (current-reference position)
		Vector3r displ() const {return pos-refPos;}
		//! Return rotation (current-reference orientation, as Vector3r)
		Vector3r rot() const { Quaternionr relRot=refOri.conjugate()*ori; AngleAxisr aa(relRot); return aa.axis()*aa.angle(); }

		// python access functions: pos and ori are references to inside Se3r and cannot be pointed to directly
		Vector3r pos_get() const {return pos;}
		void pos_set(const Vector3r p) {pos=p;}
		Quaternionr ori_get() const {return ori; }
		void ori_set(const Quaternionr o){ori=o;}

	SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(State,Serializable,"State of a body (spatial configuration, internal variables).",
		((Se3r,se3,Se3r(Vector3r::Zero(),Quaternionr::Identity()),,"Position and orientation as one object."))
		((Vector3r,vel,Vector3r::Zero(),,"Current linear velocity."))
		((Real,mass,0,,"Mass of this body"))
		((Vector3r,angVel,Vector3r::Zero(),,"Current angular velocity"))
		((Vector3r,angMom,Vector3r::Zero(),,"Current angular momentum"))
		((Vector3r,inertia,Vector3r::Zero(),,"Inertia of associated body, in local coordinate system."))
		((Vector3r,refPos,Vector3r::Zero(),,"Reference position"))
		((Quaternionr,refOri,Quaternionr::Identity(),,"Reference orientation"))
		((unsigned,blockedDOFs,,,"[Will be overridden]"))
		((bool,isDamped,true,,"Damping in :yref:`Newtonintegrator` can be deactivated for individual particles by setting this variable to FALSE. E.g. damping is inappropriate for particles in free flight under gravity but it might still be applicable to other particles in the same simulation."))
		((Real,densityScaling,1,,"|yupdate| see :yref:`GlobalStiffnessTimeStepper::targetDt`.")),
		/* additional initializers */
			((pos,se3.position))
			((ori,se3.orientation)),
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
