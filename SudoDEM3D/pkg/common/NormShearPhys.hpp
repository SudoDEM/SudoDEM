// © 2007 Janek Kozicki <cosurgi@mail.berlios.de>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include<sudodem/core/IPhys.hpp>

class NormPhys:public IPhys {
	public:
		virtual ~NormPhys() {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(NormPhys,IPhys,"Abstract class for interactions that have normal stiffness.",
		((Real,kn,0,,"Normal stiffness"))
		((Vector3r,normalForce,Vector3r::Zero(),,"Normal force after previous step (in global coordinates).")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(NormPhys,IPhys);
};
REGISTER_SERIALIZABLE(NormPhys);

class NormShearPhys: public NormPhys{
	public:
		virtual ~NormShearPhys() {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(NormShearPhys,NormPhys,
		"Abstract class for interactions that have shear stiffnesses, in addition to normal stiffness. This class is used in the PFC3d-style stiffness timestepper.",
		((Real,ks,0,,"Shear stiffness"))
		((Vector3r,shearForce,Vector3r::Zero(),,"Shear force after previous step (in global coordinates).")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(NormShearPhys,NormPhys);
};
REGISTER_SERIALIZABLE(NormShearPhys);
