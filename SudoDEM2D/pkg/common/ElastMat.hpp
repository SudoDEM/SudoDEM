// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/core/Material.hpp>

/*! Elastic material */
class ElastMat: public Material{
	public:
	virtual ~ElastMat() {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(ElastMat,Material,"Purely elastic material. The material parameters may have different meanings depending on the :yref:`IPhysFunctor` used : true Young and Poisson in :yref:`Ip2_FrictMat_FrictMat_MindlinPhys`, or contact stiffnesses in :yref:`Ip2_FrictMat_FrictMat_FrictPhys`.",
		//((Real,young,1e9,,"elastic modulus [Pa]. It has different meanings depending on the Ip functor."))
		//((Real,poisson,.25,,"Poisson's ratio or the ratio between shear and normal stiffness [-]. It has different meanings depending on the Ip functor.  "))
		((Real,Kn,1e5,,"Contact normal stiffness [N/m]."))
		((Real,Ks,0.7e5,,"Contact tangential stiffness [N/m].")),
		/*ctor*/ createIndex();
	);
	REGISTER_CLASS_INDEX(ElastMat,Material);
};
REGISTER_SERIALIZABLE(ElastMat);

/*! Granular material */
class FrictMat: public ElastMat{
	public:
	virtual ~FrictMat() {};
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(FrictMat,ElastMat,"Elastic material with contact friction. See also :yref:`ElastMat`.",
		((Real,frictionAngle,.5,,"Contact friction angle (in radians). Hint : use 'radians(degreesValue)' in python scripts.")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(FrictMat,ElastMat);
};
REGISTER_SERIALIZABLE(FrictMat);
