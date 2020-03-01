// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/core/IGeom.hpp>

/*! Abstract class that unites ScGeom and L3Geom,
	created for the purposes of GlobalStiffnessTimeStepper.
	It might be removed in the future. */
class GenericSpheresContact: public IGeom{
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(GenericSpheresContact,IGeom,
		"Class uniting :yref:`ScGeom` and :yref:`L3Geom`, for the purposes of :yref:`GlobalStiffnessTimeStepper`. (It might be removed inthe future). Do not use this class directly.",
		((Vector3r,normal,,,"Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))
		((Vector3r,contactPoint,,,"some reference point for the interaction (usually in the middle). |ycomp|"))
		((Real,refR1,,,"Reference radius of particle #1. |ycomp|"))
		((Real,refR2,,,"Reference radius of particle #2. |ycomp|")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(GenericSpheresContact,IGeom);

	virtual ~GenericSpheresContact() {};
};
REGISTER_SERIALIZABLE(GenericSpheresContact);
