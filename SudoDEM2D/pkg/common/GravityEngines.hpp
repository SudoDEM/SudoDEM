// 2004 © Janek Kozicki <cosurgi@berlios.de>
// 2007,2008 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/pkg/common/FieldApplier.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Body.hpp>

/*! Homogeneous gravity field; applies gravity×mass force on all bodies. */
class GravityEngine: public FieldApplier{
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(GravityEngine,FieldApplier,"Engine applying constant acceleration to all bodies. DEPRECATED, use :yref:`Newton::gravity` unless you need energy tracking or selective gravity application using groupMask).",
		((Vector2r,gravity,Vector2r::Zero(),,"Acceleration [kgms⁻²]"))
		((int,gravPotIx,-1,(Attr::noSave|Attr::hidden),"Index for gravPot energy"))
		((int,mask,0,,"If mask defined, only bodies with corresponding groupMask will be affected by this engine. If 0, all bodies will be affected."))
		((bool,warnOnce,true,,"For deprecation warning once."))
		,/*ctor*/,/*py*/
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(GravityEngine);


/*! Engine attracting all bodies towards a central body (doesn't depend on distance);
 *
 * @todo This code has not been yet tested at all.
 */
class CentralGravityEngine: public FieldApplier {
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(CentralGravityEngine,FieldApplier,"Engine applying acceleration to all bodies, towards a central body.",
		((Body::id_t,centralBody,Body::ID_NONE,,"The :yref:`body<Body>` towards which all other bodies are attracted."))
		((Real,accel,0,,"Acceleration magnitude [kgms⁻²]"))
		((bool,reciprocal,false,,"If true, acceleration will be applied on the central body as well."))
		((int,mask,0,,"If mask defined, only bodies with corresponding groupMask will be affected by this engine. If 0, all bodies will be affected."))
		,,
	);
};
REGISTER_SERIALIZABLE(CentralGravityEngine);

/*! Apply acceleration (independent of distance) directed towards an axis.
 *
 */
class AxialGravityEngine: public FieldApplier {
	public:
	virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(AxialGravityEngine,FieldApplier,"Apply acceleration (independent of distance) directed towards an axis.",
		((Vector2r,axisPoint,Vector2r::Zero(),,"Point through which the axis is passing."))
		((Vector2r,axisDirection,Vector2r::UnitX(),,"direction of the gravity axis (will be normalized automatically)"))
		((Real,acceleration,0,,"Acceleration magnitude [kgms⁻²]"))
		((int,mask,0,,"If mask defined, only bodies with corresponding groupMask will be affected by this engine. If 0, all bodies will be affected."))
	);
};
REGISTER_SERIALIZABLE(AxialGravityEngine);
