// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
// © 2006 Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>

#pragma once
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/IGeom.hpp>
#include<sudodem/core/State.hpp>


/*! Class representing geometry of two bodies in contact.
 *
 * The code under SCG_SHEAR is experimental and is used only if ElasticContactLaw::useShear is explicitly true
 */

#define SCG_SHEAR

class ScGeom: public IGeom {//for 2D
	private:
		//cached values
		//Vector3r twist_axis;//rotation vector around normal
		//Vector3r orthonormal_axis;//rotation vector in contact plane
		//Real rotAngle;//rotation angle
		Vector2r preNormal;
	public:
		Real &radius1, &radius2;
		virtual ~ScGeom();
		inline ScGeom& operator= (const ScGeom& source){
			normal=source.normal; contactPoint=source.contactPoint;
			//twist_axis=source.twist_axis; orthonormal_axis=source.orthonormal_axis;
			radius1=source.radius1; radius2=source.radius2;
			penetrationDepth=source.penetrationDepth; shearInc=source.shearInc;
			return *this;}

		//!precompute values of shear increment and interaction rotation data. Update contact normal to the currentNormal value. Precondition : the value of normal is not updated outside (and before) this function.
		void precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector2r& currentNormal, bool isNew, const Vector2r& shift2, bool avoidGranularRatcheting=true);

		//! Rotates a "shear" vector to keep track of contact orientation. Returns reference of the updated vector.
		Vector2r& rotate(Vector2r& tangentVector) const;
		const Vector2r& shearIncrement() const {return shearInc;}

		// Add method which returns the relative velocity (then, inside the contact law, this can be split into shear and normal component). Handle periodicity.
		Vector2r getIncidentVel(const State* rbp1, const State* rbp2, Real dt, const Vector2r& shiftVel, const Vector2r& shift2, bool avoidGranularRatcheting=true);
		// Implement another version of getIncidentVel which does not handle periodicity.
		Vector2r getIncidentVel(const State* rbp1, const State* rbp2, Real dt, bool avoidGranularRatcheting=true);
		// Add function to get the relative angular velocity (useful to determine bending moment at the contact level)
		Real getRelAngVel(const State* rbp1, const State* rbp2, Real dt);

		// convenience version to be called from python
		Vector2r getIncidentVel_py(shared_ptr<Interaction> i, bool avoidGranularRatcheting);
		Real getRelAngVel_py(shared_ptr<Interaction> i);

		SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ScGeom,IGeom,"Class representing :yref:`geometry<IGeom>` of a contact point between two :yref:`bodies<Body>`. It is more general than disk-disk contact even though it is primarily focused on disks interactions (reason for the 'Sc' namming); it is also used for representing contacts of a :yref:`Disk` with non-spherical bodies (:yref:`Facet`, :yref:`Plane`,  :yref:`Box`, :yref:`ChainedCylinder`), or between two non-spherical bodies (:yref:`ChainedCylinder`). The contact has 3 DOFs (normal and 2×shear) and uses incremental algorithm for updating shear.\n\nWe use symbols $\\vec{x}$, $\\vec{v}$, $\\vec{\\omega}$ respectively for position, linear and angular velocities (all in global coordinates) and $r$ for particles radii; subscripted with 1 or 2 to distinguish 2 disks in contact. Then we define branch length and unit contact normal\n\n.. math::\n\n\tl=||\\vec{x}_2-\\vec{x}_1||, \\vec{n}=\\frac{\\vec{x}_2-\\vec{x}_1}{||\\vec{x}_2-\\vec{x}_1||}\n\nThe relative velocity of the disks is then\n\n.. math::\n\n\t\\vec{v}_{12}=\\frac{r_1+r_2}{l}(\\vec{v}_2-\\vec{v}_1) -(r_2 \\vec{\\omega}_2 + r_1\\vec{\\omega}_1)\\times\\vec{n}\n\nwhere the fraction multplying translational velocities is to make the definition objective and avoid ratcheting effects (see :yref:`Ig2_Disk_Disk_ScGeom.avoidGranularRatcheting`). The shear component is\n\n.. math::\n\n\t\\vec{v}_{12}^s=\\vec{v}_{12}-(\\vec{n}\\cdot\\vec{v}_{12})\\vec{n}.\n\nTangential displacement increment over last step then reads\n\n.. math::\n\n\t\\Delta\\vec{x}_{12}^s=\\Delta t \\vec{v}_{12}^s.",
		((Vector2r,normal,,,"Unit vector oriented along the interaction, from particle #1, towards particle #2. |yupdate|"))
		((Vector2r,contactPoint,,,"some reference point for the interaction (usually in the middle). |ycomp|"))
		((Real,refR1,,,"Reference radius of particle #1. |ycomp|"))
		((Real,refR2,,,"Reference radius of particle #2. |ycomp|"))
		((Real,penetrationDepth,NaN,(Attr::noSave|Attr::readonly),"Penetration distance of disks (positive if overlapping)"))
		((Vector2r,shearInc,Vector2r::Zero(),(Attr::noSave|Attr::readonly),"Shear displacement increment in the last step"))
		,
		/* extra initializers */ ((radius1,refR1)) ((radius2,refR2)),
		/* ctor */ createIndex(); preNormal=Vector2r::Zero();/*twist_axis=orthonormal_axis=Vector2r::Zero();*/,
		/* py */ .def("incidentVel",&ScGeom::getIncidentVel_py,(boost::python::arg("i"),boost::python::arg("avoidGranularRatcheting")=true),"Return incident velocity of the interaction (see also :yref:`Ig2_Disk_Disk_ScGeom.avoidGranularRatcheting` for explanation of the ratcheting argument).")
		.def("relAngVel",&ScGeom::getRelAngVel_py,(boost::python::arg("i")),"Return relative angular velocity of the interaction.")
	);
	REGISTER_CLASS_INDEX(ScGeom,IGeom);
};
REGISTER_SERIALIZABLE(ScGeom);
