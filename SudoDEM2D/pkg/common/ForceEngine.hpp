// 2004 © Janek Kozicki <cosurgi@berlios.de>
// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
// 2014 © Raphael Maurin <raphael.maurin@irstea.fr>

#pragma once

#include<sudodem/core/PartialEngine.hpp>

class ForceEngine: public PartialEngine{
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(ForceEngine,PartialEngine,"Apply contact force on some particles at each step.",
		((Vector2r,force,Vector2r::Zero(),,"Force to apply."))
	);
};
REGISTER_SERIALIZABLE(ForceEngine);

/* Engine for applying force of varying magnitude but constant direction
 * on subscribed bodies. times and magnitudes must have the same length,
 * direction (normalized automatically) gives the orientation.
 *
 * As usual with interpolating engines: the first magnitude is used before the first
 * time point, last magnitude is used after the last time point. Wrap specifies whether
 * time wraps around the last time point to the first time point.
 */
class InterpolatingDirectedForceEngine: public ForceEngine{
	size_t _pos;
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(InterpolatingDirectedForceEngine,ForceEngine,"Engine for applying force of varying magnitude but constant direction on subscribed bodies. times and magnitudes must have the same length, direction (normalized automatically) gives the orientation. \n\n\
	\
	As usual with interpolating engines: the first magnitude is used before the first time point, last magnitude is used after the last time point. Wrap specifies whether time wraps around the last time point to the first time point.",
		((vector<Real>,times,,,"Time readings [s]"))
		((vector<Real>,magnitudes,,,"Force magnitudes readings [N]"))
		((Vector2r,direction,Vector2r(1.0,0),,"Contact force direction (normalized automatically)"))
		((bool,wrap,false,,"wrap to the beginning of the sequence if beyond the last time point")),
		/*ctor*/ _pos=0
	);
};
REGISTER_SERIALIZABLE(InterpolatingDirectedForceEngine);

struct RadialForceEngine: public PartialEngine{
	virtual void action();
	virtual void postLoad(RadialForceEngine&);
	SUDODEM_CLASS_BASE_DOC_ATTRS(RadialForceEngine,PartialEngine,"Apply force of given magnitude directed away from spatial axis.",
		((Vector2r,axisPt,Vector2r::Zero(),,"Point on axis"))
		((Vector2r,axisDir,Vector2r(1,0),Attr::triggerPostLoad,"Axis direction (normalized automatically)"))
		((Real,fNorm,0,,"Applied force magnitude"))
	);
};
REGISTER_SERIALIZABLE(RadialForceEngine);

class DragEngine: public PartialEngine{
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(DragEngine,PartialEngine,"Apply `drag force <http://en.wikipedia.org/wiki/Drag_equation>`__ on some particles at each step, decelerating them proportionally to their linear velocities. The applied force reads\n\n.. math:: F_{d}=-\\frac{\\vec{v}}{|\\vec{v}|}\\frac{1}{2}\\rho|\\vec{v}|^2 C_d A\n\nwhere $\\rho$ is the medium density (:yref:`density<DragEngine.Rho>`), $v$ is particle's velocity,  $A$ is particle projected area (disc), $C_d$ is the drag coefficient (0.47 for :yref:`Disk`), \n\n.. note:: Drag force is only applied to spherical particles, listed in ids.",
		((Real,Rho,1.225,,"Density of the medium (fluid or air), by default - the density of the air."))
		((Real,Cd,0.47,,"Drag coefficient <http://en.wikipedia.org/wiki/Drag_coefficient>`_."))
	);
};
REGISTER_SERIALIZABLE(DragEngine);

class LinearDragEngine: public PartialEngine{
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(LinearDragEngine,PartialEngine,"Apply `viscous resistance or linear drag <http://en.wikipedia.org/wiki/Drag_%28physics%29#Very_low_Reynolds_numbers_.E2.80.94_Stokes.27_drag>`__ on some particles at each step, decelerating them proportionally to their linear velocities. The applied force reads\n\n.. math:: F_{d}=-b{\\vec{v}} \n\nwhere $b$ is the linear drag, $\\vec{v}$ is particle's velocity. \n\n.. math:: b=6\\pi\\nu r \n\nwhere $\\nu$ is the medium viscosity, $r$ is the `Stokes radius <http://en.wikipedia.org/wiki/Stokes_radius>`__ of the particle (but in this case we accept it equal to disk radius for simplification), \n\n.. note:: linear drag is only applied to spherical particles, listed in ids.",
		((Real,nu,0.001,,"Viscosity of the medium."))
	);
};
REGISTER_SERIALIZABLE(LinearDragEngine);


class HydroForceEngine: public PartialEngine{
	private:
		vector<Vector2r> vFluct;
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(HydroForceEngine,PartialEngine,"Apply drag and lift due to a fluid flow vector (1D) to each disk + the buoyant weight.\n The applied drag force reads\n\n.. math:: F_{d}=\\frac{1}{2} C_d A\\rho^f|\\vec{v_f - v}| vec{v_f - v} \n\n where $\\rho$ is the medium density (:yref:`density<HydroForceEngine.rhoFluid>`), $v$ is particle's velocity,  $v_f$ is the velocity of the fluid at the particle center,  $A$ is particle projected area (disc), $C_d$ is the drag coefficient. The formulation of the drag coefficient depends on the local particle reynolds number and the solid volume fraction. The formulation of the drag is [Dallavalle1948]_ [RevilBaudard2013]_ with a correction of Richardson-Zaki [Richardson1954]_ to take into account the hindrance effect. This law is classical in sediment transport. It is possible to activate a fluctuation of the drag force for each particle which account for the turbulent fluctuation of the fluid velocity (:yref:`velFluct`). The model implemented for the turbulent velocity fluctuation is a simple discrete random walk which takes as input the reynolds stress tensor Re_{xz} in function of the depth and allows to recover the main property of the fluctuations by imposing <u_x'u_z'> (z) = <Re>(z)/rho^f. It requires as input <Re>(z)/rho^f called :yref:`simplifiedReynoldStresses` in the code. \n The formulation of the lift is taken from [Wiberg1985]_ and is such that : \n\n.. math:: F_{L}=\\frac{1}{2} C_L A\\rho^f((v_f - v)^2{top} - (v_f - v)^2{bottom}) \n\n Where the subscript top and bottom means evaluated at the top (respectively the bottom) of the disk considered. This formulation of the lift account for the difference of pressure at the top and the bottom of the particle inside a turbulent shear flow. As this formulation is controversial when approaching the threshold of motion [Schmeeckle2007]_ it is possible to desactivate it with the variable :yref:`lift`.\n The buoyancy is taken into account through the buoyant weight : \n\n.. math:: F_{buoyancy}= - rho^f V^p g \n\n, where g is the gravity vector along the vertical, and V^p is the volume of the particle.",
		((Real,rhoFluid,1000,,"Density of the fluid, by default - density of water"))
		((Real,viscoDyn,1e-3,,"Dynamic viscosity of the fluid, by default - viscosity of water"))
		((Real,zRef,,,"Position of the reference point which correspond to the first value of the fluid velocity"))
		((Real,nCell,,,"Size of the vector of the fluid velocity"))
		((Real,deltaZ,,,"width of the discretization cell "))
		((Real,expoRZ,3.1,,"Value of the Richardson-Zaki exponent, for the correction due to hindrance"))
                ((bool,lift,true,,"Option to activate or not the evaluation of the lift"))
		((Real,Cl,0.2,,"Value of the lift coefficient taken from [Wiberg1985]_"))
                ((Vector2r,gravity,,,"Gravity vector (may depend on the slope)."))
		((vector<Real>,vxFluid,,,"Discretized streamwise fluid velocity profile in function of the depth"))
		((vector<Real>,phiPart,,,"Discretized solid volume fraction profile in function of the depth"))
		((bool,velFluct,false,,"If true, activate the determination of turbulent fluid velocity fluctuation at the position of each particle, using a simple discrete random walk model based on the Reynolds stresses :yref:`turbStress<HydroForceEngine.squaredAverageTurbfluct>`"))
		((vector<Real>,simplifiedReynoldStresses,,,"Vector of size equal to :yref:`turbStress<HydroForceEngine.nCell>` containing the Reynolds stresses divided by the fluid density in function of the depth. simplifiedReynoldStresses(z) =  <u_x'u_z'>(z)^2 "))
		((Real,bedElevation,,,"Elevation of the bed above which the fluid flow is turbulent and the particles undergo turbulent velocity fluctuation."))
	);
};
REGISTER_SERIALIZABLE(HydroForceEngine);
