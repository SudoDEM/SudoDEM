// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#pragma once

#include<sudodem/pkg/common/BoundaryController.hpp>

class PeriIsoCompressor: public BoundaryController{
	Real avgStiffness; Real maxDisplPerStep; Vector3r sumForces, sigma;
	Real currUnbalanced;
	public:
		void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(PeriIsoCompressor,BoundaryController,"Compress/decompress cloud of spheres by controlling periodic cell size until it reaches prescribed average stress, then moving to next stress value in given stress series.",
		((vector<Real>,stresses,,,"Stresses that should be reached, one after another"))
		((Real,charLen,-1.,,"Characteristic length, should be something like mean particle diameter (default -1=invalid value))"))
		((Real,maxSpan,-1.,,"Maximum body span in terms of bbox, to prevent periodic cell getting too small. |ycomp|"))
		((Real,maxUnbalanced,1e-4,,"if actual unbalanced force is smaller than this number, the packing is considered stable,"))
		((int,globalUpdateInt,20,,"how often to recompute average stress, stiffness and unbalanced force"))
		((size_t,state,0,,"Where are we at in the stress series"))
		((string,doneHook,"",,"Python command to be run when reaching the last specified stress"))
		((bool,keepProportions,true,,"Exactly keep proportions of the cell (stress is controlled based on average, not its components")),
		/*ctor*/
			currUnbalanced=-1;
			avgStiffness=-1;
			maxDisplPerStep=-1;
			sumForces=Vector3r::Zero();
			sigma=Vector3r::Zero();
		,
		/*py*/
			.def_readonly("currUnbalanced",&PeriIsoCompressor::currUnbalanced,"Current value of unbalanced force")
			.def_readonly("sigma",&PeriIsoCompressor::sigma,"Current stress value")
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PeriIsoCompressor);
/* Engine for independently controlling stress or strain in periodic simulations.

strainStress contains absolute values for the controlled quantity, and stressMask determines
meaning of those values (0 for strain, 1 for stress): e.g. ( 1<<0 | 1<<2 ) = 1 | 4 = 5 means that
strainStress[0] and strainStress[2] are stress values, and strainStress[1] is strain.

See scripts/test/periodic-triax.py for a simple example.

*/

class PeriTriaxController: public BoundaryController{
	public:
		virtual void action();
		void strainStressStiffUpdate();
	SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(PeriTriaxController,BoundaryController,"Engine for independently controlling stress or strain in periodic simulations.\n\n``strainStress`` contains absolute values for the controlled quantity, and ``stressMask`` determines meaning of those values (0 for strain, 1 for stress): e.g. ``( 1<<0 | 1<<2 ) = 1 | 4 = 5`` means that ``strainStress[0]`` and ``strainStress[2]`` are stress values, and ``strainStress[1]`` is strain. \n\nSee scripts/test/periodic-triax.py for a simple example.",
		((bool,dynCell,false,,"Imposed stress can be controlled using the packing stiffness or by applying the laws of dynamic (dynCell=true). Don't forget to assign a :yref:`mass<PeriTriaxController.mass>` to the cell."))
		((Vector3r,goal,Vector3r::Zero(),,"Desired stress or strain values (depending on stressMask), strains defined as ``strain(i)=log(Fii)``.\n\n.. warning:: Strains are relative to the :yref:`O.cell.refSize<Cell.refSize>` (reference cell size), not the current one (e.g. at the moment when the new strain value is set)."))
		((int,stressMask,((void)"all strains",0),,"mask determining strain/stress (0/1) meaning for goal components"))
		((Vector3r,maxStrainRate,Vector3r(1,1,1),,"Maximum strain rate of the periodic cell."))
		((Real,maxUnbalanced,1e-4,,"maximum unbalanced force."))
		((Real,absStressTol,1e3,,"Absolute stress tolerance"))
		((Real,relStressTol,3e-5,,"Relative stress tolerance"))
		((Real,growDamping,.25,,"Damping of cell resizing (0=perfect control, 1=no control at all); see also ``wallDamping`` in :yref:`TriaxialStressController`."))
		((int,globUpdate,5,,"How often to recompute average stress, stiffness and unbalaced force."))
		((string,doneHook,,,"python command to be run when the desired state is reached"))
		((Vector3r,maxBodySpan,Vector3r::Zero(),,"maximum body dimension |ycomp|"))
		((Matrix3r,stressTensor,Matrix3r::Zero(),,"average stresses, updated at every step (only every globUpdate steps recomputed from interactions if !dynCell)"))
		((Vector3r,stress,Vector3r::Zero(),,"diagonal terms of the stress tensor"))
		((Vector3r,strain,Vector3r::Zero(),,"cell strain |yupdate|"))
		((Vector3r,strainRate,Vector3r::Zero(),,"cell strain rate |yupdate|"))
		((Vector3r,stiff,Vector3r::Zero(),,"average stiffness (only every globUpdate steps recomputed from interactions) |yupdate|"))
		((Real,currUnbalanced,NaN,,"current unbalanced force (updated every globUpdate) |yupdate|"))
		((Vector3r,prevGrow,Vector3r::Zero(),,"previous cell grow"))
		((Real,mass,NaN,,"mass of the cell (user set); if not set and :yref:`dynCell<PeriTriaxController.dynCell>` is used, it will be computed as sum of masses of all particles."))
		((Real,externalWork,0,,"Work input from boundary controller."))
		((int,velGradWorkIx,-1,(Attr::hidden|Attr::noSave),"Index for work done by velocity gradient, if tracking energy"))
		,,,
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PeriTriaxController);

class Peri3dController: public BoundaryController{
	public:
		Vector6r stressOld;
		Matrix3r sigma, epsilon, epsilonRate, rot, nonrot;

		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Peri3dController,BoundaryController,"Experimental controller of full strain/stress tensors on periodic cell. Detailed documentation is in py/_extraDocs.py.",
		((Vector6r,stress,Vector6r::Zero(),,"Current stress vector ($\\sigma_x$,$\\sigma_y$,$\\sigma_z$,$\\tau_{yz}$,$\\tau_{zx}$,$\\tau_{xy}$)|yupdate|."))
		((Vector6r,strain,Vector6r::Zero(),,"Current strain (deformation) vector ($\\varepsilon_x$,$\\varepsilon_y$,$\\varepsilon_z$,$\\gamma_{yz}$,$\\gamma_{zx}$,$\\gamma_{xy}$) |yupdate|."))
		((Vector6r,strainRate,Vector6r::Zero(),,"Current strain rate vector."))
		((Vector6r,stressRate,Vector6r::Zero(),,"Current stress rate vector (that is prescribed, the actual one slightly differ)."))
		((Vector6r,stressIdeal,Vector6r::Zero(),,"Ideal stress vector at current time step."))
		((Vector6r,goal,Vector6r::Zero(),,"Goal state; only the upper triangular matrix is considered; each component is either prescribed stress or strain, depending on :yref:`stressMask<Peri3dController.stressMask>`."))
		((int,stressMask,((void)"all strains",0),,"mask determining whether components of :yref:`goal<Peri3dController.goal>` are strain (0) or stress (1). The order is 00,11,22,12,02,01 from the least significant bit. (e.g. 0b000011 is stress 00 and stress 11)."))
		((int,nSteps,1000,,"Number of steps of the simulation."))
		((Real,progress,0.,,"Actual progress of the simulation with Controller."))
		((Real,mod,.1,,"Predictor modificator, by trail-and-error analysis the value 0.1 was found as the best."))
		((string,doneHook,,,"Python command (as string) to run when :yref:`nSteps<Peri3dController.nSteps>` is achieved. If empty, the engine will be set :yref:`dead<Engine.dead>`."))
		((vector<Vector2r>,xxPath,vector<Vector2r>(1,Vector2r::Ones()),,"\"Time function\" (piecewise linear) for xx direction. Sequence of couples of numbers. First number is time, second number desired value of respective quantity (stress or strain). The last couple is considered as final state (equal to (:yref:`nSteps<Peri3dController.nSteps>`, :yref:`goal<Peri3dController.goal>`)), other values are relative to this state.\n\nExample: nSteps=1000, goal[0]=300, xxPath=((2,3),(4,1),(5,2))\n\nat step 400 (=5*1000/2) the value is 450 (=3*300/2),\n\nat step 800 (=4*1000/5) the value is 150 (=1*300/2),\n\nat step 1000 (=5*1000/5=nSteps) the value is 300 (=2*300/2=goal[0]).\n\nSee example :ysrc:`scripts/test/peri3dController_example1` for illusration."))
		((vector<Vector2r>,yyPath,vector<Vector2r>(1,Vector2r::Ones()),,"Time function for yy direction, see :yref:`xxPath<Peri3dController.xxPath>`"))
		((vector<Vector2r>,zzPath,vector<Vector2r>(1,Vector2r::Ones()),,"Time function for zz direction, see :yref:`xxPath<Peri3dController.xxPath>`"))
		((vector<Vector2r>,yzPath,vector<Vector2r>(1,Vector2r::Ones()),,"Time function for yz direction, see :yref:`xxPath<Peri3dController.xxPath>`"))
		((vector<Vector2r>,zxPath,vector<Vector2r>(1,Vector2r::Ones()),,"Time function for zx direction, see :yref:`xxPath<Peri3dController.xxPath>`"))
		((vector<Vector2r>,xyPath,vector<Vector2r>(1,Vector2r::Ones()),,"Time function for xy direction, see :yref:`xxPath<Peri3dController.xxPath>`"))
		((Real,maxStrainRate,1e3,,"Maximal absolute value of strain rate (both normal and shear components of :yref:`strain<Peri3dController.strain>`)"))
		((Real,maxStrain,1e6,,"Maximal asolute value of :yref:`strain<Peri3dController.strain>` allowed in the simulation. If reached, the simulation is considered as finished"))
		((Real,youngEstimation,1e20,,"Estimation of macroscopic Young's modulus, used for the first simulation step"))
		((Real,poissonEstimation,.25,,"Estimation of macroscopic Poisson's ratio, used used for the first simulation step"))
		//
		((Vector6r,stressGoal,Vector6r::Zero(),Attr::readonly,"Peri3dController internal variable"))
		((Vector6r,strainGoal,Vector6r::Zero(),Attr::readonly,"Peri3dController internal variable"))
		((Vector6i,pe,Vector6i::Zero(),Attr::readonly,"Peri3dController internal variable"))
		((Vector6i,ps,Vector6i::Zero(),Attr::readonly,"Peri3dController internal variable"))
		((Vector6i,pathSizes,Vector6i::Zero(),Attr::readonly,"Peri3dController internal variable"))
		((Vector6i,pathsCounter,Vector6i::Zero(),Attr::readonly,"Peri3dController internal variable"))
		((int,lenPe,NaN,Attr::readonly,"Peri3dController internal variable"))
		((int,lenPs,NaN,Attr::readonly,"Peri3dController internal variable"))
		,
		/*ctor*/
		,
		/*py*/
	);
};
REGISTER_SERIALIZABLE(Peri3dController);
