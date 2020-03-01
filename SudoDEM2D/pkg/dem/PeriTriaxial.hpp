// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#pragma once

#include<sudodem/pkg/common/BoundaryController.hpp>
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
		((bool,dynCell,true,,"Imposed stress can be controlled using the packing stiffness or by applying the laws of dynamic (dynCell=true). Don't forget to assign a :yref:`mass<PeriTriaxController.mass>` to the cell."))
		((Real,z_dim,1.0,,"The extention along the z direction, like the length of cyliners in 3D. The value is involed in calculation of the stress tensor, and a recommended value is the mean particle size (FIXME)"))
		((Vector2r,goal,Vector2r::Zero(),,"Desired stress or strain values (depending on stressMask), strains defined as ``strain(i)=log(Fii)``.\n\n.. warning:: Strains are relative to the :yref:`O.cell.refSize<Cell.refSize>` (reference cell size), not the current one (e.g. at the moment when the new strain value is set)."))
		((int,stressMask,((void)"all strains",0),,"mask determining strain/stress (0/1) meaning for goal components"))
		((Vector2r,maxStrainRate,Vector2r(1,1),,"Maximum strain rate of the periodic cell."))
		((Real,maxUnbalanced,1e-4,,"maximum unbalanced force."))
		((Real,absStressTol,1e3,,"Absolute stress tolerance"))
		((Real,relStressTol,3e-5,,"Relative stress tolerance"))
		((Real,growDamping,.25,,"Damping of cell resizing (0=perfect control, 1=no control at all); see also ``wallDamping`` in :yref:`TriaxialStressController`."))
		((int,globUpdate,5,,"How often to recompute average stress, stiffness and unbalaced force."))
		((string,doneHook,,,"python command to be run when the desired state is reached"))
		((Vector2r,maxBodySpan,Vector2r::Zero(),,"maximum body dimension |ycomp|"))
		((Matrix2r,stressTensor,Matrix2r::Zero(),,"average stresses, updated at every step (only every globUpdate steps recomputed from interactions if !dynCell)"))
		((Vector2r,stress,Vector2r::Zero(),,"diagonal terms of the stress tensor"))
		((Vector2r,strain,Vector2r::Zero(),,"cell strain |yupdate|"))
		((Vector2r,strainRate,Vector2r::Zero(),,"cell strain rate |yupdate|"))
		((Vector2r,stiff,Vector2r::Zero(),,"average stiffness (only every globUpdate steps recomputed from interactions) |yupdate|"))
		((Real,currUnbalanced,NaN,,"current unbalanced force (updated every globUpdate) |yupdate|"))
		((Vector2r,prevGrow,Vector2r::Zero(),,"previous cell grow"))
		((Real,mass,NaN,,"mass of the cell (user set); if not set and :yref:`dynCell<PeriTriaxController.dynCell>` is used, it will be computed as sum of masses of all particles."))
		((Real,externalWork,0,,"Work input from boundary controller."))
		((int,velGradWorkIx,-1,(Attr::hidden|Attr::noSave),"Index for work done by velocity gradient, if tracking energy"))
		,,,
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PeriTriaxController);
