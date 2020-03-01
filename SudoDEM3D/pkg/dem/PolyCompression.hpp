/*************************************************************************
*  Copyright (C) 2016 by Zhswee                                          *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#ifdef SUDODEM_CGAL
#include"VolumeFric.hpp"
#include<sudodem/pkg/common/BoundaryController.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/dem/Polyhedra.hpp>
#include<boost/array.hpp>
#include<sudodem/core/PartialEngine.hpp>
#include<sudodem/pkg/common/Facet.hpp>
#include<fstream>
#include<string>

class State;

class Scene;
class State;


/*! \brief Controls the stress on the boundaries of a box and compute strain-like and stress-like quantities for the packing. The algorithms used have been developed initialy for simulations reported in [Chareyre2002a] and [Chareyre2005]. They have been ported to SudoDEM in a second step and used in e.g. [Kozicki2008],[Scholtes2009b],[Jerier2010b].
*/

class PolyCompressionEngine : public BoundaryController
{
	private :
		//! is this the beginning of the simulation, after reading the scene? -> it is the first time that SudoDEM passes trought the engine ThreeDTriaxialEngine
		bool firstRun;
		std::ofstream out;
		bool majorok;
		shared_ptr<Body> box [6];//pointers of facets of up box

		//! internal index values for retrieving wall,top wall and up box facets
		enum {left_wall=0,right_wall,front_wall,back_wall,bottom_wall,top_wall};

		inline const Vector3r getForce(Scene* rb, Body::id_t id){ return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */ }
	public :
		Vector3r translationAxisy;
		Vector3r translationAxisx;
		Vector3r translationAxisz;//store the rotation degrees.Not changed the variable name because of the initial test sample existing.

		//! The value of stiffness (updated according to stiffnessUpdateInterval)
		Real	stiffness;
		Real 	pos[4];
		Vector3r	strain;

		Real shearStress[2];//left and right wall
		Real left_facet_pos;
		//! Value of spheres volume (solid volume)
		Real particlesVolume;
		Real previouswallpos;
		Real previousStress [2];//shear stress at the previous step.
		//! Value of box volume
		Real boxVolume;
		//! Sample porosity
		Real porosity;
		Real fmax;//the maximum of shear f when shearing

		virtual ~PolyCompressionEngine();

		virtual void action();
		//! update the stiffness of boundary-packing interaction (sum of contacts stiffness on the boundary)
		void updateStiffness();
		//! Compute stresses on walls as "Vector3r stress[6]", compute meanStress, strain[3] and mean strain
		//void computeStressStrain();
		//! Compute the mean/max unbalanced force in the assembly (normalized by mean contact force)
    	Real ComputeUnbalancedForce(bool maxUnbalanced=false);
		Real getResultantF();//resultant force acting on the walls and particles in the down box.
		//recording
		void recordData();
		//check force on walls
		bool checkForce(int wall, Vector3r resultantForce);
		void checkMajorF(bool init);//calculate stress of bottom and top walls
		Real GetOverlpV();//calculate the overlapping volume between particles if needed.
		///Change physical propertieaabbs of interactions and/or bodies in the middle of a simulation (change only friction for the moment, complete this function to set cohesion and others before compression test)
		void setContactProperties(Real frictionDegree);

		void loadingStress();
		void getBox();
		//void getStressStrain();//
		SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(PolyCompressionEngine,BoundaryController,
		"An engine maintaining constant stresses or constant strain rates on some boundaries of a parallepipedic packing. The stress/strain control is defined for each axis using :yref:`PolyCompressionEngine::stressMask` (a bitMask) and target values are defined by goal1,goal2, and goal3. sigmaIso has to be defined during growing phases."
		"\n\n.. note::\n\t The algorithms used have been developed initialy for simulations reported in [Chareyre2002a]_ and [Chareyre2005]_. They have been ported to SudoDEM in a second step and used in e.g. [Kozicki2008]_,[Scholtes2009b]_,[Jerier2010b]."
		"The engine perform a triaxial compression with a control in direction 'i' in stress (if stressControl_i) else in strain.\n\n"
		"For a stress control the imposed stress is specified by 'sigma_i' with a 'max_veli' depending on 'strainRatei'. To obtain the same strain rate in stress control than in strain control you need to set 'wallDamping = 0.8'.\n"
		"For a strain control the imposed strain is specified by 'strainRatei'.\n"
		"With this engine you can also perform internal compaction by growing the size of particles by using ``TriaxialStressController::controlInternalStress``. For that, just switch on 'internalCompaction=1' and fix sigma_iso=value of mean pressure that you want at the end of the internal compaction.\n"
		,
   		((unsigned int,stiffnessUpdateInterval,10,,"target strain rate (./s)"))
   		((unsigned int,echo_interval,10,,"interval of iteration for printing run info."))
		((std::string,file,"recordData",,"Name of file to save to; must not be empty."))
		((Real,height,0,Attr::readonly,"size of the box z(2-axis) |yupdate|"))
		((Real,width,0,Attr::readonly,"size of the box x(0-axis) |yupdate|"))
		((Real,depth,0,Attr::readonly,"size of the box y(1-axis) |yupdate|"))
		((Real,LP_area,0.0,Attr::readonly,"area of the loading plates"))
		((Real,f_threshold,0.05,,"threshold for approaching the target fast."))
		((Real,height0,0,,"Reference size for strain definition. See :yref:`PolyCompressionEngine::height`"))
		((Real,width0,0,,"Reference size for strain definition. See :yref:`PolyCompressionEngine::width`"))
		((Real,depth0,0,,"Reference size for strain definition. See :yref:`PolyCompressionEngine::depth`"))
		((Real,wall_equal_tol,0.001,,"wall equilibrium tolerance after each increament."))
		((Real,majorF_tol,0.985,,"stress gradual reduction torlerance of top and bottom walls when compressing"))
		((Real,unbf_tol,0.5,,"tolerance of UnbalancedForce"))
		((Real,wall_max_vel,1.,,"max velocity of the top  facet when loading"))
		((int,step_interval,20,,"step interval during each displacement increament"))
		((Real,h0,0.003,,"boundary h0"))
		((Real,w0,0.003,,"boundary w0"))
		((Real,d0,0.003,,"boundary d0"))
		((bool,Sphere_on,false,,"true for spherical particles, false for polyhedral particles"))
		((bool,controlFlag,false,,"this option is used to control the stress"))
		((Real,goal,0,,"prescribed stress/strain rate on axis 1, as defined by :yref:`PolyCompressionEngine::stressMask` (see also :yref:`PolyCompressionEngine::isAxisymetric`)"))
		((unsigned int,total_incr,24000,,"the total increaments"))
		((int,outinfo_interval,5,,"the iteration interval at which displaying the info for monitoring"))
		((unsigned int,savedata_interval,12,,"the interval steps between two savedata operations"))
		((Real,max_vel,1,,"Maximum allowed walls velocity [m/s]. This value superseeds the one assigned by the stress controller if the later is higher. max_vel can be set to infinity in many cases, but sometimes helps stabilizing packings. Based on this value, different maxima are computed for each axis based on the dimensions of the sample, so that if each boundary moves at its maximum velocity, the strain rate will be isotropic (see e.g. :yref:`PolyCompressionEngine::max_vel1`)."))
		((Real,externalWork,0,Attr::readonly,"Energy provided by boundaries.|yupdate|"))
		((Real,thickness,0.,,"thickness of the board"))
		((Real, shearRate,0,,"target shear rate in direction x (./s)"))
		((Real, currentStrainRate2,0,,"current strain rate in direction 2 - converging to :yref:`ThreeDTriaxialEngine::strainRate2` (./s)"))
		((Real, UnbalancedForce,1,,"mean resultant forces divided by mean contact force"))
		((Real, frictionAngleDegree,-1,,"Value of friction used in the simulation if (updateFrictionAngle)"))
		((bool, updateFrictionAngle,false,,"Switch to activate the update of the intergranular frictionto the value :yref:`ThreeDTriaxialEngine::frictionAngleDegree`."))
		((bool, Confining,true,,"Switch to confine or shear"))
		//((Real, strainDamping,0.9997,,"factor used for smoothing changes in effective strain rate. If target rate is TR, then (1-damping)*(TR-currentRate) will be added at each iteration. With damping=0, rate=target all the time. With damping=1, it doesn't change."))
		((std::string,Key,"",,"A string appended at the end of all files, use it to name simulations."))
		,
		/* extra initializers */
		,
		//
		,
   		/* constructor */
   		firstRun = true;
		majorok = true;
		previousStress[0] = 0.;previousStress[1] = 0.;
		shearStress[0] = 0.;shearStress[1] = 0.;
		for(int j = 0; j < 4; j++) pos[j] = 0;
		previouswallpos = 0.;
		stiffness = 0.;
		left_facet_pos = 0.;
		porosity=1;
		fmax = 0;
		//
		translationAxisy=Vector3r(0,1,0);
		translationAxisx=Vector3r(1,0,0);
		translationAxisz=Vector3r(0,0,0);//CAUTION:used to store the rotation degrees along x,y and direction.

		boxVolume=0;
		//
		,
		.def_readonly("strain",&PolyCompressionEngine::strain,"Current strain in a vector (exx,eyy,ezz). The values reflect true (logarithmic) strain.")
		.def_readonly("porosity",&PolyCompressionEngine::porosity,"Porosity of the packing.")
		.def_readonly("boxVolume",&PolyCompressionEngine::boxVolume,"Total packing volume.")
		//.def("setContactProperties",&PolyCompressionEngine::setContactProperties,"Assign a new friction angle (degrees) to dynamic bodies and relative interactions")
		)
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PolyCompressionEngine);
#endif
