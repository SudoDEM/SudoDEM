/*************************************************************************
*  Copyright (C) 2017 by Sway Zhao                                       *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<sudodem/pkg/common/BoundaryController.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/pkg/dem/Superquadrics.hpp>
#include<boost/array.hpp>
#include<sudodem/core/PartialEngine.hpp>
#include<sudodem/pkg/common/Facet.hpp>

#include"Node.hpp"
#include"TriElement.hpp"

#include<fstream>
#include<string>
#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif

#define CUBIC_SAMPLE//
#define SOLE_GAIN //move two walls with their corresponding gain in the same direction
//testing contact model
#define FM//FricMat

class State;

class Scene;
class State;


/*! \brief Controls the stress on the boundaries of a box and compute strain-like and stress-like quantities for the packing. The algorithms used have been developed initialy for simulations reported in [Chareyre2002a] and [Chareyre2005]. They have been ported to SudoDEM in a second step and used in e.g. [Kozicki2008],[Scholtes2009b],[Jerier2010b].
*/

class FlexCompressionEngine : public BoundaryController
{
	private :
		//! is this the beginning of the simulation, after reading the scene? -> it is the first time that SudoDEM passes trought the engine ThreeDTriaxialEngine
		bool firstRun;
		std::ofstream out;
		bool majorok;
        int numWallFacet;//number of all walls and TriElements
        shared_ptr<Body> top_wall,bottom_wall;
        vector<shared_ptr<Node>> top_nodes;//nodes contacting the top wall
        vector<shared_ptr<Node>> bottom_nodes;//nodes contacting the bottom wall
        vector<shared_ptr<Node>> free_nodes;//other nodes without contacting the walls
        vector<shared_ptr<TriElement>> Elements;//store elements of the flexible wall
		inline const Vector3r getForce(Scene* rb, Body::id_t id){ return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */ }
	public :
		//! The value of stiffness (updated according to stiffnessUpdateInterval)
		Real wall_stiffness[6];//left_wall=0,right_wall,front_wall,back_wall,bottom_wall,top_wall

		Real 	pos[4];

		Real 	bottom_wpos[4];
		Real 	top_wpos[4];
		Real x_area; //area in the x axis
        Real loading_area;//for cylinder wall
		Vector2r wall_pos[6];//left_wall=0,right_wall,front_wall,back_wall,bottom_wall,top_wall
		//stiffness in x, y, z direction
        Real avg_rstiff,avg_zstiff;
        Real gain_r,gain_z;

        double wsxx,wsyy,wszz;//stress in the three directions
        double cylinder_force;
        double cylinder_stress;
        unsigned int num_pwall;//number of particles touching with the cylindrical wall
        unsigned int iterate_num,solve_num;
        bool Flag_ForceTarget;//forces acting on walls reach the target or not

		Real shearStress[2];//left and right wall
		Real left_facet_pos;
		//! Value of spheres volume (solid volume)
		Real particlesVolume;

		Real previousStress [2];//shear stress at the previous step.
		//! Value of box volume
		Real Init_boxVolume;
		Real Current_boxVolume;
		//! Sample porosity
		Real porosity;
		Real fmax;//the maximum of shear f when shearing
		//strain
		Real strain[3];

		virtual ~FlexCompressionEngine();

		virtual void action();
        double get_particleVolume();
        double get_boxVolume();
		Vector2r getStress();
        void initialize();
        void get_boundary_nodes();
        int get_numWallFacet();
        //for rigid cylindrical wall
        void rigid_consolidate();
        void rigid_shear();
        void fix_wall();
        void free_wall();
        void applySurfaceLoad();
        //for flexible cylindrical wall
        void flexible_consolidate();
        void flexible_shear();
		//! Compute the mean/max unbalanced force in the assembly (normalized by mean contact force)
    	Real ComputeUnbalancedForce(bool maxUnbalanced=false);
		Real getResultantF();//resultant force acting on the walls and particles in the down box.
		//recording
		void recordData();

		//check force on walls
		bool checkForce(int wall, Vector3r resultantForce);
		void checkMajorF(bool init);//calculate stress of bottom and top walls

		///Change physical propertieaabbs of interactions and/or bodies in the middle of a simulation (change only friction for the moment, complete this function to set cohesion and others before compression test)
		void setContactProperties(Real frictionDegree);
		void loadingStress();



        void get_gainz();
        void servo(double dt);
        void servo_cylinder(double dt);

        void consol_ss_cylinder();
        void checkTarget();
        void quiet_system();
        void cylinwall_go();
        void shear();
        unsigned int rampNum;//chunck number for accelerating walls
		//void getStressStrain();//
		SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(FlexCompressionEngine,BoundaryController,
		"An engine maintaining constant stresses or constant strain rates on some boundaries of a parallepipedic packing. The stress/strain control is defined for each axis using :yref:`FlexCompressionEngine::stressMask` (a bitMask) and target values are defined by goal1,goal2, and goal3. sigmaIso has to be defined during growing phases."
		"\n\n.. note::\n\t The algorithms used have been developed initialy for simulations reported in [Chareyre2002a]_ and [Chareyre2005]_. They have been ported to SudoDEM in a second step and used in e.g. [Kozicki2008]_,[Scholtes2009b]_,[Jerier2010b]."
		"The engine perform a triaxial compression with a control in direction 'i' in stress (if stressControl_i) else in strain.\n\n"
		"For a stress control the imposed stress is specified by 'sigma_i' with a 'max_veli' depending on 'strainRatei'. To obtain the same strain rate in stress control than in strain control you need to set 'wallDamping = 0.8'.\n"
		"For a strain control the imposed strain is specified by 'strainRatei'.\n"
		"With this engine you can also perform internal compaction by growing the size of particles by using ``TriaxialStressController::controlInternalStress``. For that, just switch on 'internalCompaction=1' and fix sigma_iso=value of mean pressure that you want at the end of the internal compaction.\n"
		,
		//((bool,firstRun,true,,"initialize some data at the first run"))
        ((double,wall_radius,5.0,,"the radius of the cylindrical wall"))
        ((bool,wall_fix,true,,"the radius of the cylindrical wall"))
        ((double,gain_alpha,0.5,,"alpha during geting gains"))
        ((double,target_strain,0.15,,"the target axial strain during shear"))
        ((unsigned int,iterate_num_max,100,,"interval of restart counting iterate_num"))
        ((unsigned int,solve_num_max,2000,,"interval of restart counting solve_num"))
        ((bool,z_servo,true,,"stress control or strain control in the z direction"))
        ((bool,isFlexible,false,,"use the flexible membrane?"))
        ((bool,keepRigid,true,,"keep the flexible wall rigid? Usually used during consolidation."))
   		((unsigned int,stiffnessUpdateInterval,100,,"interval of stiffness updating at walls"))
   		((unsigned int,echo_interval,10,,"interval of iteration for printing run info."))
        ((unsigned int,ramp_interval,10000,,"chunck bins for accelerating walls gradually."))
        ((unsigned int,ramp_chunks,400,,"chuncks for accelerating walls gradually."))
   		((unsigned int,stressMask,7,,"constant confining stress (non-zero) or contant volume (0)"))
		((bool,bottom_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((bool,top_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((std::string,file,"recordData",,"Name of file to save to; must not be empty."))
		((Real,height,0,Attr::readonly,"size of the box z(2-axis) |yupdate|"))
	
		((Real,f_threshold,0.05,,"threshold for approaching the target fast."))
        ((Real,fmin_threshold,0.2,,"threshold for stresses on the walls relaxation."))
		((Real,alpha,0.5,,"gain for stress control"))
        ((Real,wall_radius0,0,,"the initial radius of the cylinderical wall"))
		((Real,height0,0,,"Reference size for strain definition. See :yref:`FlexCompressionEngine::height`"))

		((Real,wall_equal_tol,0.001,,"wall equilibrium tolerance after each increament."))
		((Real,majorF_tol,0.985,,"stress gradual reduction torlerance of top and bottom walls when compressing"))
		((Real,unbf_tol,0.5,,"tolerance of UnbalancedForce"))
		((Real,wall_max_vel,1.,,"max velocity of the top  facet when loading"))
		((int,step_interval,20,,"step interval during each displacement increament"))
		((bool,Sphere_on,false,,"true for spherical particles, false for polyhedral particles"))
		((bool,hydroStrain,false,,"this option is used to perform a hydro strain rate on all directions "))
		((Real,hydroStrainRate,0.0,,"this option is used to perform a hydro strain rate on all directions "))
		((bool,debug,false,,"this option is used to output debug info"))
		((Real,goalr,0,,"prescribed stress/strain rate in the radial direction, as defined by :yref:`FlexCompressionEngine::stressMask` (see also :yref:`FlexCompressionEngine::isAxisymetric`)"))
		
		((Real,goalz,0,,"prescribed stress/strain rate on axis z, as defined by :yref:`FlexCompressionEngine::stressMask` (see also :yref:`FlexCompressionEngine::isAxisymetric`)"))
		((unsigned int,total_incr,24000,,"the total increaments"))
		((int,outinfo_interval,5,,"the iteration interval at which displaying the info for monitoring"))
		((unsigned int,savedata_interval,12,,"the interval steps between two savedata operations"))
		((Real,max_vel,1,,"Maximum allowed walls velocity [m/s]. This value superseeds the one assigned by the stress controller if the later is higher. max_vel can be set to infinity in many cases, but sometimes helps stabilizing packings. Based on this value, different maxima are computed for each axis based on the dimensions of the sample, so that if each boundary moves at its maximum velocity, the strain rate will be isotropic (see e.g. :yref:`FlexCompressionEngine::max_vel1`)."))
		((Real,externalWork,0,Attr::readonly,"Energy provided by boundaries.|yupdate|"))
		((Real, shearRate,0,,"target shear rate in direction x (./s)"))
		((Real, currentStrainRate2,0,,"current strain rate in direction 2 - converging to :yref:`ThreeDTriaxialEngine::strainRate2` (./s)"))
		((Real, UnbalancedForce,1,,"mean resultant forces divided by mean contact force"))
		((Real, frictionAngleDegree,-1,,"Value of friction used in the simulation if (updateFrictionAngle)"))
		((bool, updateFrictionAngle,false,,"Switch to activate the update of the intergranular frictionto the value :yref:`ThreeDTriaxialEngine::frictionAngleDegree`."))
		((bool, Confining,true,,"Switch to confine or shear"))
		((std::string,Key,"",,"A string appended at the end of all files, use it to name simulations."))
		,
		/* extra initializers */
		,
		//
		,
   		/* constructor */
   		firstRun =true;
		majorok = true;
		strain[0]=strain[1]=strain[2] = 0;
		x_area = 0;
		numWallFacet = 0;
		porosity=1;
		fmax = 0;
		//
		rampNum = 0;
		Init_boxVolume=0;Current_boxVolume=0;
        avg_rstiff = avg_zstiff = 0.0;
        gain_z = gain_r = 0.0;
        iterate_num = solve_num =0;
        Flag_ForceTarget = false;//forces acting on walls reach the target or not
		//
		,
		.def_readonly("strain",&FlexCompressionEngine::strain,"Current strain in a vector (exx,eyy,ezz). The values reflect true (logarithmic) strain.")
		.def_readonly("porosity",&FlexCompressionEngine::porosity,"Porosity of the packing.")
		//.def_readonly("boxVolume",&FlexCompressionEngine::boxVolume,"Total packing volume.")
		.def("getStress",&FlexCompressionEngine::getStress,"get stress within the assembly.")
        .def("getBoxVolume",&FlexCompressionEngine::get_boxVolume,"get the specimen volume.")
        .def("fix_wall",&FlexCompressionEngine::fix_wall,"set no degree of freedom for all nodes.")
        .def("free_wall",&FlexCompressionEngine::free_wall,"set XY DOF of the nodes contacting with the walls, and free other nodes.")
        .def("getParticleVolume",&FlexCompressionEngine::get_particleVolume,"get the volume of all particles within the specimen.")
		//.def("setContactProperties",&FlexCompressionEngine::setContactProperties,"Assign a new friction angle (degrees) to dynamic bodies and relative interactions")
		)
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(FlexCompressionEngine);

