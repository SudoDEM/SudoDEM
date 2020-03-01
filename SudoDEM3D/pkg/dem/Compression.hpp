/*************************************************************************
*  Copyright (C) 2016 by Zhswee                                          *
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
#include"sudodem/pkg/dem/NewtonIntegrator.hpp"
#include<fstream>
#include<string>
#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif

#define CUBIC_SAMPLE//
#define SOLE_GAIN //move two walls with their corresponding gain in the same direction
//testing contact model
#define FM//FricMat
#define NEW_SCHEME //using the new scheme for consolidation

class State;

class Scene;
class State;


/*! \brief Controls the stress on the boundaries of a box and compute strain-like and stress-like quantities for the packing. The algorithms used have been developed initialy for simulations reported in [Chareyre2002a] and [Chareyre2005]. They have been ported to SudoDEM in a second step and used in e.g. [Kozicki2008],[Scholtes2009b],[Jerier2010b].
*/

class CompressionEngine : public BoundaryController
{
	private :
		//! is this the beginning of the simulation, after reading the scene? -> it is the first time that SudoDEM passes trought the engine ThreeDTriaxialEngine
		bool firstRun;
		std::ofstream out;
		bool majorok;
		shared_ptr<Body> box [6];//pointers of facets of up box
		shared_ptr<NewtonIntegrator> newton;
		//! internal index values for retrieving wall,top wall and up box facets
		enum {left_wall=0,right_wall,front_wall,back_wall,bottom_wall,top_wall};

		inline const Vector3r getForce(Scene* rb, Body::id_t id){ return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */ }
	public :
		//Vector3r translationAxisy;
		//Vector3r translationAxisx;
		//Vector3r translationAxisz;//store the rotation degrees.Not changed the variable name because of the initial test sample existing.

		//! The value of stiffness (updated according to stiffnessUpdateInterval)
		Real wall_stiffness[6];//left_wall=0,right_wall,front_wall,back_wall,bottom_wall,top_wall

		Real 	pos[4];
		Real 	left_wpos[4];
		Real 	right_wpos[4];
		Real 	front_wpos[4];
		Real 	back_wpos[4];
		Real 	bottom_wpos[4];
		Real 	top_wpos[4];
		Real x_area; //area in the x axis
		Real y_area;
		Real z_area;
    Real loading_area;//for cylinder wall
		Vector2r wall_pos[6];//left_wall=0,right_wall,front_wall,back_wall,bottom_wall,top_wall
		//stiffness in x, y, z direction
    Real avg_xstiff;
		Real avg_ystiff;
    Real avg_zstiff;
    Real gain_x, gain_y, gain_z;
    #ifdef SOLE_GAIN
    double wsxx_left,wsxx_right,wsyy_front,wsyy_back,wszz_top,wszz_bottom;
    Real gain_x1, gain_y1, avg_xstiff1, avg_ystiff1, gain_z1, avg_zstiff1;
    #endif
    double wsxx,wsyy,wszz;//stress in the three directions
    double cylinder_force;
    double cylinder_stress;
    unsigned int num_pwall;//number of particles touching with the cylindrical wall
    unsigned int iterate_num,solve_num;
    bool Flag_ForceTarget;//forces acting on walls reach the target or not

		//Real shearStress[2];//left and right wall
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
    //wall stresses
    double left_wall_s;
    double right_wall_s;
    double front_wall_s;
    double back_wall_s;
    double bottom_wall_s;
    double top_wall_s;
		virtual ~CompressionEngine();

		virtual void action();

		Vector2r getStress();
		Matrix3r getStressTensor();
    //void action_old();
    void action_cubic();
		void action_cylindrical();
		void hydroConsolidation();
		void generalConsolidation();
		//! update the stiffness of boundary-packing interaction (sum of contacts stiffness on the boundary)
		void updateStiffness();
		//! Compute stresses on walls as "Vector3r stress[6]", compute meanStress, strain[3] and mean strain
		//void computeStressStrain();
		//! Compute the mean/max unbalanced force in the assembly (normalized by mean contact force)
  	Real ComputeUnbalancedForce(bool maxUnbalanced=false);
		//Real getResultantF();//resultant force acting on the walls and particles in the down box.
		//recording
		void recordData();

		//check force on walls
		//bool checkForce(int wall, Vector3r resultantForce);
		//void checkMajorF(bool init);//calculate stress of bottom and top walls

		//void move_wall(int wall_id, double Sign, int f_index,double (&wall_pos)[4],double goal_force, double wall_f);
		//void move_wall(int wall_id, double Sign, int f_index,double area,double goal_force, double wall_f);
		void loadingStress();
		void updateBoxSize();
		void getBox();
    void get_gain();
    void get_gainz();
    void servo(double dt);
    void servo_cylinder(double dt);
    void consol_ss();
    void consol_ss_cylinder();
    void checkTarget();
    //void quiet_system();
    void cylinwall_go();
    void shear();
    unsigned int rampNum;//chunck number for accelerating walls
		//void getStressStrain();//
		SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(CompressionEngine,BoundaryController,
		"An engine maintaining constant stresses or constant strain rates on some boundaries of a parallepipedic packing. The stress/strain control is defined for each axis using :yref:`CompressionEngine::stressMask` (a bitMask) and target values are defined by goal1,goal2, and goal3. sigmaIso has to be defined during growing phases."
		"\n\n.. note::\n\t The algorithms used have been developed initialy for simulations reported in [Chareyre2002a]_ and [Chareyre2005]_. They have been ported to SudoDEM in a second step and used in e.g. [Kozicki2008]_,[Scholtes2009b]_,[Jerier2010b]."
		"The engine perform a triaxial compression with a control in direction 'i' in stress (if stressControl_i) else in strain.\n\n"
		"For a stress control the imposed stress is specified by 'sigma_i' with a 'max_veli' depending on 'strainRatei'. To obtain the same strain rate in stress control than in strain control you need to set 'wallDamping = 0.8'.\n"
		"For a strain control the imposed strain is specified by 'strainRatei'.\n"
		"With this engine you can also perform internal compaction by growing the size of particles by using ``TriaxialStressController::controlInternalStress``. For that, just switch on 'internalCompaction=1' and fix sigma_iso=value of mean pressure that you want at the end of the internal compaction.\n"
		,
		((bool,cubic_specimen,true,,"the specimen is cubic (true) or cylindrical (false)."))
		((bool,z_servo,true,,"stress control or strain control in the z direction"))
		((bool,continueFlag,false,,"true for continuing consolidation or shear."))
		((unsigned int,shearMode,1,,"shear mode: 1 for conventional triaxial compression (CTC) (drained); 2 for CTC(undrained)"))
		//for cylinderical specimens
    ((double,wall_radius,5.0,,"the radius of the cylindrical wall"))
    ((bool,wall_fix,true,,"the radius of the cylindrical wall"))
		((Real,wall_radius0,0,,"the initial radius of the cylinderical wall"))

    ((double,gain_alpha,0.5,,"alpha during geting gains"))
    ((double,target_strain,0.15,,"the target axial strain during shear"))
    ((unsigned int,iterate_num_max,100,,"interval of restart counting iterate_num"))
    ((unsigned int,solve_num_max,2000,,"interval of restart counting solve_num"))

		//((unsigned int,stiffnessUpdateInterval,100,,"interval of stiffness updating at walls"))
		((unsigned int,echo_interval,1000,,"interval of iteration for printing run info."))
    ((unsigned int,ramp_interval,10000,,"chunck bins for accelerating walls gradually."))//deprecated
    ((unsigned int,ramp_chunks,400,,"chuncks for accelerating walls gradually."))//deprecated

		((Real,f_threshold,0.01,,"threshold for approaching the target fast."))
    ((Real,fmin_threshold,0.01,,"threshold for stresses on the walls relaxation."))
		((Real,unbf_tol,0.01,,"tolerance of UnbalancedForce"))
		((Real,wall_max_vel,1.,,"max velocity of the top  facet when loading"))
		((Real,alpha,0.5,,"gain for stress control"))//deprecated
    //for cubic specimens
		((Real,height0,0,,"Reference size for strain definition. See :yref:`CompressionEngine::height`"))
		((Real,width0,0,,"Reference size for strain definition. See :yref:`CompressionEngine::width`"))
		((Real,depth0,0,,"Reference size for strain definition. See :yref:`CompressionEngine::depth`"))
		((Real,height,0,Attr::readonly,"size of the box z(2-axis) |yupdate|"))
		((Real,width,0,Attr::readonly,"size of the box x(0-axis) |yupdate|"))
		((Real,depth,0,Attr::readonly,"size of the box y(1-axis) |yupdate|"))
		((bool,left_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((bool,right_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((bool,front_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((bool,back_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((bool,bottom_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))
		((bool,top_wall_activated,true,,"if true, this wall moves according to the target value (stress or strain rate)."))



		((bool,hydroStrain,false,,"this option is used to perform a hydro strain rate on all directions "))
		((Real,hydroStrainRate,0.0,,"this option is used to perform a hydro strain rate on all directions "))
		((bool,debug,false,,"this option is used to output debug info"))
		((Real,goalx,0,,"prescribed stress/strain rate on axis x"))
		((Real,goaly,0,,"prescribed stress/strain rate on axis y"))
		((Real,goalz,0,,"prescribed stress/strain rate on axis z"))
		((unsigned int,savedata_interval,12,,"the interval steps between two savedata operations"))
		((Real,max_vel,1,,"Maximum allowed walls velocity [m/s]. This value superseeds the one assigned by the stress controller if the later is higher. max_vel can be set to infinity in many cases, but sometimes helps stabilizing packings. Based on this value, different maxima are computed for each axis based on the dimensions of the sample, so that if each boundary moves at its maximum velocity, the strain rate will be isotropic (see e.g. :yref:`CompressionEngine::max_vel1`)."))
		((Real,externalWork,0,Attr::readonly,"Energy provided by boundaries.|yupdate|"))
		((Real, UnbalancedForce,1,,"mean resultant forces divided by mean contact force"))
		((bool, interStressControl,false,,"control walls' movement based on the principal stress within the assembly if true, otherwise based on the stresses acting on the walls"))

		((std::string,file,"recordData",,"Name of file to save to; must not be empty."))
		,
		/* extra initializers */
		,
		//
		,
   		/* constructor */
   		firstRun =true;
		majorok = true;
		strain[0]=strain[1]=strain[2] = 0;
		x_area = y_area = z_area = 0;
		for(int j = 0; j < 4; j++) pos[j] = 0;
		for(int j = 0; j < 4; j++) {left_wpos[j] = 0;right_wpos[j] = 0;front_wpos[j] = 0;back_wpos[j] = 0;bottom_wpos[j] = 0;top_wpos[j] = 0;}
		for(int j = 0; j < 6; j++) {wall_pos[j] = Vector2r::Zero();wall_stiffness[j] = 0.;}
		left_facet_pos = 0.;
		porosity=1;
		fmax = 0;
		//
		left_wall_s = 0;
    right_wall_s = 0;
    front_wall_s = 0;
    back_wall_s = 0;
    bottom_wall_s = 0;
    top_wall_s = 0;
		//
		rampNum = 0;
		Init_boxVolume=0;Current_boxVolume=0;
    avg_xstiff = avg_ystiff = avg_zstiff = 0.0;
    gain_x = gain_y = gain_z = 0.0;
    iterate_num = solve_num =0;
    Flag_ForceTarget = false;//forces acting on walls reach the target or not
		//
		,
		.def_readonly("strain",&CompressionEngine::strain,"Current strain in a vector (exx,eyy,ezz). The values reflect true (logarithmic) strain.")
		.def_readonly("porosity",&CompressionEngine::porosity,"Porosity of the packing.")
		//.def_readonly("boxVolume",&CompressionEngine::boxVolume,"Total packing volume.")
		.def("getStress",&CompressionEngine::getStress,"get stress within the assembly.")
		)
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(CompressionEngine);
