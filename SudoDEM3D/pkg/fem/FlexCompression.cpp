/*************************************************************************
*  Copyright (C) 2017 by Sway Zhao		             	        		 *
*  zhswee@gmail.com      					  							 *
*																		 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"FlexCompression.hpp"

#include<sudodem/pkg/common/Box.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/Superquadrics.hpp>
//#include<sudodem/pkg/dem/GJKParticle.hpp>
//create a new class base,'NonSphericalShape', in the furture.
#ifdef FM
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>
#endif

#include<sudodem/core/State.hpp>
#include<assert.h>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/dem/Shop.hpp>

#include<sudodem/core/Omega.hpp>
#include<sudodem/lib/base/Math.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/lambda/lambda.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>


typedef Eigen::Matrix<double, 3, 2> Matrix32d;
typedef Eigen::Quaternion<double> Quaternionr;
using namespace Eigen;

#include <sys/syscall.h>//get pid of a thread

CREATE_LOGGER(FlexCompressionEngine);
SUDODEM_PLUGIN((FlexCompressionEngine));

FlexCompressionEngine::~FlexCompressionEngine(){}

Real FlexCompressionEngine::ComputeUnbalancedForce( bool maxUnbalanced) {return Shop::unbalancedForce(maxUnbalanced,scene);}

double FlexCompressionEngine::get_particleVolume(){
    //get the volume of all solid particles
        //calculate particlesVolume
		BodyContainer::iterator bi = scene->bodies->begin();
		BodyContainer::iterator biEnd = scene->bodies->end();
		double p_v = 0;
		for ( ; bi!=biEnd; ++bi )
		{
			//if((*bi)->isClump()) continue;//no clump for Superquadrics
			const shared_ptr<Body>& b = *bi;
			//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";
            if (b->shape->getClassName()=="Sphere"){
				const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b->shape);
				p_v += 4.0/3.0*M_PI*pow(A->radius,3.0);

			}
            //#else
			if (b->shape->getClassName()=="Superquadrics"){
				const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
				p_v += A->getVolume();

			}
            /*
			if (b->shape->getClassName()=="GJKParticle"){
				const shared_ptr<GJKParticle>& A = SUDODEM_PTR_CAST<GJKParticle> (b->shape);
				p_v += A->getVolume();

			}*/

		}
    return p_v;
}
void FlexCompressionEngine::get_boundary_nodes(){//this function is just called only once when initializing the engine
    SUDODEM_PARALLEL_FOREACH_NODE_BEGIN(const shared_ptr<Node>& b, scene->nodes){
        const Node::id_t& id=b->getId();
        //to do:find the nodes contacting with the top and bottom walls, then store them to top_nodes and bottom_nodes

    } SUDODEM_PARALLEL_FOREACH_NODE_END();

}
//consolidation procedure with a rigid wall
void FlexCompressionEngine::rigid_consolidate(){
    const Real& dt = scene->dt;
    if (!Flag_ForceTarget){//forces on walls haven't reached the target, so we need to switch on servo.
        if (iterate_num < 1){solve_num = 0;quiet_system();}
        iterate_num += 1;
        servo_cylinder(dt);
        if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
            iterate_num = 1;
            get_gainz();
            checkTarget();//check the forces on the walls reach the target or not, and info is stored in Flag_ForceTarget.
        }
    }else{//keep walls fixed, then cyle to a stable state
        if (solve_num < 1){quiet_system();iterate_num = 0;}//set velocities of all walls and particles to zero.
        solve_num += 1;
        //do nothing!
        if (solve_num >= 1000){//check the stability of the sample
            UnbalancedForce = ComputeUnbalancedForce ();
            solve_num = 1;
            consol_ss_cylinder();
            if (((std::abs(cylinder_stress - goalr)/goalr) > fmin_threshold)||
                ((std::abs(wszz - goalz)/goalz) > fmin_threshold)){
                    //stresses on the walls relax too much
                    Flag_ForceTarget = false;
            }
             if (UnbalancedForce < unbf_tol){//o.5 by default
                checkTarget();


                if (Flag_ForceTarget){
                    //stop running
                    scene->stopAtIter = scene->iter + 1;
                    std::cerr<<"consolidation completed!"<<std::endl;
                }
            }
        }

    }
    //output info
    if (scene->iter % savedata_interval == 0){recordData();}
    if (scene->iter % echo_interval == 0){
	    UnbalancedForce = ComputeUnbalancedForce ();
        consol_ss_cylinder();
	    //getStressStrain();
	    std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
	    <<", conf:"<<cylinder_stress/1000.0 //average stress of the two walls, unit is kPa
	    <<", Ss_z:"<<wszz/1000.0
      <<", e:" << ((M_PI*pow(wall_radius,2.0)*height)/particlesVolume - 1.0)
	    <<"\n";
    }

}
//shear procedure with a rigid wall
void FlexCompressionEngine::rigid_shear(){
    const Real& dt = scene->dt;
    if (rampNum < ramp_interval){//accelerating the walls
        if (rampNum < 1){// shear following consolidation
            //reset the initial data
            if(Flag_ForceTarget){//consolidation completed normally
                std::cerr<<"Shear begins..."<<std::endl;
            }
        }
        rampNum ++;
        double vel= float(rampNum/(ramp_interval/ramp_chunks)+1)/float(ramp_chunks)*goalz;//velocity
        //std::cerr<<"vel of loading walls = "<<vel<<std::endl;
        bottom_wall->state->vel[2] = vel;//z
        top_wall->state->vel[2] = -vel;//z
    }
    iterate_num += 1;
    servo_cylinder(dt);
    if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
        iterate_num = 1;
        get_gainz();
    }

    //height = top_wall->state->pos.z() - bottom_wall->state->pos.z();

    if (log(height0/height) > target_strain){//
        //stop shear
        scene->stopAtIter = scene->iter + 1;
    }
    	//output info
    if (scene->iter % savedata_interval == 0){recordData();}
    if (scene->iter % echo_interval == 0){
	    UnbalancedForce = ComputeUnbalancedForce ();
        //consol_ss();
	    //getStressStrain();
	    std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
	    <<", conf:"<<cylinder_stress/1000.0 //average stress of the two walls, unit is kPa
	    <<", Str_z:"<<log(height0/height)//z
	    <<", Ss_z:"<<wszz/1000.0
      <<", e:" << ((M_PI*pow(wall_radius,2.0)*height)/particlesVolume - 1.0)
	    <<"\n";
    }

}
//consolidation procedure with a flexible wall
void FlexCompressionEngine::flexible_consolidate(){
    const Real& dt = scene->dt;
    if(keepRigid){//keepRigid 
        if (!Flag_ForceTarget){//forces on walls haven't reached the target, so we need to switch on servo.
            if (iterate_num < 1){solve_num = 0;quiet_system();}
            iterate_num += 1;
            servo_cylinder(dt);
            //Real udz = gain_z * (wszz - goalz)/dt;
            //bottom_wall->state->vel[2] = -udz;//z
            //top_wall->state->vel[2] = udz;//z
            if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
                iterate_num = 1;
                get_gainz();
                checkTarget();//check the forces on the walls reach the target or not, and info is stored in Flag_ForceTarget.
            }
        }else{//keep walls fixed, then cyle to a stable state
            if (solve_num < 1){quiet_system();iterate_num = 0;}//set velocities of all walls and particles to zero.
            solve_num += 1;
            //do nothing!
            if (solve_num >= 1000){//check the stability of the sample
                UnbalancedForce = ComputeUnbalancedForce ();
                solve_num = 1;
                consol_ss_cylinder();
                if (((std::abs(cylinder_stress - goalr)/goalr) > fmin_threshold)||
                    ((std::abs(wszz - goalz)/goalz) > fmin_threshold)){
                        //stresses on the walls relax too much
                        Flag_ForceTarget = false;
                }
                 if (UnbalancedForce < unbf_tol){//o.5 by default
                    checkTarget();


                    if (Flag_ForceTarget){
                        //stop running
                        scene->stopAtIter = scene->iter + 1;
                        std::cerr<<"consolidation completed!"<<std::endl;
                    }
                }
            }

        }
    }else{//not keepRigid begin
        if (!Flag_ForceTarget){//forces on walls haven't reached the target, so we need to switch on servo.
        if (iterate_num < 1){solve_num = 0;}
        iterate_num += 1;
        //servo_cylinder(dt);
        consol_ss_cylinder();
        //std::cout<<"wszz"<<wszz<<" goalz"<<goalz<<" gain_z"<<gain_z<<std::endl;
        Real udz = gain_z * (wszz - goalz);
        //bottom_wall->state->vel[2] = -udz;//z
        //top_wall->state->vel[2] = udz;//z
	    udz = Mathr::Sign(udz)*min(fabs(udz),max_vel*dt);
        bottom_wall->state->pos[2] += -udz;//z
        top_wall->state->pos[2] += udz;//z
        for(vector<shared_ptr<Node>>::iterator b = top_nodes.begin();b!= top_nodes.end();++b){
            (*b)->pos[2] += udz;//z 
        }
        for(vector<shared_ptr<Node>>::iterator b = bottom_nodes.begin();b!= bottom_nodes.end();++b){
            (*b)->pos[2] -= udz;//z 
        }
        if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
            iterate_num = 1;
            //get_gainz();
            avg_zstiff = 0.0;
            for(InteractionContainer::iterator ii = scene->interactions->begin();ii!=scene->interactions->end();++ii ) if ((*ii)->isReal())
            {
	            const shared_ptr<Interaction>& contact = *ii;
	            int id1 = contact->getId1(), id2 = contact->getId2();
	            int id = (id1 < id2) ? id1 : id2;
	            if (id<2){//caution: I assume that wall_id is less than 2. top and bottom walls
                    FrictPhys* currentContactPhysics =
                    static_cast<FrictPhys*> ( contact->phys.get() );

                    avg_zstiff  += currentContactPhysics->kn;
                }
            }
            gain_z = 2.0*gain_alpha*loading_area/(avg_zstiff);
            checkTarget();//check the forces on the walls reach the target or not, and info is stored in Flag_ForceTarget.
        }
        }else{//keep walls fixed, then cyle to a stable state
            if (solve_num < 1){quiet_system();iterate_num = 0;}//set velocities of all walls and particles to zero.
            solve_num += 1;
            //do nothing!
            if (solve_num >= solve_num_max){//check the stability of the sample
                UnbalancedForce = ComputeUnbalancedForce ();
                solve_num = 1;
                consol_ss_cylinder();
                if ((std::abs(wszz - goalz)/goalz) > fmin_threshold){
                        //stresses on the walls relax too much
                        Flag_ForceTarget = false;
                }
                 if (UnbalancedForce < unbf_tol){//o.5 by default
                    checkTarget();

                    if (Flag_ForceTarget){
                        //stop running
                        scene->stopAtIter = scene->iter + 1;
                        std::cerr<<"consolidation completed!"<<std::endl;
                    }
                }
            }
        }
    }//not keepRigid end
    //output info
    if (scene->iter % savedata_interval == 0){recordData();}
    if (scene->iter % echo_interval == 0){
	    UnbalancedForce = ComputeUnbalancedForce ();
        consol_ss_cylinder();
	    //getStressStrain();
	    std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
	    <<", conf:"<<cylinder_stress/1000.0 //average stress of the two walls, unit is kPa
	    <<", Ss_z:"<<wszz/1000.0
      <<", e:" << ((M_PI*pow(wall_radius,2.0)*height)/particlesVolume - 1.0)
	    <<"\n";
    }
    

}
//shear procedure with a flexible wall
void FlexCompressionEngine::flexible_shear(){
    const Real& dt = scene->dt;
    if (rampNum < ramp_interval){//accelerating the walls
        if (rampNum < 1){// shear following consolidation
            //reset the initial data
            if(Flag_ForceTarget){//consolidation completed normally
                std::cerr<<"Shear begins..."<<std::endl;
            }
        }
        rampNum ++;
        double ramp_inter = 0.5*ramp_interval;
        double vel=0.0,udz=0.0;
        if (rampNum < ramp_inter){
            vel =float(rampNum/(ramp_inter/ramp_chunks)+1)/float(ramp_chunks)*goalz*2.0;
        }else{
            vel =(float((ramp_interval - rampNum)/(ramp_inter/ramp_chunks)+1)/float(ramp_chunks)+1.0)*goalz;
        }
        udz = vel*dt*height;
        bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate
        top_wall->state->pos[2] += -udz;
        for(vector<shared_ptr<Node>>::iterator b = top_nodes.begin();b!= top_nodes.end();++b){
            (*b)->pos[2] -= udz;//z 
        }
        for(vector<shared_ptr<Node>>::iterator b = bottom_nodes.begin();b!= bottom_nodes.end();++b){
            (*b)->pos[2] += udz;//z 
        }
    }else{
        /*iterate_num += 1;
        servo_cylinder(dt);
        if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
            iterate_num = 1;
            get_gainz();
        }*/
        double udz = goalz*dt*height;
        bottom_wall->state->pos[2] += udz;//z,remain a constant loading strain rate
        top_wall->state->pos[2] += -udz;//z
        for(vector<shared_ptr<Node>>::iterator b = top_nodes.begin();b!= top_nodes.end();++b){
            (*b)->pos[2] -= udz;//z 
        }
        for(vector<shared_ptr<Node>>::iterator b = bottom_nodes.begin();b!= bottom_nodes.end();++b){
            (*b)->pos[2] += udz;//z 
        }
    }
    //height = top_wall->state->pos.z() - bottom_wall->state->pos.z();

    if (log(height0/height) > target_strain){//
        //stop shear
        scene->stopAtIter = scene->iter + 1;
    }
    	//output info
    if (scene->iter % savedata_interval == 0){recordData();}
    if (scene->iter % echo_interval == 0){
	    UnbalancedForce = ComputeUnbalancedForce ();
        //consol_ss();
        consol_ss_cylinder();
	    //getStressStrain();
	    std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
	    <<", conf:"<<cylinder_stress/1000.0 //average stress of the two walls, unit is kPa
	    <<", Str_z:"<<log(height0/height)//z
	    <<", Ss_z:"<<wszz/1000.0
      <<", e:" << ((get_boxVolume())/particlesVolume - 1.0)
	    <<"\n";
    }
}


int FlexCompressionEngine::get_numWallFacet(){
    //get the total number of walls and TriElements
    int num = 0;
    SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){//Elements should be stored in a separated list for speeding up traversing
        // clump members are handled inside clumps
        //if(b->isClumpMember()) continue;
        if ((b->shape->getClassName()=="TriElement")||(b->shape->getClassName()=="Wall")){
            num++;
        }
    } SUDODEM_PARALLEL_FOREACH_BODY_END();
    return num;
}
void FlexCompressionEngine::fix_wall(){
    if(firstRun){initialize();}
 
    for(vector<shared_ptr<Node>>::iterator b = top_nodes.begin();b!=top_nodes.end();++b){
        (*b)->blockedDOFs = Node::DOF_ALL;
    }
    for(vector<shared_ptr<Node>>::iterator b = bottom_nodes.begin();b!=bottom_nodes.end();++b){
        (*b)->blockedDOFs = Node::DOF_ALL;
    }
    for(vector<shared_ptr<Node>>::iterator b = free_nodes.begin();b!=free_nodes.end();++b){
        (*b)->blockedDOFs = Node::DOF_ALL;
    }
}

void FlexCompressionEngine::free_wall(){
    std::cout<<"is first run? "<<firstRun<<std::endl;
    if(firstRun){initialize();}

    for(vector<shared_ptr<Node>>::iterator b = top_nodes.begin();b!=top_nodes.end();++b){
        (*b)->blockedDOFs = Node::DOF_Z;//Node::DOF_XYZ;
    }
    for(vector<shared_ptr<Node>>::iterator b = bottom_nodes.begin();b!=bottom_nodes.end();++b){
        (*b)->blockedDOFs = Node::DOF_Z;//Node::DOF_XYZ;
    }
    for(vector<shared_ptr<Node>>::iterator b = free_nodes.begin();b!=free_nodes.end();++b){
        (*b)->blockedDOFs = Node::DOF_NONE;
    }
}

void FlexCompressionEngine::initialize(){
	LOG_INFO ( "First run, will initialize!" );

    bottom_wall = Body::byId(0);//bottom wall
	top_wall = Body::byId(1);//top wall

	Real top = top_wall->state->pos.z();
	Real bottom = bottom_wall->state->pos.z();
    height = height0 = top - bottom;
    //store elements of the flexible wall
    std::cout<<"Initializing FlexCompressionEngine"<<std::endl;
     
    for(BodyContainer::iterator b = scene->bodies->begin(); b!=scene->bodies->end(); ++b ){
        //const shared_ptr<Body>& b
        if ((*b)->shape->getClassName()=="TriElement"){
            const shared_ptr<TriElement>& Elm = SUDODEM_PTR_CAST<TriElement> ((*b)->shape);
            Elements.push_back(Elm);
        }

    }
  
    //group nodes
    std::cout<<"Grouping nodes..."<<std::endl;
    double element_height = Elements.size()>0 ? std::fabs(Elements[0]->nodes[2]->pos[2]-Elements[0]->nodes[0]->pos[2]):0.01*height;//FIXME
    for(NodeContainer::iterator b = scene->nodes->begin(); b!=scene->nodes->end(); ++b ){
        //group nodes
        if ((*b)->pos[2] < bottom + 0.5*element_height){//nodes contacting the bottom wall
            bottom_nodes.push_back(*b);
        }else if((*b)->pos[2] > top - 0.5*element_height){//nodes contacting the to wall
            top_nodes.push_back(*b);
        }else{//free nodes
            free_nodes.push_back(*b);
        }
    }
    applySurfaceLoad();
    //get number of all walls and TriElements
    numWallFacet = get_numWallFacet();

	
    wall_radius0 =wall_radius;
    loading_area = M_PI*pow(wall_radius,2.0);
    Init_boxVolume = height0*loading_area;
    particlesVolume = get_particleVolume();//get the volume of all particles
    rampNum = 0;
	firstRun=false;
    //
    iterate_num = 0;
    solve_num = 0;
    //get_gainz();
    //guss a gainz
    gain_z = 1e-15;
    gain_r = 1e-15;
    Flag_ForceTarget = false;
    //others
    if(isFlexible && !keepRigid){
        std::cout<<"goalr minus"<<goalr<<std::endl;
        //goalr = -std::fabs(goalr);
    }else{
        goalr = std::fabs(goalr);
        //std::cout<<"goalr plus"<<goalr<<std::endl;
    }
}
void FlexCompressionEngine::action()
{ 
    if(isFlexible){//use a flexible wall
        if ( firstRun ){initialize();}
        if (z_servo){//consolidation
            flexible_consolidate();
	    }else{//shear
            flexible_shear();
        }
    }else{//use a rigid wall
        cylinwall_go();
        if (wall_fix){return;}
        if ( firstRun ){initialize();}
        if (z_servo){//consolidation
            rigid_consolidate();
	    }else{//shear
            rigid_shear();
        }
    }
}
double FlexCompressionEngine::get_boxVolume(){
    //To do: calculate box volume with flexible membrane facets.
    //std::cout<<"watch point 1"<<std::endl;
    if(firstRun){initialize();}
    if(isFlexible){//use a flexible wall
        double box_v = 0.0;
        //std::cout<<"watch point 2"<<std::endl;
        //select a based point at the center of the specimen to calculate the box volume
        //volumes of the top and bottom parts (each one looks like a cone when the resolution is sufficient.)
        box_v += height*M_PI*pow(wall_radius,2.0)/3.0;//FIXME:the top and bottom ends are not exactly circular.
        Vector3r base_point = 0.5*(top_wall->state->pos + bottom_wall->state->pos);//FIXME:both walls are centered at x=y=0, by default.
        //std::cout<<"base_point"<<base_point<<std::endl;
        for(BodyContainer::iterator b = scene->bodies->begin(); b!=scene->bodies->end(); ++b ){//Elements should be stored in a separated list for speeding up traversing
            // clump members are handled inside clumps
            //if(b->isClumpMember()) continue;
            if ((*b)->shape->getClassName()!="TriElement") continue;
            const shared_ptr<TriElement>& Elm = SUDODEM_PTR_CAST<TriElement> ((*b)->shape);
            box_v += Elm->TetraVolume(base_point);
	    }
        return box_v;
    }else{
        return height*M_PI*pow(wall_radius,2.0);
    }
}

void FlexCompressionEngine::applySurfaceLoad(){
    for(vector<shared_ptr<TriElement>>::iterator Elm = Elements.begin();Elm != Elements.end();++Elm){
        (*Elm)->surfLoad = goalr;
    }

}
Vector2r FlexCompressionEngine::getStress()
{
// at the beginning of the file; write column titles

    consol_ss_cylinder();
    Current_boxVolume = get_boxVolume();
    //double sigma1=0.0;
    //compute the stress within the sample, sigma_ij

    InteractionContainer::iterator ii    = scene->interactions->begin();
    InteractionContainer::iterator iiEnd = scene->interactions->end();

    Matrix3d Sigma = Matrix3d::Zero();
    for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
    {	Vector3r fn,fs;
	    const shared_ptr<Interaction>& contact = *ii;
        fn = (static_cast<FrictPhys*>	(contact->phys.get()))->normalForce;
		fs = (static_cast<FrictPhys*>	(contact->phys.get()))->shearForce;

	    if (fn.norm()==0){continue;}

	    int id1 = contact->getId1(), id2 = contact->getId2();
	    int id = (id1 < id2) ? id1 : id2;
	    if (id > numWallFacet){//caution: I assume that wall_id is less than 6.
	            const shared_ptr<Body>& b1 = Body::byId(id1,scene);
                const shared_ptr<Body>& b2 = Body::byId(id2,scene);
                State* s1 = b1->state.get();
                State* s2 = b2->state.get();
                Vector3d l = s2->pos - s1->pos;
				//Vector3r f = fn+ (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->shearForce;
                Sigma += (fn+fs)*(l.transpose());
        }else{//boundary
			/*	//const shared_ptr<Body>& b1 = Body::byId(id1,scene);
                const shared_ptr<Body>& b2 = Body::byId(id2,scene);
                //State* s1 = b1->state.get();
				Vector3r contactpoint = (static_cast<SuperquadricsGeom*>(contact->geom.get()))->contactPoint;
                State* s2 = b2->state.get();
                Vector3d l = s2->pos - contactpoint;
                Sigma += fn*(l.transpose());
			*/
		}
    }
	Sigma /= Current_boxVolume;
	double meanStress = Sigma.trace()/3.0;
	Sigma -= meanStress*Matrix3d::Identity();
	double deviatorStress = Sigma.norm()*sqrt(1.5);

	return Vector2r(meanStress,deviatorStress);

}
void FlexCompressionEngine::recordData()
{
	if(!out.is_open()){
		assert(!out.is_open());

		std::string fileTemp = file;


		if(fileTemp.empty()) throw ios_base::failure(__FILE__ ": Empty filename.");
		out.open(fileTemp.c_str(), truncate ? fstream::trunc : fstream::app);
		if(!out.good()) throw ios_base::failure(__FILE__ ": I/O error opening file `"+fileTemp+"'.");
	}

	Vector2r pq = getStress();

	if(out.tellp()==0)	out<<"AxialStrain VolStrain meanStress deviatorStress"<<endl;
	//getStressStrain();
        //strain[0] = 1.0 - width/width0;//x
        //strain[1] = 1.0 - depth/depth0;//y
        //strain[2] = 1.0 - height/height0;//z
	strain[2] = log(height0/height);//z
	//Init_boxVolume = height0*M_PI*pow(wall_radius,2.0);//width0*height0*depth0;

	//out << boost::lexical_cast<string> ( scene->iter ) << " "
 	//<< boost::lexical_cast<string> ( strain[0] ) << " "  //
 	//<< boost::lexical_cast<string> ( strain[1] ) << " "  //
	out << boost::lexical_cast<string> ( strain[2] ) << " "
	<< boost::lexical_cast<string> ( log(Init_boxVolume/Current_boxVolume) ) << " " //volumetric strain
 	<< boost::lexical_cast<string> ( pq[0] ) << " "
    << boost::lexical_cast<string> ( pq[1] ) << " "
	//<< boost::lexical_cast<string> ( syscall(SYS_gettid) ) << " "
 	<< endl;
}
void FlexCompressionEngine::cylinwall_go(){//cylindrical wall
    //traverse all particles
    BodyContainer::iterator bb    = scene->bodies->begin();
	BodyContainer::iterator iiEnd = scene->bodies->end();
    cylinder_force = 0.0;
    num_pwall = 0;
    if(!Sphere_on){//non-spherical particles
	    for(  ; bb!=iiEnd ; ++bb )
	    {

            //if ((*bb)->shape->getClassName()!="Superquadrics"){continue;}
            //calculate the overlap between the particle and the cylindrical wall


            const Body::id_t& id=(*bb)->getId();
            if (id<2){continue;}//the body is a (top or bottom) wall when id < 2
            State* state=(*bb)->state.get();
            Superquadrics* B = static_cast<Superquadrics*>((*bb)->shape.get());
            Vector2r dist = Vector2r(state->pos(0),state->pos(1));//at the section plane, xoy
            Vector3r B_pos = state->pos;
            Vector3r Origin = Vector3r(0.0,0.0,B_pos[2]);//the center of the cylindrical wall is at the origin (0,0) on the section plane by defaut
            //rough contact detection
            if ((dist.norm() + B->getr_max()) < wall_radius){//not touching
                continue;
            }
            double PenetrationDepth= 0;
            Vector3r contactpoint,normal;
            normal = Origin - B_pos;
            normal.normalize();
            // is the particle spherical?
            if (B->isSphere){
                    //the particle is spherical

                    contactpoint = B_pos;


                    PenetrationDepth = dist.norm() + B->getr_max() - wall_radius;
                    contactpoint = B_pos - normal*(B->getr_max()-0.5*PenetrationDepth);
                    //std::cerr<<PenetrationDepth<<std::endl;

	        }else{
                Vector3r wall_pos = -normal*wall_radius;
                Matrix3r rot_mat1 = B->rot_mat2local;//(se32.orientation).conjugate().toRotationMatrix();//to particle's system
	            Matrix3r rot_mat2 = B->rot_mat2global;//(se32.orientation).toRotationMatrix();//to global system

	            Vector2r phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
	            Vector3r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
	            //check whether point2 is at the negative side of the wall.
	            Vector3r p1;
                    p1 = wall_pos;
                    p1[2] = point2[2];
	            Vector3r d = point2 - p1;	//the distance vector
                PenetrationDepth = d.norm();
	            //d[PA] -= wall_pos[PA];
	            //cout<<"distance v="<<d<<endl;
	            if (normal.dot(d)>=0.0){//(normal.dot(d) < 0.) {//no touching between the wall and the particle
		            //cout<<"watch dog2"<<endl;
		            PenetrationDepth= 0;

		            continue;
	            }

            }
            num_pwall ++;
            //std::cerr<<" num pwall"<<num_pwall<<std::endl;
            Vector3r f = PenetrationDepth*1e8*normal;
            cylinder_force += f.norm();
            scene->forces.addForce (id,f);
            //not considering the torque performed on the particle from the wall.
        }
    }else{
	    for(  ; bb!=iiEnd ; ++bb )
	    {

            //if ((*bb)->shape->getClassName()!="Superquadrics"){continue;}
            //calculate the overlap between the particle and the cylindrical wall


            const Body::id_t& id=(*bb)->getId();
            if (id<2){continue;}//the body is a (top or bottom) wall when id < 2
            State* state=(*bb)->state.get();
            Sphere* B = static_cast<Sphere*>((*bb)->shape.get());
            Vector2r dist = Vector2r(state->pos(0),state->pos(1));//at the section plane, xoy
            Vector3r B_pos = state->pos;
            Vector3r Origin = Vector3r(0.0,0.0,B_pos[2]);//the center of the cylindrical wall is at the origin (0,0) on the section plane by defaut
            //rough contact detection
            if ((dist.norm() + B->radius) < wall_radius){//not touching
                continue;
            }
            double PenetrationDepth= 0;
            Vector3r contactpoint,normal;
            normal = Origin - B_pos;
            normal.normalize();

            contactpoint = B_pos;
            PenetrationDepth = dist.norm() + B->radius - wall_radius;
            contactpoint = B_pos - normal*(B->radius-0.5*PenetrationDepth);
            num_pwall ++;
            //std::cerr<<" num pwall"<<num_pwall<<std::endl;
            Vector3r f = PenetrationDepth*1e8*normal;
            cylinder_force += f.norm();
            scene->forces.addForce (id,f);
            //not considering the torque performed on the particle from the wall.
        }

    }

}


void FlexCompressionEngine::get_gainz(){//determine servo gain parameters
   //cylindrical wall
    
    if(!z_servo){return;}
    avg_zstiff = 0.0;
    avg_rstiff = 0.0;//total stiffness of particles contacting with the flexible wall
    InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();

	for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	{
		const shared_ptr<Interaction>& contact = *ii;
		int id1 = contact->getId1(), id2 = contact->getId2();
		int id = (id1 < id2) ? id1 : id2;
		if (id<2){//caution: I assume that wall_id is less than 2. top and bottom walls
            FrictPhys* currentContactPhysics =
	        static_cast<FrictPhys*> ( contact->phys.get() );

	        avg_zstiff  += currentContactPhysics->kn;
	    }else if((keepRigid)&&(id < numWallFacet)){//TriElements
            FrictPhys* currentContactPhysics =
	        static_cast<FrictPhys*> ( contact->phys.get() );
            avg_rstiff  += currentContactPhysics->kn;
        }

	}
	//calc gains
    if(isFlexible){
        gain_r = 2.0*gain_alpha*height*M_PI*2.0*wall_radius/(avg_rstiff);
    }else{
        gain_r = 2.0*gain_alpha*height*M_PI*2.0*wall_radius/(num_pwall*1e8);
    }
   //std::cerr<<"gain x"<<gain_r<<" num pwall"<<num_pwall<<std::endl;
   gain_z = 2.0*gain_alpha*loading_area/(avg_zstiff);
}



void FlexCompressionEngine::servo_cylinder(double dt){
    //if (!consol_on){return;}//if the flag consol_on is false, which means not consolidating
    consol_ss_cylinder();
    Real udr, udz;//velocities of walls in the three directions
    if(!isFlexible){
        gain_r = 2.0*gain_alpha*height*M_PI*2.0*wall_radius/(num_pwall*1e8);
    }
    udr = gain_r* (cylinder_stress - goalr);
    //std::cout<<"cylinder_stress: "<<cyliner_stress<<"! gain_r "<<gain_r<<"! udr"<<std::endl;
    wall_radius += udr;
    if (z_servo){//stress control in z direction
        udz = gain_z * (wszz - goalz);
        bottom_wall->state->pos[2] -= udz;//z
        top_wall->state->pos[2] += udz;//z
    }//strain control during shear
    if(isFlexible&&keepRigid){//update nodes' positions
        SUDODEM_PARALLEL_FOREACH_NODE_BEGIN(const shared_ptr<Node>& b, scene->nodes){
            Vector2r p0 = Vector2r(b->pos[0],b->pos[1]);//x,y of the node's pos
            p0 = p0*(1.0+udr/p0.norm());
            b->pos[0]=p0[0];b->pos[1]=p0[1];
        }SUDODEM_PARALLEL_FOREACH_NODE_END();
        //top and bottom nodes
        if(z_servo){
            for(vector<shared_ptr<Node>>::iterator b = top_nodes.begin();b!= top_nodes.end();++b){
                (*b)->pos[2] += udz;//z 
            }
            for(vector<shared_ptr<Node>>::iterator b = bottom_nodes.begin();b!= bottom_nodes.end();++b){
                (*b)->pos[2] -= udz;//z 
            }
        }
    }

}


void FlexCompressionEngine::consol_ss_cylinder(){
    // sync thread storage of ForceContainer
	 scene->forces.sync();

	//calculate stress
	//get box size
	Real top = top_wall->state->pos.z();
	Real bottom = bottom_wall->state->pos.z();

	height = top - bottom;
    //if(!isFlexible){x_area = 2.0*M_PI*wall_radius * height;}
    x_area = 2.0*M_PI*wall_radius * height;
	loading_area = M_PI*pow(wall_radius,2.0);
    if(debug){
        x_area = 0.0;
        for(vector<shared_ptr<TriElement>>::iterator elm=Elements.begin();elm!=Elements.end();++elm){
            x_area += (*elm)->getArea();
        }
    }
	//wall forces
    if(isFlexible){
        if(keepRigid || debug){
                cylinder_force = 0.0;
                for(int i=2;i<numWallFacet;i++){
                    cylinder_force += getForce(scene,i).norm();//no z-component due to a frictionless wall
                }
                cylinder_stress = (cylinder_force)/x_area;
                /*
                InteractionContainer::iterator ii    = scene->interactions->begin();
	            InteractionContainer::iterator iiEnd = scene->interactions->end();

	            for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	            {
		            const shared_ptr<Interaction>& contact = *ii;
		            int id1 = contact->getId1(), id2 = contact->getId2();
		            int id = (id1 < id2) ? id1 : id2;
		            if (id<2){//caution: I assume that wall_id is less than 2. top and bottom walls
                        FrictPhys* currentContactPhysics =
	                    static_cast<FrictPhys*> ( contact->phys.get() );

	                    avg_zstiff  += currentContactPhysics->kn;
	                }else if((keepRigid)&&(id < numWallFacet)){//TriElements
                        FrictPhys* currentContactPhysics =
	                    static_cast<FrictPhys*> ( contact->phys.get() );
                        avg_rstiff  += currentContactPhysics->kn;
                    }

	            }
                */
        }
    }else{
        cylinder_stress = (cylinder_force)/x_area;
    }
    wszz = 0.5*(getForce(scene,1)[2] - getForce(scene,0)[2])/loading_area;//top_wall id1 - bottom_wall id0
    //std::cerr<<"loadingarea"<<loading_area<<" wszz "<<wszz<<std::endl;
}
void FlexCompressionEngine::checkTarget(){
    consol_ss_cylinder();
    Flag_ForceTarget = false;
    if(!isFlexible){
       if ((std::abs(cylinder_stress - goalr)/goalr) < f_threshold){

            if ((std::abs(wszz - goalz)/goalz) < f_threshold){
                Flag_ForceTarget = true;//forces acting on walls reached the target
            }
        }
    }else{
        if(keepRigid){
            if ((std::abs(cylinder_stress - goalr)/goalr) < f_threshold){

                if ((std::abs(wszz - goalz)/goalz) < f_threshold){
                    Flag_ForceTarget = true;//forces acting on walls reached the target
                }
            }
        }else{
            if ((std::abs(wszz - goalz)/goalz) < f_threshold){Flag_ForceTarget = true;}//forces acting on walls reached the target
        }
    }
}
void FlexCompressionEngine::quiet_system(){
    BodyContainer::iterator bb    = scene->bodies->begin();
	BodyContainer::iterator iiEnd = scene->bodies->end();

	for(  ; bb!=iiEnd ; ++bb )
	{
        State* state=(*bb)->state.get();
        state->vel = Vector3r(0.0,0.0,0.0);
        state->angVel = Vector3r(0.0,0.0,0.0);
    }
}


Real FlexCompressionEngine::getResultantF(){
	scene->forces.sync();
	Real f=0.;
	/*
	for(int i=0;i<downBox.size();i++){
		//f += getForce(scene,downBox[i])[0];//force along x,f=scene->forces.getForce(id);
		f += scene->forces.getForce(downBox[i])[0];
	}
	*/
	return f;
}



