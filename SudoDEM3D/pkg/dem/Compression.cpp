/*************************************************************************
*  Copyright (C) 2016 by Zhswee		        		 	        		 *
*  zhswee@gmail.com      					  							 *
*																		 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Compression.hpp"

#include<sudodem/pkg/common/Box.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/Superquadrics.hpp>
#include<sudodem/pkg/dem/PolySuperellipsoid.hpp>
//#include<sudodem/pkg/dem/GJKParticle.hpp>
//create a new class base,'NonSphericalShape', in the furture.
#ifdef FM
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>
#endif

#include<sudodem/core/State.hpp>
#include<sudodem/core/Clump.hpp>
#include<assert.h>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/dem/Shop.hpp>

#include<sudodem/core/Omega.hpp>
#include<sudodem/lib/base/Math.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/lambda/lambda.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>

//typedef Eigen::Matrix<double, 3, 2> Matrix32d;
//typedef Eigen::Quaternion<double> Quaternionr;
//using namespace Eigen;

#include <sys/syscall.h>//get pid of a thread

CREATE_LOGGER(CompressionEngine);
SUDODEM_PLUGIN((CompressionEngine));

CompressionEngine::~CompressionEngine(){}

void CompressionEngine::updateStiffness ()
{
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();

        for(int j = 0; j < 6; j++) {wall_stiffness[j] = 0.;}
	for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	{
		const shared_ptr<Interaction>& contact = *ii;
		//Real fn = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->normalForce.norm();
		//if (fn!=0)
		//{
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;
			if (id<6){//caution: I assume that wall_id is less than 6.
			        SuperquadricsPhys* currentContactPhysics =
				static_cast<SuperquadricsPhys*> ( contact->phys.get() );
		        switch(id){//
		                case 0:{
			        //left_wall
				wall_stiffness[0]  += currentContactPhysics->kn;
				break;}
			        case 1:{
			        //right_wall
				wall_stiffness[1]  += currentContactPhysics->kn;
				break;}
				case 2:{
			        //front_wall
				wall_stiffness[2]  += currentContactPhysics->kn;
				break;}
				case 3:{
			        //back_wall
				wall_stiffness[3]  += currentContactPhysics->kn;
				break;}
				case 4:{
			        //bottom_wall
				wall_stiffness[4]  += currentContactPhysics->kn;
				break;}
				case 5:{
			        //top_wall
				wall_stiffness[5]  += currentContactPhysics->kn;
				break;}

		        }
		       // }
		}
	}

}


Real CompressionEngine::ComputeUnbalancedForce( bool maxUnbalanced) {return Shop::unbalancedForce(maxUnbalanced,scene);}

void CompressionEngine::updateBoxSize(){//computing the box size
  Real left = box[left_wall]->state->pos.x();
	Real right = box[right_wall]->state->pos.x();
	Real front = box[front_wall]->state->pos.y();
	Real back = box[back_wall]->state->pos.y();
	Real top = box[top_wall]->state->pos.z();
	Real bottom = box[bottom_wall]->state->pos.z();

	width  = right - left;
	depth  = back - front;
	height = top - bottom;
	x_area = depth * height;
	y_area = width * height;
	z_area = depth * width;
}
// void CompressionEngine::action_old()
// {
//
// 	if ( firstRun )
// 	{
// 		LOG_INFO ( "First run, will initialize!" );
//
//
// 		getBox();
//
// 		Real left = box[left_wall]->state->pos.x();
// 		Real right = box[right_wall]->state->pos.x();
// 		Real front = box[front_wall]->state->pos.y();
// 		Real back = box[back_wall]->state->pos.y();
// 		Real top = box[top_wall]->state->pos.z();
// 		Real bottom = box[bottom_wall]->state->pos.z();
//
// 		width = width0 = right - left;
// 		depth = depth0 = back - front;
// 		left_facet_pos = left;
// 		height = height0 = top - bottom;
// 		x_area = depth * height;
// 	        y_area = width * height;
// 	        z_area = depth * width;
// 		//calculate particlesVolume
// 		BodyContainer::iterator bi = scene->bodies->begin();
// 		BodyContainer::iterator biEnd = scene->bodies->end();
// 		particlesVolume = 0;
// 		for ( ; bi!=biEnd; ++bi )
// 		{
// 			//if((*bi)->isClump()) continue;//no clump for Superquadrics
// 			const shared_ptr<Body>& b = *bi;
// 			//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";
//
// 			if (b->shape->getClassName()=="Superquadrics"){
// 				const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
// 				particlesVolume += A->getVolume();
//
// 			}
//
// 		}
//       rampNum = 0;
// 		firstRun=false;
// 	}
//
// 	const Real& dt = scene->dt;
// 	rampNum += 1;
// 	//stress
// 	if (scene->iter % stiffnessUpdateInterval == 0 || scene->iter<100) updateStiffness();
// 	// sync thread storage of ForceContainer
// 	 scene->forces.sync();
//
// 	//calculate stress
//
// 	//UnbalancedForce = ComputeUnbalancedForce ();
// 	//loadingStress();//add force on the top facets
// 	//get box size
// 	updateBoxSize();
// 	//wall forces
//         left_wall_s = getForce(scene,left_wall)[0]/x_area;
//         right_wall_s = getForce(scene,right_wall)[0]/x_area;
//         front_wall_s = getForce(scene,front_wall)[1]/y_area;
//         back_wall_s = getForce(scene,back_wall)[1]/y_area;
//         bottom_wall_s = getForce(scene,bottom_wall)[2]/z_area;
//         top_wall_s = getForce(scene,top_wall)[2]/z_area;
// 	if (hydroStrain)//compaction with a constant strain rate in all directions
// 	{
// 	        if ((left_wall_s > goalx) || (right_wall_s > goalx) || (front_wall_s > goaly) || (back_wall_s > goaly) || (bottom_wall_s > goalz) || (top_wall_s > goalz)){//if any wall reaches a desired force, then stop all walls
// 	        box[left_wall]->state->vel[0] = 0;
//                 //right wall
// 	        box[right_wall]->state->vel[0] = 0;
// 	        //front wall
// 	        box[front_wall]->state->vel[1] = 0;
// 	        //back wall
// 	        box[back_wall]->state->vel[1] = 0;
// 	        //bottom wall
// 	        box[bottom_wall]->state->vel[2] = 0;
// 	        //top wall
// 	        box[top_wall]->state->vel[2] = 0;
// 	        }else{
// 	        //left wall
//                 box[left_wall]->state->vel[0] = hydroStrainRate*width;
//                 //right wall
// 	        box[right_wall]->state->vel[0] = -hydroStrainRate*width;
// 	        //front wall
// 	        box[front_wall]->state->vel[1] = hydroStrainRate*depth;
// 	        //back wall
// 	        box[back_wall]->state->vel[1] = -hydroStrainRate*depth;
// 	        //bottom wall
// 	        box[bottom_wall]->state->vel[2] = hydroStrainRate*height;
// 	        //top wall
// 	        box[top_wall]->state->vel[2] = -hydroStrainRate*height;
// 	        }
//
// 	}else{
//
// 	        //moving walls
// 	         double ramp_alpha = float(rampNum)/float(ramp_interval);
// 	        //left wall
// 	        if (left_wall_activated) {
// 		        if (stressMask & 1) move_wall(left_wall,-1,0,x_area,goalx,left_wall_s);	//x
// 		        else box[left_wall]->state->vel[0] = 0.5*goalx*width;
// 	        } else box[left_wall]->state->vel=Vector3r::Zero();
//
// 	        //right wall
// 	        if (right_wall_activated) {
// 		        if (stressMask & 1)  move_wall(right_wall,1,0,x_area,goalx,right_wall_s);	//x
// 		        else box[right_wall]->state->vel[0] = -0.5*goalx*width;
// 	        } else box[right_wall]->state->vel=Vector3r::Zero();
//
// 	        //front wall
// 	        if (front_wall_activated) {
// 		        if (stressMask & 2) move_wall(front_wall,-1,1,y_area,goaly,front_wall_s);	//y
// 		        else box[front_wall]->state->vel[1] = 0.5*goaly*depth;
// 	        } else box[front_wall]->state->vel=Vector3r::Zero();
//
// 	        //back wall
// 	        if (back_wall_activated) {
// 		        if (stressMask & 2)  move_wall(back_wall,1,1,y_area,goaly,back_wall_s);	//y
// 		        else box[back_wall]->state->vel[1] = -0.5*goaly*depth;
// 	        } else box[back_wall]->state->vel=Vector3r::Zero();
//
// 	        //bottom wall
// 	        if (bottom_wall_activated) {
// 		        if (stressMask & 4)  move_wall(bottom_wall,-1,2,z_area,goalz,bottom_wall_s);	//z
// 		        else box[bottom_wall]->state->vel[2] = 0.5*goalz*height*ramp_alpha;
// 	        } else box[bottom_wall]->state->vel=Vector3r::Zero();
//
// 	        //top wall
// 	        if (top_wall_activated) {
// 		        if (stressMask & 4)  move_wall(top_wall,1,2,z_area,goalz,top_wall_s);	//z
// 		        else box[top_wall]->state->vel[2] = -0.5*goalz*height*ramp_alpha;
// 	        } else box[top_wall]->state->vel=Vector3r::Zero();
// 	}
// 	//output info
// 		//calculate stress
// 	//stress[wall_top] = getForce(scene,wall_id[wall_top])/(width*depth);
// 	if (scene->iter % savedata_interval == 0){recordData();}
// 	if (scene->iter % echo_interval == 0){
// 		UnbalancedForce = ComputeUnbalancedForce ();
// 		//getStressStrain();
// 		std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
// 		<<", Ss_x:"<<(right_wall_s - left_wall_s)/2000.0 //average stress of the two walls, unit is kPa
// 		<<", Ss_y:"<<(back_wall_s - front_wall_s)/2000.0
// 		<<", Ss_z:"<<(top_wall_s - bottom_wall_s)/2000.0
//       <<", e:" << ((width*depth*height)/particlesVolume - 1.0)
// 		<<"\n";
// 	}
//    if (rampNum >= ramp_interval){rampNum = 0;}
// }
void CompressionEngine::hydroConsolidation(){
		const Real& dt = scene->dt;
	 	scene->forces.sync();
		//get box size
		updateBoxSize();
		//wall forces
    left_wall_s = getForce(scene,left_wall)[0]/x_area;
    right_wall_s = getForce(scene,right_wall)[0]/x_area;
    front_wall_s = getForce(scene,front_wall)[1]/y_area;
    back_wall_s = getForce(scene,back_wall)[1]/y_area;
    bottom_wall_s = getForce(scene,bottom_wall)[2]/z_area;
    top_wall_s = getForce(scene,top_wall)[2]/z_area;
		Real delta_left, delta_right,delta_front,delta_back,delta_bottom,delta_top;
		delta_left = goalx - fabs(left_wall_s);
		delta_right = goalx - fabs(right_wall_s);
		delta_front = goaly - fabs(front_wall_s);
		delta_back = goaly - fabs(back_wall_s);
		delta_bottom = goalz - fabs(bottom_wall_s);
		delta_top = goalz - fabs(top_wall_s);
		if (delta_left>0){
			box[left_wall]->state->pos[0] += min(0.1*gain_x*delta_left,hydroStrainRate*width*dt);
		}//hydroStrainRate should be small enough
		if (delta_right>0){
			box[right_wall]->state->pos[0] -= min(0.1*gain_x1*delta_right,hydroStrainRate*width*dt);
		}
		if (delta_front>0){
			box[front_wall]->state->pos[1] += min(0.1*gain_y*delta_front,hydroStrainRate*depth*dt);
		}
		if (delta_back>0){
			box[back_wall]->state->pos[1] -= min(0.1*gain_y1*delta_back,hydroStrainRate*depth*dt);
		}
		if (delta_bottom>0){
			box[bottom_wall]->state->pos[2] += min(0.1*gain_z*delta_bottom,hydroStrainRate*height*dt);
		}
		if (delta_top>0){
			box[top_wall]->state->pos[2] -= min(0.1*gain_z1*delta_top,hydroStrainRate*height*dt);
		}
		//if (iterate_num < 1){quiet_system();}
		iterate_num += 1;
		if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
				iterate_num = 1;
				get_gain();
				//check stress
		}
		solve_num += 1;
		//do nothing!
		if (solve_num >= solve_num_max){
			solve_num = 1;
			Vector2r stress = getStress();
			//here we use goalx cause we regard the stress state is isotropic with the same ones along the three directions
			if ((std::abs(stress(0) - goalx)/goalx < f_threshold) && (stress(1)/stress(0) < fmin_threshold)){//we reach the target
				if (UnbalancedForce < unbf_tol){//o.5 by default
							 //stop running
							 scene->stopAtIter = scene->iter + 1;
							 std::cerr<<"consolidation completed!"<<std::endl;
			  }
	   }
	}
}

void CompressionEngine::generalConsolidation(){
	const Real& dt = scene->dt;
	//all walls remain constant stresses
	#ifdef NEW_SCHEME//new scheme
		//keep servo
		//cout<<"solve_num"<<solve_num<<"solmx="<<solve_num_max<<endl;
		//if (solve_num < 1){iterate_num = 0;}//frenquency for checking balance
		iterate_num += 1;
		solve_num += 1;
		servo(dt);
		if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles by default
				iterate_num = 0;
				get_gain();
		}

		if (solve_num >= solve_num_max){//check the stability of the sample
				solve_num = 0;
				Vector2r stress = getStress();
				//std::cerr<<"confining stress ratio= "<<stress(0)/goalx<<" deviatorStress ratio= "<<stress(1)/stress(0)<<endl;
				//here we use goalx cause we regard the stress state is isotropic with the same ones along the three directions
				if ((std::abs(stress(0) - goalx) < goalx*f_threshold) && (stress(1)/stress(0) < fmin_threshold)){//we reach the target
					UnbalancedForce = ComputeUnbalancedForce ();
					if (UnbalancedForce < unbf_tol){//o.5 by default
								 //stop running
								 scene->stopAtIter = scene->iter + 1;
								 if(newton){newton->quiet_system_flag=true;}
								 std::cerr<<"consolidation completed!"<<std::endl;
					}
			 }
		}
	#else
		if (!Flag_ForceTarget){//forces on walls haven't reached the target, so we need to switch on servo.
				if (iterate_num < 1){solve_num = 0;/*if(newton){newton->quiet_system_flag=true;}*//*quiet_system();*/}
				iterate_num += 1;
				servo(dt);
				if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
						iterate_num = 1;
						get_gain();
						checkTarget();//check the forces on the walls reach the target or not, and info is stored in Flag_ForceTarget.
						//if(newton){newton->quiet_system_flag=true;cout<<"quiet system"<<endl;}
				}
		}else{//keep walls fixed, then cyle to a stable state
				if (solve_num < 1){/*if(newton){newton->quiet_system_flag=true;}*//*quiet_system();*/iterate_num = 0;}//set velocities of all walls and particles to zero.
				solve_num += 1;
				//do nothing!
				if (solve_num >= solve_num_max){//check the stability of the sample
						UnbalancedForce = ComputeUnbalancedForce ();
						solve_num = 1;
						consol_ss();
						if ((left_wall_activated && ((std::abs(wsxx - goalx)/goalx) > fmin_threshold))||
								(front_wall_activated && ((std::abs(wsyy - goaly)/goaly) > fmin_threshold))||
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
	#endif
}

void CompressionEngine::action_cubic(){

		if ( firstRun )
		{
			LOG_INFO ( "First run, will initialize!" );


			getBox();
			FOREACH(shared_ptr<Engine>& e, scene->engines){ newton=SUDODEM_PTR_DYN_CAST<NewtonIntegrator>(e);if(newton) {break;} }
			Real left = box[left_wall]->state->pos.x();
			Real right = box[right_wall]->state->pos.x();
			Real front = box[front_wall]->state->pos.y();
			Real back = box[back_wall]->state->pos.y();
			Real top = box[top_wall]->state->pos.z();
			Real bottom = box[bottom_wall]->state->pos.z();
	    //continue running program
	    //LOG_WARN("Please set continueFlag to True for continuing consolidating or shearing after loading data.");
			if (continueFlag){//cintinue consolidating or shearing after loading data
	        width  = right - left;
			    depth  = back - front;
			    height = top - bottom;
	        }else{//new running
	            width = width0 = right - left;
			    depth = depth0 = back - front;
			    height = height0 = top - bottom;
	        }
					x_area = depth * height;
	        y_area = width * height;
	        z_area = depth * width;
	        Init_boxVolume = height0*width0*depth0;
			//calculate particlesVolume
			BodyContainer::iterator bi = scene->bodies->begin();
			BodyContainer::iterator biEnd = scene->bodies->end();
			particlesVolume = 0;
			for ( ; bi!=biEnd; ++bi )
			{
				//if((*bi)->isClump()) continue;//no clump for Superquadrics
				const shared_ptr<Body>& b = *bi;
				//std::cerr << "watch point at particle volume calc " << "\n";
				//#ifdef FM
				if(b->isClump()){
				    //std::cerr << "density = "<<b->material->density << "\n";
	                const shared_ptr<Clump>& clump=SUDODEM_PTR_CAST<Clump>(b->shape);
				    if(clump->members.size()>0){
	                    //MemberMap::iterator I=clump->members.begin();
			            //shared_ptr<Body> subBody=Body::byId(I->first);
				        shared_ptr<Body> subBody=Body::byId(clump->members.begin()->first);//use the material of the first member
				        //const shared_ptr<FrictMat>& BodyMat = SUDODEM_PTR_CAST<FrictMat>(b->material);
				        particlesVolume += b->state->mass/subBody->material->density;//this should be a general method to access a particle volume.
				    }
				    continue;
				}else if(b->isClumpMember()){
				    continue;
				}
				//std::cerr << "watch point at particle volume calc--end " << "\n";
	     if (b->shape->getClassName()=="Sphere"){
					const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b->shape);
					particlesVolume += 4.0/3.0*M_PI*pow(A->radius,3.0);

				}
	            //#else
				if (b->shape->getClassName()=="Superquadrics"){
					const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
					particlesVolume += A->getVolume();

				}
				if (b->shape->getClassName()=="PolySuperellipsoid"){
					const shared_ptr<PolySuperellipsoid>& A = SUDODEM_PTR_CAST<PolySuperellipsoid> (b->shape);
					particlesVolume += A->getVolume();

				}
	            /*
				if (b->shape->getClassName()=="GJKParticle"){
					const shared_ptr<GJKParticle>& A = SUDODEM_PTR_CAST<GJKParticle> (b->shape);
					particlesVolume += A->getVolume();

				}*/
				//#endif
			}
    	rampNum = 0;
			firstRun=false;
      //
      iterate_num = 0;
      solve_num = 0;
      get_gain();
      Flag_ForceTarget = false;

		}

		const Real& dt = scene->dt;
	  if (z_servo){//consolidation

			if (hydroStrain)//compaction with a constant strain rate in all directions
			{
				hydroConsolidation();
			}else{
				generalConsolidation();
			}

	        //output info
		    if (scene->iter % savedata_interval == 0){recordData();}
		    if (scene->iter % echo_interval == 0){
			    UnbalancedForce = ComputeUnbalancedForce ();
	            consol_ss();
			    //getStressStrain();
			    std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
			    <<", Ss_x:"<<wsxx/1000.0 //average stress of the two walls, unit is kPa
			    <<", Ss_y:"<<wsyy/1000.0
			    <<", Ss_z:"<<wszz/1000.0
	          <<", e:" << ((width*depth*height)/particlesVolume - 1.0)
			    <<"\n";
		    }
		}else{//shear
	        if ((!continueFlag) && (rampNum < ramp_interval)){//accelerating the walls
	            //Caution:continuing shearing after shear step excceeds ramp_interval is allowed.
	            if (rampNum < 1){// shear following consolidation
	                //reset the initial data
	                if(Flag_ForceTarget){//consolidation completed normally
	                    std::cerr<<"Shear begins..."<<std::endl;
	                    updateBoxSize();
	                    width0 = width;
	                    height0 = height;
	                    depth0 = depth;
	                }
	            }
	            rampNum ++;

	            double ramp_inter = 0.5*ramp_interval;
	            //double vel= float((ramp_inter - std::abs(rampNum-ramp_inter))/(ramp_inter/ramp_chunks)+1)/float(ramp_chunks)*goalz;//velocity
	            //double vel= float(rampNum/(ramp_interval/ramp_chunks)+1)/float(ramp_chunks)*goalz;//velocity
	            //std::cerr<<"vel of loading walls = "<<vel<<std::endl;
	            double vel=0.0;
	            if (rampNum < ramp_inter){
	                vel =float(rampNum/(ramp_inter/ramp_chunks)+1)/float(ramp_chunks)*goalz*2.0;
	            }else{
	                vel =(float((ramp_interval - rampNum)/(ramp_inter/ramp_chunks)+1)/float(ramp_chunks)+1.0)*goalz;
	            }
	            box[bottom_wall]->state->pos[2] += vel*dt*height;//z,remain a constant loading strain rate
	            box[top_wall]->state->pos[2] += -vel*dt*height;
	            /*
	            double vel= float(rampNum/(ramp_interval/ramp_chunks)+1)/float(ramp_chunks)*goalz;//velocity
	            box[bottom_wall]->state->vel[2] = vel;//z
	            box[top_wall]->state->vel[2] = -vel;//z
	            */
	        }else{
	            //confining stress control
	            if(shearMode & 1){
		            iterate_num += 1;
		            if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
		            iterate_num = 1;
		            get_gain();

		            }
	            }
	            box[bottom_wall]->state->pos[2] += goalz*dt*height;//z,remain a constant loading strain rate
	            box[top_wall]->state->pos[2] += -goalz*dt*height;//z
	        }
					//control the positions of walls
					if(shearMode & 1){//conventional triaxial compression (drained)
						servo(dt);//conventional triaxial shear with constant confining stress applied
					}else{
						//constant volume, modeling undrained triaxial test
	          Real left = box[left_wall]->state->pos.x();
		        Real right = box[right_wall]->state->pos.x();
		        Real front = box[front_wall]->state->pos.y();
		        Real back = box[back_wall]->state->pos.y();
		        Real top = box[top_wall]->state->pos.z();
		        Real bottom = box[bottom_wall]->state->pos.z();

		        width  = right - left;
		        depth  = back - front;
		        height = top - bottom;
	            double S = Init_boxVolume/height;
	            double delta_d,delta_w;
	            delta_d = 0.5*(S/width - depth);
	            delta_w = S/(depth+delta_d) - width;

	            box[left_wall]->state->pos[0] += -0.5*delta_w;//x
	            box[right_wall]->state->pos[0] += 0.5*delta_w;//x
	            box[front_wall]->state->pos[1] += -0.5*delta_d;//y
	            box[back_wall]->state->pos[1] += 0.5*delta_d;//y
	        }

	        //height = box[top_wall]->state->pos.z() - box[bottom_wall]->state->pos.z();

	        if (log(height0/height) > target_strain){//
	            //stop shear
	            std::cerr<<"Shear ends!"<<std::endl;
	            scene->stopAtIter = scene->iter + 1;
	        }
	        	//output info
		    if (scene->iter % savedata_interval == 0){recordData();}
		    if (scene->iter % echo_interval == 0){
			    UnbalancedForce = ComputeUnbalancedForce ();
	            //consol_ss();
			    //getStressStrain();
					if(shearMode & 2){//constant volume
						consol_ss();
					}

						std::cerr<<"Iter "<<scene->iter<<" Ubf "<<UnbalancedForce
						<<", conf:"<<(wsxx+wsyy)/2000.0 //average stress of the two walls, unit is kPa
						<<", Str_z:"<<log(height0/height)//z
						<<", Ss_z:"<<wszz/1000.0
							<<", e:" << ((width*depth*height)/particlesVolume - 1.0)
						<<"\n";
		    }
	    }

}

void CompressionEngine::action_cylindrical(){
	cylinwall_go();
 if (wall_fix){return;}
 if ( firstRun )
 {
	 LOG_INFO ( "First run, will initialize!" );


	 //getBox();

			 box[4] = Body::byId(0);//bottom wall
	 box[5] = Body::byId(1);//top wall
	 Real top = box[top_wall]->state->pos.z();
	 Real bottom = box[bottom_wall]->state->pos.z();


	 height = height0 = top - bottom;
			 wall_radius0 =wall_radius;
		 loading_area = M_PI*pow(wall_radius,2.0);
			 Init_boxVolume = height0*loading_area;
	 //calculate particlesVolume
	 BodyContainer::iterator bi = scene->bodies->begin();
	 BodyContainer::iterator biEnd = scene->bodies->end();
	 particlesVolume = 0;
	 for ( ; bi!=biEnd; ++bi )
	 {
		 //if((*bi)->isClump()) continue;//no clump for Superquadrics
		 const shared_ptr<Body>& b = *bi;
		 //std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		 if (b->shape->getClassName()=="Superquadrics"){
			 const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
			 particlesVolume += A->getVolume();

		 }
		 if (b->shape->getClassName()=="PolySuperellipsoid"){
			 const shared_ptr<PolySuperellipsoid>& A = SUDODEM_PTR_CAST<PolySuperellipsoid> (b->shape);
			 particlesVolume += A->getVolume();

		 }

	 }
		 rampNum = 0;
	 firstRun=false;
			 //
			 iterate_num = 0;
			 solve_num = 0;
			 get_gainz();
			 Flag_ForceTarget = false;

 }


	 const Real& dt = scene->dt;
	 if (z_servo){//consolidation
			 if (!Flag_ForceTarget){//forces on walls haven't reached the target, so we need to switch on servo.
					 if (iterate_num < 1){solve_num = 0;if(newton){newton->quiet_system_flag=true;}/*quiet_system();*/}
					 iterate_num += 1;
					 servo_cylinder(dt);
					 if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
							 iterate_num = 1;
							 get_gainz();
							 checkTarget();//check the forces on the walls reach the target or not, and info is stored in Flag_ForceTarget.
					 }
			 }else{//keep walls fixed, then cyle to a stable state
					 if (solve_num < 1){if(newton){newton->quiet_system_flag=true;}/*quiet_system();*/iterate_num = 0;}//set velocities of all walls and particles to zero.
					 solve_num += 1;
					 //do nothing!
					 if (solve_num >= 1000){//check the stability of the sample
							 UnbalancedForce = ComputeUnbalancedForce ();
							 solve_num = 1;
							 consol_ss_cylinder();
							 if (((std::abs(cylinder_stress - goalx)/goalx) > fmin_threshold)||
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
 }else{//shear
			 if (rampNum < ramp_interval){//accelerating the walls
					 if (rampNum < 1){// shear following consolidation
							 //reset the initial data
							 if(Flag_ForceTarget){//consolidation completed normally
									 std::cerr<<"Shear begins..."<<std::endl;
									 //updateBoxSize();
									 //width0 = width;
									 //height0 = height;
									 //depth0 = depth;
							 }
					 }
					 rampNum ++;
					 double vel= float(rampNum/(ramp_interval/ramp_chunks)+1)/float(ramp_chunks)*goalz;//velocity
					 //std::cerr<<"vel of loading walls = "<<vel<<std::endl;
					 box[bottom_wall]->state->vel[2] = vel;//z
					 box[top_wall]->state->vel[2] = -vel;//z
			 }
			 iterate_num += 1;
			 servo_cylinder(dt);
			 if (iterate_num >= iterate_num_max){//update the servo gains every 100 cycles
					 iterate_num = 1;
					 get_gainz();
			 }

			 //height = box[top_wall]->state->pos.z() - box[bottom_wall]->state->pos.z();

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
}

void CompressionEngine::action()
{
	if(cubic_specimen){
		action_cubic();
	}else{//cylindrical wall
		action_cylindrical();
	}
}

Matrix3r CompressionEngine::getStressTensor(){
	//scene->forces.sync();
	if(cubic_specimen){
	//get box size
  updateBoxSize();
	Current_boxVolume = height*z_area;//width*depth*height;
  }else{
	//consol_ss_cylinder();
	Current_boxVolume = height*M_PI*pow(wall_radius,2.0);//width*depth*height;
  }
	//double sigma1=0.0;
	//compute the stress within the sample, sigma_ij

	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();

	Matrix3r Sigma = Matrix3r::Zero();
	for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	{	Vector3r fn,fs;
		const shared_ptr<Interaction>& contact = *ii;
			//#ifdef FM
			fn = (static_cast<FrictPhys*>	(contact->phys.get()))->normalForce;
	fs = (static_cast<FrictPhys*>	(contact->phys.get()))->shearForce;
			//#else
	/*if(isHertzMindlin){
		SuperquadricsHertzMindlinPhys* currentContactPhysics =
					static_cast<SuperquadricsHertzMindlinPhys*> ( contact->phys.get() );
		fn = currentContactPhysics->normalForce;
		fs = currentContactPhysics->shearForce;
	}else{*/
		//fn = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->normalForce;
		//fs = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->shearForce;
	//}
			//#endif
		if (fn.norm()==0){continue;}

		int id1 = contact->getId1(), id2 = contact->getId2();
		int id = (id1 < id2) ? id1 : id2;
		Vector3r pos1,pos2,l;
		{//if (id>5){//caution: I assume that wall_id is less than 6.
				const shared_ptr<Body>& b1 = Body::byId(id1,scene);
							const shared_ptr<Body>& b2 = Body::byId(id2,scene);
	//if(!b1->isClump()){//clump does not have interactions with others
			if(b1->isClumpMember()){//clump member
		id1 = b1->clumpId;//asign the clump id
		pos1 = Body::byId(id1,scene)->state->pos;
			}else{//general particle
		pos1 = b1->state->pos;
			}
	//}
	//if(!b2->isClump()){
			if(b2->isClumpMember()){//clump member
		id2 = b2->clumpId;//asign the clump id
		pos2 = Body::byId(id2,scene)->state->pos;
			}else{//general particle
		pos2 = b2->state->pos;
			}
	//}
	if(id1!=id2){//the two clump members if possible are from different clumps
			//l = pos2 - pos1;
			Sigma += (fn+fs)*((pos2 - pos1).transpose());
	}
			}
	}
	Sigma /= Current_boxVolume;
	return	Sigma;
}
Vector2r CompressionEngine::getStress()
{

	Matrix3r Sigma = getStressTensor();
	//std::cout<<"stress"<<Sigma<<std::endl;
	double meanStress = Sigma.trace()/3.0;
	Sigma -= meanStress*Matrix3r::Identity();
	double deviatorStress = Sigma.norm()*sqrt(1.5);

	return Vector2r(meanStress,deviatorStress);

}
void CompressionEngine::recordData()
{
	if(!out.is_open()){
		assert(!out.is_open());

		std::string fileTemp = file;


		if(fileTemp.empty()) throw ios_base::failure(__FILE__ ": Empty filename.");
		out.open(fileTemp.c_str(), truncate ? fstream::trunc : fstream::app);
		if(!out.good()) throw ios_base::failure(__FILE__ ": I/O error opening file `"+fileTemp+"'.");
	}
	/*
	// at the beginning of the file; write column titles
    #ifdef CUBIC_SAMPLE
    consol_ss();
    double volume = height*z_area;//width*depth*height;
    #else
    consol_ss_cylinder();
    double volume = height*M_PI*pow(wall_radius,2.0);//width*depth*height;
    #endif
    double sigma1=0.0;

    if (1){//compute the stress within the sample, sigma_ij

        InteractionContainer::iterator ii    = scene->interactions->begin();
	    InteractionContainer::iterator iiEnd = scene->interactions->end();

        Matrix3d Sigma = Matrix3d::Zero();
	    for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	    {
		    const shared_ptr<Interaction>& contact = *ii;
            #ifdef FM
            Vector3r fn = (static_cast<FrictPhys*>	(contact->phys.get()))->normalForce;
            #else
		    Vector3r fn = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->normalForce;
            #endif
		    if (fn.norm()==0){continue;}

		    int id1 = contact->getId1(), id2 = contact->getId2();
		    int id = (id1 < id2) ? id1 : id2;
		    if (id>6){//caution: I assume that wall_id is less than 6.
		            const shared_ptr<Body>& b1 = Body::byId(id1,scene);
                    const shared_ptr<Body>& b2 = Body::byId(id2,scene);
                    State* s1 = b1->state.get();
                    State* s2 = b2->state.get();
                    Vector3r l = s2->pos - s1->pos;
                    Sigma += fn*(l.transpose());
            }
        }
        EigenSolver<Matrix3d> es;
        es.compute(Sigma/volume,false);
        std::complex<double> E1, E2,E3;
        E1 = es.eigenvalues().col(0)[0];
        E2 = es.eigenvalues().col(0)[1];
        E3 = es.eigenvalues().col(0)[2];
        sigma1 = E1.real()>E2.real()?E1.real():E2.real();
        sigma1 = sigma1 > E3.real()?sigma1:E3.real();

    }
	*/
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
 	//<< boost::lexical_cast<string> ( translationAxisz[0] ) << " "
	//<< boost::lexical_cast<string> ( translationAxisz[1] ) << " "
 	//<< boost::lexical_cast<string> ( translationAxisz[2] ) << " "
	//<< boost::lexical_cast<string> ( syscall(SYS_gettid) ) << " "
 	<< endl;
}
void CompressionEngine::cylinwall_go(){//cylindrical wall
    //traverse all particles
    BodyContainer::iterator bb    = scene->bodies->begin();
	BodyContainer::iterator iiEnd = scene->bodies->end();
    cylinder_force = 0.0;
    num_pwall = 0;
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


}
// void CompressionEngine::move_wall(int wall_id, double Sign, int f_index,double area,double goal_stress, double wall_stress){
//         //wall_id
//         //Sign
//         //f_index:0,1,2 for forces in the x, y and z directions.
//         //wall_pos: an array storing the wall_position
// 	//scene->forces.sync();
//
// 	//Real wall_f = getForce(scene,wall_id)[f_index];
// 	Real translation= wall_stress - goal_stress*Sign;//force in f_index direction
// 	Real threshold = f_threshold * goal_stress;
//    if (std::abs(translation) < std::abs(threshold)){
//       box[wall_id]->state->vel = Vector3r::Zero();
//       return;
//    }
// 	Real trans=0.;
// 	//
// 	Real wpos = box[wall_id]->state->pos[f_index];//
// 	//FIXME:if force on the top wall is equal to zero
// 	//}
// 	//moving the wall with the maximum velocity if the current fouce on the wall is less than 80% of the target
// 	//this case including that the wall is not touching with the particles
// 	//if (std::abs(wall_f) < 0.8*goal*area){//
// 	Vector3r v = Vector3r::Zero();
// 	v[f_index] = 1.0;
// 	if (std::abs(wall_stress) < 1.0){//wall stress is very small, less than 1 kPa
// 		 box[wall_id]->state->vel = wall_max_vel * Mathr::Sign(translation) * v;//Vector3r(0.,0.0,);
// 		 if (debug) std::cerr<<"1"<<"\n";
// 		 //Caution:the wall may move too fast to make the most-top particles escape through the wall
// 	}else{//
//
//         	if (wall_stiffness[wall_id]!=0){
// 			if (debug) std::cerr<<"4"<<"\n";
// 			trans = alpha*area*translation / wall_stiffness[wall_id];
// 			//std::cerr << "translation:" <<translation<<"wall_max_vel:"<<wall_max_vel<< "\n";
//
// 			//translation = std::abs(translation) * Mathr::Sign(translation);
//
// 		}
// 		else {
// 		        trans = wall_max_vel * Mathr::Sign(translation)*scene->dt;
// 		        if (debug) std::cerr<<"5"<<"\n";
// 		}
//
// 		trans = std::min( std::abs(trans), wall_max_vel*scene->dt ) * Mathr::Sign(translation);
// 		box[wall_id]->state->vel = trans/scene->dt *v;//Vector3r(0.,0.,trans/scene->dt);
// 	}
// }
/*
void CompressionEngine::move_wall(int wall_id, double Sign, int f_index,double (&wall_pos)[4],double goal_stress, double wall_stress){
        //wall_id
        //Sign
        //f_index:0,1,2 for forces in the x, y and z directions.
        //wall_pos: an array storing the wall_position
	//scene->forces.sync();

	//Real wall_f = getForce(scene,wall_id)[f_index];
	Real translation= wall_stress - goal_stress*Sign;//force in f_index direction
	Real threshold = f_threshold * goal_stress;
	Real trans=0.;
	//
	Real wpos = box[wall_id]->state->pos[f_index];//
	//FIXME:if force on the top wall is equal to zero
	//}
	//moving the wall with the maximum velocity if the current fouce on the wall is less than 80% of the target
	//this case including that the wall is not touching with the particles
	//if (std::abs(wall_f) < 0.8*goal*area){//
	Vector3r v = Vector3r::Zero();
	v[f_index] = 1.0;
	if (std::abs(wall_stress) < 1.0){//wall stress is very small, less than 1 kPa
		 box[wall_id]->state->vel = wall_max_vel * Mathr::Sign(translation) * v;//Vector3r(0.,0.0,);
		 if (debug) std::cerr<<"1"<<"\n";
		 //Caution:the wall may move too fast to make the most-top particles escape through the wall
	}else{//
	        //correct position of wall
	        if (scene->iter % stiffnessUpdateInterval == 0){
	        	if (wall_stiffness[wall_id]!=0){
				if (debug) std::cerr<<"4"<<"\n";
				trans = translation / wall_stiffness[wall_id];
				//std::cerr << "translation:" <<translation<<"wall_max_vel:"<<wall_max_vel<< "\n";

				//translation = std::abs(translation) * Mathr::Sign(translation);

			}
			else {
			        trans = wall_max_vel * Mathr::Sign(translation)*scene->dt;
			        if (debug) std::cerr<<"5"<<"\n";
			}
	        }else{
		//numerically find the target velocity of the top wall
		if (translation != wall_pos[3] && wall_pos[3]!=0) {
			/*if (std::abs(translation) > threshold && translation*wall_pos[3] >0 && std::abs(translation) > std::abs(wall_pos[3])){
				trans = 1.1*(wall_pos[1] - wpos);
				if (debug) std::cerr<<"2"<<"\n";
			}
			//else{
				Real x3 = wpos - translation*(wpos-wall_pos[1])/(translation -wall_pos[3]);
				trans = x3 - wpos;
				if (debug) std::cerr<<"3"<<"\n";
			//}
		}
		else{
			if (wall_stiffness[wall_id]!=0){
				if (debug) std::cerr<<"4"<<"\n";
				trans = translation / wall_stiffness[wall_id];
				//std::cerr << "translation:" <<translation<<"wall_max_vel:"<<wall_max_vel<< "\n";

				//translation = std::abs(translation) * Mathr::Sign(translation);

			}
			else {
			  trans = wall_max_vel * Mathr::Sign(translation)*scene->dt;
			  if (debug) std::cerr<<"5"<<"\n";
			}

		}
		}
		trans = std::min( std::abs(trans), wall_max_vel*scene->dt ) * Mathr::Sign(translation);
		box[wall_id]->state->vel = trans/scene->dt *v;//Vector3r(0.,0.,trans/scene->dt);

		wall_pos[0] = wall_pos[1];
		wall_pos[1] = wpos;
		wall_pos[2] = wall_pos[3];
		wall_pos[3] = translation ;
	}
}
*/
void CompressionEngine::get_gain(){//determine servo gain parameters

   avg_xstiff = 0.0;
   avg_ystiff = 0.0;
   avg_zstiff = 0.0;
    #ifdef SOLE_GAIN
    avg_xstiff1 = 0.0;
    avg_ystiff1 = 0.0;
    avg_zstiff1 = 0.0;
    #endif
   InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();

	for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	{
		const shared_ptr<Interaction>& contact = *ii;
		//Real fn = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->normalForce.norm();
		//if (fn!=0)
		//{
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;

			if (id<6){//caution: I assume that wall_id is less than 6.
					Real kn;
            //#ifdef FM
          FrictPhys* currentContactPhysics = static_cast<FrictPhys*> ( contact->phys.get() );
					kn = currentContactPhysics->kn;
                    //#else
                    /*if(isHertzMindlin){
						SuperquadricsHertzMindlinPhys* currentContactPhysics =
				    static_cast<SuperquadricsHertzMindlinPhys*> ( contact->phys.get() );
						kn = currentContactPhysics->kn;
					}else{*/
			        //SuperquadricsPhys* currentContactPhysics =
				    //static_cast<SuperquadricsPhys*> ( contact->phys.get() );
					//	kn = currentContactPhysics->kn;
					//}
                    //#endif
		        switch(id){//
		                case 0:{
			        //left_wall
				      avg_xstiff  += kn;
				break;}
			        case 1:{
			        //right_wall
                      #ifdef SOLE_GAIN
                      avg_xstiff1 += kn;
                      #else
				      avg_xstiff  += kn;
                      #endif
				break;}
				case 2:{
			        //front_wall
				   avg_ystiff  += kn;
				break;}
				case 3:{
			        //back_wall
                    #ifdef SOLE_GAIN
				    avg_ystiff1  += kn;
                    #else
                    avg_ystiff  += kn;
                    #endif
				break;}
				case 4:{
			        //bottom_wall
				   avg_zstiff  += kn;
				break;}
				case 5:{
			        //top_wall
                    #ifdef SOLE_GAIN
                        avg_zstiff1  += kn;
                    #else
				        avg_zstiff  += kn;
                    #endif
				break;}

		        }
		       // }
		}
	}
	//calc gains

    #ifdef SOLE_GAIN
    gain_x = gain_alpha*x_area/(avg_xstiff);
    gain_y = gain_alpha*y_area/(avg_ystiff);
    gain_z = gain_alpha*z_area/(avg_zstiff);
    gain_x1 = gain_alpha*x_area/(avg_xstiff1);
    gain_y1 = gain_alpha*y_area/(avg_ystiff1);
    gain_z1 = gain_alpha*z_area/(avg_zstiff1);
    #else
    gain_x = 2.0*gain_alpha*x_area/(avg_xstiff);
    gain_y = 2.0*gain_alpha*y_area/(avg_ystiff);
    gain_z = 2.0*gain_alpha*z_area/(avg_zstiff);
    #endif
    //std::cout<<"avg_xstiffness"<<avg_xstiff<<"avg_ystiffness"<<avg_ystiff <<"avg_zstiffness"<<avg_zstiff<<std::endl;
}

void CompressionEngine::get_gainz(){//determine servo gain parameters
   //cylindrical wall
    gain_x = 2.0*gain_alpha*height*M_PI*2.0*wall_radius/(num_pwall*1e8);
    if(!z_servo){return;}
   avg_zstiff = 0.0;
   InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();

	for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	{
		const shared_ptr<Interaction>& contact = *ii;
		//Real fn = (static_cast<SuperquadricsPhys*>	(contact->phys.get()))->normalForce.norm();
		//if (fn!=0)
		//{
			int id1 = contact->getId1(), id2 = contact->getId2();
			int id = (id1 < id2) ? id1 : id2;
			if (id<2){//caution: I assume that wall_id is less than 2. top and bottom walls
				FrictPhys* currentContactPhysics =
				static_cast<FrictPhys*> ( contact->phys.get() );
				avg_zstiff  +=  currentContactPhysics->kn;
		    }
		       // }

	}
	//calc gains

   //std::cerr<<"gain x"<<gain_x<<" num pwall"<<num_pwall<<std::endl;
   gain_z = 2.0*gain_alpha*loading_area/(avg_zstiff);
}

void CompressionEngine::servo(double dt){
    //if (!consol_on){return;}//if the flag consol_on is false, which means not consolidating
    consol_ss();

    #ifdef SOLE_GAIN
    Real udx, udy, udx1, udy1;//velocities of walls in the three directions
		if(left_wall_activated){
	    udx = gain_x * (wsxx_left - goalx);///dt;
		  udx = Mathr::Sign(udx)*min(fabs(udx),max_vel*dt);//constrain the maximum velocity of the wall
	    box[left_wall]->state->pos[0] += -udx;//x
	  }
		if(right_wall_activated){
	    udx1 = gain_x1 * (wsxx_right - goalx);///dt;
		  udx1 = Mathr::Sign(udx1)*min(fabs(udx1),max_vel*dt);//constrain the maximum velocity of the wall
	    box[right_wall]->state->pos[0] += udx1;//x
		}
		if(front_wall_activated){
	    udy = gain_y * (wsyy_front - goaly);///dt;
		  udy = Mathr::Sign(udy)*min(fabs(udy),max_vel*dt);//constrain the maximum velocity of the wall
	    box[front_wall]->state->pos[1] += -udy;//y
	  }
		if(back_wall_activated){
	    udy1 = gain_y1 * (wsyy_back - goaly);///dt;
		  udy1 = Mathr::Sign(udy1)*min(fabs(udy1),max_vel*dt);//constrain the maximum velocity of the wall
	    box[back_wall]->state->pos[1] += udy1;//y
	  }
    #else

    Real udx, udy;//velocities of walls in the three directions
		if(left_wall_activated && right_wall_activated){
	    udx = gain_x * (wsxx - goalx)/dt;
		  udx = Mathr::Sign(udx)*min(fabs(udx),max_vel);//constrain the maximum velocity of the wall
	    box[left_wall]->state->vel[0] = -udx;//x
	    box[right_wall]->state->vel[0] = udx;//x
		}
		if(front_wall_activated && back_wall_activated){
	    udy = gain_y * (wsyy - goaly)/dt;
		  udy = Mathr::Sign(udy)*min(fabs(udy),max_vel);
	    box[front_wall]->state->vel[1] = -udy;//y
	    box[back_wall]->state->vel[1] = udy;//y
	  }
    #endif
    //std::cout<<"servo--moving wall:"<<udx<<" "<<udy<<std::endl;
    if (z_servo){//stress control in z direction
        #ifdef SOLE_GAIN
            Real udz,udz1;
            udz = gain_z * (wszz_bottom - goalz);///dt;
			udz = Mathr::Sign(udz)*min(fabs(udz),max_vel*dt);//constrain the maximum velocity of the wall
            box[bottom_wall]->state->pos[2] += -udz;//z
            udz1 = gain_z1 * (wszz_top - goalz);///dt;
			udz1 = Mathr::Sign(udz1)*min(fabs(udz1),max_vel*dt);//constrain the maximum velocity of the wall
            box[top_wall]->state->pos[2] += udz1;//z
        #else
            Real udz;
            udz = gain_z * (wszz - goalz);///dt;
			udz = Mathr::Sign(udz)*min(fabs(udz),max_vel*dt);
            box[bottom_wall]->state->pos[2] += -udz;//z
            box[top_wall]->state->pos[2] += udz;//z
        #endif
    }//else{//strain control during shear
     //    box[bottom_wall]->state->pos[2] += goalz*dt;//z
     //    box[top_wall]->state->pos[2] += -goalz*dt;//z
    //}
}

void CompressionEngine::servo_cylinder(double dt){
    //if (!consol_on){return;}//if the flag consol_on is false, which means not consolidating
    consol_ss_cylinder();
    Real udr, udz;//velocities of walls in the three directions
    gain_x = 2.0*gain_alpha*height*M_PI*2.0*wall_radius/(num_pwall*1e8);
    udr = gain_x* (cylinder_stress - goalx);
    wall_radius += udr;
    if (z_servo){//stress control in z direction
        udz = gain_z * (wszz - goalz)/dt;
        box[bottom_wall]->state->vel[2] = -udz;//z
        box[top_wall]->state->vel[2] = udz;//z
    }//strain control during shear
}
void CompressionEngine::consol_ss(){
    // sync thread storage of ForceContainer
	 scene->forces.sync();
	//calculate stress
	//wall forces
    #ifdef SOLE_GAIN
		//get box size
		updateBoxSize();
    wsxx_right = getForce(scene,right_wall)[0]/x_area;
    wsxx_left = - getForce(scene,left_wall)[0]/x_area;
    wsxx = 0.5*(wsxx_left + wsxx_right);
    wsyy_back = getForce(scene,back_wall)[1]/y_area;
    wsyy_front = - getForce(scene,front_wall)[1]/y_area;
    wsyy = 0.5*(wsyy_back + wsyy_front);
    wszz_top = getForce(scene,top_wall)[2]/z_area;
    wszz_bottom = - getForce(scene,bottom_wall)[2]/z_area;
    wszz = 0.5*(wszz_top + wszz_bottom);
    #else
		if(interStressControl){//using the stress within the assembly
			Matrix3r Sigma = getStressTensor();
			wsxx = Sigma(0,0);
			wsyy = Sigma(1,1);
			wszz = Sigma(2,2);
		}else{
			//get box size
			updateBoxSize();
			wsxx = 0.5*(getForce(scene,right_wall)[0]- getForce(scene,left_wall)[0])/x_area ;
			wsyy = 0.5*(getForce(scene,back_wall)[1] - getForce(scene,front_wall)[1])/y_area;
			wszz = 0.5*(getForce(scene,top_wall)[2] - getForce(scene,bottom_wall)[2])/z_area;
		}
    #endif


}

void CompressionEngine::consol_ss_cylinder(){
    // sync thread storage of ForceContainer
	 scene->forces.sync();

	//calculate stress
	//get box size
	Real top = box[top_wall]->state->pos.z();
	Real bottom = box[bottom_wall]->state->pos.z();

	height = top - bottom;
	x_area = 2.0*M_PI*wall_radius * height;

	loading_area = M_PI*pow(wall_radius,2.0);
	//wall forces
    cylinder_stress = (cylinder_force)/x_area;

    wszz = 0.5*(getForce(scene,1)[2] - getForce(scene,0)[2])/loading_area;//top_wall id1 - bottom_wall id0
    //std::cerr<<"loadingarea"<<loading_area<<" wszz "<<wszz<<std::endl;
}
void CompressionEngine::checkTarget(){
	if(cubic_specimen){
    consol_ss();
    Flag_ForceTarget = false;
       if ((std::abs(wsxx - goalx)/goalx) < f_threshold){
        if ((std::abs(wsyy - goaly)/goaly) < f_threshold){
            if ((std::abs(wszz - goalz)/goalz) < f_threshold){
                Flag_ForceTarget = true;//forces acting on walls reached the target
            }
        }
    }
	}else{
    consol_ss_cylinder();
    Flag_ForceTarget = false;
       if ((std::abs(cylinder_stress - goalx)/goalx) < f_threshold){

            if ((std::abs(wszz - goalz)/goalz) < f_threshold){
                Flag_ForceTarget = true;//forces acting on walls reached the target
            }
        }
  }
}
// void CompressionEngine::quiet_system(){
//     BodyContainer::iterator bb    = scene->bodies->begin();
// 	BodyContainer::iterator iiEnd = scene->bodies->end();
//
// 	for(  ; bb!=iiEnd ; ++bb )
// 	{
//         State* state=(*bb)->state.get();
//         state->vel = Vector3r(0.0,0.0,0.0);
//         state->angVel = Vector3r(0.0,0.0,0.0);
//     }
// }

void CompressionEngine::getBox(){
	for(int i = 0; i < 6; i++){
		 box[i] = Body::byId(i);
	}
	//find particles at the boundary

}

// Real CompressionEngine::getResultantF(){
// 	scene->forces.sync();
// 	Real f=0.;
// 	/*
// 	for(int i=0;i<downBox.size();i++){
// 		//f += getForce(scene,downBox[i])[0];//force along x,f=scene->forces.getForce(id);
// 		f += scene->forces.getForce(downBox[i])[0];
// 	}
// 	*/
// 	return f;
// }
/*
void CompressionEngine::getStressStrain(){
	scene->forces.sync();
	//get forces on the left and right facets of the up box
	//Real f = getForce(scene,facet_uleft1_id)[0] + getForce(scene,facet_uleft2_id)[0]+getForce(scene,facet_uright1_id)[0] + getForce(scene,facet_uright2_id)[0] ;
	//Real f = getForce(scene,facet_uright1_id)[0] + getForce(scene,facet_uright2_id)[0] ;

	Real left = box[facet_uleft1_id]->state->pos.x();
	width = width0 - std::abs(left - left_facet_pos);
	strain[0] = std::abs(left - left_facet_pos)/width0;//shear strain
	Real top = box[top_wall]->state->pos.y();
	Real bottom = box[bottom_wall]->state->pos.y();
	height = top - bottom - thickness;
	strain[1] = (height - height0)/height0;//volumetric strain; positive sign for expanding
	//shearStress[0] = (getForce(scene,facet_dleft1_id)[0]) / (depth0*width);
	//shearStress[1] = (getForce(scene,facet_dright1_id)[0]) / (depth0*width);
	Real f=getResultantF()/ (depth0*width);
	shearStress[0] = -f;
	shearStress[1] = (getForce(scene,facet_dright1_id)[0]) / (depth0*width)+(getForce(scene,facet_dleft1_id)[0]) / (depth0*width);
}
*/
