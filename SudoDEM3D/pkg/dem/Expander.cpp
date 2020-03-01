/*************************************************************************
*  Copyright (C) 2016 by Zhswee		        		 	        		 *
*  zhswee@gmail.com      					  							 *
*																		 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Expander.hpp"

#include<sudodem/pkg/common/Box.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/Superquadrics.hpp>
#include<sudodem/pkg/dem/PolySuperellipsoid.hpp>
//#include<sudodem/pkg/dem/GJKParticle.hpp>
//create a new class base,'NonSphericalShape', in the furture.
#include<sudodem/pkg/common/Sphere.hpp>
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

#include <sys/syscall.h>//get pid of a thread

CREATE_LOGGER(Expander);
SUDODEM_PLUGIN((Expander));

Expander::~Expander(){}

void Expander::initial(){
	expandAlphaList.clear();
	double delta = (1.0 - initExpandAlpha)/expandSlice;
	for(short int i=0;i<expandSlice;i++){
		expandAlphaList.push_back(initExpandAlpha+delta*i);
	}
	expandAlphaList.push_back(1.0);//put 1. to avoid computing round errors.
	curExpandSlice = 0;
	preExpandAlpha = 1.0;//initExpandAlpha;
}
void Expander::expand(){
	//expand a partile by (1+expandAlpha) times
	double expandAlpha = expandAlphaList[curExpandSlice];
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	for (unsigned int i=0 ; bi!=biEnd; ++bi,++i )
	{
		//if((*bi)->isClump()) continue;//no clump for Superquadrics
		const shared_ptr<Body>& b = *bi;

		//std::cerr << "watch point at particle volume calc--end " << "\n";
		if (b->shape->getClassName()=="Sphere"){
			const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b->shape);
			A->radius *= expandAlpha/preExpandAlpha;
			//Aabb* aabb=static_cast<Aabb*>(b->bound->get());
			Vector3r pos = b->state->pos;
			Vector3r mincoords = b->bound->min - pos;
			Vector3r maxcoords = b->bound->max - pos;
			mincoords *= expandAlpha/preExpandAlpha;
			maxcoords *= expandAlpha/preExpandAlpha;
			b->bound->min = pos + mincoords;
			b->bound->max = pos + maxcoords;
		}
		if (b->shape->getClassName()=="PolySuperellipsoid"){
			const shared_ptr<PolySuperellipsoid>& A = SUDODEM_PTR_CAST<PolySuperellipsoid> (b->shape);
			//we expand a particle without change the properties, mass , momonent of inertia, etc, whose values are the target.
			Vector6r rxyz = A->getrxyz_ref();
			Vector3r mc = A->getMassCenter_ref();
			Vector3r mc1 = A->getMassCenter();
			double r_max = A->getr_max_ref();
			A->setRxyz(rxyz*expandAlpha);
			A->setMassCenter(mc*expandAlpha);
			A->setRmax(r_max*expandAlpha);
			//necessary to change AABB manually?
			//Aabb* aabb=static_cast<Aabb*>(b->bound->get());
			Vector3r pos = b->state->pos;
			Vector3r mincoords = b->bound->min - pos + mc1;
			Vector3r maxcoords = b->bound->max - pos + mc1;
			mincoords *= expandAlpha/preExpandAlpha;
			maxcoords *= expandAlpha/preExpandAlpha;
			b->bound->min = pos - mc1 + mincoords;
			b->bound->max = pos - mc1 + maxcoords;
		}
	}
	preExpandAlpha = expandAlpha;
	if(curExpandSlice==expandSlice){dead = true;}
	curExpandSlice++;
}

void Expander::action()
{
	if(firstRun){initial();firstRun=false;}
	if((curExpandSlice < expandSlice+1) && (scene->iter % cycleInterval == 0)){
		expand();
	}
}
