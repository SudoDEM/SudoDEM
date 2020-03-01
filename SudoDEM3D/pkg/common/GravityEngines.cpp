/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<sudodem/pkg/common/GravityEngines.hpp>
#include<sudodem/pkg/common/PeriodicEngines.hpp>
#include<sudodem/core/BodyContainer.hpp>
#include<sudodem/core/Scene.hpp>
#include<boost/regex.hpp>

SUDODEM_PLUGIN((GravityEngine)(CentralGravityEngine)(AxialGravityEngine)(HdapsGravityEngine));
CREATE_LOGGER(GravityEngine);

void GravityEngine::action(){
	if (warnOnce) {warnOnce=false; LOG_ERROR("GravityEngine is deprecated, consider using Newton::gravity instead (unless gravitational energy has to be tracked - not implemented with the newton attribute).")}
	const bool trackEnergy(scene->trackEnergy);
	const Real dt(scene->dt);
	SUDODEM_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
		// skip clumps, only apply forces on their constituents
		if(b->isClump()) continue;
		if(mask!=0 && !b->maskCompatible(mask)) continue;
		scene->forces.addForce(b->getId(),gravity*b->state->mass);
		// work done by gravity is "negative", since the energy appears in the system from outside
		if(trackEnergy) scene->energy->add(-gravity.dot(b->state->vel)*b->state->mass*dt,"gravWork",fieldWorkIx,/*non-incremental*/false);
	} SUDODEM_PARALLEL_FOREACH_BODY_END();
}

void CentralGravityEngine::action(){
	const Vector3r& centralPos=Body::byId(centralBody)->state->pos;
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(b->isClump() || b->getId()==centralBody) continue; // skip clumps and central body
		if(mask!=0 && !b->maskCompatible(mask)) continue;
		Real F=accel*b->state->mass;
		Vector3r toCenter=centralPos-b->state->pos; toCenter.normalize();
		scene->forces.addForce(b->getId(),F*toCenter);
		if(reciprocal) scene->forces.addForce(centralBody,-F*toCenter);
	}
}

void AxialGravityEngine::action(){
	FOREACH(const shared_ptr<Body>&b, *scene->bodies){
		if(!b || b->isClump()) continue;
		if(mask!=0 && !b->maskCompatible(mask)) continue;
		/* http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html */
		const Vector3r& x0=b->state->pos;
		const Vector3r& x1=axisPoint;
		const Vector3r x2=axisPoint+axisDirection;
		Vector3r closestAxisPoint=x1+(x2-x1) * /* t */ (-(x1-x0).dot(x2-x1))/((x2-x1).squaredNorm());
		Vector3r toAxis=closestAxisPoint-x0; toAxis.normalize();
		if(toAxis.squaredNorm()==0) continue;
		scene->forces.addForce(b->getId(),acceleration*b->state->mass*toAxis);
	}
}


Vector2i HdapsGravityEngine::readSysfsFile(const string& name){
	char buf[256];
	ifstream f(name.c_str());
	if(!f.is_open()) throw std::runtime_error(("HdapsGravityEngine: unable to open file "+name).c_str());
	f.read(buf,256);f.close();
	const boost::regex re("\\(([0-9+-]+),([0-9+-]+)\\).*");
   boost::cmatch matches;
	if(!boost::regex_match(buf,matches,re)) throw std::runtime_error(("HdapsGravityEngine: error parsing data from "+name).c_str());
	//cerr<<matches[1]<<","<<matches[2]<<endl;
	return Vector2i(boost::lexical_cast<int>(matches[1]),boost::lexical_cast<int>(matches[2]));

}

void HdapsGravityEngine::action(){
	if(!calibrated) { calibrate=readSysfsFile(hdapsDir+"/calibrate"); calibrated=true; }
	Real now=PeriodicEngine::getClock();
	if(now-lastReading>1e-3*msecUpdate){
		Vector2i a=readSysfsFile(hdapsDir+"/position");
		lastReading=now;
		a-=calibrate;
		if(std::abs(a[0]-accel[0])>updateThreshold) accel[0]=a[0];
		if(std::abs(a[1]-accel[1])>updateThreshold) accel[1]=a[1];
		Quaternionr trsf(AngleAxisr(.5*accel[0]*M_PI/180.,-Vector3r::UnitY())*AngleAxisr(.5*accel[1]*M_PI/180.,-Vector3r::UnitX()));
		gravity=trsf*zeroGravity;
	}
	GravityEngine::action();
}
