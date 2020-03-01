/*************************************************************************
*  Copyright (C) 2018 by Shiwei Zhao                                     *
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*  Copyright (C) 2006 by Bruno Chareyre                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include<sudodem/lib/base/Math.hpp>
#include"Ig2_Basic_ScGeom.hpp"
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/common/Disk.hpp>
//#include<sudodem/pkg/common/Box.hpp>
//#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/core/Scene.hpp>

#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/common/InteractionLoop.hpp>


SUDODEM_PLUGIN((Ig2_Disk_Disk_ScGeom)
							 (Ig2_Fwall_Disk_ScGeom)
							 //(Ig2_Box_Disk_ScGeom)(Ig2_Box_Disk_ScGeom6D)
							 (Ig2_Wall_Disk_ScGeom));

CREATE_LOGGER(Ig2_Fwall_Disk_ScGeom);


bool Ig2_Disk_Disk_ScGeom::go(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	TIMING_DELTAS_START();
	const Se2r& se21=state1.se2; const Se2r& se22=state2.se2;
	const Disk *s1=static_cast<Disk*>(cm1.get()), *s2=static_cast<Disk*>(cm2.get());
	Vector2r normal=(se22.position+shift2)-se21.position;
	if (!c->isReal() && !force) {//don't fast-check distance if geometry will be updated anyway
		Real penetrationDepthSq=pow(interactionDetectionFactor*(s1->radius+s2->radius),2) - normal.squaredNorm();
		if (penetrationDepthSq<0) {
			TIMING_DELTAS_CHECKPOINT("Ig2_Disk_Disk_ScGeom");
			return false;
		}
	}
	shared_ptr<ScGeom> scm;
	bool isNew = !c->geom;
	if(!isNew) scm=SUDODEM_PTR_CAST<ScGeom>(c->geom);
	else { scm=shared_ptr<ScGeom>(new ScGeom()); c->geom=scm; }
	Real norm=normal.norm(); normal/=norm; // normal is unit vector now
#ifdef SUDODEM_DEBUG
	if(norm==0) throw runtime_error(("Zero distance between disks #"+boost::lexical_cast<string>(c->getId1())+" and #"+boost::lexical_cast<string>(c->getId2())+".").c_str());
#endif
	Real penetrationDepth=s1->radius+s2->radius-norm;
	scm->contactPoint=se21.position+(s1->radius-0.5*penetrationDepth)*normal;//0.5*(pt1+pt2);
	scm->penetrationDepth=penetrationDepth;
	scm->radius1=s1->radius;
	scm->radius2=s2->radius;
	scm->precompute(state1,state2,scene,c,normal,isNew,shift2,avoidGranularRatcheting);
	TIMING_DELTAS_CHECKPOINT("Ig2_Disk_Disk_ScGeom");
	return true;
}

bool Ig2_Disk_Disk_ScGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector2r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	return go(cm1,cm2,state2,state1,-shift2,force,c);
}

/********* Wall + Disk **********/

bool Ig2_Wall_Disk_ScGeom::go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& c){
	Wall* wall=static_cast<Wall*>(cm1.get());
	const Real radius=static_cast<Disk*>(cm2.get())->radius;
	const int& ax(wall->axis);
	Real dist=(state2.pos)[ax]+shift2[ax]-state1.pos[ax]; // signed "distance" between centers
	if(!c->isReal() && std::abs(dist)>radius && !force) { return false; }// wall and disk too far from each other

	// contact point is disk center projected onto the wall
	Vector2r contPt=state2.pos+shift2; contPt[ax]=state1.pos[ax];
	Vector2r normal(0.,0.);
	// wall interacting from both sides: normal depends on disk's position
	assert(wall->sense==-1 || wall->sense==0 || wall->sense==1);
	if(wall->sense==0) normal[ax]=dist>0?1.:-1.;
	else normal[ax]=wall->sense==1?1.:-1;

	bool isNew=!c->geom;
	if(isNew) c->geom=shared_ptr<ScGeom>(new ScGeom());
	const shared_ptr<ScGeom>& ws=SUDODEM_PTR_CAST<ScGeom>(c->geom);
	ws->radius1=ws->radius2=radius; // do the same as for facet-disk: wall's "radius" is the same as the disk's radius
	ws->contactPoint=contPt;
	ws->penetrationDepth=-(std::abs(dist)-radius);
	// ws->normal is assigned by precompute
	ws->precompute(state1,state2,scene,c,normal,isNew,shift2,noRatch);
	return true;
}


bool Ig2_Fwall_Disk_ScGeom::go(const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const State& state1,
							const State& state2,
							const Vector2r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
	//TIMING_DELTAS_START();
	const Se2r& se21=state1.se2; const Se2r& se22=state2.se2;
	Fwall*  fwall = static_cast<Fwall*>(cm1.get());
	bool isNew = !c->geom;
 	if (isNew) c->geom = shared_ptr<ScGeom>(new ScGeom());
	const shared_ptr<ScGeom>& scm=SUDODEM_PTR_CAST<ScGeom>(c->geom);
	/* could be written as (needs to be tested):
	 * Vector3r cl=se31.orientation.Conjugate()*(se32.position-se31.position);
	 */
	//Matrix2r fw_rot = se21.rotation.toRotationMatrix();
	//Matrix2r fw = fw_rot.transpose();
	// local orientation
	//no need to rotate as we do not rotate fwall during simulations
	//Vector2r cl = fw*(se22.position + shift2 - se21.position);  // "contact line" in Fwall-local coords
	Vector2r cl = se22.position + shift2 - se21.position;

	Vector2r normal = fwall->normal;
	Real L = normal.dot(cl);
	if (L<0) {normal=-normal; L=-L; }

	Real diskRadius = static_cast<Disk*>(cm2.get())->radius;
	Real penetrationDepth = diskRadius - L;
	if (penetrationDepth < 1E-18 && !c->isReal() && !force) { // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false
	//	TIMING_DELTAS_CHECKPOINT("Ig2_Fwall_Disk_ScGeom");
		return false;
	}
  //projection of cl along vu of Fwall
	Vector2r cp = cl - L*normal;
	Real cp_l = cp.norm();
	//cout<<"ffff"<<(cp_l > fwall->vl*0.5 +diskRadius)<<endl;
	/*if (cp_l > fwall->vl*0.5 +diskRadius){
		scm->penetrationDepth = 0;
		return true;
	} //not contact
	*/
	if (cp_l > fwall->vl*0.5){
		penetrationDepth *= 0.5;
	}
  //cout<<"c-phys="<<c->phys<<"c-geom"<<c->geom<<endl;


	//normal = fw_rot*normal; // in global orientation
	scm->contactPoint = se22.position + shift2 - (diskRadius-0.5*penetrationDepth)*normal;
	scm->penetrationDepth = penetrationDepth;
	scm->radius1 = 2*diskRadius;
	scm->radius2 = diskRadius;

	//cout<<"contactpoint"<<scm->contactPoint<<"isnew"<<isNew<<endl;
	scm->precompute(state1,state2,scene,c,normal,isNew,shift2,false/*avoidGranularRatcheting only for disk-disk*/);
	return true;
}


bool Ig2_Fwall_Disk_ScGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector2r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	//LOG_WARN("Swapped interaction order for "<<c->getId2()<<"&"<<c->getId1());
	return go(cm2,cm1,state2,state1,-shift2,force,c);
}
