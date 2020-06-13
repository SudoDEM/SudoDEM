/*************************************************************************
*  Copyright (C) 2020 by Sway Zhao                                       *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include "TriElement_Ig2.hpp"
#include <sudodem/pkg/dem/ScGeom.hpp>

SUDODEM_PLUGIN((Ig2_TriElement_Sphere_ScGeom));

CREATE_LOGGER(Ig2_TriElement_Sphere_ScGeom);

bool Ig2_TriElement_Sphere_ScGeom::go(const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const State& state1,
							const State& state2,
							const Vector3r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
    //FIXME:updateNode() should be called before contact detection when nodes' positions change.
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	TriElement*   f = static_cast<TriElement*>(cm1.get());
    Vector3r cl = (se32.position + shift2 - f->pos);  // "contact line" in facet-local coords
    //std::cout<<"segeom trielement sphere"<<std::endl;
	// BEGIN everything in facet-local coordinates
	//
	Vector3r normal = f->normal;
	Real L = normal.dot(cl);
	if (L<0) {normal=-normal; L=-L; }

	Real sphereRadius = static_cast<Sphere*>(cm2.get())->radius;
	if (L>sphereRadius && !c->isReal() && !force) { // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false
		return false;
	}

	Vector3r cp = cl - L*normal;
	const Vector3r* ne = f->ne;

	Real penetrationDepth=0;

	Real bm = ne[0].dot(cp);
	int m=0;
	for (int i=1; i<3; ++i)
	{
		Real b=ne[i].dot(cp);
		if (bm<b) {bm=b; m=i;}
	}

	Real sh = sphereRadius*shrinkFactor;
	Real icr = f->icr-sh;
    
	if (icr<0)
	{
		LOG_WARN("a radius of a facet's inscribed circle less than zero! So, shrinkFactor is too large and would be reduced to zero.");
		shrinkFactor=0;
		icr =f->icr;
		sh = 0;
	}


	if (bm<icr) // contact with facet's surface
	{
		penetrationDepth = sphereRadius - L;
		normal.normalize();
	}
	else
	{
		cp = cp + ne[m]*(icr-bm);
		if (cp.dot(ne[(m-1<0)?2:m-1])>icr) // contact with vertex m
//			cp = facet->vertices[m];
			cp = f->vu[m]*(f->vl[m]-sh);
		else if (cp.dot(ne[m=(m+1>2)?0:m+1])>icr) // contact with vertex m+1
//			cp = facet->vertices[(m+1>2)?0:m+1];
			cp = f->vu[m]*(f->vl[m]-sh);
		normal = cl-cp;
		Real norm=normal.norm(); normal/=norm;
		penetrationDepth = sphereRadius - norm;
	}
	//
	// END everything in facet-local coordinates
	//

	if (penetrationDepth>0 || c->isReal())
	{
		shared_ptr<ScGeom> scm;
		bool isNew = !c->geom;
		if (c->geom)
			scm = SUDODEM_PTR_CAST<ScGeom>(c->geom);
		else
			scm = shared_ptr<ScGeom>(new ScGeom());

		//normal = facetAxisT*normal; // in global orientation
		scm->contactPoint = se32.position + shift2 - (sphereRadius-0.5*penetrationDepth)*normal;
		scm->penetrationDepth = penetrationDepth;
		scm->radius1 = 2*sphereRadius;
		scm->radius2 = sphereRadius;
		if (isNew) c->geom = scm;
		scm->precompute(state1,state2,scene,c,normal,isNew,shift2,false/*avoidGranularRatcheting only for sphere-sphere*/);
		return true;
	}
	return false;
}

bool Ig2_TriElement_Sphere_ScGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector3r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	//LOG_WARN("Swapped interaction order for "<<c->getId2()<<"&"<<c->getId1());
	return go(cm2,cm1,state2,state1,-shift2,force,c);
}

//PolySuperellipsoid and TriElement
SUDODEM_PLUGIN((Ig2_TriElement_PolySuperellipsoid_ScGeom));

CREATE_LOGGER(Ig2_TriElement_PolySuperellipsoid_ScGeom);

bool Ig2_TriElement_PolySuperellipsoid_ScGeom::go(const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const State& state1,
							const State& state2,
							const Vector3r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
    //FIXME:updateNode() should be called before contact detection when nodes' positions change.
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	TriElement*   f = static_cast<TriElement*>(cm1.get());//particle A is a TriElement by default.
    PolySuperellipsoid* B = static_cast<PolySuperellipsoid*>(cm2.get());
    Vector3r B_pos = se32.position + shift2;
    Vector3r cl = (B_pos - f->pos);  // "contact line" in facet-local coords
    //std::cout<<"segeom trielement sphere"<<std::endl;
	// BEGIN everything in facet-local coordinates
     

	Vector3r normal = f->normal;
	Real L = normal.dot(cl);
	if (L<0) {normal=-normal; L=-L; }	//check whether partcle B is at the negative side of the facet.

	Real sphereRadius = B->getr_max();
	if (L>sphereRadius && !c->isReal() && !force) { // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false
		return false;
	}
    //get the support point
    Matrix3r rot_mat1 = B->rot_mat2local;//to particle's system
	Matrix3r rot_mat2 = B->rot_mat2global;//to global system

	Vector2r phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
	Vector3r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
    //check whether point2 is at the negative side of the facet.
    Vector3r cl2 = point2 - f->pos;
    L = normal.dot(cl2);
    Real penetrationDepth = -L;
	bool isNew = !c->geom;
	shared_ptr<ScGeom> scm;   
	if (isNew) {
		// new interaction
		scm = shared_ptr<ScGeom>(new ScGeom());
		c->geom = scm;
	}else{
		// use data from old interaction
		scm = SUDODEM_PTR_CAST<ScGeom>(c->geom);
	}
	scm->penetrationDepth = penetrationDepth;
	if (L > 0){//no penetration
        return true;
    }
    //check if the contact point is inside the facet region.
    if(f->isPointInTriangle(point2)){

            //normal = facetAxisT*normal; // in global orientation
            scm->contactPoint = point2;//here we simply use point2 as the contact point
            //scm->penetrationDepth = penetrationDepth;
            scm->radius1 = 2*sphereRadius;
            scm->radius2 = (point2 - B_pos).norm();
            scm->precompute(state1,state2,scene,c,normal,isNew,shift2,false/*avoidGranularRatcheting only for sphere-sphere*/);
            return true;
    }else{//outside, no contact
		scm->penetrationDepth = -penetrationDepth;
        return true;
    }
	//return true;
}

SUDODEM_PLUGIN((Ip2_RolFrictMat_PolySuperellipsoidMat_FrictPhys));

void Ip2_RolFrictMat_PolySuperellipsoidMat_FrictPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;

	const shared_ptr<RolFrictMat>& mat1 = SUDODEM_PTR_CAST<RolFrictMat>(b1);
	const shared_ptr<PolySuperellipsoidMat>& mat2 = SUDODEM_PTR_CAST<PolySuperellipsoidMat>(b2);

	interaction->phys = shared_ptr<FrictPhys>(new FrictPhys());
	//const shared_ptr<FrictPhys>& contactPhysics = SUDODEM_PTR_CAST<FrictPhys>(interaction->phys);
	//ScGeom* geom = SUDODEM_CAST<ScGeom*>(interaction->geom.get());

	FrictPhys* contactPhysics = SUDODEM_CAST<FrictPhys*>(interaction->phys.get());

    Real frictionAngle = std::min(mat1->frictionAngle,mat2->frictionAngle);
    contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
    contactPhysics->kn = 2.0*mat1->Kn*mat2->Kn/(mat1->Kn+mat2->Kn);
    contactPhysics->ks = 2.0*mat1->Ks*mat2->Ks/(mat1->Ks+mat2->Ks);
};


//GJKParticle and TriElement
SUDODEM_PLUGIN((Ig2_TriElement_GJKParticle_ScGeom));

CREATE_LOGGER(Ig2_TriElement_GJKParticle_ScGeom);

bool Ig2_TriElement_GJKParticle_ScGeom::go(const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const State& state1,
							const State& state2,
							const Vector3r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
    //FIXME:updateNode() should be called before contact detection when nodes' positions change.
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	TriElement*   f = static_cast<TriElement*>(cm1.get());//particle A is a TriElement by default.
    GJKParticle* B = static_cast<GJKParticle*>(cm2.get());
    Vector3r B_pos = se32.position + shift2;
    Vector3r cl = (B_pos - f->pos);  // "contact line" in facet-local coords
    //std::cout<<"segeom trielement sphere"<<std::endl;
	// BEGIN everything in facet-local coordinates
     

	Vector3r normal = f->normal;
	Real L = normal.dot(cl);
	if (L<0) {normal=-normal; L=-L; }	//check whether partcle B is at the negative side of the facet.
	/*
	Real sphereRadius = B->getr_max();
	if (L>sphereRadius && !c->isReal() && !force) { // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false
		return false;
	}
	*/
    //get the support point
	Vector3r point2 = B->support(-normal) + B_pos;	//the closest point on particle B to the wall
    //check whether point2 is at the negative side of the facet.
    Vector3r cl2 = point2 - f->pos;
    L = normal.dot(cl2);
    Real penetrationDepth = -L;
	bool isNew = !c->geom;
	shared_ptr<ScGeom> scm;   
	if (isNew) {
		// new interaction
		scm = shared_ptr<ScGeom>(new ScGeom());
		c->geom = scm;
	}else{
		// use data from old interaction
		scm = SUDODEM_PTR_CAST<ScGeom>(c->geom);
	}
	scm->penetrationDepth = penetrationDepth;
	if (L > 0){//no penetration
        return true;
    }
    //check if the contact point is inside the facet region.
    if(f->isPointInTriangle(point2)){

            //normal = facetAxisT*normal; // in global orientation
            scm->contactPoint = point2;//here we simply use point2 as the contact point
            //scm->penetrationDepth = penetrationDepth;
			Real radius = (point2 - B_pos).norm();
            scm->radius1 = radius;//2*sphereRadius;
            scm->radius2 = radius;
            scm->precompute(state1,state2,scene,c,normal,isNew,shift2,false/*avoidGranularRatcheting only for sphere-sphere*/);
            return true;
    }else{//outside, no contact
		scm->penetrationDepth = -penetrationDepth;
        return true;
    }
	//return true;
}

SUDODEM_PLUGIN((Ip2_RolFrictMat_GJKParticleMat_FrictPhys));

void Ip2_RolFrictMat_GJKParticleMat_FrictPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;

	const shared_ptr<RolFrictMat>& mat1 = SUDODEM_PTR_CAST<RolFrictMat>(b1);
	const shared_ptr<GJKParticleMat>& mat2 = SUDODEM_PTR_CAST<GJKParticleMat>(b2);

	interaction->phys = shared_ptr<FrictPhys>(new FrictPhys());
	//const shared_ptr<FrictPhys>& contactPhysics = SUDODEM_PTR_CAST<FrictPhys>(interaction->phys);
	//ScGeom* geom = SUDODEM_CAST<ScGeom*>(interaction->geom.get());

	FrictPhys* contactPhysics = SUDODEM_CAST<FrictPhys*>(interaction->phys.get());

    Real frictionAngle = std::min(mat1->frictionAngle,mat2->frictionAngle);
    contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
    contactPhysics->kn = 2.0*mat1->Kn*mat2->Kn/(mat1->Kn+mat2->Kn);
    contactPhysics->ks = 2.0*mat1->Ks*mat2->Ks/(mat1->Ks+mat2->Ks);
};