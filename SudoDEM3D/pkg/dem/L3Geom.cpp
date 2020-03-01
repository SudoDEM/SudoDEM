
#include<sudodem/pkg/dem/L3Geom.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/Wall.hpp>
#include<sudodem/pkg/common/Facet.hpp>

#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	#include<sudodem/lib/opengl/GLUtils.hpp>
	#include<GL/glu.h>
#endif

SUDODEM_PLUGIN((L3Geom)(L6Geom)(Ig2_Sphere_Sphere_L3Geom)(Ig2_Wall_Sphere_L3Geom)(Ig2_Facet_Sphere_L3Geom)(Ig2_Sphere_Sphere_L6Geom)(Law2_L3Geom_FrictPhys_ElPerfPl)(Law2_L6Geom_FrictPhys_Linear)
	#ifdef SUDODEM_OPENGL
		(Gl1_L3Geom)(Gl1_L6Geom)
	#endif

);

L3Geom::~L3Geom(){}
void L3Geom::applyLocalForceTorque(const Vector3r& localF, const Vector3r& localT, const Interaction* I, Scene* scene, NormShearPhys* nsp) const {

	Vector2r foo; // avoid undefined ~Vector2r with clang?
	#ifdef L3_TRSF_QUATERNION
		Vector3r globF=trsf.conjugate()*localF;
	#else
		Vector3r globF=trsf.transpose()*localF; // trsf is orthonormal, therefore inverse==transpose
	#endif
	Vector3r x1c(normal*(refR1+.5*u[0])), x2c(-normal*(refR2+.5*u[0]));
	if(nsp){ nsp->normalForce=normal*globF.dot(normal); nsp->shearForce=globF-nsp->normalForce; }
	Vector3r globT=Vector3r::Zero();
	// add torque, if any
	#ifdef L3_TRSF_QUATERNION
		if(localT!=Vector3r::Zero()){	globT=trsf.conjugate()*localT; }
	#else
		if(localT!=Vector3r::Zero()){	globT=trsf.transpose()*localT; }
	#endif
	// apply force and torque
	scene->forces.addForce(I->getId1(), globF); scene->forces.addTorque(I->getId1(),x1c.cross( globF)+globT);
	scene->forces.addForce(I->getId2(),-globF); scene->forces.addTorque(I->getId2(),x2c.cross(-globF)-globT);
}

void L3Geom::applyLocalForce(const Vector3r& localF, const Interaction* I, Scene* scene, NormShearPhys* nsp) const {
	applyLocalForceTorque(localF,Vector3r::Zero(),I,scene,nsp);
}


L6Geom::~L6Geom(){}

bool Ig2_Sphere_Sphere_L3Geom::go(const shared_ptr<Shape>& s1, const shared_ptr<Shape>& s2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& I){
	return genericGo(/*is6Dof*/false,s1,s2,state1,state2,shift2,force,I);
};

bool Ig2_Sphere_Sphere_L6Geom::go(const shared_ptr<Shape>& s1, const shared_ptr<Shape>& s2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& I){
	return genericGo(/*is6Dof*/true,s1,s2,state1,state2,shift2,force,I);
};


CREATE_LOGGER(Ig2_Sphere_Sphere_L3Geom);

bool Ig2_Sphere_Sphere_L3Geom::genericGo(bool is6Dof, const shared_ptr<Shape>& s1, const shared_ptr<Shape>& s2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& I){
	// temporary hack only, to not have elastic potential energy in rigid packings with overlapping spheres
	//if(state1.blockedDOFs==State::DOF_ALL && state2.blockedDOFs==State::DOF_ALL) return false;

	const Real& r1=s1->cast<Sphere>().radius; const Real& r2=s2->cast<Sphere>().radius;
	Vector3r relPos=state2.pos+shift2-state1.pos;
	Real unDistSq=relPos.squaredNorm()-pow(std::abs(distFactor)*(r1+r2),2);
	if (unDistSq>0 && !I->isReal() && !force) return false;

	// contact exists, go ahead

	Real dist=relPos.norm();
	Real uN=dist-(r1+r2);
	Vector3r normal=relPos/dist;
	Vector3r contPt=state1.pos+(r1+0.5*uN)*normal;

	handleSpheresLikeContact(I,state1,state2,shift2,is6Dof,normal,contPt,uN,r1,r2);

	return true;

};

/*
Generic function to compute L3Geom (with colinear points), used for {sphere,facet,wall}+sphere contacts now
*/
void Ig2_Sphere_Sphere_L3Geom::handleSpheresLikeContact(const shared_ptr<Interaction>& I, const State& state1, const State& state2, const Vector3r& shift2, bool is6Dof, const Vector3r& normal, const Vector3r& contPt, Real uN, Real r1, Real r2){
	// create geometry
	if(!I->geom){
		if(is6Dof) I->geom=shared_ptr<L6Geom>(new L6Geom);
		else       I->geom=shared_ptr<L3Geom>(new L3Geom);
		L3Geom& g(I->geom->cast<L3Geom>());
		g.contactPoint=contPt;
		g.refR1=r1; g.refR2=r2;
		g.normal=normal;
		// g.trsf.setFromTwoVectors(Vector3r::UnitX(),g.normal); // quaternion just from the X-axis; does not seem to work for some reason?!
		const Vector3r& locX(g.normal);
		// initial local y-axis orientation, in the xz or xy plane, depending on which component is larger to avoid singularities
		Vector3r locY=normal.cross(std::abs(normal[1])<std::abs(normal[2])?Vector3r::UnitY():Vector3r::UnitZ()); locY-=locX*locY.dot(locX); locY.normalize();
		Vector3r locZ=normal.cross(locY);
		#ifdef L3_TRSF_QUATERNION
			Matrix3r trsf; trsf.row(0)=locX; trsf.row(1)=locY; trsf.row(2)=locZ;
			g.trsf=Quaternionr(trsf); // from transformation matrix
		#else
			g.trsf.row(0)=locX; g.trsf.row(1)=locY; g.trsf.row(2)=locZ;
		#endif
		g.u=Vector3r(uN,0,0); // zero shear displacement
		if(distFactor<0) g.u0[0]=uN;
		// L6Geom::phi is initialized to Vector3r::Zero() automatically
		//cerr<<"Init trsf=\n"<<g.trsf<<endl<<"locX="<<locX<<", locY="<<locY<<", locZ="<<locZ<<endl;
		return;
	}

	// update geometry

	/* motion of the conctact consists in rigid motion (normRotVec, normTwistVec) and mutual motion (relShearDu);
	   they are used to update trsf and u
	*/

	L3Geom& g(I->geom->cast<L3Geom>());
	const Vector3r& currNormal(normal); const Vector3r& prevNormal(g.normal);
	// normal rotation vector, between last steps
	Vector3r normRotVec=prevNormal.cross(currNormal);
	// contrary to what ScGeom::precompute does now (r2486), we take average normal, i.e. .5*(prevNormal+currNormal),
	// so that all terms in the equation are in the previous mid-step
	// the re-normalization might not be necessary for very small increments, but better do it
	Vector3r avgNormal=(approxMask&APPROX_NO_MID_NORMAL) ? prevNormal : .5*(prevNormal+currNormal);
	if(!(approxMask&APPROX_NO_RENORM_MID_NORMAL) && !(approxMask&APPROX_NO_MID_NORMAL)) avgNormal.normalize(); // normalize only if used and if requested via approxMask
	// twist vector of the normal from the last step
	Vector3r normTwistVec=avgNormal*scene->dt*.5*avgNormal.dot(state1.angVel+state2.angVel);
	// compute relative velocity
	// noRatch: take radius or current distance as the branch vector; see discussion in ScGeom::precompute (avoidGranularRatcheting)
	Vector3r c1x=((noRatch && !(r1>0)) ? ( r1*normal).eval() : (contPt-state1.pos).eval()); // used only for sphere-sphere
	Vector3r c2x=(noRatch ? (-r2*normal).eval() : (contPt-state2.pos+shift2).eval());
	//Vector3r state2velCorrected=state2.vel+(scene->isPeriodic?scene->cell->intrShiftVel(I->cellDist):Vector3r::Zero()); // velocity of the second particle, corrected with meanfield velocity if necessary
	//cerr<<"correction "<<(scene->isPeriodic?scene->cell->intrShiftVel(I->cellDist):Vector3r::Zero())<<endl;
	Vector3r relShearVel=(state2.vel+state2.angVel.cross(c2x))-(state1.vel+state1.angVel.cross(c1x));
	// account for relative velocity of particles in different cell periods
	if(scene->isPeriodic) relShearVel+=scene->cell->intrShiftVel(I->cellDist);
	relShearVel-=avgNormal.dot(relShearVel)*avgNormal;
	Vector3r relShearDu=relShearVel*scene->dt;

	/* Update of quantities in global coords consists in adding 3 increments we have computed; in global coords (a is any vector)

		1. +relShearVel*scene->dt;      // mutual motion of the contact
		2. -a.cross(normRotVec);   // rigid rotation perpendicular to the normal
		3. -a.cross(normTwistVec); // rigid rotation parallel to the normal

	*/

	// compute current transformation, by updating previous axes
	// the X axis can be prescribed directly (copy of normal)
	// the mutual motion on the contact does not change transformation
	#ifdef L3_TRSF_QUATERNION
		const Matrix3r prevTrsf(g.trsf.toRotationMatrix());
		Quaternionr prevTrsfQ(g.trsf);
	#else
		const Matrix3r prevTrsf(g.trsf); // could be reference perhaps, but we need it to compute midTrsf (if applicable)
	#endif
	Matrix3r currTrsf; currTrsf.row(0)=currNormal;
	for(int i=1; i<3; i++){
		currTrsf.row(i)=prevTrsf.row(i)-prevTrsf.row(i).cross(normRotVec)-prevTrsf.row(i).cross(normTwistVec);
	}
	#ifdef L3_TRSF_QUATERNION
		Quaternionr currTrsfQ(currTrsf);
		if((scene->iter % trsfRenorm)==0 && trsfRenorm>0) currTrsfQ.normalize();
	#else
		if((scene->iter % trsfRenorm)==0 && trsfRenorm>0){
			#if 1
				currTrsf.row(0).normalize();
				currTrsf.row(1)-=currTrsf.row(0)*currTrsf.row(1).dot(currTrsf.row(0)); // take away y projected on x, to stabilize numerically
				currTrsf.row(1).normalize();
				currTrsf.row(2)=currTrsf.row(0).cross(currTrsf.row(1));
				currTrsf.row(2).normalize();
			#else
				currTrsf=Matrix3r(Quaternionr(currTrsf).normalized());
			#endif
			#ifdef SUDODEM_DEBUG
				if(std::abs(currTrsf.determinant()-1)>.05){
					LOG_ERROR("##"<<I->getId1()<<"+"<<I->getId2()<<", |trsf|="<<currTrsf.determinant());
					g.trsf=currTrsf;
					throw runtime_error("Transformation matrix far from orthonormal.");
				}
			#endif
		}
	#endif

	/* Previous local trsf u'⁻ must be updated to current u'⁰. We have transformation T⁻ and T⁰,
		δ(a) denotes increment of a as defined above.  Two possibilities:

		1. convert to global, update, convert back: T⁰(T⁻*(u'⁻)+δ(T⁻*(u'⁻))). Quite heavy.
		2. update u'⁻ straight, using relShearVel in local coords; since relShearVel is computed
			at (t-Δt/2), we would have to find intermediary transformation (same axis, half angle;
			the same as slerp at t=.5 between the two).

			This could be perhaps simplified by using T⁰ or T⁻ since they will not differ much,
			but it would have to be verified somehow.
	*/
	// if requested via approxMask, just use prevTrsf
	#ifdef L3_TRSF_QUATERNION
		Quaternionr midTrsf=(approxMask&APPROX_NO_MID_TRSF) ? prevTrsfQ : prevTrsfQ.slerp(.5,currTrsfQ);
	#else
		Quaternionr midTrsf=(approxMask&APPROX_NO_MID_TRSF) ? Quaternionr(prevTrsf) : Quaternionr(prevTrsf).slerp(.5,Quaternionr(currTrsf));
	#endif
	//cerr<<"prevTrsf=\n"<<prevTrsf<<", currTrsf=\n"<<currTrsf<<", midTrsf=\n"<<Matrix3r(midTrsf)<<endl;

	// updates of geom here

	// midTrsf*relShearVel should have the 0-th component (approximately) zero -- to be checked
	g.u+=midTrsf*relShearDu;
	//cerr<<"midTrsf=\n"<<midTrsf<<",relShearDu="<<relShearDu<<", transformed "<<midTrsf*relShearDu<<endl;
	g.u[0]=uN; // this does not have to be computed incrementally
	#ifdef L3_TRSF_QUATERNION
		g.trsf=currTrsfQ;
	#else
		g.trsf=currTrsf;
	#endif

	// GenericSpheresContact
	g.refR1=r1; g.refR2=r2;
	g.normal=currNormal;
	g.contactPoint=contPt;

	if(is6Dof){
		// update phi, from the difference of angular velocities
		// the difference is transformed to local coord using the midTrsf transformation
		// perhaps not consistent when spheres have different radii (depends how bending moment is computed)
		I->geom->cast<L6Geom>().phi+=midTrsf*(scene->dt*(state2.angVel-state1.angVel));
	}
};

bool Ig2_Wall_Sphere_L3Geom::go(const shared_ptr<Shape>& s1, const shared_ptr<Shape>& s2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& I){
	if(scene->isPeriodic) throw std::logic_error("Ig2_Wall_Sphere_L3Geom does not handle periodic boundary conditions.");
	const Real& radius=s2->cast<Sphere>().radius; const int& ax(s1->cast<Wall>().axis); const int& sense(s1->cast<Wall>().sense);
	Real dist=state2.pos[ax]+shift2[ax]-state1.pos[ax]; // signed "distance" between centers
	if(!I->isReal() && std::abs(dist)>radius && !force) { return false; }// wall and sphere too far from each other
	// contact point is sphere center projected onto the wall
	Vector3r contPt=state2.pos+shift2; contPt[ax]=state1.pos[ax];
	Vector3r normal=Vector3r::Zero();
	// wall interacting from both sides: normal depends on sphere's position
	assert(sense==-1 || sense==0 || sense==1);
	if(sense==0) normal[ax]=dist>0?1.:-1.;
	else normal[ax]=(sense==1?1.:-1);
	Real uN=normal[ax]*dist-radius; // takes in account sense, radius and distance

	// check that the normal did not change orientation (would be abrup here)
	if(I->geom && I->geom->cast<L3Geom>().normal!=normal){
		ostringstream oss; oss<<"Ig2_Wall_Sphere_L3Geom: normal changed from ("<<I->geom->cast<L3Geom>().normal<<" to "<<normal<<" in Wall+Sphere ##"<<I->getId1()<<"+"<<I->getId2()<<" (with Wall.sense=0, a particle might cross the Wall plane, if Δt is too high)"; throw std::logic_error(oss.str().c_str());
	}
	handleSpheresLikeContact(I,state1,state2,shift2,/*is6Dof*/false,normal,contPt,uN,/*r1*/0,/*r2*/radius);
	return true;
};

bool Ig2_Facet_Sphere_L3Geom::go(const shared_ptr<Shape>& s1, const shared_ptr<Shape>& s2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& I){
	const Facet& facet(s1->cast<Facet>());
	Real radius=s2->cast<Sphere>().radius;
	// begin facet-local coordinates
		Vector3r cogLine=state1.ori.conjugate()*(state2.pos+shift2-state1.pos); // connect centers of gravity
		Vector3r normal=facet.normal; // trial contact normal
		Real planeDist=normal.dot(cogLine);
		if(std::abs(planeDist)>radius && !I->isReal() && !force) return false; // sphere too far
		if(planeDist<0){normal*=-1; planeDist*=-1; }
		Vector3r planarPt=cogLine-planeDist*normal; // project sphere center to the facet plane
		Vector3r contactPt; // facet's point closes to the sphere
		Real normDotPt[3];  // edge outer normals dot products
		for(int i=0; i<3; i++) normDotPt[i]=facet.ne[i].dot(planarPt-facet.vertices[i]);
		short w=(normDotPt[0]>0?1:0)+(normDotPt[1]>0?2:0)+(normDotPt[2]>0?4:0); // bitmask whether the closest point is outside (1,2,4 for respective edges)
		switch(w){
			case 0: contactPt=planarPt; break; // ---: inside triangle
			case 1: contactPt=getClosestSegmentPt(planarPt,facet.vertices[0],facet.vertices[1]); break; // +-- (n1)
			case 2: contactPt=getClosestSegmentPt(planarPt,facet.vertices[1],facet.vertices[2]); break; // -+- (n2)
			case 4: contactPt=getClosestSegmentPt(planarPt,facet.vertices[2],facet.vertices[0]); break; // --+ (n3)
			case 3: contactPt=facet.vertices[1]; break; // ++- (v1)
			case 5: contactPt=facet.vertices[0]; break; // +-+ (v0)
			case 6: contactPt=facet.vertices[2]; break; // -++ (v2)
			case 7: throw logic_error("Ig2_Facet_Sphere_L3Geom: Impossible sphere-facet intersection (all points are outside the edges). (please report bug)"); // +++ (impossible)
			default: throw logic_error("Ig2_Facet_Sphere_L3Geom: Nonsense intersection value. (please report bug)");
		}
		normal=cogLine-contactPt; // normal is now the contact normal, still in local coords
		if(!I->isReal() && normal.squaredNorm()>radius*radius && !force) { return false; } // fast test before sqrt
		Real dist=normal.norm(); normal/=dist; // normal is unit vector now
	// end facet-local coordinates
	normal=state1.ori*normal; // normal is in global coords now
	handleSpheresLikeContact(I,state1,state2,shift2,/*is6Dof*/false,normal,/*contact pt*/state2.pos+shift2-normal*dist,dist-radius,0,radius);
	return true;
}

bool Law2_L3Geom_FrictPhys_ElPerfPl::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){
	L3Geom* geom=static_cast<L3Geom*>(ig.get()); FrictPhys* phys=static_cast<FrictPhys*>(ip.get());

	// compute force
	Vector3r& localF(geom->F);
	localF=geom->relU().cwiseProduct(Vector3r(phys->kn,phys->ks,phys->ks));
	// break if necessary
	if(localF[0]>0 && !noBreak) return false;

	if(!noSlip){
		// plastic slip, if necessary; non-zero elastic limit only for compression
		Real maxFs=-min((Real)0.,localF[0]*phys->tangensOfFrictionAngle); Eigen::Map<Vector2r> Fs(&localF[1]);
		//cerr<<"u="<<geom->relU()<<", maxFs="<<maxFs<<", Fn="<<localF[0]<<", |Fs|="<<Fs.norm()<<", Fs="<<Fs<<endl;
		if(Fs.squaredNorm()>maxFs*maxFs){
			Real ratio=sqrt(maxFs*maxFs/Fs.squaredNorm());
			Vector3r u0slip=(1-ratio)*Vector3r(/*no slip in the normal sense*/0,geom->relU()[1],geom->relU()[2]);
			geom->u0+=u0slip; // increment plastic displacement
			Fs*=ratio; // decrement shear force value;
			if(scene->trackEnergy){ Real dissip=Fs.norm()*u0slip.norm(); if(dissip>0) scene->energy->add(dissip,"plastDissip",plastDissipIx,/*reset*/false); }
		}
	}
	if(scene->trackEnergy)	{ scene->energy->add(0.5*(pow(geom->relU()[0],2)*phys->kn+(pow(geom->relU()[1],2)+pow(geom->relU()[2],2))*phys->ks),"elastPotential",elastPotentialIx,/*reset at every timestep*/true); }
	// apply force: this converts the force to global space, updates NormShearPhys::{normal,shear}Force, applies to particles
	geom->applyLocalForce(localF,I,scene,phys);
	return true;
}


bool Law2_L6Geom_FrictPhys_Linear::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){
	L6Geom& geom=ig->cast<L6Geom>(); FrictPhys& phys=ip->cast<FrictPhys>();

	// simple linear relationships
	Vector3r localF=geom.relU().cwiseProduct(Vector3r(phys.kn,phys.ks,phys.ks));
	Vector3r localT=charLen*(geom.relPhi().cwiseProduct(Vector3r(phys.kn,phys.ks,phys.ks)));

	geom.applyLocalForceTorque(localF,localT,I,scene,static_cast<NormShearPhys*>(ip.get()));
	return true;
}

#ifdef SUDODEM_OPENGL
bool Gl1_L3Geom::axesLabels;
Real Gl1_L3Geom::axesWd;
Real Gl1_L3Geom::axesScale;
Real Gl1_L3Geom::uPhiWd;
Real Gl1_L3Geom::uScale;
Real Gl1_L6Geom::phiScale;

void Gl1_L3Geom::go(const shared_ptr<IGeom>& ig, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool){ draw(ig); }
void Gl1_L6Geom::go(const shared_ptr<IGeom>& ig, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool){ draw(ig,true,phiScale); }

void Gl1_L3Geom::draw(const shared_ptr<IGeom>& ig, bool isL6Geom, const Real& phiScale){
	const L3Geom& g(ig->cast<L3Geom>());
	glTranslatev(g.contactPoint);
	#ifdef L3_TRSF_QUATERNION
		//glMultMatrixd(Eigen::Affine3d(Matrix3r(g.trsf).transpose()).data());
		glMultMatrix(Eigen::Transform<Real,3,Eigen::Affine>(Matrix3r(g.trsf).transpose()).data());
	#else
		//glMultMatrixd(Eigen::Affine3d(g.trsf.transpose()).data());
		glMultMatrix(Eigen::Transform<Real,3,Eigen::Affine>(g.trsf.transpose()).data());
	#endif
	Real rMin=g.refR1<=0?g.refR2:(g.refR2<=0?g.refR1:min(g.refR1,g.refR2));
	if(axesWd>0){
		glLineWidth(axesWd);
		for(int i=0; i<3; i++){
			Vector3r pt=Vector3r::Zero(); pt[i]=.5*rMin*axesScale; Vector3r color=.3*Vector3r::Ones(); color[i]=1;
			GLUtils::GLDrawLine(Vector3r::Zero(),pt,color);
			if(axesLabels) GLUtils::GLDrawText(string(i==0?"x":(i==1?"y":"z")),pt,color);
		}
	}
	if(uPhiWd>0){
		glLineWidth(uPhiWd);
		if(uScale!=0) GLUtils::GLDrawLine(Vector3r::Zero(),uScale*g.relU(),Vector3r(0,1,.5));
		if(isL6Geom && phiScale>0) GLUtils::GLDrawLine(Vector3r::Zero(),ig->cast<L6Geom>().relPhi()/Mathr::PI*rMin*phiScale,Vector3r(.8,0,1));
	}
	glLineWidth(1.);
};

#endif
