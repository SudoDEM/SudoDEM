#include<sudodem/pkg/dem/DomainLimiter.hpp>
#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/pkg/dem/Shop.hpp>

SUDODEM_PLUGIN((DomainLimiter)(LawTester)
	#ifdef SUDODEM_OPENGL
		(GlExtra_LawTester)(GlExtra_OctreeCubes)
	#endif
);

void DomainLimiter::action(){
	std::list<Body::id_t> out;
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if((!b) or ((mask>0) and ((b->groupMask & mask)==0))) continue;
		const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		if (sphere){ //Delete only spheres
			const Vector3r& p(b->state->pos);
			if(p[0]<lo[0] || p[0]>hi[0] || p[1]<lo[1] || p[1]>hi[1] || p[2]<lo[2] || p[2]>hi[2]) {
				out.push_back(b->id);
				nDeleted++;
				mDeleted+=b->state->mass;
				Real r = sphere->radius; vDeleted+=(4/3.)*Mathr::PI*pow(r,3);
			}
		}
	}
	FOREACH(Body::id_t id, out){
		scene->bodies->erase(id,false);
	}
}

#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/L3Geom.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
#include<sudodem/lib/smoothing/LinearInterpolate.hpp>
#include<sudodem/lib/pyutil/gil.hpp>

CREATE_LOGGER(LawTester);

void LawTester::postLoad(LawTester&){
	if(ids.size()==0) return; // uninitialized object, don't do nothing at all
	if(ids.size()!=2) throw std::invalid_argument("LawTester.ids: exactly two values must be given.");
	if(disPath.empty() && rotPath.empty()) throw invalid_argument("LawTester.{disPath,rotPath}: at least one point must be given.");
	if(pathSteps.empty()) throw invalid_argument("LawTester.pathSteps: at least one value must be given.");
	size_t pathSize=max(disPath.size(),rotPath.size());
	// update path points
	_path.clear(); _path.push_back(Vector6r::Zero());
	for(size_t i=0; i<pathSize; i++) {
		Vector6r pt;
		pt.head<3>()=Vector3r(i<disPath.size()?disPath[i]:(disPath.empty()?Vector3r::Zero():*(disPath.rbegin())));
		pt.tail<3>()=Vector3r(i<rotPath.size()?rotPath[i]:(rotPath.empty()?Vector3r::Zero():*(rotPath.rbegin())));
		_path.push_back(pt);

	}
	// update time points from distances, repeat last distance if shorter than path
	_pathT.clear(); _pathT.push_back(0);
	for(size_t i=0; i<pathSteps.size(); i++) _pathT.push_back(_pathT[i]+pathSteps[i]);
	int lastDist=pathSteps[pathSteps.size()-1];
	for(size_t i=pathSteps.size(); i<pathSize; i++) _pathT.push_back(*(_pathT.rbegin())+lastDist);
}

void LawTester::action(){
	Vector2r foo; // avoid undefined ~Vector2r with clang?
	if(ids.size()!=2) throw std::invalid_argument("LawTester.ids: exactly two values must be given.");
	LOG_DEBUG("=================== LawTester step "<<step<<" ========================");
	const shared_ptr<Interaction> Inew=scene->interactions->find(ids[0],ids[1]);
	string strIds("##"+boost::lexical_cast<string>(ids[0])+"+"+boost::lexical_cast<string>(ids[1]));
	// interaction not found at initialization
	if(!I && (!Inew || !Inew->isReal())){
		LOG_WARN("Interaction "<<strIds<<" does not exist (yet?), no-op."); return;
		//throw std::runtime_error("LawTester: interaction "+strIds+" does not exist"+(Inew?" (to be honest, it does exist, but it is not real).":"."));
	}
	// interaction was deleted meanwhile
	if(I && (!Inew || !Inew->isReal())) throw std::runtime_error("LawTester: interaction "+strIds+" was deleted"+(Inew?" (is not real anymore).":"."));
	// different interaction object
	if(I && Inew && I!=Inew) throw std::logic_error("LawTester: interacion "+strIds+" is a different object now?!");
	assert(Inew);
	bool doInit=(!I);
	if(doInit) I=Inew;

	id1=I->getId1(); id2=I->getId2();
	// test object types
	GenericSpheresContact* gsc=dynamic_cast<GenericSpheresContact*>(I->geom.get());
	ScGeom* scGeom=dynamic_cast<ScGeom*>(I->geom.get());
	L3Geom* l3Geom=dynamic_cast<L3Geom*>(I->geom.get());
	L6Geom* l6Geom=dynamic_cast<L6Geom*>(I->geom.get());
	ScGeom6D* scGeom6d=dynamic_cast<ScGeom6D*>(I->geom.get());
	bool hasRot=(l6Geom || scGeom6d);
	//NormShearPhys* phys=dynamic_cast<NormShearPhys*>(I->phys.get());			//Disabled because of warning
	if(!gsc) throw std::invalid_argument("LawTester: IGeom of "+strIds+" not a GenericSpheresContact.");
	if(!scGeom && !l3Geom) throw std::invalid_argument("LawTester: IGeom of "+strIds+" is neither ScGeom, nor L3Geom (or L6Geom).");
	assert(!((bool)scGeom && (bool)l3Geom)); // nonsense
	// get body objects
	State *state1=Body::byId(id1,scene)->state.get(), *state2=Body::byId(id2,scene)->state.get();
	scene->forces.sync();
	if(state1->blockedDOFs!=State::DOF_ALL) { LOG_INFO("Blocking all DOFs for #"<<id1); state1->blockedDOFs=State::DOF_ALL;}
	if(state2->blockedDOFs!=State::DOF_ALL) { LOG_INFO("Blocking all DOFs for #"<<id2); state2->blockedDOFs=State::DOF_ALL;}


	if(step-1>*(_pathT.rbegin())){
		LOG_INFO("Last step done, setting zero velocities on #"<<id1<<", #"<<id2<<".");
		state1->vel=state1->angVel=state2->vel=state2->angVel=Vector3r::Zero();
		uTest=uTestNext;
		if(doneHook.empty()){ LOG_INFO("No doneHook set, dying."); dead=true; }
		else{ LOG_INFO("Running doneHook: "<<doneHook);	pyRunString(doneHook);}
		return;
	}
	/* initialize or update local axes and trsf */
	uGeom.tail<3>()=Vector3r(NaN,NaN,NaN);
	if(!l3Geom){ // IGeom's that don't have local axes
		axX=gsc->normal; /* just in case */ axX.normalize();
		if(doInit){ // initialization of the new interaction -- define local axes
			// take vector in the y or z direction, depending on its length; arbitrary, but one of them is sure to be non-zero
			axY=axX.cross(std::abs(axX[1])<std::abs(axX[2])?Vector3r::UnitY():Vector3r::UnitZ());
			axY.normalize();
			axZ=axX.cross(axY);
			LOG_DEBUG("Initial axes x="<<axX<<", y="<<axY<<", z="<<axZ);
			if(scGeom6d) uGeom.tail<3>()=Vector3r::Zero();
		} else { // udpate of an existing interaction
			if(scGeom){
				scGeom->rotate(axY); scGeom->rotate(axZ);
				scGeom->rotate(shearTot);
				shearTot+=scGeom->shearIncrement();
				uGeom.head<3>()=Vector3r(-scGeom->penetrationDepth,shearTot.dot(axY),shearTot.dot(axZ));
				if(scGeom6d) uGeom.tail<3>()=-1.*Vector3r(scGeom6d->getTwist(),scGeom6d->getBending().dot(axY),scGeom6d->getBending().dot(axZ));
			}
			else{ // d3dGeom
				throw runtime_error("Geom type not yet supported.");
			}
		}
		// update the transformation
		// the matrix is orthonormal, since axX, axY are normalized and and axZ is their cross-product
		trsf.row(0)=axX; trsf.row(1)=axY; trsf.row(2)=axZ;
	} else {
		trsf=Matrix3r(l3Geom->trsf);
		axX=trsf.row(0); axY=trsf.row(1); axZ=trsf.row(2);
		uGeom.head<3>()=l3Geom->u;
		if(l6Geom) uGeom.tail<3>()=l6Geom->phi;
	}
	// perform all shearing by translation, as it does not induce bending
	if(hasRot && rotWeight!=0){ LOG_INFO("LawTester.rotWeight set to 0 (was "<<rotWeight<<"), since rotational DoFs are in use."); rotWeight=0; }
	contPt=gsc->contactPoint;
	refLength=gsc->refR1+gsc->refR2;
	renderLength=.5*refLength;

	// here we go ahead, finally
	Vector6r uu=linearInterpolate<Vector6r,int>(step,_pathT,_path,_interpPos);
	Vector6r dUU=uu-uuPrev; uuPrev=uu;
	Vector3r dU(dUU.head<3>()), dPhi(dUU.tail<3>());
	//Vector3r dU=u-uPrev.head<3>(); uPrev.head<3>()=u;
	//Vector3r dPhi=phi-uPrev.tail<3>(); uPrev.tail<3>()=phi;
	if(displIsRel){
		LOG_DEBUG("Relative displacement diff is "<<dU<<" (will be normalized by "<<gsc->refR1+gsc->refR2<<")");
		dU*=refLength;
	}
	LOG_DEBUG("Absolute diff is: displacement "<<dU<<", rotation "<<dPhi);
	uTest=uTestNext; // the value that was next in the previous step is the current one now
	uTestNext.head<3>()+=dU; uTestNext.tail<3>()+=dPhi;

	// reset velocities where displacement is controlled
	//for(int i=0; i<3; i++){ if(forceControl[i]==0){ state1.vel[i]=0; state2.vel[i]=0; }

	// shear is applied as rotation of id2: dε=r₁dθ → dθ=dε/r₁;
	Vector3r vel[2],angVel[2];
	//State* states[]={state1,state2};
	for(int i=0; i<2; i++){
		int sign=(i==0?-1:1);
		Real weight=(i==0?1-idWeight:idWeight);
		// FIXME: this should not use refR1, but real CP-particle distance perhaps?
		Real radius=(i==0?gsc->refR1:gsc->refR2);
		Real relRad=radius/refLength;
		// signed and weighted displacement/rotation to be applied on this sphere (reversed for #0)
		// some rotations must cancel the sign, by multiplying by sign again
		Vector3r ddU=sign*dU*weight;

		// twist can be still distributed with idWeight (!)
		Vector3r ddPhi=sign*dPhi*(1-relRad); /* shear angles must distribute to both, otherwise it would induce shear */
		ddPhi[0]=sign*dPhi[0]*weight; // twist can be still distributed with idWeight
		vel[i]=angVel[i]=Vector3r::Zero();

		// normal displacement

		vel[i]+=axX*ddU[0]/scene->dt;

		// shear rotation

		//   multiplication by sign cancels sign in ddU, since rotation is non-symmetric (to increase shear, both spheres have the same rotation)
		//   (unlike shear displacement, which is symmetric)
		// rotation around Z (which gives y-shear) must be inverted: +ry gives +εzm while -rz gives +εy
		Real rotZ=-sign*rotWeight*ddU[1]/radius, rotY=sign*rotWeight*ddU[2]/radius;
		angVel[i]+=(rotY*axY+rotZ*axZ)/scene->dt;

		// shear displacement

		// angle that is traversed by a sphere in order to give desired ddU when displaced on the branch of r1+r2
		// FIXME: is the branch value correct here?!
		Real arcAngleY=atan((1-rotWeight)*ddU[1]/radius), arcAngleZ=atan((1-rotWeight)*ddU[2]/radius);
		// same, but without the atan, which can be disregarded for small increments:
		//    Real arcAngleY=(1-rotWeight)*ddU[1]/radius, arcAngleZ=(1-rotWeight)*ddU[2]/radius;
		vel[i]+=axY*radius*sin(arcAngleY)/scene->dt; vel[i]+=axZ*radius*sin(arcAngleZ)/scene->dt;

		// compensate distance increase caused by motion along the perpendicular axis
		// cos(argAngle*) is always positive, regardless of the orientation
		// and the compensation is always in the -εx sense (-sign → +1 for #0, -1 for #1)
		vel[i]+=-sign*axX*radius*((1-cos(arcAngleY))+(1-cos(arcAngleZ)))/scene->dt;

		// rotation, convert from local to global
		angVel[i]+=trsf.transpose()*ddPhi;

		LOG_DEBUG("vel="<<vel[i]<<", angVel="<<angVel[i]<<", rotY,rotZ="<<rotY<<","<<rotZ<<", arcAngle="<<arcAngleY<<","<<arcAngleZ<<", sign="<<sign<<", weight="<<weight);

	}
	state1->vel=vel[0]; state1->angVel=angVel[0];
	state2->vel=vel[1]; state2->angVel=angVel[1];
	LOG_DEBUG("Body #"<<id1<<", setting vel="<<vel[0]<<", angVel="<<angVel[0]);
	LOG_DEBUG("Body #"<<id2<<", setting vel="<<vel[1]<<", angVel="<<angVel[1]);

	/* find out where are we at in the path, run hooks if approriate */
	// _pathT has the first (zero) value added by us, so we skip it
	int nPathT=_pathT.size();
	for(int i=1; i<nPathT; i++){
		// i-th point on _pathT is (i-1)th on path; run corresponding hook, if it exists
		if(step==_pathT[i] && ((int) hooks.size())>(i-1) && !hooks[i-1].empty()) pyRunString(hooks[i-1]);
	}
	step++;
}

#ifdef SUDODEM_OPENGL
#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/lib/opengl/GLUtils.hpp>
#include<sudodem/pkg/common/GLDrawFunctors.hpp>
#include<sudodem/pkg/common/OpenGLRenderer.hpp>
#include<GL/glu.h>

CREATE_LOGGER(GlExtra_LawTester);

void GlExtra_LawTester::render(){
	// scene object changed (after reload, for instance), for re-initialization
	if(tester && tester->scene!=scene) tester=shared_ptr<LawTester>();

	if(!tester){ FOREACH(shared_ptr<Engine> e, scene->engines){ tester=SUDODEM_PTR_DYN_CAST<LawTester>(e); if(tester) break; } }
	if(!tester){ LOG_ERROR("No LawTester in O.engines, killing myself."); dead=true; return; }

	//if(tester->renderLength<=0) return;
	glColor3v(Vector3r(1,0,1));

	// switch to local coordinates
	glTranslatev(tester->contPt);
	//glMultMatrixd(Eigen::Affine3d(tester->trsf.transpose()).data());
	glMultMatrix(Eigen::Transform<Real,3,Eigen::Affine>(tester->trsf.transpose()).data());


	glDisable(GL_LIGHTING);
	//glColor3v(Vector3r(1,0,1));
	//glBegin(GL_LINES); glVertex3v(Vector3r::Zero()); glVertex3v(.1*Vector3r::Ones()); glEnd();
	//GLUtils::GLDrawText(string("This is the contact point!"),Vector3r::Zero(),Vector3r(1,0,1));

	// local axes
	glLineWidth(2.);
	for(int i=0; i<3; i++){
		Vector3r pt=Vector3r::Zero(); pt[i]=.5*tester->renderLength; Vector3r color=.3*Vector3r::Ones(); color[i]=1;
		GLUtils::GLDrawLine(Vector3r::Zero(),pt,color);
		GLUtils::GLDrawText(string(i==0?"x":(i==1?"y":"z")),pt,color);
	}

	// put the origin to the initial (no-shear) point, so that the current point appears at the contact point
	glTranslatev(Vector3r(0,tester->uTestNext[1],tester->uTestNext[2]));


	const int t(tester->step); const vector<int>& TT(tester->_pathT); const vector<Vector6r>& VV(tester->_path);
	size_t numSegments=TT.size();
	const Vector3r colorBefore=Vector3r(.7,1,.7), colorAfter=Vector3r(1,.7,.7);

	// scale displacement, if they have the strain meaning
	Real scale=1;
	if(tester->displIsRel) scale=tester->refLength;

	// find maximum displacement, draw axes in the shear plane
	Real displMax=0;
	FOREACH(const Vector6r& v, VV) displMax=max(v.head<3>().squaredNorm(),displMax);
	displMax=1.2*scale*sqrt(displMax);

	glLineWidth(1.);
	GLUtils::GLDrawLine(Vector3r(0,-displMax,0),Vector3r(0,displMax,0),Vector3r(.5,0,0));
	GLUtils::GLDrawLine(Vector3r(0,0,-displMax),Vector3r(0,0,displMax),Vector3r(.5,0,0));

	// draw displacement path
	glLineWidth(4.);
	for(size_t segment=0; segment<numSegments-1; segment++){
		// different colors before and after the current point
		Real t0=TT[segment],t1=TT[segment+1];
		const Vector3r &from=-VV[segment].head<3>()*scale, &to=-VV[segment+1].head<3>()*scale;
		// current segment
		if(t>t0 && t<t1){
			Real norm=(t-t0)/(t1-t0);
			GLUtils::GLDrawLine(from,from+(to-from)*norm,colorBefore);
			GLUtils::GLDrawLine(from+(to-from)*norm,to,colorAfter);
		} else {	// other segment
			GLUtils::GLDrawLine(from,to,t<t0?colorAfter:colorBefore);
		}
	}

	glLineWidth(1.);
}


void GlExtra_OctreeCubes::postLoad(GlExtra_OctreeCubes&){
	if(boxesFile.empty()) return;
	boxes.clear();
	ifstream txt(boxesFile.c_str());
	while(!txt.eof()){
		Real data[8];
		for(int i=0; i<8; i++){ if(i<7 && txt.eof()) goto done; txt>>data[i]; }
		OctreeBox ob; Vector3r mn(data[0],data[1],data[2]), mx(data[3],data[4],data[5]);
		ob.center=.5*(mn+mx); ob.extents=(.5*(mx-mn)); ob.level=(int)data[6]; ob.fill=(int)data[7];
		// for(int i=0; i<=ob.level; i++) cerr<<"\t"; cerr<<ob.level<<": "<<mn<<"; "<<mx<<"; "<<ob.center<<"; "<<ob.extents<<"; "<<ob.fill<<endl;
		boxes.push_back(ob);
	}
	done:
	std::cerr<<"GlExtra_OctreeCubes::postLoad: loaded "<<boxes.size()<<" boxes."<<std::endl;
}

void GlExtra_OctreeCubes::render(){
	FOREACH(const OctreeBox& ob, boxes){
		if(ob.fill<fillRangeDraw[0] || ob.fill>fillRangeDraw[1]) continue;
		if(ob.level<levelRangeDraw[0] || ob.level>levelRangeDraw[1]) continue;
		bool doFill=(ob.fill>=fillRangeFill[0] && ob.fill<=fillRangeFill[1] && (ob.fill!=0 || !noFillZero));
		// -2: empty
		// -1: recursion limit, empty
		// 0: subdivided
		// 1: recursion limit, full
		// 2: full
		Vector3r color=(ob.fill==-2?Vector3r(1,0,0):(ob.fill==-1?Vector3r(1,1,0):(ob.fill==0?Vector3r(0,0,1):(ob.fill==1)?Vector3r(0,1,0):(ob.fill==2)?Vector3r(0,1,1):Vector3r(1,1,1))));
		glColor3v(color);
		glPushMatrix();
			glTranslatev(ob.center);
			glScalef(2*ob.extents[0],2*ob.extents[1],2*ob.extents[2]);
		 	if (doFill) glutSolidCube(1);
			else glutWireCube(1);
		glPopMatrix();
	}
}

#endif /* SUDODEM_OPENGL */
