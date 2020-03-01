// 2007 © Václav Šmilauer <eudoxos@arcig.cz>
#include"Shop.hpp"

#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Body.hpp>
#include<sudodem/core/Interaction.hpp>

#include<sudodem/pkg/common/Aabb.hpp>
#include<sudodem/core/Clump.hpp>
#include<sudodem/pkg/common/InsertionSortCollider.hpp>

#include<sudodem/pkg/common/Box.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/dem/ViscoelasticPM.hpp>


//#include<sudodem/pkg/common/Bo1_Aabb.hpp>
#include<sudodem/pkg/dem/NewtonIntegrator.hpp>
#include<sudodem/pkg/dem/Ig2_Basic_ScGeom.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>

#include<sudodem/pkg/common/ForceResetter.hpp>

#include<sudodem/pkg/common/Dispatching.hpp>
#include<sudodem/pkg/common/InteractionLoop.hpp>
#include<sudodem/pkg/common/GravityEngines.hpp>

#include<sudodem/pkg/dem/GlobalStiffnessTimeStepper.hpp>
#include<sudodem/pkg/dem/ElasticContactLaw.hpp>

#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/pkg/dem/FrictPhys.hpp>





#ifdef SUDODEM_OPENGL
	#include<sudodem/pkg/common/Gl1_NormPhys.hpp>
#endif

CREATE_LOGGER(Shop);

/*! Flip periodic cell by given number of cells.

Still broken, some interactions are missed. Should be checked.
*/

/* Detremination of time step as according to Rayleigh wave speed of force propagation (see Thornton 2000, ref. MillerPursey_1955) */
Real Shop::RayleighWaveTimeStep(const shared_ptr<Scene> _rb){
	shared_ptr<Scene> rb=(_rb?_rb:Omega::instance().getScene());
	Real dt=std::numeric_limits<Real>::infinity();
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b || !b->material || !b->shape) continue;

		shared_ptr<ElastMat> ebp=SUDODEM_PTR_DYN_CAST<ElastMat>(b->material);
		shared_ptr<Sphere> s=SUDODEM_PTR_DYN_CAST<Sphere>(b->shape);
		if(!ebp || !s) continue;

		Real density=b->state->mass/((4/3.)*Mathr::PI*pow(s->radius,3));
		Real ShearModulus=ebp->young/(2.*(1+ebp->poisson));
		Real lambda=0.1631*ebp->poisson+0.876605;
		dt=min(dt,Mathr::PI*s->radius/lambda*sqrt(density/ShearModulus));
	}
	return dt;
}

/* Project 3d point into 2d using spiral projection along given axis;
 * the returned tuple is
 *
 *  (height relative to the spiral, distance from axis, theta )
 *
 * dH_dTheta is the inclination of the spiral (height increase per radian),
 * theta0 is the angle for zero height (by given axis).
 */
boost::tuple<Real,Real,Real> Shop::spiralProject(const Vector3r& pt, Real dH_dTheta, int axis, Real periodStart, Real theta0){
	int ax1=(axis+1)%3,ax2=(axis+2)%3;
	Real r=sqrt(pow(pt[ax1],2)+pow(pt[ax2],2));
	Real theta;
	if(r>Mathr::ZERO_TOLERANCE){
		theta=acos(pt[ax1]/r);
		if(pt[ax2]<0) theta=Mathr::TWO_PI-theta;
	}
	else theta=0;
	Real hRef=dH_dTheta*(theta-theta0);
	long period;
	if(isnan(periodStart)){
		Real h=Shop::periodicWrap(pt[axis]-hRef,hRef-Mathr::PI*dH_dTheta,hRef+Mathr::PI*dH_dTheta,&period);
		return boost::make_tuple(r,h,theta);
	}
	else{
		// Real hPeriodStart=(periodStart-theta0)*dH_dTheta;
		//TRVAR4(hPeriodStart,periodStart,theta0,theta);
		//Real h=Shop::periodicWrap(pt[axis]-hRef,hPeriodStart,hPeriodStart+2*Mathr::PI*dH_dTheta,&period);
		theta=Shop::periodicWrap(theta,periodStart,periodStart+2*Mathr::PI,&period);
		Real h=pt[axis]-hRef+period*2*Mathr::PI*dH_dTheta;
		//TRVAR3(pt[axis],pt[axis]-hRef,period);
		//TRVAR2(h,theta);
		return boost::make_tuple(r,h,theta);
	}
}

shared_ptr<Interaction> Shop::createExplicitInteraction(Body::id_t id1, Body::id_t id2, bool force){
	IGeomDispatcher* geomMeta=NULL;
	IPhysDispatcher* physMeta=NULL;
	shared_ptr<Scene> rb=Omega::instance().getScene();
	if(rb->interactions->find(Body::id_t(id1),Body::id_t(id2))!=0) throw runtime_error(string("Interaction #")+boost::lexical_cast<string>(id1)+"+#"+boost::lexical_cast<string>(id2)+" already exists.");
	FOREACH(const shared_ptr<Engine>& e, rb->engines){
		if(!geomMeta) { geomMeta=dynamic_cast<IGeomDispatcher*>(e.get()); if(geomMeta) continue; }
		if(!physMeta) { physMeta=dynamic_cast<IPhysDispatcher*>(e.get()); if(physMeta) continue; }
		InteractionLoop* id(dynamic_cast<InteractionLoop*>(e.get()));
		if(id){ geomMeta=id->geomDispatcher.get(); physMeta=id->physDispatcher.get(); }
		if(geomMeta&&physMeta){break;}
	}
	if(!geomMeta) throw runtime_error("No IGeomDispatcher in engines or inside InteractionLoop.");
	if(!physMeta) throw runtime_error("No IPhysDispatcher in engines or inside InteractionLoop.");
	shared_ptr<Body> b1=Body::byId(id1,rb), b2=Body::byId(id2,rb);
	if(!b1) throw runtime_error(("No body #"+boost::lexical_cast<string>(id1)).c_str());
	if(!b2) throw runtime_error(("No body #"+boost::lexical_cast<string>(id2)).c_str());
	shared_ptr<Interaction> i=geomMeta->explicitAction(b1,b2,/*force*/force);
	assert(force && i);
	if(!i) return i;
	physMeta->explicitAction(b1->material,b2->material,i);
	i->iterMadeReal=rb->iter;
	rb->interactions->insert(i);
	return i;
}

Vector3r Shop::inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2)
{
	return v0+((v2-v0)*(v1-v0).norm()+(v1-v0)*(v2-v0).norm())/((v1-v0).norm()+(v2-v1).norm()+(v0-v2).norm());
}

void Shop::getViscoelasticFromSpheresInteraction( Real tc, Real en, Real es, shared_ptr<ViscElMat> b)
{
    throw runtime_error("Setting parameters in ViscoElastic model is changed. You do not need to use getViscoelasticFromSpheresInteraction function any more, because this functino is deprecated. You need to set the parameters tc, en and es directly in material properties. Please, update your scripts. How to do it you can see in the following commit https://github.com/sudodem/trunk/commit/1987c2febdb8a6ce2d27f2dc1bb29df0dc5f686e");
}

/* This function is copied almost verbatim from scientific python, module Visualization, class ColorScale
 *
 */
Vector3r Shop::scalarOnColorScale(Real x, Real xmin, Real xmax){
	Real xnorm=min((Real)1.,max((x-xmin)/(xmax-xmin),(Real)0.));
	if(xnorm<.25) return Vector3r(0,4.*xnorm,1);
	if(xnorm<.5)  return Vector3r(0,1,1.-4.*(xnorm-.25));
	if(xnorm<.75) return Vector3r(4*(xnorm-.5),1.,0);
	return Vector3r(1,1-4*(xnorm-.75),0);
}

/* Wrap floating point number into interval (x0,x1〉such that it is shifted
 * by integral number of the interval range. If given, *period will hold
 * this number. The wrapped value is returned.
 */
Real Shop::periodicWrap(Real x, Real x0, Real x1, long* period){
	Real xNorm=(x-x0)/(x1-x0);
	Real xxNorm=xNorm-floor(xNorm);
	if(period) *period=(long)floor(xNorm);
	return x0+xxNorm*(x1-x0);
}

void Shop::getStressForEachBody(vector<Shop::bodyState>& bodyStates){
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	bodyStates.resize(scene->bodies->size());
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		Vector3r normalStress,shearStress;
		if(!I->isReal()) continue;

		const FrictPhys* physFP = SUDODEM_CAST<FrictPhys*>(I->phys.get());
		ScGeom* geomScG=SUDODEM_CAST<ScGeom*>(I->geom.get());

		const Body::id_t id1=I->getId1(), id2=I->getId2();

		if((physFP) and (geomScG)){
			Real minRad=(geomScG->radius1<=0?geomScG->radius2:(geomScG->radius2<=0?geomScG->radius1:min(geomScG->radius1,geomScG->radius2)));
			Real crossSection=Mathr::PI*pow(minRad,2);

			normalStress=((1./crossSection)*geomScG->normal.dot(physFP->normalForce))*geomScG->normal;
			for(int i=0; i<3; i++){
				int ix1=(i+1)%3,ix2=(i+2)%3;
				shearStress[i]=geomScG->normal[ix1]*physFP->shearForce[ix1]+geomScG->normal[ix2]*physFP->shearForce[ix2];
				shearStress[i]/=crossSection;
			}
			bodyStates[id1].normStress+=normalStress;
			bodyStates[id2].normStress+=normalStress;
			bodyStates[id1].shearStress+=shearStress;
			bodyStates[id2].shearStress+=shearStress;
		}
	}
}

/* Return the stress tensor decomposed in 2 contributions, from normal and shear forces.
The formulation follows the [Thornton2000]_ article
"Numerical simulations of deviatoric shear deformation of granular media", eq (3) and (4)
 */
py::tuple Shop::normalShearStressTensors(bool compressionPositive, bool splitNormalTensor, Real thresholdForce){

	//*** Stress tensor split into shear and normal contribution ***/
	Scene* scene=Omega::instance().getScene().get();
	if (!scene->isPeriodic){ throw runtime_error("Can't compute stress of periodic cell in aperiodic simulation."); }
	Matrix3r sigN(Matrix3r::Zero()),sigT(Matrix3r::Zero());
	//const Matrix3r& cellHsize(scene->cell->Hsize);   //Disabled because of warning.
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(I->geom.get());
		NormShearPhys* phys=SUDODEM_CAST<NormShearPhys*>(I->phys.get());
		const Vector3r& n=geom->normal;
		// if compression has positive sign, we need to change sign for both normal and shear force
		// make copy to Fs since shearForce is used at multiple places (less efficient, but avoids confusion)
		Vector3r Fs=(compressionPositive?-1:1)*phys->shearForce;
		Real N=(compressionPositive?-1:1)*phys->normalForce.dot(n), T=Fs.norm();
		bool hasShear=(T>0);
		Vector3r t; if(hasShear) t=Fs/T;
		// Real R=(Body::byId(I->getId2(),scene)->state->pos+cellHsize*I->cellDist.cast<Real>()-Body::byId(I->getId1(),scene)->state->pos).norm();
		Real R=.5*(geom->refR1+geom->refR2);
		for(int i=0; i<3; i++) for(int j=i; j<3; j++){
			sigN(i,j)+=R*N*n[i]*n[j];
			if(hasShear) sigT(i,j)+=R*T*n[i]*t[j];
		}
	}
	Real vol=scene->cell->getVolume();
	sigN*=2/vol; sigT*=2/vol;
	// fill terms under the diagonal
	sigN(1,0)=sigN(0,1); sigN(2,0)=sigN(0,2); sigN(2,1)=sigN(1,2);
	sigT(1,0)=sigT(0,1); sigT(2,0)=sigT(0,2); sigT(2,1)=sigT(1,2);

	// *** Normal stress tensor split into two parts according to subnetworks of strong and weak forces (or other distinction if a threshold value for the force is assigned) ***/
	Real Fmean(0); Matrix3r f, fs, fw;
	fabricTensor(Fmean,f,fs,fw,false,compressionPositive,NaN);
	Matrix3r sigNStrong(Matrix3r::Zero()), sigNWeak(Matrix3r::Zero());
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(I->geom.get());
		NormShearPhys* phys=SUDODEM_CAST<NormShearPhys*>(I->phys.get());
		const Vector3r& n=geom->normal;
		Real N=(compressionPositive?-1:1)*phys->normalForce.dot(n);
		// Real R=(Body::byId(I->getId2(),scene)->state->pos+cellHsize*I->cellDist.cast<Real>()-Body::byId(I->getId1(),scene)->state->pos).norm();
		Real R=.5*(geom->refR1+geom->refR2);
		Real Fsplit=(!isnan(thresholdForce))?thresholdForce:Fmean;
		if (compressionPositive?(N<Fsplit):(N>Fsplit)){
			for(int i=0; i<3; i++) for(int j=i; j<3; j++){
				sigNStrong(i,j)+=R*N*n[i]*n[j];}
		}
		else{
			for(int i=0; i<3; i++) for(int j=i; j<3; j++){
				sigNWeak(i,j)+=R*N*n[i]*n[j];}
		}
	}
	sigNStrong*=2/vol; sigNWeak*=2/vol;
	// fill terms under the diagonal
	sigNStrong(1,0)=sigNStrong(0,1); sigNStrong(2,0)=sigNStrong(0,2); sigNStrong(2,1)=sigNStrong(1,2);
	sigNWeak(1,0)=sigNWeak(0,1); sigNWeak(2,0)=sigNWeak(0,2); sigNWeak(2,1)=sigNWeak(1,2);

	/// tensile forces are taken as positive!
	if (splitNormalTensor){return py::make_tuple(sigNStrong,sigNWeak);} // return strong-weak or tensile-compressive parts of the stress tensor (only normal part)
	return py::make_tuple(sigN,sigT); // return normal and shear components
}

/* Return the fabric tensor as according to [Satake1982]. */
/* as side-effect, set Gl1_NormShear::strongWeakThresholdForce */
void Shop::fabricTensor(Real& Fmean, Matrix3r& fabric, Matrix3r& fabricStrong, Matrix3r& fabricWeak, bool splitTensor, bool revertSign, Real thresholdForce){
	Scene* scene=Omega::instance().getScene().get();
	if (!scene->isPeriodic){ throw runtime_error("Can't compute fabric tensor of periodic cell in aperiodic simulation."); }

	// *** Fabric tensor ***/
	fabric=Matrix3r::Zero();
	int count=0; // number of interactions
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(I->geom.get());
		const Vector3r& n=geom->normal;
		for(int i=0; i<3; i++) for(int j=i; j<3; j++){
			fabric(i,j)+=n[i]*n[j];
		}
		count++;
	}
	// fill terms under the diagonal
	fabric(1,0)=fabric(0,1); fabric(2,0)=fabric(0,2); fabric(2,1)=fabric(1,2);
	fabric/=count;

	// *** Average contact force ***/
	// calculate average contact force
	Fmean=0; // initialize
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(I->geom.get());
		NormShearPhys* phys=SUDODEM_CAST<NormShearPhys*>(I->phys.get());
		const Vector3r& n=geom->normal;
		Real f=(revertSign?-1:1)*phys->normalForce.dot(n);
		//Real f=phys->normalForce.norm();
		Fmean+=f;
	}
	Fmean/=count;

	#ifdef SUDODEM_OPENGL
		Gl1_NormPhys::maxWeakFn=Fmean;
	#endif

	// *** Weak and strong fabric tensors ***/
	// evaluate two different parts of the fabric tensor
	// make distinction between strong and weak network of contact forces
	fabricStrong=Matrix3r::Zero();
	fabricWeak=Matrix3r::Zero();
	int nStrong(0), nWeak(0); // number of strong and weak contacts respectively
	if (!splitTensor & !isnan(thresholdForce)) {LOG_WARN("The bool splitTensor should be set to True if you specified a threshold value for the contact force, otherwise the function will return only the fabric tensor and not the two separate contributions.");}
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(I->geom.get());
		NormShearPhys* phys=SUDODEM_CAST<NormShearPhys*>(I->phys.get());
		const Vector3r& n=geom->normal;
		Real  f=(revertSign?-1:1)*phys->normalForce.dot(n);
		// slipt the tensor according to the mean contact force or a threshold value if this is given
		Real Fsplit=(!isnan(thresholdForce))?thresholdForce:Fmean;
		if (revertSign?(f<Fsplit):(f>Fsplit)){ // reminder: forces are compared with their sign
			for(int i=0; i<3; i++) for(int j=i; j<3; j++){
				fabricStrong(i,j)+=n[i]*n[j];
			}
			nStrong++;
		}
		else{
			for(int i=0; i<3; i++) for(int j=i; j<3; j++){
				fabricWeak(i,j)+=n[i]*n[j];
			}
			nWeak++;
		}
	}
	// fill terms under the diagonal
	fabricStrong(1,0)=fabricStrong(0,1); fabricStrong(2,0)=fabricStrong(0,2); fabricStrong(2,1)=fabricStrong(1,2);
	fabricWeak(1,0)=fabricWeak(0,1); fabricWeak(2,0)=fabricWeak(0,2); fabricWeak(2,1)=fabricWeak(1,2);
	fabricStrong/=nStrong;
	fabricWeak/=nWeak;

	// *** Compute total fabric tensor from the two tensors above ***/
	Matrix3r fabricTot(Matrix3r::Zero());
	int q(0);
	if(!count){ // compute only if there are some interactions
		q=nStrong*1./count;
		fabricTot=(1-q)*fabricWeak+q*fabricStrong;
	}
}

py::tuple Shop::fabricTensor(bool splitTensor, bool revertSign, Real thresholdForce){
	Real Fmean; Matrix3r fabric, fabricStrong, fabricWeak;
	fabricTensor(Fmean,fabric,fabricStrong,fabricWeak,splitTensor,revertSign,thresholdForce);
	// returns fabric tensor or alternatively the two distinct contributions according to strong and weak subnetworks (or, if thresholdForce is specified, the distinction is made according to that value and not the mean one)
	if (!splitTensor){return py::make_tuple(fabric);}
	else{return py::make_tuple(fabricStrong,fabricWeak);}
}

Matrix3r Shop::getStress(Real volume){
	Scene* scene=Omega::instance().getScene().get();
	if (volume==0) volume = scene->isPeriodic?scene->cell->hSize.determinant():1;
	Matrix3r stressTensor = Matrix3r::Zero();
	const bool isPeriodic = scene->isPeriodic;
	FOREACH(const shared_ptr<Interaction>&I, *scene->interactions){
		if (!I->isReal()) continue;
		shared_ptr<Body> b1 = Body::byId(I->getId1(),scene);
		shared_ptr<Body> b2 = Body::byId(I->getId2(),scene);
		Vector3r pos1 = b1->state->pos;
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );
		/*if (b1->shape->getClassName()=="Wall"){
			if(b2->shape->getClassName()=="Sphere"){//for ball-wall contact
				const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b2->shape);

				pos1 = b2->state->pos-A->radius*nsi->normalForce.normalized();
		}else{//other

		}
	}*/

		Vector3r branch= pos1 - b2->state->pos;
		if (isPeriodic) branch-= scene->cell->hSize*I->cellDist.cast<Real>();
		stressTensor += (nsi->normalForce+nsi->shearForce)*branch.transpose();
	}
	return stressTensor/volume;
}


void Shop::getStressLWForEachBody(vector<Matrix3r>& bStresses){
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	bStresses.resize(scene->bodies->size());
	for (size_t k=0;k<scene->bodies->size();k++) bStresses[k]=Matrix3r::Zero();
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		GenericSpheresContact* geom=SUDODEM_CAST<GenericSpheresContact*>(I->geom.get());
		NormShearPhys* phys=SUDODEM_CAST<NormShearPhys*>(I->phys.get());
		Vector3r f=phys->normalForce+phys->shearForce;
		//Sum f_i*l_j/V for each contact of each particle
		bStresses[I->getId1()]-=(3.0/(4.0*Mathr::PI*pow(geom->refR1,3)))*f*((geom->contactPoint-Body::byId(I->getId1(),scene)->state->pos).transpose());
		if (!scene->isPeriodic)
			bStresses[I->getId2()]+=(3.0/(4.0*Mathr::PI*pow(geom->refR2,3)))*f*((geom->contactPoint- (Body::byId(I->getId2(),scene)->state->pos)).transpose());
		else
			bStresses[I->getId2()]+=(3.0/(4.0*Mathr::PI*pow(geom->refR2,3)))*f* ((geom->contactPoint- (Body::byId(I->getId2(),scene)->state->pos + (scene->cell->hSize*I->cellDist.cast<Real>()))).transpose());
	}
}

py::list Shop::getStressLWForEachBody(){
	py::list ret;
	vector<Matrix3r> bStresses;
	getStressLWForEachBody(bStresses);
	FOREACH(const Matrix3r& m, bStresses) ret.append(m);
	return ret;
// 	return py::make_tuple(bStresses);
}

void Shop::calm(const shared_ptr<Scene>& _scene, int mask){
	const shared_ptr<Scene> scene=(_scene?_scene:Omega::instance().getScene());
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b || !b->isDynamic()) continue;
		if(((mask>0) and ((b->groupMask & mask)==0))) continue;
		b->state->vel=Vector3r::Zero();
		b->state->angVel=Vector3r::Zero();
		b->state->angMom=Vector3r::Zero();
	}
}

py::list Shop::getBodyIdsContacts(Body::id_t bodyID) {
	py::list ret;
	if (bodyID < 0) {
		throw std::logic_error("BodyID should be a positive value!");
	}

	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	const shared_ptr<Body>& b = Body::byId(bodyID,scene);

	for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {
		ret.append((*it).first);
	}
	return ret;
}

void Shop::setContactFriction(Real angleRad){
	Scene* scene = Omega::instance().getScene().get();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	FOREACH(const shared_ptr<Body>& b,*scene->bodies){
		if(b->isClump()) continue;
		if (b->isDynamic())
		SUDODEM_PTR_CAST<FrictMat> (b->material)->frictionAngle = angleRad;
	}
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions){
		if (!ii->isReal()) continue;
		const shared_ptr<FrictMat>& sdec1 = SUDODEM_PTR_CAST<FrictMat>((*bodies)[(Body::id_t) ((ii)->getId1())]->material);
		const shared_ptr<FrictMat>& sdec2 = SUDODEM_PTR_CAST<FrictMat>((*bodies)[(Body::id_t) ((ii)->getId2())]->material);
		//FIXME - why dynamic_cast fails here?
		FrictPhys* contactPhysics = SUDODEM_CAST<FrictPhys*>((ii)->phys.get());
		const Real& fa = sdec1->frictionAngle;
		const Real& fb = sdec2->frictionAngle;
		contactPhysics->tangensOfFrictionAngle = std::tan(std::min(fa,fb));
	}
}

void Shop::growParticles(Real multiplier, bool updateMass, bool dynamicOnly)
{
	Scene* scene = Omega::instance().getScene().get();
	const int sphereIdx = Sphere::getClassIndexStatic();
	FOREACH(const shared_ptr<Body>& b,*scene->bodies){
		if (dynamicOnly && !b->isDynamic()) continue;
		//We grow only spheres and clumps
		if(b->isClump() or sphereIdx == b->shape->getClassIndex()){
			if (updateMass) {b->state->mass*=pow(multiplier,3); b->state->inertia*=pow(multiplier,5);}
			// for clumps we updated inertia, nothing else to do
			if (b->isClump()) continue;
			// for spheres, we update radius
			(SUDODEM_CAST<Sphere*> (b->shape.get()))->radius *= multiplier;
			// and if they are clump members,clump volume variation with homothetic displacement of all members
			if (b->isClumpMember()) b->state->pos += (multiplier-1) * (b->state->pos - Body::byId(b->clumpId, scene)->state->pos);
		}
	}
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions){
		int ci1=(*(scene->bodies))[ii->getId1()]->shape->getClassIndex();
		int ci2=(*(scene->bodies))[ii->getId2()]->shape->getClassIndex();
		if (ii->isReal()) {
			GenericSpheresContact* contact = SUDODEM_CAST<GenericSpheresContact*>(ii->geom.get());
			if ((!dynamicOnly || (*(scene->bodies))[ii->getId1()]->isDynamic()) && ci1==sphereIdx)
				contact->refR1 = SUDODEM_CAST<Sphere*>((* (scene->bodies))[ii->getId1()]->shape.get())->radius;
			if ((!dynamicOnly || (*(scene->bodies))[ii->getId2()]->isDynamic()) && ci2==sphereIdx)
				contact->refR2 = SUDODEM_CAST<Sphere*>((* (scene->bodies))[ii->getId2()]->shape.get())->radius;
			const shared_ptr<FrictPhys>& contactPhysics = SUDODEM_PTR_CAST<FrictPhys>(ii->phys);
			contactPhysics->kn*=multiplier; contactPhysics->ks*=multiplier;
		}

	}
}

void Shop::growParticle(Body::id_t bodyID, Real multiplier, bool updateMass)
{
	const shared_ptr<Body>& b = Body::byId(bodyID);
	Real& rad = SUDODEM_CAST<Sphere*>(b->shape.get())->radius;
	rad *= multiplier;
	if (updateMass) {b->state->mass*=pow(multiplier,3); b->state->inertia*=pow(multiplier,5);}
	for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		if(!(*it).second->isReal()) continue;
		GenericSpheresContact* contact = SUDODEM_CAST<GenericSpheresContact*>((*it).second->geom.get());
		if (bodyID==it->second->getId1()) contact->refR1 = rad;
		else contact->refR2 = rad;
	}
}
