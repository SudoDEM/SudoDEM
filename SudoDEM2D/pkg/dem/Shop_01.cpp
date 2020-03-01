// 2007 © Václav Šmilauer <eudoxos@arcig.cz>
#include"Shop.hpp"
#include<boost/tokenizer.hpp>

#include"sudodem/core/Scene.hpp"
#include"sudodem/core/Body.hpp"
#include"sudodem/core/Interaction.hpp"

#include"sudodem/pkg/common/Aabb.hpp"
#include"sudodem/core/Clump.hpp"
#include"sudodem/pkg/common/InsertionSortCollider.hpp"

//#include"sudodem/pkg/common/Box.hpp"
#include"sudodem/pkg/common/Disk.hpp"
#include"sudodem/pkg/common/ElastMat.hpp"
//#include"sudodem/pkg/dem/ViscoelasticPM.hpp"


//#include"sudodem/pkg/common/Bo1_Aabb.hpp"
#include"sudodem/pkg/dem/NewtonIntegrator.hpp"
#include"sudodem/pkg/dem/Ig2_Basic_ScGeom.hpp"
#include"sudodem/pkg/dem/FrictPhys.hpp"

#include"sudodem/pkg/common/ForceResetter.hpp"

#include"sudodem/pkg/common/Dispatching.hpp"
#include"sudodem/pkg/common/InteractionLoop.hpp"
#include"sudodem/pkg/common/GravityEngines.hpp"

//#include"sudodem/pkg/dem/GlobalStiffnessTimeStepper.hpp"
#include"sudodem/pkg/dem/ElasticContactLaw.hpp"

#include"sudodem/pkg/dem/ScGeom.hpp"
#include"sudodem/pkg/dem/FrictPhys.hpp"


CREATE_LOGGER(Shop);

Real Shop::unbalancedForce(bool useMaxForce, Scene* _rb){
	Scene* rb=_rb ? _rb : Omega::instance().getScene().get();
	rb->forces.sync();
	shared_ptr<NewtonIntegrator> newton;
	Vector2r gravity = Vector2r::Zero();
	FOREACH(shared_ptr<Engine>& e, rb->engines){ newton=SUDODEM_PTR_DYN_CAST<NewtonIntegrator>(e); if(newton) {gravity=newton->gravity; break;} }
	// get maximum force on a body and sum of all forces (for averaging)
	Real sumF=0,maxF=0,currF; int nb=0;
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b || b->isClumpMember() || !b->isDynamic()) continue;
		currF=(rb->forces.getForce(b->id)+b->state->mass*gravity).norm();
		if(b->isClump() && currF==0){ // this should not happen unless the function is called by an engine whose position in the loop is before Newton (with the exception of bodies which really have null force), because clumps forces are updated in Newton. Typical triaxial loops are using such ordering unfortunately (triaxEngine before Newton). So, here we make sure that they will get correct unbalance. In the future, it is better for optimality to check unbalancedF inside scripts at the end of loops, so that this "if" is never active.
			Vector2r f(rb->forces.getForce(b->id));
			Real m(0.0);
			b->shape->cast<Clump>().addForceTorqueFromMembers(b->state.get(),rb,f,m);
			currF=(f+b->state->mass*gravity).norm();
		}
		maxF=max(currF,maxF); sumF+=currF; nb++;
	}
	Real meanF=sumF/nb;
	// get mean force on interactions
	sumF=0; nb=0;
	FOREACH(const shared_ptr<Interaction>& I, *rb->interactions){
		if(!I->isReal()) continue;
		shared_ptr<NormShearPhys> nsi=SUDODEM_PTR_CAST<NormShearPhys>(I->phys); assert(nsi);
		sumF+=(nsi->normalForce+nsi->shearForce).norm(); nb++;
	}
	sumF/=nb;
	return (useMaxForce?maxF:meanF)/(sumF);
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

/*
Real Shop::PWaveTimeStep(const shared_ptr<Scene> _rb){
	shared_ptr<Scene> rb=(_rb?_rb:Omega::instance().getScene());
	Real dt=std::numeric_limits<Real>::infinity();
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b || !b->material || !b->shape) continue;
		shared_ptr<ElastMat> ebp=SUDODEM_PTR_DYN_CAST<ElastMat>(b->material);
		shared_ptr<Disk> s=SUDODEM_PTR_DYN_CAST<Disk>(b->shape);
		if(!ebp || !s) continue;
		Real density=b->state->mass/((4/3.)*Mathr::PI*pow(s->radius,3));
		dt=min(dt,s->radius/sqrt(ebp->young/density));
	}
	if (dt==std::numeric_limits<Real>::infinity()) {
		dt = 1.0;
		LOG_WARN("PWaveTimeStep has not found any suitable spherical body to calculate dt. dt is set to 1.0");
	}
	return dt;
}
*/
//get stress tensor and tangent operator tensor.
py::tuple Shop::getStressAndTangent2D(Real z_dim, bool symmetry){
	Scene* scene=Omega::instance().getScene().get();
	if(z_dim == 0) z_dim = 1;
	Real volume = scene->isPeriodic?scene->cell->hSize.determinant()*z_dim:1;
	Matrix2r stressTensor = Matrix2r::Zero();
	Vector6r tangent = Vector6r::Zero();//coresponding to 6 elements in the 6x6 matrix in 3D, i.e., sequentially, (0,0),(0,1),(0,5),(1,1),(1,5),(5,5)
	const bool isPeriodic = scene->isPeriodic;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		//not considering clumped particles. CAUTION:SudoDEM does not prefer clumped particles.
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );

		//Contact force
		Vector2r fn = nsi->normalForce;
		if (fn.norm() == 0)continue;
		Vector2r ft = nsi->shearForce;
		Vector2r f= fn + ft ;//here we constructe a positive stress tensor
		Real kN=nsi->kn;
		Real kT=nsi->ks;
		Vector2r n = fn.normalized();
		//cout<<"n"<<n<<endl;
		Vector2r t = ft.normalized();

		Vector2r branch=Body::byId(I->getId1(),scene)->state->pos - Body::byId(I->getId2(),scene)->state->pos;
		if (isPeriodic) branch -= scene->cell->hSize*I->cellDist.cast<Real>();

		stressTensor+=f*branch.transpose();
		Real n0b0 = n[0]*branch[0]; Real t0b0 = t[0]*branch[0];
		Real n1b0 = n[1]*branch[0]; Real t1b0 = t[1]*branch[0];
		Real n0b1 = n[0]*branch[1]; Real t0b1 = t[0]*branch[1];
		Real n1b1 = n[1]*branch[1]; Real t1b1 = t[1]*branch[1];
		tangent(0) += kN*n0b0*n0b0+kT*t0b0*t0b0;
		tangent(1) += kN*n0b0*n1b1+kT*t0b0*t1b1;
		tangent(2) += kN*n0b0*(n0b1+n1b0)*0.5+kT*t0b0*(t0b1+t1b0)*0.5;
		tangent(3) += kN*n1b1*n1b1+kT*t1b1*t1b1;
		tangent(4) += kN*n1b1*(n0b1+n1b0)*0.5+kT*t1b1*(t0b1+t1b0)*0.5;
		tangent(5) += kN*(n0b1*n0b1+n0b1*n1b0*2+n1b0*n1b0)*0.25+kT*(t0b1*t0b1+t0b1*t1b0*2+t1b0*t1b0)*0.25;


		/*

		tangent(0,0) += kN*n[0]*branch[0]*n[0]*branch[0]+kT*t[0]*branch[0]*t[0]*branch[0];
		tangent(0,1) += kN*n[0]*branch[0]*n[1]*branch[1]+kT*t[0]*branch[0]*t[1]*branch[1];
		//tangent(0,2) += kN*n[0]*branch[0]*n[2]*branch[2]+kT*t[0]*branch[0]*t[2]*branch[2];
    //tangent(0,3) += kN*n[0]*branch[0]*(n[1]*branch[2]+n[2]*branch[1])*0.5+kT*t[0]*branch[0]*(t[1]*branch[2]+t[2]*branch[1])*0.5;
    //tangent(0,4) += kN*n[0]*branch[0]*(n[0]*branch[2]+n[2]*branch[0])*0.5+kT*t[0]*branch[0]*(t[0]*branch[2]+t[2]*branch[0])*0.5;
    tangent(0,5) += kN*n[0]*branch[0]*(n[0]*branch[1]+n[1]*branch[0])*0.5+kT*t[0]*branch[0]*(t[0]*branch[1]+t[1]*branch[0])*0.5;
		tangent(1,1) += kN*n[1]*branch[1]*n[1]*branch[1]+kT*t[1]*branch[1]*t[1]*branch[1];
		//tangent(1,2) += kN*n[1]*branch[1]*n[2]*branch[2]+kT*t[1]*branch[1]*t[2]*branch[2];
    //tangent(1,3) += kN*n[1]*branch[1]*(n[1]*branch[2]+n[2]*branch[1])*0.5+kT*t[1]*branch[1]*(t[1]*branch[2]+t[2]*branch[1])*0.5;
    //tangent(1,4) += kN*n[1]*branch[1]*(n[0]*branch[2]+n[2]*branch[0])*0.5+kT*t[1]*branch[1]*(t[0]*branch[2]+t[2]*branch[0])*0.5;
    tangent(1,5) += kN*n[1]*branch[1]*(n[0]*branch[1]+n[1]*branch[0])*0.5+kT*t[1]*branch[1]*(t[0]*branch[1]+t[1]*branch[0])*0.5;
		//tangent(2,2) += kN*n[2]*branch[2]*n[2]*branch[2]+kT*t[2]*branch[2]*t[2]*branch[2];
    //tangent(2,3) += kN*n[2]*branch[2]*(n[1]*branch[2]+n[2]*branch[1])*0.5+kT*t[2]*branch[2]*(t[1]*branch[2]+t[2]*branch[1])*0.5;
    //tangent(2,4) += kN*n[2]*branch[2]*(n[0]*branch[2]+n[2]*branch[0])*0.5+kT*t[2]*branch[2]*(t[0]*branch[2]+t[2]*branch[0])*0.5;
    //tangent(2,5) += kN*n[2]*branch[2]*(n[0]*branch[1]+n[1]*branch[0])*0.5+kT*t[2]*branch[2]*(t[0]*branch[1]+t[1]*branch[0])*0.5;
    //tangent(3,3) += kN*(n[1]*branch[2]*n[1]*branch[2]+n[1]*branch[2]*n[2]*branch[1]*2+n[2]*branch[1]*n[2]*branch[1])*0.25+kT*(t[1]*branch[2]*t[1]*branch[2]+t[1]*branch[2]*t[2]*branch[1]*2+t[2]*branch[1]*t[2]*branch[1])*0.25;
    //tangent(3,4) += kN*(n[1]*branch[2]*n[0]*branch[2]+n[1]*branch[2]*n[2]*branch[0]+n[2]*branch[1]*n[0]*branch[2]+n[2]*branch[1]*n[2]*branch[0])*0.25+kT*(t[1]*branch[2]*t[0]*branch[2]+t[1]*branch[2]*t[2]*branch[0]+t[2]*branch[1]*t[0]*branch[2]+t[2]*branch[1]*t[2]*branch[0])*0.25;
    //tangent(3,5) += kN*(n[1]*branch[2]*n[0]*branch[1]+n[1]*branch[2]*n[1]*branch[0]+n[2]*branch[1]*n[0]*branch[1]+n[2]*branch[1]*n[1]*branch[0])*0.25+kT*(t[1]*branch[2]*t[0]*branch[1]+t[1]*branch[2]*t[1]*branch[0]+t[2]*branch[1]*t[0]*branch[1]+t[2]*branch[1]*t[1]*branch[0])*0.25;
    //tangent(4,4) += kN*(n[0]*branch[2]*n[0]*branch[2]+n[0]*branch[2]*n[2]*branch[0]*2+n[2]*branch[0]*n[2]*branch[0])*0.25+kT*(t[0]*branch[2]*t[0]*branch[2]+t[0]*branch[2]*t[2]*branch[0]*2+t[2]*branch[0]*t[2]*branch[0])*0.25;
    //tangent(4,5) += kN*(n[0]*branch[2]*n[0]*branch[1]+n[0]*branch[2]*n[1]*branch[0]+n[2]*branch[0]*n[0]*branch[1]+n[2]*branch[0]*n[1]*branch[0])*0.25+kT*(t[0]*branch[2]*t[0]*branch[1]+t[0]*branch[2]*t[1]*branch[0]+t[2]*branch[0]*t[0]*branch[1]+t[2]*branch[0]*t[1]*branch[0])*0.25;
    tangent(5,5) += kN*(n[0]*branch[1]*n[0]*branch[1]+n[0]*branch[1]*n[1]*branch[0]*2+n[1]*branch[0]*n[1]*branch[0])*0.25+kT*(t[0]*branch[1]*t[0]*branch[1]+t[0]*branch[1]*t[1]*branch[0]*2+t[1]*branch[0]*t[1]*branch[0])*0.25;
*/

	}
	stressTensor/=volume;
	tangent/=volume;
	return py::make_tuple(stressTensor,tangent);
}

//get stress tensor, tangent operator tensor and thermal conductivity tensor.
py::tuple Shop::getStressTangentThermal2D(Real z_dim, bool symmetry){
	Scene* scene=Omega::instance().getScene().get();
	if(z_dim == 0) z_dim = 1;
	Real volume = scene->isPeriodic?scene->cell->hSize.determinant()*z_dim:1;
	Matrix2r stressTensor = Matrix2r::Zero();
	Matrix2r thermalTensor = Matrix2r::Zero();//effective thermal conductivity based on the formula in PFC.
	Vector6r tangent = Vector6r::Zero();//coresponding to 6 elements in the 6x6 matrix in 3D, i.e., sequentially, (0,0),(0,1),(0,5),(1,1),(1,5),(5,5)
	const bool isPeriodic = scene->isPeriodic;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		//not considering clumped particles. CAUTION:SudoDEM does not prefer clumped particles.
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );

		//Contact force
		Vector2r fn = nsi->normalForce;
		if (fn.norm() == 0)continue;
		Vector2r ft = nsi->shearForce;
		Vector2r f= fn + ft ;//here we constructe a positive stress tensor
		Real kN=nsi->kn;
		Real kT=nsi->ks;
		Vector2r n = fn.normalized();
		//cout<<"n"<<n<<endl;
		Vector2r t = ft.normalized();
		Vector2r pos1 = Body::byId(I->getId1(),scene)->state->pos;
		Vector2r pos2 = Body::byId(I->getId2(),scene)->state->pos;
		Vector2r branch=pos1 - pos2;
		if (isPeriodic) branch -= scene->cell->hSize*I->cellDist.cast<Real>();

		stressTensor+=f*branch.transpose();
		//thermal conductivity tensor
		Vector2r cp = SUDODEM_PTR_CAST<SuperellipseGeom>(I->geom)->contactPoint;
		thermalTensor += (cp-pos1)*(cp-pos1).transpose()/(cp-pos1).norm();
		thermalTensor += (cp-pos2)*(cp-pos2).transpose()/(cp-pos2).norm();
		//thermalTensor += branch.norm()*n*n.transpose();
		Real n0b0 = n[0]*branch[0]; Real t0b0 = t[0]*branch[0];
		Real n1b0 = n[1]*branch[0]; Real t1b0 = t[1]*branch[0];
		Real n0b1 = n[0]*branch[1]; Real t0b1 = t[0]*branch[1];
		Real n1b1 = n[1]*branch[1]; Real t1b1 = t[1]*branch[1];
		tangent(0) += kN*n0b0*n0b0+kT*t0b0*t0b0;
		tangent(1) += kN*n0b0*n1b1+kT*t0b0*t1b1;
		tangent(2) += kN*n0b0*(n0b1+n1b0)*0.5+kT*t0b0*(t0b1+t1b0)*0.5;
		tangent(3) += kN*n1b1*n1b1+kT*t1b1*t1b1;
		tangent(4) += kN*n1b1*(n0b1+n1b0)*0.5+kT*t1b1*(t0b1+t1b0)*0.5;
		tangent(5) += kN*(n0b1*n0b1+n0b1*n1b0*2+n1b0*n1b0)*0.25+kT*(t0b1*t0b1+t0b1*t1b0*2+t1b0*t1b0)*0.25;

	}
	stressTensor/=volume;
	thermalTensor/=volume;//note: the thermal resistance per unit length has not been divided yet.
	tangent/=volume;
	return py::make_tuple(stressTensor,tangent,thermalTensor);
}
