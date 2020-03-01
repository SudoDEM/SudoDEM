// 2007 © Václav Šmilauer <eudoxos@arcig.cz>
#include"Shop.hpp"
#include<boost/tokenizer.hpp>

#include"sudodem/core/Scene.hpp"
#include"sudodem/core/Body.hpp"
#include"sudodem/core/Interaction.hpp"

#include"sudodem/pkg/common/Aabb.hpp"
#include"sudodem/core/Clump.hpp"
#include"sudodem/pkg/common/InsertionSortCollider.hpp"

#include"sudodem/pkg/common/Box.hpp"
#include"sudodem/pkg/common/Sphere.hpp"
#include"sudodem/pkg/common/Wall.hpp"
#include"sudodem/pkg/common/ElastMat.hpp"
#include"sudodem/pkg/dem/ViscoelasticPM.hpp"


//#include"sudodem/pkg/common/Bo1_Aabb.hpp"
#include"sudodem/pkg/dem/NewtonIntegrator.hpp"
#include"sudodem/pkg/dem/Ig2_Basic_ScGeom.hpp"
#include"sudodem/pkg/dem/FrictPhys.hpp"

#include"sudodem/pkg/common/ForceResetter.hpp"

#include"sudodem/pkg/common/Dispatching.hpp"
#include"sudodem/pkg/common/InteractionLoop.hpp"
#include"sudodem/pkg/common/GravityEngines.hpp"

#include"sudodem/pkg/dem/GlobalStiffnessTimeStepper.hpp"
#include"sudodem/pkg/dem/ElasticContactLaw.hpp"

#include"sudodem/pkg/dem/ScGeom.hpp"
#include"sudodem/pkg/dem/FrictPhys.hpp"





#ifdef SUDODEM_OPENGL
	#include"sudodem/pkg/common/Gl1_NormPhys.hpp"
#endif

CREATE_LOGGER(Shop);

/*! Flip periodic cell by given number of cells.

Still broken, some interactions are missed. Should be checked.
*/

Matrix3r Shop::flipCell(const Matrix3r& _flip){
	Scene* scene=Omega::instance().getScene().get(); const shared_ptr<Cell>& cell(scene->cell); const Matrix3r& trsf(cell->trsf);
	Vector3r size=cell->getSize();
	Matrix3r flip;
	if(_flip==Matrix3r::Zero()){
		bool hasNonzero=false;
		for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
			if(i==j){ flip(i,j)=0; continue; }
			flip(i,j)=-floor(.5+trsf(i,j)/(size[j]/size[i]));
			if(flip(i,j)!=0) hasNonzero=true;
		}
		if(!hasNonzero) {LOG_TRACE("No flip necessary."); return Matrix3r::Zero();}
		LOG_DEBUG("Computed flip matrix: upper "<<flip(0,1)<<","<<flip(0,2)<<","<<flip(1,2)<<"; lower "<<flip(1,0)<<","<<flip(2,0)<<","<<flip(2,1));
	} else {
		flip=_flip;
	}

	// current cell coords of bodies
	vector<Vector3i > oldCells; oldCells.resize(scene->bodies->size());
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue; cell->wrapShearedPt(b->state->pos,oldCells[b->getId()]);
	}

	// change cell trsf here
	Matrix3r trsfInc;
	for(int i=0; i<3; i++) for(int j=0; j<3; j++){
		if(i==j) { if(flip(i,j)!=0) LOG_WARN("Non-zero diagonal term at ["<<i<<","<<j<<"] is meaningless and will be ignored."); trsfInc(i,j)=0; continue; }
		// make sure non-diagonal entries are "integers"
		if(flip(i,j)!=double(int(flip(i,j)))) LOG_WARN("Flip matrix entry "<<flip(i,j)<<" at ["<<i<<","<<j<<"] not integer?! (will be rounded)");
		trsfInc(i,j)=int(flip(i,j))*size[j]/size[i];
	}
	TRWM3MAT(cell->trsf);
	TRWM3MAT(trsfInc);
	cell->trsf+=trsfInc;
	cell->postLoad(*cell);

	// new cell coords of bodies
	vector<Vector3i > newCells; newCells.resize(scene->bodies->size());
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		cell->wrapShearedPt(b->state->pos,newCells[b->getId()]);
	}

	// remove all potential interactions
	scene->interactions->eraseNonReal();
	// adjust Interaction::cellDist for real interactions;
	FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
		Body::id_t id1=i->getId1(),id2=i->getId2();
		// this must be the same for both old and new interaction: cell2-cell1+cellDist
		// c₂-c₁+c₁₂=k; c'₂+c₁'+c₁₂'=k   (new cell coords have ')
		// c₁₂'=(c₂-c₁+c₁₂)-(c₂'-c₁')
		i->cellDist=(oldCells[id2]-oldCells[id1]+i->cellDist)-(newCells[id2]-newCells[id1]);
	}


	// force reinitialization of the collider
	bool colliderFound=false;
	FOREACH(const shared_ptr<Engine>& e, scene->engines){
		Collider* c=dynamic_cast<Collider*>(e.get());
		if(c){ colliderFound=true; c->invalidatePersistentData(); }
	}
	if(!colliderFound) LOG_WARN("No collider found while flipping cell; continuing simulation might give garbage results.");
	return flip;
}

/* Apply force on contact point to 2 bodies; the force is oriented as it applies on the first body and is reversed on the second.
 */
void Shop::applyForceAtContactPoint(const Vector3r& force, const Vector3r& contPt, Body::id_t id1, const Vector3r& pos1, Body::id_t id2, const Vector3r& pos2, Scene* scene){
	scene->forces.addForce(id1,force);
	scene->forces.addForce(id2,-force);
	scene->forces.addTorque(id1,(contPt-pos1).cross(force));
	scene->forces.addTorque(id2,-(contPt-pos2).cross(force));
}


/*! Compute sum of forces in the whole simulation and averages stiffness.

Designed for being used with periodic cell, where diving the resulting components by
areas of the cell will give average stress in that direction.

Requires all .isReal() interaction to have phys deriving from NormShearPhys.
*/
Vector3r Shop::totalForceInVolume(Real& avgIsoStiffness, Scene* _rb){
	Scene* rb=_rb ? _rb : Omega::instance().getScene().get();
	Vector3r force(Vector3r::Zero()); Real stiff=0; long n=0;
	FOREACH(const shared_ptr<Interaction>&I, *rb->interactions){
		if(!I->isReal()) continue;
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*>(I->phys.get());
		force+=Vector3r(std::abs(nsi->normalForce[0]+nsi->shearForce[0]),std::abs(nsi->normalForce[1]+nsi->shearForce[1]),std::abs(nsi->normalForce[2]+nsi->shearForce[2]));
		stiff+=(1/3.)*nsi->kn+(2/3.)*nsi->ks; // count kn in one direction and ks in the other two
		n++;
	}
	avgIsoStiffness= n>0 ? (1./n)*stiff : -1;
	return force;
}

Real Shop::unbalancedForce(bool useMaxForce, Scene* _rb){
	Scene* rb=_rb ? _rb : Omega::instance().getScene().get();
	rb->forces.sync();
	shared_ptr<NewtonIntegrator> newton;
	Vector3r gravity = Vector3r::Zero();
	FOREACH(shared_ptr<Engine>& e, rb->engines){ newton=SUDODEM_PTR_DYN_CAST<NewtonIntegrator>(e); if(newton) {gravity=newton->gravity; break;} }
	// get maximum force on a body and sum of all forces (for averaging)
	Real sumF=0,maxF=0,currF; int nb=0;int id_max;
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b || b->isClumpMember() || !b->isDynamic()) continue;
		currF=(rb->forces.getForce(b->id)+b->state->mass*gravity).norm();
		if(b->isClump() && currF==0){ // this should not happen unless the function is called by an engine whose position in the loop is before Newton (with the exception of bodies which really have null force), because clumps forces are updated in Newton. Typical triaxial loops are using such ordering unfortunately (triaxEngine before Newton). So, here we make sure that they will get correct unbalance. In the future, it is better for optimality to check unbalancedF inside scripts at the end of loops, so that this "if" is never active.
			Vector3r f(rb->forces.getForce(b->id)),m(Vector3r::Zero());
			b->shape->cast<Clump>().addForceTorqueFromMembers(b->state.get(),rb,f,m);
			currF=(f+b->state->mass*gravity).norm();
		}
		//maxF=max(currF,maxF);
		if(currF>maxF){
			maxF = currF;
			id_max = b->id;
		}
		sumF+=currF; nb++;
	}
	Real meanF=sumF;///nb;
	// get mean force on interactions
	sumF=0; //nb=0;
	Real maxCF=0;
	FOREACH(const shared_ptr<Interaction>& I, *rb->interactions){
		if(!I->isReal()) continue;
		shared_ptr<NormShearPhys> nsi=SUDODEM_PTR_CAST<NormShearPhys>(I->phys); assert(nsi);
		if(nsi->normalForce.norm()==0){continue;}//virtual contact
		Real f = (nsi->normalForce+nsi->shearForce).norm();
		sumF += f; //nb++;
		if((I->id1 == id_max) || (I->id2 == id_max)){
			maxCF += f;
		}
	}
	//sumF/=nb;
	sumF *= 2.0;
	//return (useMaxForce?maxF:meanF)/(sumF);
	return useMaxForce?maxF/maxCF:meanF/sumF;
}

Real Shop::kineticEnergy(Scene* _scene, Body::id_t* maxId){
	Scene* scene=_scene ? _scene : Omega::instance().getScene().get();
	Real ret=0.;
	Real maxE=0; if(maxId) *maxId=Body::ID_NONE;
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->isDynamic()) continue;
		const State* state(b->state.get());
		// ½(mv²+ωIω)
		Real E=0;
		if(scene->isPeriodic){
			/* Only take in account the fluctuation velocity, not the mean velocity of homothetic resize. */
			E=.5*state->mass*scene->cell->bodyFluctuationVel(state->pos,state->vel,scene->cell->velGrad).squaredNorm();
		} else {
			E=.5*(state->mass*state->vel.squaredNorm());
		}
		if(b->isAspherical()){
			Matrix3r T(state->ori);
			// the tensorial expression http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor
			// inertia tensor rotation from http://www.kwon3d.com/theory/moi/triten.html
			Matrix3r mI; mI<<state->inertia[0],0,0, 0,state->inertia[1],0, 0,0,state->inertia[2];
			//E+=.5*state->angVel.transpose().dot((T.transpose()*mI*T)*state->angVel);
			E+=.5*state->angVel.transpose().dot((T*mI*T.transpose())*state->angVel);
		}
		else { E+=0.5*state->angVel.dot(state->inertia.cwiseProduct(state->angVel));}
		if(maxId && E>maxE) { *maxId=b->getId(); maxE=E; }
		ret+=E;
	}
	return ret;
}

Vector3r Shop::momentum() {
	Vector3r ret = Vector3r::Zero();
	Scene* scene=Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Body> b, *scene->bodies){
		ret += b->state->mass * b->state->vel;
	}
	return ret;
}

Vector3r Shop::angularMomentum(Vector3r origin) {
	Vector3r ret = Vector3r::Zero();
	Scene* scene=Omega::instance().getScene().get();
	Matrix3r T,Iloc;
	FOREACH(const shared_ptr<Body> b, *scene->bodies){
		ret += (b->state->pos-origin).cross(b->state->mass * b->state->vel);
		ret += b->state->angMom;
	}
	return ret;
}


shared_ptr<FrictMat> Shop::defaultGranularMat(){
	shared_ptr<FrictMat> mat(new FrictMat);
	mat->density=2e3;
	mat->young=30e9;
	mat->poisson=.3;
	mat->frictionAngle=.5236; //30˚
	return mat;
}

/*! Create body - sphere. */
shared_ptr<Body> Shop::sphere(Vector3r center, Real radius, shared_ptr<Material> mat){
	shared_ptr<Body> body(new Body);
	body->material=mat ? mat : boost::static_pointer_cast<Material>(defaultGranularMat());
	body->state->pos=center;
	body->state->mass=4.0/3.0*Mathr::PI*radius*radius*radius*body->material->density;
	body->state->inertia=Vector3r(2.0/5.0*body->state->mass*radius*radius,2.0/5.0*body->state->mass*radius*radius,2.0/5.0*body->state->mass*radius*radius);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->shape=shared_ptr<Sphere>(new Sphere(radius));
	return body;
}

/*! Create body - box. */
shared_ptr<Body> Shop::box(Vector3r center, Vector3r extents, shared_ptr<Material> mat){
	shared_ptr<Body> body(new Body);
	body->material=mat ? mat : boost::static_pointer_cast<Material>(defaultGranularMat());
	body->state->pos=center;
	Real mass=8.0*extents[0]*extents[1]*extents[2]*body->material->density;
	body->state->mass=mass;
	body->state->inertia=Vector3r(mass*(4*extents[1]*extents[1]+4*extents[2]*extents[2])/12.,mass*(4*extents[0]*extents[0]+4*extents[2]*extents[2])/12.,mass*(4*extents[0]*extents[0]+4*extents[1]*extents[1])/12.);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->shape=shared_ptr<Box>(new Box(extents));
	return body;
}



void Shop::saveSpheresToFile(string fname){
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	ofstream f(fname.c_str());
	if(!f.good()) throw runtime_error("Unable to open file `"+fname+"'");

	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b->isDynamic()) continue;
		shared_ptr<Sphere>	intSph=SUDODEM_PTR_DYN_CAST<Sphere>(b->shape);
		if(!intSph) continue;
		const Vector3r& pos=b->state->pos;
		f<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<intSph->radius<<endl; // <<" "<<1<<" "<<1<<endl;
	}
	f.close();
}

Real Shop::getSpheresVolume(const shared_ptr<Scene>& _scene, int mask){
	const shared_ptr<Scene> scene=(_scene?_scene:Omega::instance().getScene());
	Real vol=0;
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b || !b->isDynamic()) continue;
		Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
		if((!s) or ((mask>0) and ((b->groupMask & mask)==0))) continue;
		vol += (4/3.)*Mathr::PI*pow(s->radius,3);
	}
	return vol;
}

Real Shop::getSpheresMass(const shared_ptr<Scene>& _scene, int mask){
	const shared_ptr<Scene> scene=(_scene?_scene:Omega::instance().getScene());
	Real mass=0;
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b || !b->isDynamic()) continue;
		Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
		if((!s) or ((mask>0) and ((b->groupMask & mask)==0))) continue;
		mass += b->state->mass;
	}
	return mass;
}

Real Shop::getPorosity(const shared_ptr<Scene>& _scene, Real _volume){
	const shared_ptr<Scene> scene=(_scene?_scene:Omega::instance().getScene());
	Real V;
	if(!scene->isPeriodic){
		if(_volume<=0) throw std::invalid_argument("utils.porosity must be given (positive) *volume* for aperiodic simulations.");
		V=_volume;
	} else {
		V=scene->cell->getVolume();
	}
	Real Vs=Shop::getSpheresVolume();
	return (V-Vs)/V;
}

Real Shop::getVoxelPorosity(const shared_ptr<Scene>& _scene, int _resolution, Vector3r _start,Vector3r _end){
	const shared_ptr<Scene> scene=(_scene?_scene:Omega::instance().getScene());
	if(_start==_end) throw std::invalid_argument("utils.voxelPorosity: cannot calculate porosity when start==end of the volume box.");
	if(_resolution<50)     throw std::invalid_argument("utils.voxelPorosity: it doesn't make sense to calculate porosity with voxel resolution below 50.");

	// prepare the gird, it eats a lot of memory.
	// I am not optimizing for using bits. A single byte for each cell is used.
	std::vector<std::vector<std::vector<unsigned char> > > grid;
	int S(_resolution);
	grid.resize(S);
	for(int i=0 ; i<S ; ++i) {
		grid[i].resize(S);
		for(int j=0 ; j<S ; ++j) {
			grid[i][j].resize(S,0);
		}
	}

	Vector3r start(_start), size(_end - _start);

	FOREACH(shared_ptr<Body> bi, *scene->bodies){
		if((bi)->isClump()) continue;
		const shared_ptr<Body>& b = bi;
		if ( b->isDynamic() || b->isClumpMember() ) {
			const shared_ptr<Sphere>& sphere = SUDODEM_PTR_CAST<Sphere> ( b->shape );
			Real r       = sphere->radius;
			Real rr      = r*r;
			Vector3r pos = b->state->se3.position;
			// we got sphere with radius r, at position pos.
			// and a box of size S, scaled to 'size'
			// mark cells that are iniside a sphere
			int ii(0),II(S),jj(0),JJ(S),kk(0),KK(S);
			// make sure to loop only in AABB of that sphere. No need to waste cycles outside of it.
			ii = std::max((int)((Real)(pos[0]-start[0]-r)*(Real)(S)/(Real)(size[0]))-1,0); II = std::min(ii+(int)((Real)(2.0*r)*(Real)(S)/(Real)(size[0]))+3,S);
			jj = std::max((int)((Real)(pos[1]-start[1]-r)*(Real)(S)/(Real)(size[1]))-1,0); JJ = std::min(jj+(int)((Real)(2.0*r)*(Real)(S)/(Real)(size[1]))+3,S);
			kk = std::max((int)((Real)(pos[2]-start[2]-r)*(Real)(S)/(Real)(size[2]))-1,0); KK = std::min(kk+(int)((Real)(2.0*r)*(Real)(S)/(Real)(size[2]))+3,S);
			for(int i=ii ; i<II ; ++i) {
			for(int j=jj ; j<JJ ; ++j) {
			for(int k=kk ; k<KK ; ++k) {
				Vector3r a(i,j,k);
				a = a/(Real)(S);
				Vector3r b( a[0]*size[0] , a[1]*size[1] , a[2]*size[2] );
				b = b+start;
				Vector3r c(0,0,0);
				c = pos-b;
				Real x = c[0];
				Real y = c[1];
				Real z = c[2];
				if(x*x+y*y+z*z < rr)
					grid[i][j][k]=1;
			}}}
		}
	}

	Real Vv=0;
	for(int i=0 ; i<S ; ++i ) {
		for(int j=0 ; j<S ; ++j ) {
			for(int k=0 ; k<S ; ++k ) {
				if(grid[i][j][k]==1)
					Vv += 1.0;
			}
		}
	}

	return ( std::pow(S,3) - Vv ) / std::pow(S,3);
};

vector<boost::tuple<Vector3r,Real,int> > Shop::loadSpheresFromFile(const string& fname, Vector3r& minXYZ, Vector3r& maxXYZ, Vector3r* cellSize){
	if(!boost::filesystem::exists(fname)) {
		throw std::invalid_argument(string("File with spheres `")+fname+"' doesn't exist.");
	}
	vector<boost::tuple<Vector3r,Real,int> > spheres;
	ifstream sphereFile(fname.c_str());
	if(!sphereFile.good()) throw std::runtime_error("File with spheres `"+fname+"' couldn't be opened.");
	Vector3r C;
	Real r=0;
	int clumpId=-1;
	string line;
	size_t lineNo=0;
	while(std::getline(sphereFile, line, '\n')){
		lineNo++;
		boost::tokenizer<boost::char_separator<char> > toks(line,boost::char_separator<char>(" \t"));
		vector<string> tokens; FOREACH(const string& s, toks) tokens.push_back(s);
		if(tokens.empty()) continue;
		if(tokens[0]=="##PERIODIC::"){
			if(tokens.size()!=4) throw std::invalid_argument(("Spheres file "+fname+":"+boost::lexical_cast<string>(lineNo)+" contains ##PERIODIC::, but the line is malformed.").c_str());
			if(cellSize){ *cellSize=Vector3r(boost::lexical_cast<Real>(tokens[1]),boost::lexical_cast<Real>(tokens[2]),boost::lexical_cast<Real>(tokens[3])); }
			continue;
		}
		if(tokens.size()!=5 and tokens.size()!=4) throw std::invalid_argument(("Line "+boost::lexical_cast<string>(lineNo)+" in the spheres file "+fname+" has "+boost::lexical_cast<string>(tokens.size())+" columns (must be 4 or 5).").c_str());
		C=Vector3r(boost::lexical_cast<Real>(tokens[0]),boost::lexical_cast<Real>(tokens[1]),boost::lexical_cast<Real>(tokens[2]));
		r=boost::lexical_cast<Real>(tokens[3]);
		for(int j=0; j<3; j++) { minXYZ[j]=(spheres.size()>0?min(C[j]-r,minXYZ[j]):C[j]-r); maxXYZ[j]=(spheres.size()>0?max(C[j]+r,maxXYZ[j]):C[j]+r);}
		if(tokens.size()==5)clumpId=boost::lexical_cast<int>(tokens[4]);
		spheres.push_back(boost::tuple<Vector3r,Real,int>(C,r,clumpId));
	}
	return spheres;
}

Real Shop::PWaveTimeStep(const shared_ptr<Scene> _rb){
	shared_ptr<Scene> rb=(_rb?_rb:Omega::instance().getScene());
	Real dt=std::numeric_limits<Real>::infinity();
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b || !b->material || !b->shape) continue;
		shared_ptr<ElastMat> ebp=SUDODEM_PTR_DYN_CAST<ElastMat>(b->material);
		shared_ptr<Sphere> s=SUDODEM_PTR_DYN_CAST<Sphere>(b->shape);
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
