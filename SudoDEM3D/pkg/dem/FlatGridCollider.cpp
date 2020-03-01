// 2010 © Václav Šmilauer <eudoxos@arcig.cz>
#include<sudodem/pkg/dem/FlatGridCollider.hpp>
#include<sudodem/core/Scene.hpp>

#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/dem/NewtonIntegrator.hpp>
//#include<sudodem/pkg/common/Facet.hpp>

SUDODEM_PLUGIN((FlatGridCollider));
CREATE_LOGGER(FlatGridCollider);

bool FlatGridCollider::isActivated(){
	// keep interactions trequested for deletion as potential (forget removal requests)
// 	scene->interactions->clearPendingErase();
	if(!newton) return true;
	// handle verlet distance
	fastestBodyMaxDist+=sqrt(newton->maxVelocitySq)*scene->dt;
	if(fastestBodyMaxDist>verletDist) return true;
	return false;
}

void FlatGridCollider::action(){
	if(!newton){
		FOREACH(const shared_ptr<Engine>& e, scene->engines){ newton=SUDODEM_PTR_DYN_CAST<NewtonIntegrator>(e); if(newton) break; }
		if(!newton){ throw runtime_error("FlatGridCollider: Unable to find NewtonIntegrator in engines."); }
	}
	fastestBodyMaxDist=0;
	// make interaction loop delete unseen potential interactions
	scene->interactions->iterColliderLastRun=scene->iter;
	// adjust grid if necessary
	updateGrid();
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue; // deleted bodies
		updateBodyCells(b);
	}
	updateCollisions();
}


// called from the ctor
void FlatGridCollider::initIndices(){
	sphereIdx=facetIdx=wallIdx=boxIdx=-1;
	sphereIdx=Sphere::getClassIndexStatic();
	LOG_DEBUG("sphereIdx="<<sphereIdx);
}

void FlatGridCollider::updateGrid(){
	// checks
	if(step<=0) throw std::runtime_error("FlatGridCollider::step must be positive.");
	Vector3r aabbSize=aabbMax-aabbMin; if(aabbSize[0]<=0 || aabbSize[1]<=0 || aabbSize[2]<=0) throw std::runtime_error("FlatGridCollider::{aabbMin,aabbMax} must give positive volume.");
	// compute new size
	grid.step=step;
	grid.mn=aabbMin;
	for(int i=0;i<3;i++)grid.size[i]=(int)ceil((aabbMax[i]-aabbMin[i])/step);
	grid.mx=aabbMin+Vector3r(grid.size[0]*step,grid.size[1]*step,grid.size[2]*step);
	size_t len=grid.size[0]*grid.size[1]*grid.size[2];
	// reset grid data
	grid.data.clear(); grid.data.resize(len); // delete and recreate; could be optimized somehow
	LOG_DEBUG("New grid "<<grid.size[0]<<"×"<<grid.size[1]<<"×"<<grid.size[2]<<"="<<len<<" cells, step "<<step<<"; spans ("<<grid.mn<<")--("<<grid.mx<<").");
}

void FlatGridCollider::updateBodyCells(const shared_ptr<Body>& b){
	if(!b->shape) return; const Shape* shape(b->shape.get());
	// Sphere
	if(shape->getClassIndex()==sphereIdx){
		Real r=static_cast<const Sphere*>(shape)->radius+verletDist;
		const Vector3r& C=b->state->pos;
		// create "bounding cells" and traverse them one by one; integrized coords can be _outside_ grid, they will be forced to closest cells before insertion
		Vector3i cMn=grid.pt2int(C-Vector3r(r,r,r)), cMx=grid.pt2int(C+Vector3r(r,r,r)), cC=grid.pt2int(C), cPt;
		LOG_TRACE("Sphere #"<<b->id<<", C="<<C<<", cMn="<<cMn<<", cC="<<cC<<", cMx="<<cMx);
		for(cPt[0]=cMn[0]; cPt[0]<=cMx[0]; cPt[0]++) for(cPt[1]=cMn[1]; cPt[1]<=cMx[1]; cPt[1]++) for(cPt[2]=cMn[2]; cPt[2]<=cMx[2]; cPt[2]++){
			// find closest cell point (to cC); keep coordinate in same line (cPt[i]==cC[i]); take upper/lower point in lower/upper lines (cPt[i]<cC[i])
			Vector3r ccp;
			for(int i=0;i<3;i++) ccp[i]=(cPt[i]==cC[i] ? C[i] : (grid.mn[i]+grid.step*(cPt[i]+(cPt[i]<cC[i] ? 1 : 0))));
			if((C-ccp).squaredNorm()<=r*r){ // closest cell point it inside the spehre; add the sphere to cell at cell position
				Vector3i cPtIn(grid.fitGrid(cPt));
				// perhaps slower, but inserts each body only once into the cell (meaningful if outside grid: multiple integer coords coolapse in one cell)
				{ Grid::idVector& vv=grid(cPtIn); if(vv.size()==0 || *(vv.rbegin())!=b->id) vv.push_back(b->id);}
				//grid(cPtIn).push_back(b->id);
				LOG_TRACE("Added sphere #"<<b->id<<" to cell ("<<cPtIn<<")←["<<cPt<<"]; center ("<<C<<"), closest ("<<ccp<<"), dist "<<(C-ccp).norm());
			}
		}
		return;
	}
	throw std::runtime_error("FlatGridCollider::updateBodyCells does not handle Shape of type "+shape->getClassName()+"!");
}

void FlatGridCollider::updateCollisions(){
	shared_ptr<InteractionContainer>& intrs=scene->interactions;
	const long& iter=scene->iter;
	// create interactions for all combinations of bodies within one cell
	FOREACH(const Grid::idVector& v, grid.data){
		size_t sz=v.size();
		for(size_t i=0; i<sz; i++) for(size_t j=i+1; j<sz; j++){
			Body::id_t id1=v[i], id2=v[j];
			if(id1==id2) continue; // at grid boundary, it is possible to have one body more times in one cell
			const shared_ptr<Interaction>& I=intrs->find(id1,id2);
			if(I){ I->iterLastSeen=iter; continue; }
			// no interaction yet
			if(!Collider::mayCollide(Body::byId(id1,scene).get(),Body::byId(id2,scene).get())) continue;
			intrs->insert(shared_ptr<Interaction>(new Interaction(id1,id2)));
			LOG_TRACE("Created new interaction #"<<id1<<"+#"<<id2);
		}
	}
}
