#include<sudodem/pkg/dem/UnbalancedForceCallbacks.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Body.hpp>
#include<sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((SumIntrForcesCb)
#ifdef SUDODEM_BODY_CALLBACK
	(SumBodyForcesCb)
#endif
);

IntrCallback::FuncPtr SumIntrForcesCb::stepInit(){
	// if(scene->iter%100 != 0) return NULL;

	cerr<<"("<<(Real)force<<","<<(int)numIntr<<")";
	// reset accumulators
	force.reset(); numIntr.reset();
	// return function pointer
	return &SumIntrForcesCb::go;
}

void SumIntrForcesCb::go(IntrCallback* _self, Interaction* i){
	SumIntrForcesCb* self=static_cast<SumIntrForcesCb*>(_self);
	NormShearPhys* nsp=SUDODEM_CAST<NormShearPhys*>(i->phys.get());
	assert(nsp!=NULL); // only effective in debug mode
	Vector3r f=nsp->normalForce+nsp->shearForce;
	if(f==Vector3r::Zero()) return;
	self->numIntr+=1;
	self->force+=f.norm();
	//cerr<<"[cb#"<<i->getId1()<<"+"<<i->getId2()<<"]";
}

#ifdef SUDODEM_BODY_CALLBACK
	BodyCallback::FuncPtr SumBodyForcesCb::stepInit(){
		cerr<<"{"<<(Real)force<<","<<(int)numBodies<<",this="<<this<<",scene="<<scene<<",forces="<<&(scene->forces)<<"}";
		force.reset(); numBodies.reset(); // reset accumulators
		return &SumBodyForcesCb::go;
	}
	void SumBodyForcesCb::go(BodyCallback* _self, Body* b){
		if(b->state->blockedDOFs==State::DOF_ALL) return;
		SumBodyForcesCb* self=static_cast<SumBodyForcesCb*>(_self);
	#ifdef SUDODEM_OPENMP
		cerr<<"["<<omp_get_thread_num()<<",#"<<b->id<<",scene="<<self->scene<<"]";
	#endif
		cerr<<"[force="<<self->scene->forces.getForce(b->id)<<"]";
		self->numBodies+=1;
		//self->scene->forces.sync();
		self->force+=self->scene->forces.getForce(b->id).norm();
	}
#endif
