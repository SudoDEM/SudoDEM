// 2010 © Václav Šmilauer <eudoxos@arcig.cz>

#include "Scene.hpp"
#include "Body.hpp"
#include "BodyContainer.hpp"
#include "Clump.hpp"
#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif


CREATE_LOGGER(BodyContainer);

void BodyContainer::clear(){
	body.clear();
}

Body::id_t BodyContainer::insert(shared_ptr<Body>& b){
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	b->iterBorn=scene->iter;
	b->timeBorn=scene->time;
	b->id=body.size();
	scene->doSort = true;
	body.push_back(b);
	// Notify ForceContainer about new id
	scene->forces.addMaxId(b->id);
	return b->id;
}

bool BodyContainer::erase(Body::id_t id, bool eraseClumpMembers){//default is false (as before)
	if(!body[id]) return false;
	const shared_ptr<Body>& b=Body::byId(id);
	if ((b) and (b->isClumpMember())) {
		const shared_ptr<Body> clumpBody=Body::byId(b->clumpId);
		const shared_ptr<Clump> clump=SUDODEM_PTR_CAST<Clump>(clumpBody->shape);
		Clump::del(clumpBody, b);
		if (clump->members.size()==0) this->erase(clumpBody->id,false);	//Clump has no members any more. Remove it
	}

	if ((b) and (b->isClump())){
		//erase all members if eraseClumpMembers is true:
		const shared_ptr<Clump>& clump=SUDODEM_PTR_CAST<Clump>(b->shape);
		std::map<Body::id_t,Se3r>& members = clump->members;
		FOREACH(MemberMap::value_type& mm, members){
			const Body::id_t& memberId=mm.first;
			if (eraseClumpMembers) {
				this->erase(memberId,false);	// erase members
			} else {
				//when the last members is erased, the clump will be erased automatically, see above
				Body::byId(memberId)->clumpId=Body::ID_NONE; // make members standalones
			}
		}
		body[id].reset();
		return true;
	}
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {  //Iterate over all body's interactions
		scene->interactions->requestErase((*it).second);
	}
	body[id].reset();
	return true;
}
