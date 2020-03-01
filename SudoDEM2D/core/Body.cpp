
#include<sudodem/core/Body.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/InteractionContainer.hpp>

//! This could be -1 if id_t is re-typedef'ed as `int'
const Body::id_t Body::ID_NONE=Body::id_t(-1);

const shared_ptr<Body>& Body::byId(Body::id_t _id, Scene* rb){return (*((rb?rb:Omega::instance().getScene().get())->bodies))[_id];}
const shared_ptr<Body>& Body::byId(Body::id_t _id, shared_ptr<Scene> rb){return (*(rb->bodies))[_id];}

// return list of interactions of this particle
boost::python::list Body::py_intrs(){
  boost::python::list ret;
	for(Body::MapId2IntrT::iterator it=this->intrs.begin(),end=this->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		if(!(*it).second->isReal()) continue;
		ret.append((*it).second);
	}
	return ret;
}

// return list of interactions of this particle
unsigned int Body::coordNumber(){
	unsigned int intrSize = 0;
	for(Body::MapId2IntrT::iterator it=this->intrs.begin(),end=this->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		if(!(*it).second->isReal()) continue;
		intrSize++;
	}
	return intrSize;
}


bool Body::maskOk(int mask) const { return (mask==0 || ((groupMask & mask) != 0)); }
bool Body::maskCompatible(int mask) const { return (groupMask & mask) != 0; }
#ifdef SUDODEM_MASK_ARBITRARY
bool Body::maskOk(const mask_t& mask) const { return (mask==0 || ((groupMask & mask) != 0)); }
bool Body::maskCompatible(const mask_t& mask) const { return (groupMask & mask) != 0; }
#endif



