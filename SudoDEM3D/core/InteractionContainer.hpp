// 2004 © Olivier Galizzi <olivier.galizzi@imag.fr>
// 2004 © Janek Kozicki <cosurgi@berlios.de>
// 2010 © Václav Šmilauer <eudoxos@arcig.cz>

#pragma once

#include<sudodem/lib/serialization/Serializable.hpp>
#include<boost/thread/mutex.hpp>

#ifdef SUDODEM_OPENMP
	#include<omp.h>
#endif

#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/BodyContainer.hpp>

/* This InteractionContainer implementation has reference to the body container and
stores interactions in 2 places:

* Internally in a std::vector; that allows for const-time linear traversal.
  Each interaction internally holds back-reference to the position in this container in Interaction::linIx.
* Inside Body::intrs (in the body with min(id1,id2)).

Both must be kep in sync, which is handled by insert & erase methods.

It was originally written by 2008 © Sergei Dorofeenko <sega@users.berlios.de>,
later devirtualized and put here.

Alternative implementations of InteractionContainer should implement the same API. Due to performance
reasons, no base class with virtual methods defining such API programatically is defined (it could
be possible to create class template for this, though).

Future (?):

* the linear vector might be removed; in favor of linear traversal of bodies by their subdomains,
  then traversing the map in each body. If the previous point would come to realization, half of the
  interactions would have to be skipped explicitly in such a case.

*/
class InteractionContainer: public Serializable{
	private:
		typedef vector<shared_ptr<Interaction> > ContainerT;
		// linear array of container interactions
		ContainerT linIntrs;
		// pointer to body container, since each body holds (some) interactions
		// this must always point to scene->bodies->body
		const BodyContainer::ContainerT* bodies;
		// always in sync with intrs.size(), to avoid that function call
		size_t currSize;
		shared_ptr<Interaction> empty;
		// used only during serialization/deserialization
		vector<shared_ptr<Interaction> > interaction;
	public:
		// flag for notifying the collider that persistent data should be invalidated
		bool dirty;
		// required by the class factory... :-|
		InteractionContainer(): currSize(0),dirty(false),serializeSorted(false),iterColliderLastRun(-1){
			bodies=NULL;
// 			#ifdef SUDODEM_OPENMP
// 				threadsPendingErase.resize(omp_get_max_threads());
// 			#endif
		}
		void clear();
		// iterators
		typedef ContainerT::iterator iterator;
		typedef ContainerT::const_iterator const_iterator;
		iterator begin(){return linIntrs.begin();}
		iterator end()  {return linIntrs.end();}
		const_iterator begin() const {return linIntrs.begin();}
		const_iterator end()   const {return linIntrs.end();}
		// insertion/deletion
		bool insert(Body::id_t id1,Body::id_t id2);
		bool insert(const shared_ptr<Interaction>& i);
		//3rd parameter is used to remove I from linIntrs (in conditionalyEraseNonReal()) when body b1 has been removed
		bool erase(Body::id_t id1,Body::id_t id2,int linPos=-1);
		const shared_ptr<Interaction>& find(Body::id_t id1,Body::id_t id2);
// 		bool found(Body::id_t id1,Body::id_t id2);
		inline bool found(const Body::id_t& id1,const Body::id_t& id2){
			assert(bodies); return (id1>id2)?(*bodies)[id2]->intrs.count(id1):(*bodies)[id1]->intrs.count(id2);}
		// index access
		shared_ptr<Interaction>& operator[](size_t id){return linIntrs[id];}
		const shared_ptr<Interaction>& operator[](size_t id) const { return linIntrs[id];}
		size_t size(){ return currSize; }
		// simulation API

		//! Erase all non-real (in term of Interaction::isReal()) interactions
		void eraseNonReal();

		// mutual exclusion to avoid crashes in the rendering loop
		boost::mutex drawloopmutex;
		// sort interactions before serializations; useful if comparing XML files from different runs (false by default)
		bool serializeSorted;
		// iteration number when the collider was last run; set by the collider, if it wants interactions that were not encoutered in that step to be deleted by InteractionLoop (such as SpatialQuickSortCollider). Other colliders (such as InsertionSortCollider) set it it -1, which is the default
		long iterColliderLastRun;
		//! Ask for erasing the interaction given (from the constitutive law); this resets the interaction (to the initial=potential state) and collider should traverse potential interactions to decide whether to delete them completely or keep them potential
		void requestErase(Body::id_t id1, Body::id_t id2);
		void requestErase(const shared_ptr<Interaction>& I);
		void requestErase(Interaction* I);

		/*! Traverse all interactions and erase them if they are not real and the (T*)->shouldBeErased(id1,id2) return true, or if body(id1) has been deleted
			Class using this interface (which is presumably a collider) must define the
				bool shouldBeErased(Body::id_t, Body::id_t) const
		*/
		template<class T> size_t conditionalyEraseNonReal(const T& t, Scene* rb){
			// beware iterators here, since erase is invalidating them. We need to iterate carefully, and keep in mind that erasing one interaction is moving the last one to the current position.
			size_t initSize=currSize;
		 	for (size_t linPos=0; linPos<currSize;){
				const shared_ptr<Interaction>& i=linIntrs[linPos];
				if(!i->isReal() && t.shouldBeErased(i->getId1(),i->getId2(),rb)) erase(i->getId1(),i->getId2(),linPos);
				else linPos++;}
			return initSize-currSize;
		}

	// we must call Scene's ctor (and from Scene::postLoad), since we depend on the existing BodyContainer at that point.
	void postLoad__calledFromScene(const shared_ptr<BodyContainer>&);
	void preLoad(InteractionContainer&);
	void preSave(InteractionContainer&);
	void postSave(InteractionContainer&);


	REGISTER_ATTRIBUTES(Serializable,(interaction)(serializeSorted)(dirty));
	REGISTER_CLASS_AND_BASE(InteractionContainer,Serializable);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(InteractionContainer);
