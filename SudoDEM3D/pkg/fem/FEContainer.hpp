/*************************************************************************
*  Copyright (C) 20017 by Sway Zhao                                      *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include<sudodem/lib/base/Math.hpp>//fix _POXIC_C_SOURCE warning
#include<sudodem/core/Shape.hpp>

//#include<sudodem/core/State.hpp>
//#include<sudodem/core/Material.hpp>


#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/lib/multimethods/Indexable.hpp>
#include<boost/tuple/tuple.hpp>



class Scene;
class Interaction;

class Element: public Serializable{
	public:
		// numerical types for storing ids
		typedef int id_t;
		// internal structure to hold some interaction of a Element; used by InteractionContainer;
		//typedef std::map<Element::id_t, shared_ptr<Interaction> > MapId2IntrT;
		// groupMask type
		//! mutex for updating the parameters from within the interaction loop (only used rarely)
		boost::mutex updateMutex;
		//! symbolic constant for Element that doesn't exist.
		static const Element::id_t ID_NONE;
		//! get Element pointer given its id.
		static const shared_ptr<Element>& byId(Element::id_t _id,Scene* rb=NULL);
		static const shared_ptr<Element>& byId(Element::id_t _id,shared_ptr<Scene> rb);

		boost::python::list py_intrs();

		Element::id_t getId() const {return id;};

		// only ElementContainer can set the id of a Element
		friend class ElementContainer;

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Element,Serializable,"A general FE element.",
		((Element::id_t,id,Element::ID_NONE,Attr::readonly,"Unique id of this Element."))
        ((Vector3r,pos,Vector3r::Zero(),,"Position."))
		((Quaternionr,ori,Quaternionr::Identity(),,"Orientation :math:`q` of this element."))      
		,
		/* ctor */,
		/* py */
		//
		// references must be set using wrapper funcs
		//.add_property("pos",&Element::pos_get,&Element::pos_set,"Current position.")
	);
};

REGISTER_SERIALIZABLE(Element);

#if SUDODEM_OPENMP
	#define SUDODEM_PARALLEL_FOREACH_ELEMENT_BEGIN(b_,elements) const Element::id_t _sz(elements->size()); _Pragma("omp parallel for") for(Element::id_t _id=0; _id<_sz; _id++){ if(!(*elements)[_id])  continue; b_((*elements)[_id]);
	#define SUDODEM_PARALLEL_FOREACH_ELEMENT_END() }
#else
	#define SUDODEM_PARALLEL_FOREACH_ELEMENT_BEGIN(b,elements) FOREACH(b,*(elements)){
	#define SUDODEM_PARALLEL_FOREACH_ELEMENT_END() }
#endif

/*
Container of Elements implemented as flat std::vector. It handles Element removal and
intelligently reallocates free ids for newly added ones.
The nested iterators and the specialized FOREACH_Element macros above will silently skip null Element pointers which may exist after removal. The null pointers can still be accessed via the [] operator.

Any alternative implementation should use the same API.
*/
class FEContainer: public Serializable{
	private:
		typedef std::vector<shared_ptr<Element> > ContainerT;
		typedef std::map<Element::id_t,Se3r> MemberMap;
		ContainerT element;
	public:
		//friend class InteractionContainer;  // accesses the Element vector directly

		//An iterator that will automatically jump slots with null Elements
		class smart_iterator : public ContainerT::iterator {
			public:
			ContainerT::iterator end;
			smart_iterator& operator++() {
				ContainerT::iterator::operator++();
				while (!(this->operator*()) && end!=(*this)) ContainerT::iterator::operator++();
				return *this;}
			smart_iterator operator++(int) {smart_iterator temp(*this); operator++(); return temp;}
			smart_iterator& operator=(const ContainerT::iterator& rhs) {ContainerT::iterator::operator=(rhs); return *this;}
			smart_iterator& operator=(const smart_iterator& rhs) {ContainerT::iterator::operator=(rhs); end=rhs.end; return *this;}
			smart_iterator() {}
			smart_iterator(const ContainerT::iterator& source) {(*this)=source;}
			smart_iterator(const smart_iterator& source) {(*this)=source; end=source.end;}
		};
		typedef smart_iterator iterator;
		typedef const smart_iterator const_iterator;

		FEContainer() {};
		virtual ~FEContainer() {};
		Element::id_t insert(shared_ptr<Element>&);
		void clear();
		iterator begin() {
			iterator temp(element.begin()); temp.end=element.end();
			return (element.begin()==element.end() || *temp)?temp:++temp;}
		iterator end() { iterator temp(element.end()); temp.end=element.end(); return temp;}
		const_iterator begin() const { return begin();}
		const_iterator end() const { return end();}

		size_t size() const { return element.size(); }
		shared_ptr<Element>& operator[](unsigned int id){ return element[id];}
		const shared_ptr<Element>& operator[](unsigned int id) const { return element[id]; }

		bool exists(Element::id_t id) const { return (id>=0) && ((size_t)id<element.size()) && ((bool)element[id]); }
		bool erase(Element::id_t id);

		REGISTER_CLASS_AND_BASE(FEContainer,Serializable);
		REGISTER_ATTRIBUTES(Serializable,(element));
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(FEContainer);
