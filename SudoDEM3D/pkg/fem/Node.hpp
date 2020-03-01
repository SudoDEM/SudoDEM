/*************************************************************************
*  Copyright (C) 20017 by Sway Zhao                                      *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include<sudodem/lib/base/Math.hpp>//Sway, fix _POXIC_C_SOURCE warning
#include<sudodem/core/Shape.hpp>
#include<sudodem/core/Bound.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/core/Material.hpp>


#include<sudodem/lib/serialization/Serializable.hpp>
#include<sudodem/lib/multimethods/Indexable.hpp>
#include<boost/tuple/tuple.hpp>



class Scene;
class Interaction;

class Node: public Serializable{
	public:
		// numerical types for storing ids
		typedef int id_t;
		// internal structure to hold some interaction of a Node; used by InteractionContainer;
		typedef std::map<Node::id_t, shared_ptr<Interaction> > MapId2IntrT;
		// groupMask type
		//! mutex for updating the parameters from within the interaction loop (only used rarely)
		boost::mutex updateMutex;
		// bits for Node::flags
		enum { FLAG_BOUNDED=1, FLAG_ASPHERICAL=2 }; /* add powers of 2 as needed */
		//! symbolic constant for Node that doesn't exist.
		static const Node::id_t ID_NONE;
		//! get Node pointer given its id.
		static const shared_ptr<Node>& byId(Node::id_t _id,Scene* rb=NULL);
		static const shared_ptr<Node>& byId(Node::id_t _id,shared_ptr<Scene> rb);
        // bits for blockedDOFs
		enum {DOF_NONE=0,DOF_X=1,DOF_Y=2,DOF_Z=4,DOF_RX=8,DOF_RY=16,DOF_RZ=32};
		//! shorthand for all DOFs blocked
		static const unsigned DOF_ALL=DOF_X|DOF_Y|DOF_Z|DOF_RX|DOF_RY|DOF_RZ;
		//! shorthand for all displacements blocked
		static const unsigned DOF_XYZ=DOF_X|DOF_Y|DOF_Z;
		//! shorthand for all rotations blocked
		static const unsigned DOF_RXRYRZ=DOF_RX|DOF_RY|DOF_RZ;

		//! Return DOF_* constant for given axis∈{0,1,2} and rotationalDOF∈{false(default),true}; e.g. axisDOF(0,true)==DOF_RX
		static unsigned axisDOF(int axis, bool rotationalDOF=false){return 1<<(axis+(rotationalDOF?3:0));}
		//! set DOFs according to two Vector3r arguments (blocked is when disp[i]==1.0 or rot[i]==1.0)
		void setDOFfromVector3r(Vector3r disp,Vector3r rot=Vector3r::Zero());
		//! Getter of blockedDOFs for list of strings (e.g. DOF_X | DOR_RX | DOF_RZ → 'xXZ')
		std::string blockedDOFs_vec_get() const;
		//! Setter of blockedDOFs from string ('xXZ' → DOF_X | DOR_RX | DOF_RZ)
		#ifdef SUDODEM_DEPREC_DOF_LIST
			void blockedDOFs_vec_set(const python::object&);
		#else
			void blockedDOFs_vec_set(const std::string& dofs);
		#endif

		//! Return displacement (current-reference position)
		Vector3r displ() const {return pos-refPos;}

		// python access functions: pos and ori are references to inside Se3r and cannot be pointed to directly
		Vector3r pos_get() const {return pos;}
		void pos_set(const Vector3r p) {pos=p;}

		/*! Hook for clump to update position of members when user-forced reposition and redraw (through GUI) occurs.
		 * This is useful only in cases when engines that do that in every iteration are not active - i.e. when the simulation is paused.
		 * (otherwise, GLViewer would depend on Clump and therefore Clump would have to go to core...) */
		virtual void userForcedDisplacementRedrawHook(){return;}

		boost::python::list py_intrs();

		Node::id_t getId() const {return id;};
		unsigned int coordNumber();  // Number of neighboring particles

		mask_t getGroupMask() const {return groupMask; };
		bool maskOk(int mask) const;
		bool maskCompatible(int mask) const;
#ifdef SUDODEM_MASK_ARBITRARY
		bool maskOk(const mask_t& mask) const;
		bool maskCompatible(const mask_t& mask) const;
#endif

		// only NodeContainer can set the id of a Node
		friend class NodeContainer;

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Node,Serializable,"A node for fem elements; no interaction with other Nodes.",
		((Node::id_t,id,Node::ID_NONE,Attr::readonly,"Unique id of this Node."))
		((mask_t,groupMask,1,,"Bitmask for determining interactions."))
		((int,flags,FLAG_BOUNDED,Attr::readonly,"Bits of various Node-related flags. *Do not access directly*. In c++, use isDynamic/setDynamic, isBounded/setBounded, isAspherical/setAspherical. In python, use :yref:`Node.dynamic`, :yref:`Node.bounded`, :yref:`Node.aspherical`."))
        ((Vector3r,pos,Vector3r::Zero(),,"Position."))
		((Vector3r,vel,Vector3r::Zero(),,"Current linear velocity."))
        ((Vector3r,angVel,Vector3r::Zero(),,"Current angular velocity."))
        ((Vector3r,inertia,Vector3r::Zero(),,"Inertia of associated body, in local coordinate system."))
		((Real,mass,0,,"Mass of this Node"))
		((Vector3r,refPos,Vector3r::Zero(),,"Reference position"))
		((Quaternionr,ori,Quaternionr::Identity(),,"Orientation :math:`q` of this node.")) 
		((unsigned,blockedDOFs,,,"[Will be overridden]"))
		((bool,isDamped,true,,"Damping in :yref:`Newtonintegrator` can be deactivated for individual particles by setting this variable to FALSE. E.g. damping is inappropriate for particles in free flight under gravity but it might still be applicable to other particles in the same simulation."))        

		//((MapId2IntrT,intrs,,Attr::hidden,"Map from otherId to Interaction with otherId, managed by InteractionContainer. NOTE: (currently) does not contain all interactions with this Node (only those where otherId>id), since performance issues with such data duplication have not yet been investigated."))
		((long,iterBorn,-1,,"Step number at which the Node was added to simulation."))
		((Real,timeBorn,-1,,"Time at which the Node was added to simulation."))
		,
		/* ctor */,
		/* py */
		//
		.add_property("mask",boost::python::make_getter(&Node::groupMask,boost::python::return_value_policy<boost::python::return_by_value>()),boost::python::make_setter(&Node::groupMask,boost::python::return_value_policy<boost::python::return_by_value>()),"Shorthand for :yref:`Node::groupMask`")
		.add_property("iterBorn",&Node::iterBorn,"Returns step number at which the Node was added to simulation.")
		.add_property("timeBorn",&Node::timeBorn,"Returns time at which the Node was added to simulation.")
		//.def("intrs",&Node::py_intrs,"Return all interactions in which this Node participates.")
        .add_property("blockedDOFs",&Node::blockedDOFs_vec_get,&Node::blockedDOFs_vec_set,"Degress of freedom where linear/angular velocity will be always constant (equal to zero, or to an user-defined value), regardless of applied force/torque. String that may contain 'xyzXYZ' (translations and rotations).")
		// references must be set using wrapper funcs
		.add_property("pos",&Node::pos_get,&Node::pos_set,"Current position.")
		.def("displ",&Node::displ,"Displacement from :yref:`reference position<Node.refPos>` (:yref:`pos<Node.pos>` - :yref:`refPos<State.refPos>`)")
	);
};

REGISTER_SERIALIZABLE(Node);

#if SUDODEM_OPENMP
	#define SUDODEM_PARALLEL_FOREACH_NODE_BEGIN(b_,nodes) const Node::id_t _sz(nodes->size()); _Pragma("omp parallel for") for(Node::id_t _id=0; _id<_sz; _id++){ if(!(*nodes)[_id])  continue; b_((*nodes)[_id]);
	#define SUDODEM_PARALLEL_FOREACH_NODE_END() }
#else
	#define SUDODEM_PARALLEL_FOREACH_NODE_BEGIN(b,nodes) FOREACH(b,*(nodes)){
	#define SUDODEM_PARALLEL_FOREACH_NODE_END() }
#endif

/*
Container of Nodes implemented as flat std::vector. It handles Node removal and
intelligently reallocates free ids for newly added ones.
The nested iterators and the specialized FOREACH_Node macros above will silently skip null Node pointers which may exist after removal. The null pointers can still be accessed via the [] operator.

Any alternative implementation should use the same API.
*/
class NodeContainer: public Serializable{
	private:
		typedef std::vector<shared_ptr<Node> > ContainerT;
		typedef std::map<Node::id_t,Se3r> MemberMap;
		ContainerT node;
	public:
		//friend class InteractionContainer;  // accesses the Node vector directly

		//An iterator that will automatically jump slots with null Nodes
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

		NodeContainer() {};
		virtual ~NodeContainer() {};
		Node::id_t insert(shared_ptr<Node>&);
		void clear();
		iterator begin() {
			iterator temp(node.begin()); temp.end=node.end();
			return (node.begin()==node.end() || *temp)?temp:++temp;}
		iterator end() { iterator temp(node.end()); temp.end=node.end(); return temp;}
		const_iterator begin() const { return begin();}
		const_iterator end() const { return end();}

		size_t size() const { return node.size(); }
		shared_ptr<Node>& operator[](unsigned int id){ return node[id];}
		const shared_ptr<Node>& operator[](unsigned int id) const { return node[id]; }

		bool exists(Node::id_t id) const { return (id>=0) && ((size_t)id<node.size()) && ((bool)node[id]); }
		bool erase(Node::id_t id);

		REGISTER_CLASS_AND_BASE(NodeContainer,Serializable);
		REGISTER_ATTRIBUTES(Serializable,(node));
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(NodeContainer);
