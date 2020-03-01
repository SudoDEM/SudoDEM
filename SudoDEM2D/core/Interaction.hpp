#pragma once
#include<sudodem/lib/serialization/Serializable.hpp>
// keep those two here, template instantiation & boost::python gets broken otherwise, e.g. https://bugs.launchpad.net/bugs/618766
#include<sudodem/core/IGeom.hpp>
#include<sudodem/core/IPhys.hpp>
#include<sudodem/core/Body.hpp>


class IGeomFunctor;
class IPhysFunctor;
class LawFunctor;
class Scene;

class Interaction: public Serializable{
	private:
		friend class IPhysDispatcher;
		friend class InteractionLoop;
	public:
		bool isReal() const {return (bool)geom && (bool)phys;}
		//! If this interaction was just created in this step (for the constitutive law, to know that it is the first time there)
		bool isFresh(Scene* rb);
		bool isActive;

		Interaction(Body::id_t newId1,Body::id_t newId2);

		const Body::id_t& getId1() const {return id1;};
		const Body::id_t& getId2() const {return id2;};

		//! swaps order of bodies within the interaction
		void swapOrder();

		bool operator<(const Interaction& other) const { return getId1()<other.getId1() || (getId1()==other.getId1() && getId2()<other.getId2()); }

		//! cache functors that are called for this interaction. Currently used by InteractionLoop.
		struct {
			// Whether geometry dispatcher exists at all; this is different from !geom, since that can mean we haven't populated the cache yet.
			// Therefore, geomExists must be initialized to true first (done in Interaction::reset() called from ctor).
			bool geomExists;
			// shared_ptr's are initialized to NULLs automagically
			shared_ptr<IGeomFunctor> geom;
			shared_ptr<IPhysFunctor> phys;
			shared_ptr<LawFunctor> constLaw;
		} functorCache;

		//! Reset interaction to the intial state (keep only body ids)
		void reset();
		//! common initialization called from both constructor and reset()
		void init();

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Interaction,Serializable,"Interaction between pair of bodies.",
		((Body::id_t,id1,0,Attr::readonly,":yref:`Id<Body::id>` of the first body in this interaction."))
		((Body::id_t,id2,0,Attr::readonly,":yref:`Id<Body::id>` of the second body in this interaction."))
		((long,iterMadeReal,-1,,"Step number at which the interaction was fully (in the sense of geom and phys) created. (Should be touched only by :yref:`IPhysDispatcher` and :yref:`InteractionLoop`, therefore they are made friends of Interaction"))
		((long,iterLastSeen,-1,(Attr::noSave|Attr::hidden),"At which step this interaction was last detected by the collider. InteractionLoop will remove it if InteractionContainer::iterColliderLastRun==scene->iter, InteractionContainer::iterColliderLastRun is positive (some colliders manage interaction deletion themselves, such as :yref:`InsertionSortCollider`) and iterLastSeen<scene->iter."))
		((shared_ptr<IGeom>,geom,,,"Geometry part of the interaction."))
		((shared_ptr<IPhys>,phys,,,"Physical (material) part of the interaction."))
		((Vector2i,cellDist,Vector2i(0,0),,"Distance of bodies in cell size units, if using periodic boundary conditions; id2 is shifted by this number of cells from its :yref:`State::pos` coordinates for this interaction to exist. Assigned by the collider.\n\n.. warning::\n\t(internal)  cellDist must survive Interaction::reset(), it is only initialized in ctor. Interaction that was cancelled by the constitutive law, was reset() and became only potential must have thepriod information if the geometric functor again makes it real. Good to know after few days of debugging that :-)"))
		((int,linIx,-1,(Attr::noSave|Attr::hidden),"Index in the linear interaction container. For internal use by InteractionContainer only."))
		((long,iterBorn,-1,,"Step number at which the interaction was added to simulation."))
		,
		/* ctor */ init(),
		/*py*/
		.add_property("isReal",&Interaction::isReal,"True if this interaction has both geom and phys; False otherwise.")
		.def_readwrite("isActive",&Interaction::isActive,"True if this interaction is active. Otherwise the forces from this interaction will not be taken into account. True by default.")
	);
};

REGISTER_SERIALIZABLE(Interaction);
