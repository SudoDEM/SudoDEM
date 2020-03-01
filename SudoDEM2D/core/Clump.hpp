// (c) 2007 Vaclav Smilauer <eudoxos@arcig.cz>

#pragma once

#include<sudodem/core/Body.hpp>
#include<sudodem/lib/base/Logging.hpp>
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/core/PartialEngine.hpp>


/*! Body representing clump (rigid aggregate) composed by other existing bodies.

	Clump is one of bodies that reside in scene->bodies.
	When an existing body is added to ::Clump, it's ::Body::dynamic flag is set to false
	(it is still subscribed to all its engines, to make it possible to remove it from the clump again).
	All forces acting on Clump::members are made to act on the clump itself, which will ensure that they
	influence all Clump::members as if the clump were a rigid particle.

	What are clump requirements so that they function?
	-# Given any body, tell
		- if it is a clump member: Body::isClumpMember()
	 	- if it is a clump: Body:: isClump(). (Correct result is assured at each call to Clump::add).
		 (we could use RTTI instead? Would that be more reliable?)
		- if it is a standalone Body: Body::isStandalone()
		- what is it's clump id (Body::clumpId)
	-# given the root body, tell
		- what clumps it contains (enumerate all bodies and filter clumps, see above)
	-#	given a clump, tell
		- what bodies it contains (keys of ::Clump::members)
		- what are se3 of these bodies (values of ::Clump::members)
	-# add/delete bodies from/to clump (::Clump::add, ::Clump::del)
		- This includes saving se3 of the subBody: it \em must be in clump's local coordinates so that it is constant. The transformation from global to local is given by clump's se3 at the moment of addition. Clump's se3 is initially (origin,identity)
	-# Update clump's physical properties (Clump::updateProperties)
		- This \em must reposition members so that they have the same se3 globally
	-# Apply forces acting on members to the clump instead (done in NewtonsForceLaw, NewtonsMomentumLaw) - uses world coordinates to calculate effect on the clump's centroid
	-# Integrate position and orientation of the clump
		- LeapFrogPositionIntegrator and LeapFrogOrientationIntegrator move clump as whole
			- clump members are skipped, since they have Body::dynamic==false.
		- ClumpMemberMover is an engine that updates positions of the clump memebers in each timestep (calls Clump::moveSubBodies internally)

	Some more information can be found http://beta.arcig.cz/~eudoxos/phd/index.cgi/YaDe/HighLevelClumps

	For an example how to generate a clump, see ClumpTestGen::createOneClump.

	@todo GravityEngine should be applied to members, not to clump as such?! Still not sure. Perhaps Clumps should have mass and inertia set to zeros so that engines unaware of clumps do not act on it. It would have some private mass and insertia that would be used in NewtonsForceLaw etc for clumps specially...

	@note Collider::mayCollide (should be used by all colliders) bypasses Clumps explicitly. This no longer depends on the absence of bound.
	@note Clump relies on its id being assigned (as well as id of its components); therefore, only bodies that have already been inserted to the container may be added to Clump which has been itself already added to the container. We further requier that clump id is greater than ids of clumped bodies

 */

class NewtonIntegrator;

class Clump: public Shape {
	public:
		typedef std::map<Body::id_t,Se2r> MemberMap;

		static void add(const shared_ptr<Body>& clump, const shared_ptr<Body>& subBody);
		static void del(const shared_ptr<Body>& clump, const shared_ptr<Body>& subBody);
		//! Recalculate physical properties of Clump.
		static void updateProperties(const shared_ptr<Body>& clump, unsigned int discretization);
		//! Calculate positions and orientations of members based on relative Se3; newton pointer (if non-NULL) calls NewtonIntegrator::saveMaximaVelocity
		// done as template to avoid cross-dependency between clump and newton (not necessary if all plugins are linked together)
		template<class IntegratorT>
		static void moveMembers(const shared_ptr<Body>& clumpBody, Scene* scene, IntegratorT* integrator=NULL){
			const shared_ptr<Clump>& clump=SUDODEM_PTR_CAST<Clump>(clumpBody->shape);
			const shared_ptr<State>& clumpState=clumpBody->state;
			FOREACH(MemberMap::value_type& B, clump->members){
				// B.first is Body::id_t, B.second is local Se3r of that body in the clump
				const shared_ptr<Body>& b = Body::byId(B.first,scene);
				const shared_ptr<State>& subState=b->state; const Vector2r& subPos(B.second.position); const Rotationr& subOri(B.second.rotation);
				// position update
				subState->pos=clumpState->pos+clumpState->ori*subPos;
				subState->ori=clumpState->ori*subOri;
				// velocity update
				subState->vel=clumpState->vel+clumpState->angVel*(subState->pos-clumpState->pos);
				subState->angVel=clumpState->angVel;
				if(integrator) integrator->saveMaximaDisplacement(b);
			}
		}


		//! update member positions after clump being moved by mouse (in case simulation is paused and engines will not do that).
		void userForcedDisplacementRedrawHook(){ throw runtime_error("Clump::userForcedDisplacementRedrawHook not yet implemented (with Clump as subclass of Shape).");}

		//! get force and torque on the clump itself, from forces/torques on members; does not include force on clump itself
		void addForceTorqueFromMembers(const State* clumpState, Scene* scene, Vector2r& F, Real& T);


		//! Recalculates inertia tensor of a body after translation away from (default) or towards its centroid.
		//static Real inertiaTensorTranslate(const Matrix3r& I,const Real m, const Vector2r& off);
		//! Recalculate body's inertia tensor in rotated coordinates.
		//static Real inertiaTensorRotate(const Matrix3r& I, const Matrix3r& T);
		//! Recalculate body's inertia tensor in rotated coordinates.
		//static Real inertiaTensorRotate(const Matrix3r& I, const Quaternionr& rot);

    boost::python::dict members_get();

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Clump,Shape,"Rigid aggregate of bodies",
		((MemberMap,members,,Attr::hidden,"Ids and relative positions+orientations of members of the clump (should not be accessed directly)"))
		// ((vector<int>,ids,,Attr::readonly,"Ids of constituent particles (only informative; direct modifications will have no effect)."))
		,/*ctor*/ createIndex();
		,/*py*/ .add_property("members",&Clump::members_get,"Return clump members as {'id1':(relPos,relOri),...}")
	);
	DECLARE_LOGGER;
	REGISTER_CLASS_INDEX(Clump,Shape);
};
REGISTER_SERIALIZABLE(Clump);
