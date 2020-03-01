// (c) 2007-2010 Vaclav Smilauer <eudoxos@arcig.cz>

#include"Clump.hpp"
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/BodyContainer.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/Sphere.hpp>

SUDODEM_PLUGIN((Clump));
CREATE_LOGGER(Clump);

boost::python::dict Clump::members_get(){
  boost::python::dict ret;
	FOREACH(MemberMap::value_type& b, members){
		ret[b.first]=boost::python::make_tuple(b.second.position,b.second.orientation);
	}
	return ret;
}

void Clump::add(const shared_ptr<Body>& clumpBody, const shared_ptr<Body>& subBody){
	Body::id_t subId=subBody->getId();
	const shared_ptr<Clump> clump=SUDODEM_PTR_CAST<Clump>(clumpBody->shape);
	if(clump->members.count(subId)!=0) throw std::invalid_argument(("Body #"+boost::lexical_cast<string>(subId)+" is already part of this clump #"+boost::lexical_cast<string>(clumpBody->id)).c_str());
	if(subBody->isClumpMember()) throw std::invalid_argument(("Body #"+boost::lexical_cast<string>(subId)+" is already a clump member of #"+boost::lexical_cast<string>(subBody->clumpId)).c_str());
	else if(subBody->isClump()){
		const shared_ptr<Clump> subClump=SUDODEM_PTR_CAST<Clump>(subBody->shape);
		FOREACH(const MemberMap::value_type& mm, subClump->members){
			const Body::id_t& memberId=mm.first;
			Scene* scene(Omega::instance().getScene().get());	// get scene
			const shared_ptr<Body>& member=Body::byId(memberId,scene);
			assert(member->isClumpMember());
			member->clumpId=clumpBody->id;
			clump->members[memberId]=Se3r();// meaningful values will be put in by Clump::updateProperties
			//LOG_DEBUG("Added body #"<<memberId->id<<" to clump #"<<clumpBody->id);
		}
		//LOG_DEBUG("Clump #"<<subClump->id<<" will be erased.");// see addToClump() in sudodemWrapper.cpp
	}
	else{	// subBody must be a standalone!
		clump->members[subId]=Se3r();// meaningful values will be put in by Clump::updateProperties
		subBody->clumpId=clumpBody->id;
	}
	clumpBody->clumpId=clumpBody->id; // just to make sure
	clumpBody->setBounded(false); // disallow collisions with the clump itself
	if(subBody->isStandalone()){LOG_DEBUG("Added body #"<<subBody->id<<" to clump #"<<clumpBody->id);}
}

void Clump::del(const shared_ptr<Body>& clumpBody, const shared_ptr<Body>& subBody){
	// erase the subBody; removing body that is not part of the clump throws
	const shared_ptr<Clump> clump=SUDODEM_PTR_CAST<Clump>(clumpBody->shape);
	if(clump->members.erase(subBody->id)!=1) throw std::invalid_argument(("Body #"+boost::lexical_cast<string>(subBody->id)+" not part of clump #"+boost::lexical_cast<string>(clumpBody->id)+"; not removing.").c_str());
	subBody->clumpId=Body::ID_NONE;
	LOG_DEBUG("Removed body #"<<subBody->id<<" from clump #"<<clumpBody->id);
}

void Clump::addForceTorqueFromMembers(const State* clumpState, Scene* scene, Vector3r& F, Vector3r& T){
	FOREACH(const MemberMap::value_type& mm, members){
		const Body::id_t& memberId=mm.first;
		const shared_ptr<Body>& member=Body::byId(memberId,scene);
		assert(member->isClumpMember());
		State* memberState=member->state.get();
		const Vector3r& f=scene->forces.getForce(memberId);
		const Vector3r& t=scene->forces.getTorque(memberId);
		F+=f;
		T+=t+(memberState->pos-clumpState->pos).cross(f);
	}
}

/*! Clump's se3 will be updated (origin at centroid and axes coincident with principal inertia axes) and subSe3 modified in such a way that members positions in world coordinates will not change.

	Note: velocities and angularVelocities of constituents are zeroed.

	OLD DOCS (will be cleaned up):

	-# Clump::members values and Clump::physicalParameters::se3 are invalid from this point
	-# M=0; S=vector3r(0,0,0); I=zero tensor; (ALL calculations are in world coordinates!)
	-# loop over Clump::members (position x_i, mass m_i, inertia at subBody's centroid I_i) [this loop will be replaced by numerical integration (rasterization) for the intersecting case; the rest will be the same]
		- M+=m_i
		- S+=m_i*x_i (local static moments are zero (centroid)
		- get inertia tensor of subBody in world coordinates, by rotating the principal (local) tensor against subBody->se3->orientation; then translate it to world origin (parallel axes theorem), then I+=I_i_world
	-# clumpPos=S/M
	-# translate aggregate's inertia tensor; parallel axes on I (R=clumpPos): I^c_jk=I'_jk-M*(delta_jk R.R - R_j*R_k) [http://en.wikipedia.org/wiki/Moments_of_inertia#Parallel_axes_theorem]
	-# eigen decomposition of I, get principal inertia and rotation matrix of the clump
	-# se3->orientation=quaternion(rotation_matrix); se3->position=clumpPos
	-#	update subSe3s

*/

void Clump::updateProperties(const shared_ptr<Body>& clumpBody, unsigned int discretization){
	LOG_DEBUG("Updating clump #"<<clumpBody->id<<" parameters");
	const shared_ptr<State> state(clumpBody->state);
	const shared_ptr<Clump> clump(SUDODEM_PTR_CAST<Clump>(clumpBody->shape));

	if(clump->members.empty()){ throw std::runtime_error("Clump::updateProperties: clump has zero members."); }
	// trivial case
	if(clump->members.size()==1){
		LOG_DEBUG("Clump of size one will be treated specially.")
		MemberMap::iterator I=clump->members.begin();
		shared_ptr<Body> subBody=Body::byId(I->first);
		//const shared_ptr<RigidBodyParameters>& subRBP(SUDODEM_PTR_CAST<RigidBodyParameters>(subBody->physicalParameters));
		State* subState=subBody->state.get();
		// se3 of the clump as whole is the same as the member's se3
		state->pos=subState->pos;
		state->ori=subState->ori;
		// relative member's se3 is identity
		I->second.position=Vector3r::Zero(); I->second.orientation=Quaternionr::Identity();
		state->inertia=subState->inertia;
		state->mass=subState->mass;
		state->vel=Vector3r::Zero();
		state->angVel=Vector3r::Zero();
		return;
	}
	//check for intersections:
	bool intersecting = false;
	shared_ptr<Sphere> sph (new Sphere);
	int Sph_Index = sph->getClassIndexStatic();		// get sphere index for checking if bodies are spheres
	if (discretization>0){
		FOREACH(MemberMap::value_type& mm, clump->members){
			const shared_ptr<Body> subBody1=Body::byId(mm.first);
			FOREACH(MemberMap::value_type& mm, clump->members){
				const shared_ptr<Body> subBody2=Body::byId(mm.first);
				if ((subBody1->shape->getClassIndex() ==  Sph_Index) && (subBody2->shape->getClassIndex() ==  Sph_Index) && (subBody1!=subBody2)){//clump members should be spheres
					Vector3r dist = subBody1->state->pos - subBody2->state->pos;
					const Sphere* sphere1 = SUDODEM_CAST<Sphere*> (subBody1->shape.get());
					const Sphere* sphere2 = SUDODEM_CAST<Sphere*> (subBody2->shape.get());
					Real un = (sphere1->radius+sphere2->radius) - dist.norm();
					if (un > 0.001*min(sphere1->radius,sphere2->radius)) {intersecting = true; break;}
				}
			}
			if (intersecting) break;
		}
	}
	/* quantities suffixed by
		g: global (world) coordinates
		s: local subBody's coordinates
		c: local clump coordinates
	*/
	Real M=0; // mass
	Real dens=0;//density
	Vector3r Sg(0,0,0); // static moment, for getting clump's centroid
	Matrix3r Ig(Matrix3r::Zero()), Ic(Matrix3r::Zero()); // tensors of inertia; is upper triangular, zeros instead of symmetric elements

	/**
	algorithm for estimation of volumes and inertia tensor from clumps using summation/integration scheme with regular grid spacing
	(some parts copied from woo: http://bazaar.launchpad.net/~eudoxos/woo/trunk/view/head:/pkg/dem/Clump.cpp)
	*/
	if(intersecting){
		//get boundaries and minimum radius of clump:
		Real rMin=1./0.; AlignedBox3r aabb;
		FOREACH(MemberMap::value_type& mm, clump->members){
			const shared_ptr<Body> subBody = Body::byId(mm.first);
			if (subBody->shape->getClassIndex() == Sph_Index){//clump member should be a sphere
				const Sphere* sphere = SUDODEM_CAST<Sphere*> (subBody->shape.get());
				aabb.extend(subBody->state->pos + Vector3r::Constant(sphere->radius));
				aabb.extend(subBody->state->pos - Vector3r::Constant(sphere->radius));
				rMin=min(rMin,sphere->radius);
			}
		}
		//get volume and inertia tensor using regular cubic cell array inside bounding box of the clump:
		Real dx = rMin/discretization; 	//edge length of cell
		//Real aabbMax = max(max(aabb.max().x()-aabb.min().x(),aabb.max().y()-aabb.min().y()),aabb.max().z()-aabb.min().z());
		//if (aabbMax/dx > 150) dx = aabbMax/150;//limit dx: leads to bug, when long chain of clump members is created (https://bugs.launchpad.net/sudodem/+bug/1273172)
		Real dv = pow(dx,3);		//volume of cell
		long nCells=(aabb.sizes()/dx).prod();
		if(nCells>1e7) LOG_WARN("Clump::updateProperties: Cell array has "<<nCells<<" cells. Integrate inertia may take a while ...");
		Vector3r x;			//position vector (center) of cell
		for(x.x()=aabb.min().x()+dx/2.; x.x()<aabb.max().x(); x.x()+=dx){
			for(x.y()=aabb.min().y()+dx/2.; x.y()<aabb.max().y(); x.y()+=dx){
				for(x.z()=aabb.min().z()+dx/2.; x.z()<aabb.max().z(); x.z()+=dx){
					FOREACH(MemberMap::value_type& mm, clump->members){
						const shared_ptr<Body> subBody = Body::byId(mm.first);
						if (subBody->shape->getClassIndex() == Sph_Index){//clump member should be a sphere
							dens = subBody->material->density;
							const Sphere* sphere = SUDODEM_CAST<Sphere*> (subBody->shape.get());
							if((x-subBody->state->pos).squaredNorm() < pow(sphere->radius,2)){
								M += dv;
								Sg += dv*x;
								//inertia I = sum_i( mass_i*dist^2 + I_s) )	//steiners theorem
								Ig += dv*( x.dot(x)*Matrix3r::Identity()-x*x.transpose())/*dist^2*/+Matrix3r(Vector3r::Constant(dv*pow(dx,2)/6.).asDiagonal())/*I_s/m = d^2: along princial axes of dv; perhaps negligible?*/;
								break;
							}
						}
					}
				}
			}
		}
	}
	if(!intersecting){
		FOREACH(MemberMap::value_type& mm, clump->members){
			// mm.first is Body::id_t, mm.second is Se3r of that body
			const shared_ptr<Body> subBody=Body::byId(mm.first);
			dens = subBody->material->density;
			if (subBody->shape->getClassIndex() ==  Sph_Index){//clump member should be a sphere
				State* subState=subBody->state.get();
				const Sphere* sphere = SUDODEM_CAST<Sphere*> (subBody->shape.get());
				Real vol = (4./3.)*Mathr::PI*pow(sphere->radius,3.);
				M+=vol;
				Sg+=vol*subState->pos;
				Ig+=Clump::inertiaTensorTranslate(Vector3r::Constant((2/5.)*vol*pow(sphere->radius,2)).asDiagonal(),vol,-1.*subState->pos);
			}
		}
	}
	assert(M>0); LOG_TRACE("M=\n"<<M<<"\nIg=\n"<<Ig<<"\nSg=\n"<<Sg);
	// clump's centroid
	state->pos=Sg/M;
	// this will calculate translation only, since rotation is zero
	Matrix3r Ic_orientG=inertiaTensorTranslate(Ig, -M /* negative mass means towards centroid */, state->pos); // inertia at clump's centroid but with world orientation
	LOG_TRACE("Ic_orientG=\n"<<Ic_orientG);
	Ic_orientG(1,0)=Ic_orientG(0,1); Ic_orientG(2,0)=Ic_orientG(0,2); Ic_orientG(2,1)=Ic_orientG(1,2); // symmetrize
	Eigen::SelfAdjointEigenSolver<Matrix3r> decomposed(Ic_orientG);
	const Matrix3r& R_g2c(decomposed.eigenvectors());
	// has NaNs for identity matrix??
	LOG_TRACE("R_g2c=\n"<<R_g2c);
	// set quaternion from rotation matrix
	state->ori=Quaternionr(R_g2c); state->ori.normalize();
	state->inertia=dens*decomposed.eigenvalues();
	state->mass=dens*M;

	// TODO: these might be calculated from members... but complicated... - someone needs that?!
	state->vel=state->angVel=Vector3r::Zero();
	clumpBody->setAspherical(state->inertia[0]!=state->inertia[1] || state->inertia[0]!=state->inertia[2]);

	// update subBodySe3s; subtract clump orientation (=apply its inverse first) to subBody's orientation
	FOREACH(MemberMap::value_type& I, clump->members){
		shared_ptr<Body> subBody=Body::byId(I.first);
		State* subState=subBody->state.get();
		I.second.orientation=state->ori.conjugate()*subState->ori;
		I.second.position=state->ori.conjugate()*(subState->pos-state->pos);
	}
}

/*! @brief Recalculates inertia tensor of a body after translation away from (default) or towards its centroid.
 *
 * @param I inertia tensor in the original coordinates; it is assumed to be upper-triangular (elements below the diagonal are ignored).
 * @param m mass of the body; if positive, translation is away from the centroid; if negative, towards centroid.
 * @param off offset of the new origin from the original origin
 * @return inertia tensor in the new coordinate system; the matrix is symmetric.
 */
Matrix3r Clump::inertiaTensorTranslate(const Matrix3r& I,const Real m, const Vector3r& off){
	return I+m*(off.dot(off)*Matrix3r::Identity()-off*off.transpose());
}

/*! @brief Recalculate body's inertia tensor in rotated coordinates.
 *
 * @param I inertia tensor in old coordinates
 * @param T rotation matrix from old to new coordinates
 * @return inertia tensor in new coordinates
 */
Matrix3r Clump::inertiaTensorRotate(const Matrix3r& I,const Matrix3r& T){
	/* [http://www.kwon3d.com/theory/moi/triten.html] */
	return T.transpose()*I*T;
}

/*! @brief Recalculate body's inertia tensor in rotated coordinates.
 *
 * @param I inertia tensor in old coordinates
 * @param rot quaternion that describes rotation from old to new coordinates
 * @return inertia tensor in new coordinates
 */
Matrix3r Clump::inertiaTensorRotate(const Matrix3r& I, const Quaternionr& rot){
	Matrix3r T=rot.toRotationMatrix();
	return inertiaTensorRotate(I,T);
}
