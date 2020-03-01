// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2004 Janek Kozicki <cosurgi@berlios.de>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
// © 2006 Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>

#include<sudodem/pkg/dem/ScGeom.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((ScGeom));
ScGeom::~ScGeom(){};

Vector2r& ScGeom::rotate(Vector2r& shearForce) const {
	// approximated rotations
	/*double sin_theta = sin(rotAngle);
	double cos_theta = cos(rotAngle);
	Matrix2r rot;
	rot << cos_theta, -sin_theta,
				 sin_theta, cos_theta;
	shearForce = rot*shearForce;
	*/
	//shearForce -= shearForce.cross(orthonormal_axis);
	//shearForce -= shearForce.cross(twist_axis);
	//NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	//shearForce -= normal.dot(shearForce)*normal;
	Vector2r tangent = Vector2r(-normal[1],normal[0]);
	Vector2r preTangent = Vector2r(-preNormal[1],preNormal[0]);//tangent is left-hand side of normal
  shearForce = shearForce.dot(preTangent)*tangent;
	return shearForce;
}

//!Precompute data needed for rotating tangent vectors attached to the interaction
void ScGeom::precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector2r& currentNormal, bool isNew, const Vector2r& shift2, bool avoidGranularRatcheting){
	if(!isNew) {
		preNormal = normal;//rotAngle = scene->dt*0.5*(rbp1.angVel + rbp2.angVel);
	}
	else preNormal = currentNormal;
	//Update contact normal
	normal=currentNormal;
	//Precompute shear increment
	Vector2r relativeVelocity=getIncidentVel(&rbp1,&rbp2,scene->dt,shift2,scene->isPeriodic ? scene->cell->intrShiftVel(c->cellDist) : Vector2r::Zero(),avoidGranularRatcheting);
	//keep the shear part only
	relativeVelocity = relativeVelocity-normal.dot(relativeVelocity)*normal;
	shearInc = relativeVelocity*scene->dt;
}

Vector2r ScGeom::getIncidentVel(const State* rbp1, const State* rbp2, Real dt, const Vector2r& shift2, const Vector2r& shiftVel, bool avoidGranularRatcheting){
	//left-hand vector perpendicular to a given vector (x,y) is (-y,x)
	Vector2r del = Vector2r(-normal[1],normal[0]);
	if(avoidGranularRatcheting){
		/* B.C. Comment :
		Short explanation of what we want to avoid :
		Numerical ratcheting is best understood considering a small elastic cycle at a contact between two grains : assuming b1 is fixed, impose this displacement to b2 :
		1. translation "dx" in the normal direction
		2. rotation "a"
		3. translation "-dx" (back to initial position)
		4. rotation "-a" (back to initial orientation)
		If the branch vector used to define the relative shear in rotation×branch is not constant (typically if it is defined from the vector center→contactPoint), then the shear displacement at the end of this cycle is not zero: rotations *a* and *-a* are multiplied by branches of different lengths.
		It results in a finite contact force at the end of the cycle even though the positions and orientations are unchanged, in total contradiction with the elastic nature of the problem. It could also be seen as an *inconsistent energy creation or loss*. Given that DEM simulations tend to generate oscillations around equilibrium (damped mass-spring), it can have a significant impact on the evolution of the packings, resulting for instance in slow creep in iterations under constant load.
		The solution adopted here to avoid ratcheting is as proposed by McNamara and co-workers.
		They analyzed the ratcheting problem in detail - even though they comment on the basis of a cycle that differs from the one shown above. One will find interesting discussions in e.g. DOI 10.1103/PhysRevE.77.031304, even though solution it suggests is not fully applied here (equations of motion are not incorporating alpha, in contradiction with what is suggested by McNamara et al.).
		 */
		// For disk-facet contact this will give an erroneous value of relative velocity...
		Real alpha = (radius1+radius2)/(radius1+radius2-penetrationDepth);

		Vector2r relativeVelocity = (rbp2->vel-rbp1->vel)*alpha - rbp2->angVel*radius2*del - rbp1->angVel*radius1*del;
		relativeVelocity+=alpha*shiftVel;
		return relativeVelocity;
	} else {
		// It is correct for disk-disk and disk-facet contact
		Vector2r c1x = (contactPoint - rbp1->pos);
		Vector2r c2x = (contactPoint - rbp2->pos - shift2);
		Vector2r relativeVelocity = (rbp2->vel - rbp2->angVel*c2x.norm()*del) - (rbp1->vel+rbp1->angVel*c1x.norm()*del);
		relativeVelocity+=shiftVel;
		return relativeVelocity;}
}

Vector2r ScGeom::getIncidentVel(const State* rbp1, const State* rbp2, Real dt, bool avoidGranularRatcheting){
	//Just pass null shift to the periodic version
	return getIncidentVel(rbp1,rbp2,dt,Vector2r::Zero(),Vector2r::Zero(),avoidGranularRatcheting);
}

Vector2r ScGeom::getIncidentVel_py(shared_ptr<Interaction> i, bool avoidGranularRatcheting){
	if(i->geom.get()!=this) throw invalid_argument("ScGeom object is not the same as Interaction.geom.");
	Scene* scene=Omega::instance().getScene().get();
	return getIncidentVel(Body::byId(i->getId1(),scene)->state.get(),Body::byId(i->getId2(),scene)->state.get(),
		scene->dt,
		scene->isPeriodic ? scene->cell->intrShiftPos(i->cellDist): Vector2r::Zero(), // shift2
		scene->isPeriodic ? scene->cell->intrShiftVel(i->cellDist): Vector2r::Zero(), // shiftVel
		avoidGranularRatcheting);
}

Real ScGeom::getRelAngVel(const State* rbp1, const State* rbp2, Real dt){
  	return (rbp2->angVel-rbp1->angVel);
}

Real ScGeom::getRelAngVel_py(shared_ptr<Interaction> i){
	if(i->geom.get()!=this) throw invalid_argument("ScGeom object is not the same as Interaction.geom.");
	Scene* scene=Omega::instance().getScene().get();
	return getRelAngVel(Body::byId(i->getId1(),scene)->state.get(),Body::byId(i->getId2(),scene)->state.get(),scene->dt);
}
