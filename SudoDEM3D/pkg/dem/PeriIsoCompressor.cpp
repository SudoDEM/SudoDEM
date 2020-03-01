
// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#include<sudodem/pkg/dem/PeriIsoCompressor.hpp>
#include<sudodem/pkg/dem/Shop.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/lib/pyutil/gil.hpp>

SUDODEM_PLUGIN((PeriIsoCompressor)(PeriTriaxController)(Peri3dController))

CREATE_LOGGER(PeriIsoCompressor);
void PeriIsoCompressor::action(){
	if(!scene->isPeriodic){ LOG_FATAL("Being used on non-periodic simulation!"); throw; }
	if(state>=stresses.size()) return;
	// initialize values
	if(charLen<=0){
		Bound* bv=Body::byId(0,scene)->bound.get();
		if(!bv){ LOG_FATAL("No charLen defined and body #0 has no bound"); throw; }
		const Vector3r sz=bv->max-bv->min;
		charLen=(sz[0]+sz[1]+sz[2])/3.;
		LOG_INFO("No charLen defined, taking avg bbox size of body #0 = "<<charLen);
	}
	if(maxSpan<=0){
		FOREACH(const shared_ptr<Body>& b, *scene->bodies){
			if(!b || !b->bound) continue;
			for(int i=0; i<3; i++) maxSpan=max(maxSpan,b->bound->max[i]-b->bound->min[i]);
		}
	}
	if(maxDisplPerStep<0) maxDisplPerStep=1e-2*charLen; // this should be tuned somehow…
	const long& step=scene->iter;
	Vector3r cellSize=scene->cell->getSize(); //unused: Real cellVolume=cellSize[0]*cellSize[1]*cellSize[2];
	Vector3r cellArea=Vector3r(cellSize[1]*cellSize[2],cellSize[0]*cellSize[2],cellSize[0]*cellSize[1]);
	Real minSize=min(cellSize[0],min(cellSize[1],cellSize[2])), maxSize=max(cellSize[0],max(cellSize[1],cellSize[2]));
	if(minSize<2.1*maxSpan){ throw runtime_error("Minimum cell size is smaller than 2.1*span_of_the_biggest_body! (periodic collider requirement)"); }
	if(((step%globalUpdateInt)==0) || avgStiffness<0 || sigma[0]<0 || sigma[1]<0 || sigma[2]<0){
		Vector3r sumForces=Shop::totalForceInVolume(avgStiffness,scene);
		sigma=-Vector3r(sumForces[0]/cellArea[0],sumForces[1]/cellArea[1],sumForces[2]/cellArea[2]);
		LOG_TRACE("Updated sigma="<<sigma<<", avgStiffness="<<avgStiffness);
	}
	Real sigmaGoal=stresses[state]; assert(sigmaGoal<0);
	// expansion of cell in this step (absolute length)
	Vector3r cellGrow(Vector3r::Zero());
	// is the stress condition satisfied in all directions?
	bool allStressesOK=true;
	if(keepProportions){ // the same algo as below, but operating on quantitites averaged over all dimensions
		Real sigAvg=(sigma[0]+sigma[1]+sigma[2])/3., avgArea=(cellArea[0]+cellArea[1]+cellArea[2])/3., avgSize=(cellSize[0]+cellSize[1]+cellSize[2])/3.;
		Real avgGrow=1e-4*(sigmaGoal-sigAvg)*avgArea/(avgStiffness>0?avgStiffness:1);
		Real maxToAvg=maxSize/avgSize;
		if(std::abs(maxToAvg*avgGrow)>maxDisplPerStep) avgGrow=Mathr::Sign(avgGrow)*maxDisplPerStep/maxToAvg;
		Real okGrow=-(minSize-2.1*maxSpan)/maxToAvg;
		if(avgGrow<okGrow) throw runtime_error("Unable to shring cell due to maximum body size (although required by stress condition). Increase particle rigidity, increase total sample dimensions, or decrease goal stress.");
		// avgGrow=max(avgGrow,-(minSize-2.1*maxSpan)/maxToAvg);
		if(avgStiffness>0) { sigma+=(avgGrow*avgStiffness)*Vector3r::Ones(); sigAvg+=avgGrow*avgStiffness; }
		if(std::abs((sigAvg-sigmaGoal)/sigmaGoal)>5e-3) allStressesOK=false;
		cellGrow=(avgGrow/avgSize)*cellSize;
	}
	else{ // handle each dimension separately
		for(int axis=0; axis<3; axis++){
			// Δσ=ΔεE=(Δl/l)×(l×K/A) ↔ Δl=Δσ×A/K
			// FIXME: either NormShearPhys::{kn,ks} is computed wrong or we have dimensionality problem here
			// FIXME: that is why the fixup 1e-4 is needed here
			// FIXME: or perhaps maxDisplaPerStep=1e-2*charLen is too big??
			cellGrow[axis]=1e-4*(sigmaGoal-sigma[axis])*cellArea[axis]/(avgStiffness>0?avgStiffness:1);  // FIXME: wrong dimensions? See PeriTriaxController
			if(std::abs(cellGrow[axis])>maxDisplPerStep) cellGrow[axis]=Mathr::Sign(cellGrow[axis])*maxDisplPerStep;
			cellGrow[axis]=max(cellGrow[axis],-(cellSize[axis]-2.1*maxSpan));
			// crude way of predicting sigma, for steps when it is not computed from intrs
			if(avgStiffness>0) sigma[axis]+=cellGrow[axis]*avgStiffness; // FIXME: dimensions
			if(std::abs((sigma[axis]-sigmaGoal)/sigmaGoal)>5e-3) allStressesOK=false;
		}
	}
	TRVAR4(cellGrow,sigma,sigmaGoal,avgStiffness);
	assert(scene->dt>0);
	for(int axis=0; axis<3; axis++){ scene->cell->velGrad(axis,axis)=cellGrow[axis]/(scene->dt*scene->cell->getSize()[axis]); }

	// handle state transitions
	if(allStressesOK){
		if((step%globalUpdateInt)==0) currUnbalanced=Shop::unbalancedForce(/*useMaxForce=*/false,scene);
		if(currUnbalanced<maxUnbalanced){
			state+=1;
			// sigmaGoal reached and packing stable
			if(state==stresses.size()){ // no next stress to go for
				LOG_INFO("Finished");
				if(!doneHook.empty()){ LOG_DEBUG("Running doneHook: "<<doneHook); pyRunString(doneHook); }
			} else { LOG_INFO("Loaded to "<<sigmaGoal<<" done, going to "<<stresses[state]<<" now"); }
		} else {
			if((step%globalUpdateInt)==0) LOG_DEBUG("Stress="<<sigma<<", goal="<<sigmaGoal<<", unbalanced="<<currUnbalanced);
		}
	}
}

void PeriTriaxController::strainStressStiffUpdate(){
	//"Natural" strain, still correct for large deformations, used for comparison with goals
	for (int i=0;i<3;i++) strain[i]=log(scene->cell->trsf(i,i));

	//Compute volume of the deformed cell
	Real volume=scene->cell->hSize.determinant();

	//Compute sum(fi*lj) and stiffness
	stressTensor = Matrix3r::Zero();
	Vector3r sumStiff(Vector3r::Zero());
	int n=0;
	// NOTE : This sort of loops on interactions could be removed if we had callbacks in e.g. constitutive laws
	// → very likely performance hit; do you have some concrete design in mind?
	// → a vector with functors so we can law->functs->pushback(myThing), and access to the fundamental members (forces, stiffness, normal, etc.). Implementing the second part is less clear in my mind. Inheriting from law::funct(force&, stiffness&, ...)?
	FOREACH(const shared_ptr<Interaction>&I, *scene->interactions){
		if ( !I->isReal() ) continue;
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );
		GenericSpheresContact* gsc=SUDODEM_CAST<GenericSpheresContact*> ( I->geom.get() );
		//Contact force
		Vector3r f= (-1.)*( nsi->normalForce+nsi->shearForce );
		Vector3r branch=Body::byId(I->getId2(),scene)->state->pos + scene->cell->hSize*I->cellDist.cast<Real>() -Body::byId(I->getId1(),scene)->state->pos;
		stressTensor+=f*branch.transpose();
		if( !dynCell ){
			for ( int i=0; i<3; i++ ) sumStiff[i]+=std::abs ( gsc->normal[i] ) *nsi->kn+ ( 1-std::abs ( gsc->normal[i] ) ) *nsi->ks;
			n++;}
	}
	// Divide by volume as in stressTensor=sum(fi*lj)/Volume (Love equation)
	stressTensor /= volume;
	for(int axis=0; axis<3; axis++) stress[axis]=stressTensor(axis,axis);
	LOG_DEBUG ( "stressTensor : "<<endl
				<<stressTensor(0,0)<<" "<<stressTensor(0,1)<<" "<<stressTensor(0,2)<<endl
				<<stressTensor(1,0)<<" "<<stressTensor(1,1)<<" "<<stressTensor(1,2)<<endl
				<<stressTensor(2,0)<<" "<<stressTensor(2,1)<<" "<<stressTensor(2,2)<<endl
				<<"unbalanced = "<<Shop::unbalancedForce ( /*useMaxForce=*/false,scene ) );

	if (n>0) stiff=(1./n)*sumStiff;
	else stiff=Vector3r::Zero();
}


CREATE_LOGGER(PeriTriaxController);



void PeriTriaxController::action()
{
	if (!scene->isPeriodic){ throw runtime_error("PeriTriaxController run on aperiodic simulation."); }
	const Vector3r& cellSize=scene->cell->getSize();
	//FIXME : this is wrong I think (almost sure, B.)
	Vector3r cellArea=Vector3r(cellSize[1]*cellSize[2],cellSize[0]*cellSize[2],cellSize[0]*cellSize[1]);
	// initial updates
	if (maxBodySpan[0]<=0){
		FOREACH(const shared_ptr<Body>& b,*scene->bodies){
			if(!b || !b->bound) continue;
			for(int i=0; i<3; i++) if ((b->bound->max[i]-b->bound->min[i])<cellSize[i]) maxBodySpan[i]=max(maxBodySpan[i],b->bound->max[i]-b->bound->min[i]);}
	}
	// check current size
	if(2.1*maxBodySpan[0]>cellSize[0] || 2.1*maxBodySpan[1]>cellSize[1] || 2.1*maxBodySpan[2]>cellSize[2]){
		LOG_DEBUG("cellSize="<<cellSize<<", maxBodySpan="<<maxBodySpan);
		throw runtime_error("Minimum cell size is smaller than 2.1*maxBodySpan (periodic collider requirement)");}
	bool doUpdate((scene->iter%globUpdate)==0);
	if(doUpdate || min(stiff[0],min(stiff[1],stiff[2])) <=0 || dynCell){ strainStressStiffUpdate(); }

	// set mass to be sum of masses, if not set by the user
	if(dynCell && isnan(mass)){
		mass=0; FOREACH(const shared_ptr<Body>& b, *scene->bodies){ if(b && b->state) mass+=b->state->mass; }
		LOG_INFO("Setting cell mass to "<<mass<<" automatically.");}
	bool allOk=true;
	// apply condition along each axis separately (stress or strain)
	assert(scene->dt>0.);
	for(int axis=0; axis<3; axis++){
 		Real& strain_rate = scene->cell->velGrad(axis,axis);//strain rate on axis
		if(stressMask & (1<<axis)){   // control stress
			if(!dynCell){
				// stiffness K=EA; σ₁=goal stress; Δσ wanted stress difference to apply
				// ΔεE=(Δl/l)(K/A) - Grow is Δε, obtained by imposing the strain rate Δε/dt
				strain_rate=1/scene->dt*(goal[axis]-stress[axis])*cellArea[axis]/(stiff[axis]>0?stiff[axis]:1.);
				LOG_TRACE(axis<<": stress="<<stress[axis]<<", goal="<<goal[axis]<<", cellGrow="<<strain_rate*scene->dt);
			} else {  //accelerate the deformation using the density of the period, includes Cundall's damping
				assert( mass>0 );//user set
				Real dampFactor = 1 - growDamping*Mathr::Sign ( strain_rate * ( goal[axis]-stress[axis] ) );
				strain_rate+=dampFactor*scene->dt* ( goal[axis]-stress[axis] ) /mass;
				LOG_TRACE ( axis<<": stress="<<stress[axis]<<", goal="<<goal[axis]<<", velGrad="<<strain_rate );}

		} else {    // control strain, see "true strain" definition here http://en.wikipedia.org/wiki/Finite_strain_theory
			///NOTE : everything could be generalized to 9 independant components by comparing F[i,i] vs. Matrix3r goal[i,i], but it would be simpler in that case to let the user set the prescribed loading rates velGrad[i,i] when [i,i] is not stress-controlled. This "else" would disappear.
			strain_rate = (exp ( goal[axis]-strain[axis] ) -1)/scene->dt;
			LOG_TRACE ( axis<<": strain="<<strain[axis]<<", goal="<<goal[axis]<<", cellGrow="<<strain_rate*scene->dt);
		}
		// steady evolution with fluctuations; see TriaxialStressController
		if (!dynCell) strain_rate=(1-growDamping)*strain_rate+.8*prevGrow[axis];
		// limit maximum strain rate
		if (std::abs(strain_rate)>maxStrainRate[axis]) strain_rate = Mathr::Sign(strain_rate)*maxStrainRate[axis];
		// do not shrink below minimum cell size (periodic collider condition), although it is suboptimal WRT resulting stress
		strain_rate=max(strain_rate,-(cellSize[axis]-2.1*maxBodySpan[axis])/scene->dt);

		// crude way of predicting stress, for steps when it is not computed from intrs
		if(doUpdate) LOG_DEBUG(axis<<": cellGrow="<<strain_rate*scene->dt<<", new stress="<<stress[axis]<<", new strain="<<strain[axis]);
		// used only for temporary goal comparisons. The exact value is assigned in strainStressStiffUpdate
		strain[axis]+=strain_rate*scene->dt;
		// signal if condition not satisfied
		if(stressMask&(1<<axis)){
			Real curr=stress[axis];
			if((goal[axis]!=0 && std::abs((curr-goal[axis])/goal[axis])>relStressTol) || std::abs(curr-goal[axis])>absStressTol) allOk=false;
		}else{
			Real curr=strain[axis];
			// since strain is prescribed exactly, tolerances need just to accomodate rounding issues
			if((goal[axis]!=0 && std::abs((curr-goal[axis])/goal[axis])>1e-6) || std::abs(curr-goal[axis])>1e-6){
				allOk=false;
				if(doUpdate) LOG_DEBUG("Strain not OK; "<<std::abs(curr-goal[axis])<<">1e-6");}
		}
	}
	// update stress and strain
	if (!dynCell) for ( int axis=0; axis<3; axis++ ){
		// take in account something like poisson's effect here…
		//Real bogusPoisson=0.25; int ax1=(axis+1)%3,ax2=(axis+2)%3;
		//don't modify stress if dynCell, testing only stiff[axis]>0 would not allow switching the control mode in simulations,
		if (stiff[axis]>0) stress[axis]+=(scene->cell->velGrad(axis,axis)*scene->dt/cellSize[axis])*(stiff[axis]/cellArea[axis]);
		//-bogusPoisson*(cellGrow[ax1]/refSize[ax1])*(stiff[ax1]/cellArea[ax1])-bogusPoisson*(cellGrow[ax2]/refSize[ax2])*(stiff[ax2]/cellArea[ax2]);
	}
 	for (int k=0;k<3;k++) strainRate[k]=scene->cell->velGrad(k,k);
	//Update energy input FIXME: replace trace by norm, so it works for any kind of deformation
	Real dW=(0.5*(scene->cell->prevVelGrad + scene->cell->velGrad)*stressTensor).trace()*scene->dt*(scene->cell->hSize.determinant());
	externalWork+=dW;
	if(scene->trackEnergy) scene->energy->add(-dW,"velGradWork",velGradWorkIx,/*non-incremental*/false);
	prevGrow = strainRate;

	if(allOk){
		if(doUpdate || currUnbalanced<0){
			currUnbalanced=Shop::unbalancedForce(/*useMaxForce=*/false,scene);
			LOG_DEBUG("Stress/strain="<< (stressMask&1?stress[0]:strain[0]) <<"," <<(stressMask&2?stress[1]:strain[1])<<"," <<(stressMask&4?stress[2]:strain[2]) <<", goal="<<goal<<", unbalanced="<<currUnbalanced );}
		if(currUnbalanced<maxUnbalanced){
			// LOG_INFO("Goal reached, packing stable.");
			if (!doneHook.empty()){
				LOG_DEBUG ( "Running doneHook: "<<doneHook );
				pyRunString(doneHook);}
// 			else { Omega::instance().pause(); }
		}
	}
}

CREATE_LOGGER(Peri3dController);
void Peri3dController::action(){
	if(!scene->isPeriodic){ LOG_FATAL("Being used on non-periodic simulation!"); throw; }
	const Real& dt=scene->dt;
	assert(dt>0);

	/* "Constructor" (if (step==0) )
			ps is the vector of indices, where stress is prescribed
			pe is the vector of indices, where strain is prescribed
		 	example: goal = 0b000110 : ps=(1,2,0,0,0,0), pe=(0,3,4,5,0,0)
			lenPs (lenPe) is the meaningful length of ps (pe) (the zeros at the end of ps and pe has no meaning),
				i.e. the number of indices with prescribed stress (strain)
	*/
	bool stressBasedSimulation=false; // true when all stresses are prescribed or if all prescribed strains equal zero
	if (progress==0) {
		lenPs=0; lenPe=0;
		ps=Vector6i::Zero();
		pe=Vector6i::Zero();
		stressGoal = Vector6r::Zero();
		strainGoal = Vector6r::Zero();
		for (int i=0; i<6; i++){
			if (stressMask&(1<<i)){ // if stress is prescribed at direction i, add this direction to ps and increase lenPs by one
				ps(lenPs++)=i;
				stressGoal(i) = goal(i);
			} else{ // if strain is prescribed at direction i, add this direction to pe and increase lenPe by one
				pe(lenPe++)=i;
				strainGoal(i) = goal(i);
			}
		}

		// variables used in evaluation of ideal stress and ideal strain for each part defined by ##Path
		//paths[0]=&xxPath; paths[1]=&yyPath; paths[2]=&zzPath; paths[3]=&yzPath; paths[4]=&zxPath; paths[5]=&xyPath; // pointers to the Paths
		pathSizes[0]=xxPath.size(); pathSizes[1]=yyPath.size(); pathSizes[2]=zzPath.size();
		pathSizes[3]=yzPath.size(); pathSizes[4]=zxPath.size(); pathSizes[5]=xyPath.size();
		for (int i=0; i<6; i++) {pathsCounter[i] = 0;} // inidicator in which part of the path we are

		// abcPath[j] is j-th Vector2r in path[0]
		// PATH_OP_OP(0,j,k) = path[0]->operator[](j).operator(k) is k-th element of j-th Vector2r of xxPath
		//#define PATH_OP_OP(pi,op1i,op2i) paths[pi]->operator[](op1i).operator()(op2i)
		//#define PATH_OP_OP(pi,op1i,op2i) (pi==0?xxPath:pi==1?yyPath:pi==2?zzPath:pi==3?yzPath:pi==4?zxPath:pi==5?xyPath:NULL)->operator[](op1i).operator()(op2i)
		#define PATH_OP_OP(pi,op1i,op2i) (pi==0?xxPath:pi==1?yyPath:pi==2?zzPath:pi==3?yzPath:pi==4?zxPath:xyPath)[op1i].operator()(op2i)

		for (int i=0; i<6; i++) {
			for (int j=1; j<pathSizes[i]; j++) {
				// check if the user defined time axis is monothonically increasing
				{ if ( PATH_OP_OP(i,j-1,0) >= PATH_OP_OP(i,j,0) ) {
					throw runtime_error("Peri3dCoontroller.##Path: Time in ##Path must be monothonically increasing");
				}}
			}
			for (int j=0; j<pathSizes[i]; j++) {
				// convert relative progress values of ##Path to absolute values
				PATH_OP_OP(i,j,0) *= 1./PATH_OP_OP(i,pathSizes[i]-1,0);
				// convert relative stress/strain values of ##Path to absolute stress strain values
				if (std::abs(PATH_OP_OP(i,pathSizes[i]-1,1)) >= 1e-9) { // the final value is not 0 (otherwise always absolute values are considered)
					PATH_OP_OP(i,j,1) *= goal(i)/PATH_OP_OP(i,pathSizes[i]-1,1);
				}
			}
		}

		// set weather the simulation is "stress based" (all stress components prescribed or all prescribed strains equal zero)
		if (lenPe == 0) { stressBasedSimulation = true; }
		else {
			stressBasedSimulation = true;
			for (int i=0; i<lenPe; i++) { stressBasedSimulation = stressBasedSimulation && PATH_OP_OP(pe(i),1,1)<1e9; }
		}
	}

	// increase the pathCounter by one if we cross to the next part of path
	for (int i=0; i<6; i++) {
		if (progress >= PATH_OP_OP(i,pathsCounter[i],0)) { pathsCounter[i]++; }
	}

	/* values of prescribed stress (strain) rate in respect to prescribed path.
	   The strain indices where stress is prescribed will be overwritten by predictor */
	for (int i=0; i<lenPe; i++){
		int j = pe(i);
		if (pathSizes[j] == 1) { // path has only one part (only final values are prescribed)
			strainRate(j) = strainGoal(j)/(nSteps*dt); // ideal strain rate in respect of dSteps and dValue
		}
		else if (pathsCounter[j] == 0) { // path has more parts, but we are still at the first one
			const Real& dProgress = PATH_OP_OP(j,0,0); // progress difference of respective part of the path
			const Real& dValue = PATH_OP_OP(j,0,1); // strain difference at the respective part of the path
			strainRate(j) = dValue/(dProgress*nSteps*dt); // ideal strain rate in respect of dSteps and dValue
		}
		else if (progress < 1.) {
			const Real dProgress = PATH_OP_OP(j,pathsCounter[j],0) - PATH_OP_OP(j,pathsCounter[j]-1,0); // progress difference of respective part of the path
			const Real dValue = PATH_OP_OP(j,pathsCounter[j],1) - PATH_OP_OP(j,pathsCounter[j]-1,1); // strain difference at the respective part of the path
			strainRate(j) = dValue/(dProgress*nSteps*dt); // ideal strain rate in respect of dSteps and dValue
		}
		else { strainRate(j) = 0; }
	}
	for (int i=0; i<lenPs; i++){
		int j = ps(i);
		if (pathSizes[j] == 1) { // path has only one part (only final values are prescribed)
			stressRate(j) = stressGoal(j)/(nSteps*dt); // ideal stress rate in respect of dSteps and dValue
		}
		else if (pathsCounter[j] == 0) { // path has more parts, but we are still at the first one
			const Real& dProgress = PATH_OP_OP(j,0,0); // progress difference of respective part of the path
			const Real& dValue = PATH_OP_OP(j,0,1); // stress difference at the respective part of the path
			stressRate(j) = dValue/(dProgress*nSteps*dt); // ideal stress rate in respect of dSteps and dValue
		}
		else if (progress < 1.) {
			const Real dProgress = PATH_OP_OP(j,pathsCounter[j],0) - PATH_OP_OP(j,pathsCounter[j]-1,0); //  progress difference of respective part of the path
			const Real dValue = PATH_OP_OP(j,pathsCounter[j],1) - PATH_OP_OP(j,pathsCounter[j]-1,1); // stress difference at the respective part of the path
			stressRate(j) = dValue/(dProgress*nSteps*dt); // ideal stress rate in respect of dSteps and dValue
		}
		else { stressRate(j) = 0; }
	}

	// Update - update values from previous step to current step
	stressOld = stress; // stresssOld = stress at previous step
	//sigma = Shop::stressTensorOfPeriodicCell(/*smallStrains=*/true); // current stress tensor
	sigma = Shop::getStress(); // current stress tensor
	stress = tensor_toVoigt(sigma); // current stress vector
	stressIdeal += stressRate*dt; // stress that would be obtained if the predictor would be perfect
	strain += strainRate*dt; // current strain vector
	epsilon = voigt_toSymmTensor(strain,/*strain=*/true); // current strain tensor

	/* StrainPredictor
			extremely primitive predictor, but roboust enough and works fine :-) could be replaced by some more rigorous in future.
			In the direction with prescribed strain rate this prescribed strain rate is applied.
			In direction with prescribed stress: from values of stress and strain in previous two steps the value of strain rate
			is predicted so as the stress in the next step would be as close as possible to the ideal stress value,
			see the documentation for more info
	*/
	if (lenPs > 0){ // if at least 1 stress component is prescribed (otherwise prescribed strain is applied in all 6 directions
		if (progress == 0 && stressBasedSimulation) { // the very first step, use compliance estimation (compliance matrix for elastic isotropic material)
			Real complianceEstimation[6][6] = {
				{1/youngEstimation, -poissonEstimation/youngEstimation, -poissonEstimation/youngEstimation, 0,0,0},
				{-poissonEstimation/youngEstimation, 1/youngEstimation, -poissonEstimation/youngEstimation, 0,0,0},
				{-poissonEstimation/youngEstimation, -poissonEstimation/youngEstimation, 1/youngEstimation, 0,0,0},
				{0,0,0, (1+poissonEstimation)/youngEstimation,0,0},
				{0,0,0, 0,(1+poissonEstimation)/youngEstimation,0},
				{0,0,0, 0,0,(1+poissonEstimation)/youngEstimation}};
			for (int i=0; i<lenPs; i++) {
				strainRate(ps(i)) = 0;
				for (int j=0; j<lenPs; j++) { strainRate(ps(i)) += complianceEstimation[ps(i)][ps(j)]*stressRate[ps(j)]; }
				for (int j=0; j<lenPe; j++) { strainRate(ps(i)) += complianceEstimation[ps(i)][ps(j)]*stressRate[ps(j)]; }
			}
			//for (int i=0; i<lenPs; i++) { strainRate(ps(i)) -= maxStrainRate; }
		}
		else { // actual predictor
			Real sr=strainRate.cwiseAbs().maxCoeff();
			for (int i=0; i<lenPs; i++) {
				int j=ps(i);
				// linear extrapolation of stress error (difference between actual and ideal stress)
				Real linPred = 2*(stress(j)-stressIdeal(j)) - (stressOld(j)-(stressIdeal(j)-stressRate(j)*dt));
				// correction of strain in respect to the extrapolated stress error
				if (linPred>0){strainRate(j) -= sr*mod;}
				else {strainRate(j) += sr*mod;}
			}
		}
	}

	// correction coefficient ix strainRate.maxstd::abs() > maxStrainRate
	Real srCorr = (strainRate.cwiseAbs().maxCoeff() > maxStrainRate)? (maxStrainRate/strainRate.cwiseAbs().maxCoeff()):1.;
	strainRate *= srCorr;

	// Actual action (see the documentation for more info)
	const Matrix3r& trsf=scene->cell->trsf;
	// compute rotational and nonrotational (strain in local coordinates) part of trsf
	Matrix_computeUnitaryPositive(trsf,&rot,&nonrot);

	// prescribed velocity gradient (strain tensor rate) in global coordinates
	epsilonRate = voigt_toSymmTensor(strainRate,/*strain=*/true);
	/* transformation of prescribed strain rate (computed by predictor) into local cell coordinates,
	   multiplying by time to obtain strain increment and adding it to nonrot (current strain in local coordinates)*/
	nonrot += rot.transpose()*(epsilonRate*dt)*rot;
	Matrix3r& velGrad=scene->cell->velGrad;
	// compute velGrad:
	//   trsf = rot*nonrot
	//   dTrsf = dt*velGrad
	//   trsfNew = trsf + dTrsf*trsf = (I+dTrsf)*trsf = (I+dt*velGrad)*trsf   ->   velGrad
	velGrad = ((rot*nonrot)*trsf.inverse()- Matrix3r::Identity()) / dt ;
	progress += srCorr/nSteps;

	if (progress >= 1. || strain.cwiseAbs().maxCoeff() > maxStrain) {
		if(doneHook.empty()){ LOG_INFO("No doneHook set, dying."); dead=true; Omega::instance().pause(); }
		else{ LOG_INFO("Running doneHook: "<<doneHook);	pyRunString(doneHook);}
	}
}
