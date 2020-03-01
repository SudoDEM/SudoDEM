
// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#include<sudodem/pkg/dem/PeriTriaxial.hpp>
#include<sudodem/pkg/dem/Shop.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/NormShearPhys.hpp>
//#include<sudodem/pkg/dem/DemXDofGeom.hpp>
#include<sudodem/lib/pyutil/gil.hpp>

SUDODEM_PLUGIN((PeriTriaxController))


void PeriTriaxController::strainStressStiffUpdate(){
	//"Natural" strain, still correct for large deformations, used for comparison with goals
	for (int i=0;i<2;i++) strain[i]=log(scene->cell->trsf(i,i));

	//Compute volume of the deformed cell
	Real volume=scene->cell->hSize.determinant()*z_dim;//uinit one in the z direction

	//Compute sum(fi*lj) and stiffness
	stressTensor = Matrix2r::Zero();
	Vector2r sumStiff(Vector2r::Zero());
	int n=0;
	// NOTE : This sort of loops on interactions could be removed if we had callbacks in e.g. constitutive laws
	// → very likely performance hit; do you have some concrete design in mind?
	// → a vector with functors so we can law->functs->pushback(myThing), and access to the fundamental members (forces, stiffness, normal, etc.). Implementing the second part is less clear in my mind. Inheriting from law::funct(force&, stiffness&, ...)?
	FOREACH(const shared_ptr<Interaction>&I, *scene->interactions){
		if ( !I->isReal() ) continue;
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );
		//Contact force
		Vector2r f= (-1.)*( nsi->normalForce+nsi->shearForce );
		Vector2r branch=Body::byId(I->getId2(),scene)->state->pos + scene->cell->hSize*I->cellDist.cast<Real>() -Body::byId(I->getId1(),scene)->state->pos;
		stressTensor+=f*branch.transpose();
		if( !dynCell ){//not use stiffness for 2d temporarily
			//GenericDisksContact* gsc=SUDODEM_CAST<GenericDisksContact*> ( I->geom.get() );
			//for ( int i=0; i<2; i++ ) sumStiff[i]+=std::abs ( gsc->normal[i] ) *nsi->kn+ ( 1-std::abs ( gsc->normal[i] ) ) *nsi->ks;
			//n++;
		  }
	}
	// Divide by volume as in stressTensor=sum(fi*lj)/Volume (Love equation)
	stressTensor /= volume;
	for(int axis=0; axis<2; axis++) stress[axis]=stressTensor(axis,axis);
	LOG_DEBUG ( "stressTensor : "<<endl
				<<stressTensor(0,0)<<" "<<stressTensor(0,1)<<endl
				<<stressTensor(1,0)<<" "<<stressTensor(1,1)<<endl
				<<"unbalanced = "<<Shop::unbalancedForce ( /*useMaxForce=*/false,scene ) );

	if (n>0) stiff=(1./n)*sumStiff;
	else stiff=Vector2r::Zero();
}


CREATE_LOGGER(PeriTriaxController);



void PeriTriaxController::action()
{
	if (!scene->isPeriodic){ throw runtime_error("PeriTriaxController run on aperiodic simulation."); }
	const Vector2r& cellSize=scene->cell->getSize();
	Vector2r cellArea=Vector2r(cellSize[1],cellSize[0]);//unit one along the z direction
	// initial updates
	if (maxBodySpan[0]<=0){
		FOREACH(const shared_ptr<Body>& b,*scene->bodies){
			if(!b || !b->bound) continue;
			for(int i=0; i<2; i++) if ((b->bound->max[i]-b->bound->min[i])<cellSize[i]) maxBodySpan[i]=max(maxBodySpan[i],b->bound->max[i]-b->bound->min[i]);}
	}
	// check current size
	if(2.1*maxBodySpan[0]>cellSize[0] || 2.1*maxBodySpan[1]>cellSize[1]){
		LOG_DEBUG("cellSize="<<cellSize<<", maxBodySpan="<<maxBodySpan);
		throw runtime_error("Minimum cell size is smaller than 2.1*maxBodySpan (periodic collider requirement)");}
	bool doUpdate((scene->iter%globUpdate)==0);
	if(doUpdate || min(stiff[0],stiff[1]) <=0 || dynCell){ strainStressStiffUpdate(); }

	// set mass to be sum of masses, if not set by the user
	if(dynCell && isnan(mass)){
		mass=0; FOREACH(const shared_ptr<Body>& b, *scene->bodies){ if(b && b->state) mass+=b->state->mass; }
		LOG_INFO("Setting cell mass to "<<mass<<" automatically.");}
	bool allOk=true;
	// apply condition along each axis separately (stress or strain)
	assert(scene->dt>0.);
	for(int axis=0; axis<2; axis++){
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
			///NOTE : everything could be generalized to 9 independant components by comparing F[i,i] vs. Matrix2r goal[i,i], but it would be simpler in that case to let the user set the prescribed loading rates velGrad[i,i] when [i,i] is not stress-controlled. This "else" would disappear.
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
			if((goal[axis]!=0 && std::abs((curr-goal[axis])/goal[axis])>relStressTol) || std::abs(curr-goal[axis])>absStressTol){
				//cout<<std::abs((curr-goal[axis])/goal[axis])<<", "<< std::abs(curr-goal[axis])<<endl;
				allOk=false;
			}
		}else{
			Real curr=strain[axis];
			// since strain is prescribed exactly, tolerances need just to accomodate rounding issues
			if((goal[axis]!=0 && std::abs((curr-goal[axis])/goal[axis])>1e-6) || std::abs(curr-goal[axis])>1e-6){
				allOk=false;
				if(doUpdate) LOG_DEBUG("Strain not OK; "<<std::abs(curr-goal[axis])<<">1e-6");}
		}
	}
	// update stress and strain
	if (!dynCell) for ( int axis=0; axis<2; axis++ ){
		// take in account something like poisson's effect here…
		//Real bogusPoisson=0.25; int ax1=(axis+1)%3,ax2=(axis+2)%3;
		//don't modify stress if dynCell, testing only stiff[axis]>0 would not allow switching the control mode in simulations,
		if (stiff[axis]>0) stress[axis]+=(scene->cell->velGrad(axis,axis)*scene->dt/cellSize[axis])*(stiff[axis]/cellArea[axis]);
		//-bogusPoisson*(cellGrow[ax1]/refSize[ax1])*(stiff[ax1]/cellArea[ax1])-bogusPoisson*(cellGrow[ax2]/refSize[ax2])*(stiff[ax2]/cellArea[ax2]);
	}
 	for (int k=0;k<2;k++) strainRate[k]=scene->cell->velGrad(k,k);
	//Update energy input FIXME: replace trace by norm, so it works for any kind of deformation
	Real dW=(0.5*(scene->cell->prevVelGrad + scene->cell->velGrad)*stressTensor).trace()*scene->dt*(scene->cell->hSize.determinant());
	externalWork+=dW;
	if(scene->trackEnergy) scene->energy->add(-dW,"velGradWork",velGradWorkIx,/*non-incremental*/false);
	prevGrow = strainRate;
	//cout<<"allok?"<<allOk<<endl;
	if(allOk){
		if(doUpdate || currUnbalanced<0){
			currUnbalanced=Shop::unbalancedForce(/*useMaxForce=*/false,scene);
			cout<<"Stress/strain="<< (stressMask&1?stress[0]:strain[0]) <<"," <<(stressMask&2?stress[1]:strain[1])<<"," <<", goal="<<goal<<", unbalanced="<<currUnbalanced<<endl;;}
		if(currUnbalanced<maxUnbalanced){
			// LOG_INFO("Goal reached, packing stable.");
			if (!doneHook.empty()){
				LOG_DEBUG ( "Running doneHook: "<<doneHook );
				pyRunString(doneHook);}
			else {
				cout<<"Completed!"<<endl;
				Omega::instance().pause(); 
			 }
		}
	}
}
