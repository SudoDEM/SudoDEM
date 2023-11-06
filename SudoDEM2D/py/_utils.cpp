#include<sudodem/pkg/dem/Shop.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/dem/Superellipse.hpp>
#include<sudodem/pkg/common/Disk.hpp>
//#include<sudodem/pkg/dem/ScGeom.hpp>
//#include<sudodem/pkg/dem/DemXDofGeom.hpp>
//#include<sudodem/pkg/common/Facet.hpp>


//#include<sudodem/pkg/common/NormShearPhys.hpp>

#include<sudodem/lib/pyutil/doc_opts.hpp>
//#include<sudodem/pkg/dem/ViscoelasticPM.hpp>

#include<numpy/ndarrayobject.h>


namespace py = boost::python;
Real Shop__unbalancedForce(bool useMaxForce /*false by default*/){return Shop::unbalancedForce(useMaxForce);}
py::tuple Shop__getStressAndTangent2D(Real z_dim=1, bool symmetry=true){return Shop::getStressAndTangent2D(z_dim,symmetry);}
py::tuple Shop__getStressTangentThermal2D(Real z_dim=1, bool symmetry=true){return Shop::getStressTangentThermal2D(z_dim,symmetry);}
//mean coordination number
Real Shop__meanCoordinationNumber(){
	const shared_ptr<Scene> scene=Omega::instance().getScene();
	int num_contacts = 0;
	int num_particles = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );
		//Contact force
		if(nsi->normalForce.norm() == 0) continue;//make sure it is a valid contact//FIXME:using euqality test is not a good choice
		num_contacts += 1;
	}
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b) continue;
		num_particles += 1;
	}
	return 2.0*num_contacts/num_particles;
}

bool Shop__changeParticleSize2D(double alpha){
	const shared_ptr<Scene> scene=Omega::instance().getScene();
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b) continue;
		if(b->shape->getClassName() == "Superellipse") {
				const shared_ptr<Superellipse>& A = SUDODEM_PTR_CAST<Superellipse> (b->shape);
				if(!A) break;
				//vol += A->getArea();
				A->rxy = A->ref_rxy*alpha;
				A->rx = A->ref_rxy(0)*alpha;
				A->ry = A->ref_rxy(1)*alpha;

		}else if(b->shape->getClassName() == "Disk"){
				const shared_ptr<Disk>& A = SUDODEM_PTR_CAST<Disk> (b->shape);
				if(!A) break;
				A->radius = A->ref_radius*alpha;
		}
		}
	return true;
}

//actually particle area for 2D
Real Shop__getParticlesVolume2D(){
	const shared_ptr<Scene> scene=Omega::instance().getScene();
	Real vol=0;
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b) continue;
		if(b->shape->getClassName() == "Superellipse") {
				const shared_ptr<Superellipse>& A = SUDODEM_PTR_CAST<Superellipse> (b->shape);
				if(!A) break;
				vol += A->getArea();

		}else if(b->shape->getClassName() == "Disk"){
				const shared_ptr<Disk>& A = SUDODEM_PTR_CAST<Disk> (b->shape);
				if(!A) break;
				vol += Mathr::PI*pow(A->radius,2);
		}
		}
	return vol;
}

//void ratio
Real Shop__getVoidRatio2D( Real cellArea = 1){
	const shared_ptr<Scene> scene= Omega::instance().getScene();
	Real V;
	if(!scene->isPeriodic){
		V = cellArea;
	} else {
		V=scene->cell->getVolume();//it is area actually for 2D //dterminant of hSize of cell may be negative, but here it is impossible.
	}

	Real Vs=Shop__getParticlesVolume2D();
	return (V-Vs)/Vs;
}

//get stress tensor and tangent operator tensor.
Matrix2r Shop__getStress2D(Real z_dim){
	Scene* scene=Omega::instance().getScene().get();
	if(z_dim == 0) z_dim = 1;
	Real volume = scene->isPeriodic?scene->cell->hSize.determinant()*z_dim:1;
	Matrix2r stressTensor = Matrix2r::Zero();
	const bool isPeriodic = scene->isPeriodic;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		//not considering clumped particles. CAUTION:SudoDEM does not prefer clumped particles.
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );

		//Contact force
		Vector2r fn = nsi->normalForce;
		if (fn.norm() == 0)continue;
		Vector2r ft = nsi->shearForce;
		Vector2r f= fn + ft ;//here we constructe a positive stress tensor
		Real kN=nsi->kn;
		Real kT=nsi->ks;
		Vector2r n = fn.normalized();
		//cout<<"n"<<n<<endl;
		Vector2r t = ft.normalized();

		Vector2r branch=Body::byId(I->getId1(),scene)->state->pos - Body::byId(I->getId2(),scene)->state->pos;
		if (isPeriodic) branch -= scene->cell->hSize*I->cellDist.cast<Real>();//shifting on particle 2

		stressTensor+=f*branch.transpose();
	}
	stressTensor/=volume;
	return stressTensor;
}
//fabric tensor based on contact normal
Matrix2r Shop__getFabricTensorCN2D(){
	Scene* scene=Omega::instance().getScene().get();
	//Real volume = scene->isPeriodic?scene->cell->hSize.determinant()*z_dim:1;
	Matrix2r fabricTensor = Matrix2r::Zero();
	const bool isPeriodic = scene->isPeriodic;
	int num = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		//not considering clumped particles. CAUTION:SudoDEM does not prefer clumped particles.
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );

		//Contact force
		Vector2r fn = nsi->normalForce;
		if (fn.norm() == 0)continue;
		Vector2r n = fn.normalized();
		fabricTensor+=n*n.transpose();
		num++;
	}
	fabricTensor/=num;
	return fabricTensor;
}

//fabric tensor based on particle orientation
Matrix2r Shop__getFabricTensorPO2D(){
	Scene* scene=Omega::instance().getScene().get();
	//Real volume = scene->isPeriodic?scene->cell->hSize.determinant()*z_dim:1;
	Matrix2r fabricTensor = Matrix2r::Zero();
	const bool isPeriodic = scene->isPeriodic;
	int num = 0;
	FOREACH(shared_ptr<Body> b, *scene->bodies){
		if (!b) continue;
		Vector2r n = b->state->ori*Vector2r(1.0,0.0);
		fabricTensor+=n*n.transpose();
		num++;
		}
	fabricTensor/=num;
	return fabricTensor;
}


#ifdef SUDODEM_OPENGL
/*#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sudodem/gui/qt4/OpenGLManager.hpp>
void SaveSnapshot(){
	if(!OpenGLManager::self) throw logic_error("No OpenGLManager instance?!");
	const shared_ptr<GLViewer>& glv=OpenGLManager::self->views[0];
	glv->saveSnapshot(false,false);
}*/
#endif /* SUDODEM_OPENGL */

BOOST_PYTHON_MODULE(_utils){
	// http://numpy.scipy.org/numpydoc/numpy-13.html mentions this must be done in module init, otherwise we will crash
	//import_array();

	SUDODEM_SET_DOCSTRING_OPTS;
	py::def("unbalancedForce",&Shop__unbalancedForce,(py::args("useMaxForce")=false),"Compute the ratio of mean (or maximum, if *useMaxForce*) summary force on bodies and mean force magnitude on interactions. For perfectly static equilibrium, summary force on all bodies is zero (since forces from interactions cancel out and induce no acceleration of particles); this ratio will tend to zero as simulation stabilizes, though zero is never reached because of finite precision computation. Sufficiently small value can be e.g. 1e-2 or smaller, depending on how much equilibrium it should be.");
	py::def("getParticleVolume2D",Shop__getParticlesVolume2D,"Compute the total volume (area for 2D) of particles in the scene.");
	py::def("getMeanCN",Shop__meanCoordinationNumber,"Get the mean coordination number of a packing.");
	py::def("changeParticleSize2D",Shop__changeParticleSize2D,"expand a particle by a given coefficient");
	py::def("getVoidRatio2D", Shop__getVoidRatio2D, (py::args("cellArea")=1),"Compute 2D void ratio. Keyword cellArea is effective only for aperioidc cell.");
	py::def("getStress2D",Shop__getStress2D,(py::args("z_dim")=1),"Compute overall stress of periodic cell.");
	py::def("getFabricTensorCN2D",Shop__getFabricTensorCN2D,"Fabric tensor of contact normal.");
	py::def("getFabricTensorPO2D",Shop__getFabricTensorPO2D,"Fabric tensor of particle orientation along the rx axis.");
	py::def("getStressAndTangent2D",Shop__getStressAndTangent2D,(py::args("z_dim")=1,py::args("symmetry")=true),"Compute overall stress of periodic cell using the same equation as function getStress. In addition, the tangent operator is calculated using the equation published in [Kruyt and Rothenburg1998]_:\n\n.. math:: S_{ijkl}=\\frac{1}{V}\\sum_{c}(k_n n_i l_j n_k l_l + k_t t_i l_j t_k l_l)\n\n:param float volume: same as in function getStress\n:param bool symmetry: make the tensors symmetric.\n\n:return: macroscopic stress tensor and tangent operator as py::tuple");
	py::def("getStressTangentThermal2D",Shop__getStressTangentThermal2D,(py::args("z_dim")=1,py::args("symmetry")=true),"Compute overall stress of periodic cell using the same equation as function getStress. In addition, the tangent operator is calculated using the equation published in [Kruyt and Rothenburg1998]_:\n\n.. math:: S_{ijkl}=\\frac{1}{V}\\sum_{c}(k_n n_i l_j n_k l_l + k_t t_i l_j t_k l_l)\n\n:param float volume: same as in function getStress\n:param bool symmetry: make the tensors symmetric. Finally, the thermal conductivity tensor is calculated based on the formular in PFC.\n\n:return: macroscopic stress tensor and tangent operator as py::tuple");
	#ifdef SUDODEM_OPENGL
	//py::def("SaveSnapshot",SaveSnapshot,"TODO");
	#endif
}
