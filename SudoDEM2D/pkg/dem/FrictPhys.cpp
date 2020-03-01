#include "FrictPhys.hpp"
#include <sudodem/pkg/dem/ScGeom.hpp>
SUDODEM_PLUGIN((FrictPhys)(ViscoFrictPhys)(Ip2_FrictMat_FrictMat_ViscoFrictPhys)(Ip2_FrictMat_FrictMat_FrictPhys));

// The following code was moved from Ip2_FrictMat_FrictMat_FrictPhys.hpp

void Ip2_FrictMat_FrictMat_FrictPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<FrictMat>& mat1 = SUDODEM_PTR_CAST<FrictMat>(b1);
	const shared_ptr<FrictMat>& mat2 = SUDODEM_PTR_CAST<FrictMat>(b2);
	interaction->phys = shared_ptr<FrictPhys>(new FrictPhys());
	const shared_ptr<FrictPhys>& contactPhysics = SUDODEM_PTR_CAST<FrictPhys>(interaction->phys);
	/*
	Real Ra,Rb;//Vector2r normal;
	assert(dynamic_cast<ScGeom*>(interaction->geom.get()));//only in debug mode
	ScGeom* sphCont=SUDODEM_CAST<ScGeom*>(interaction->geom.get());
	Ra=sphCont->refR1>0?sphCont->refR1:sphCont->refR2;
	Rb=sphCont->refR2>0?sphCont->refR2:sphCont->refR1;

	interaction->phys = shared_ptr<FrictPhys>(new FrictPhys());
	const shared_ptr<FrictPhys>& contactPhysics = SUDODEM_PTR_CAST<FrictPhys>(interaction->phys);
	Real Ea 	= mat1->young;
	Real Eb 	= mat2->young;
	Real Va 	= mat1->poisson;
	Real Vb 	= mat2->poisson;

	//harmonic average of the two stiffnesses when (2*Ri*Ei) is the stiffness of a contact point on disk "i"
	Real Kn = 2*Ea*Ra*Eb*Rb/(Ea*Ra+Eb*Rb);
	//same for shear stiffness
	Real Ks = 2*Ea*Ra*Va*Eb*Rb*Vb/(Ea*Ra*Va+Eb*Rb*Vb);
	*/
	//harmonic average of the two stiffnesses
	Real Kn = 2.0*mat1->Kn*mat2->Kn/(mat1->Kn + mat2->Kn);
	Real Ks = 2.0*mat1->Ks*mat2->Ks/(mat1->Ks + mat2->Ks);
	Real frictionAngle = (!frictAngle) ? std::min(mat1->frictionAngle,mat2->frictionAngle) : (*frictAngle)(mat1->id,mat2->id,mat1->frictionAngle,mat2->frictionAngle);
	contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
	contactPhysics->kn = Kn;
	contactPhysics->ks = Ks;
};

void Ip2_FrictMat_FrictMat_ViscoFrictPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<FrictMat>& mat1 = SUDODEM_PTR_CAST<FrictMat>(b1);
	const shared_ptr<FrictMat>& mat2 = SUDODEM_PTR_CAST<FrictMat>(b2);
	interaction->phys = shared_ptr<ViscoFrictPhys>(new ViscoFrictPhys());
	const shared_ptr<ViscoFrictPhys>& contactPhysics = SUDODEM_PTR_CAST<ViscoFrictPhys>(interaction->phys);
/*
	Real Ea 	= mat1->young;
	Real Eb 	= mat2->young;
	Real Va 	= mat1->poisson;
	Real Vb 	= mat2->poisson;

	Real Ra,Rb;//Vector2r normal;
	assert(dynamic_cast<ScGeom*>(interaction->geom.get()));//only in debug mode
	ScGeom* sphCont=SUDODEM_CAST<ScGeom*>(interaction->geom.get());
	Ra=sphCont->refR1>0?sphCont->refR1:sphCont->refR2;
	Rb=sphCont->refR2>0?sphCont->refR2:sphCont->refR1;

	//harmonic average of the two stiffnesses when (Ri.Ei/2) is the stiffness of a contact point on disk "i"
	Real Kn = 2*Ea*Ra*Eb*Rb/(Ea*Ra+Eb*Rb);
	//same for shear stiffness
	Real Ks = 2*Ea*Ra*Va*Eb*Rb*Vb/(Ea*Ra*Va+Eb*Rb*Vb);
*/
	Real Kn = 2.0*mat1->Kn*mat2->Kn/(mat1->Kn + mat2->Kn);
	Real Ks = 2.0*mat1->Ks*mat2->Ks/(mat1->Ks + mat2->Ks);
	Real frictionAngle = (!frictAngle) ? std::min(mat1->frictionAngle,mat2->frictionAngle) : (*frictAngle)(mat1->id,mat2->id,mat1->frictionAngle,mat2->frictionAngle);
	contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
	contactPhysics->kn = Kn;
	contactPhysics->ks = Ks;
};
