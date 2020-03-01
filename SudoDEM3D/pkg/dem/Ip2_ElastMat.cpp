#include "Ip2_ElastMat.hpp"

#include <sudodem/pkg/common/NormShearPhys.hpp>
#include <sudodem/pkg/dem/DemXDofGeom.hpp>

SUDODEM_PLUGIN((Ip2_ElastMat_ElastMat_NormPhys)(Ip2_ElastMat_ElastMat_NormShearPhys));


void Ip2_ElastMat_ElastMat_NormPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<ElastMat>& mat1 = SUDODEM_PTR_CAST<ElastMat>(b1);
	const shared_ptr<ElastMat>& mat2 = SUDODEM_PTR_CAST<ElastMat>(b2);
	Real Ea 	= mat1->young;
	Real Eb 	= mat2->young;
	interaction->phys = shared_ptr<NormPhys>(new NormPhys());
	const shared_ptr<NormPhys>& phys = SUDODEM_PTR_CAST<NormPhys>(interaction->phys);
	Real Kn;
	const GenericSpheresContact* geom=dynamic_cast<GenericSpheresContact*>(interaction->geom.get());
	if (geom) {
		Real Ra,Rb;//Vector3r normal;
		Ra=geom->refR1>0?geom->refR1:geom->refR2;
		Rb=geom->refR2>0?geom->refR2:geom->refR1;
		//harmonic average of the two stiffnesses when (Ri.Ei/2) is the stiffness of a contact point on sphere "i"
		Kn = 2*Ea*Ra*Eb*Rb/(Ea*Ra+Eb*Rb);
	} else {
		Kn = 2*Ea*Eb/(Ea+Eb);
	}
	phys->kn = Kn;
};


void Ip2_ElastMat_ElastMat_NormShearPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<ElastMat>& mat1 = SUDODEM_PTR_CAST<ElastMat>(b1);
	const shared_ptr<ElastMat>& mat2 = SUDODEM_PTR_CAST<ElastMat>(b2);
	Real Ea 	= mat1->young;
	Real Eb 	= mat2->young;
	Real Va 	= mat1->poisson;
	Real Vb 	= mat2->poisson;
	interaction->phys = shared_ptr<NormShearPhys>(new NormShearPhys());
	const shared_ptr<NormShearPhys>& phys = SUDODEM_PTR_CAST<NormShearPhys>(interaction->phys);
	Real Kn=0.0, Ks=0.0;
	GenericSpheresContact* geom=dynamic_cast<GenericSpheresContact*>(interaction->geom.get());
	if (geom) {
		Real Ra,Rb;//Vector3r normal;
		Ra=geom->refR1>0?geom->refR1:geom->refR2;
		Rb=geom->refR2>0?geom->refR2:geom->refR1;
		//harmonic average of the two stiffnesses when (Ri.Ei/2) is the stiffness of a contact point on sphere "i"
		Kn = 2*Ea*Ra*Eb*Rb/(Ea*Ra+Eb*Rb);
		Ks = 2*Ea*Ra*Va*Eb*Rb*Vb/(Ea*Ra*Va+Eb*Rb*Vb);
	} else {
		Kn = 2*Ea*Eb/(Ea+Eb);
		Kn = 2*Ea*Va*Eb*Vb/(Ea*Va+Eb*Vb);
	}
	phys->kn = Kn;
	phys->ks = Ks;
};
