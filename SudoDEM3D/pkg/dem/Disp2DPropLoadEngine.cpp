/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<sudodem/lib/base/Math.hpp>
#include "Disp2DPropLoadEngine.hpp"
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/Box.hpp>
#include<sudodem/core/Scene.hpp>
//#include<sudodem/lib/base/Math.hpp>


SUDODEM_PLUGIN((Disp2DPropLoadEngine));

void Disp2DPropLoadEngine::postLoad(Disp2DPropLoadEngine&)
{
	std::string outputFile="DirSearch" + Key + "SudoDEM";
	bool file_exists = (bool)std::ifstream (outputFile.c_str()); //if file does not exist, we will write colums titles
	ofile.open(outputFile.c_str(), std::ios::app);
	if (!file_exists) ofile<<"theta (!angle in plane (gamma,-du) ) dtau (kPa) dsigma (kPa) dgamma (m) du (m) tau0 (kPa) sigma0 (kPa) d2W coordSs0 coordTot0 coordSsF coordTotF (SudoDEM)" << endl;
}




void Disp2DPropLoadEngine::action()
{
	if(LOG) cerr << "debut applyCondi !!" << endl;
	leftbox = Body::byId(id_boxleft);
	rightbox = Body::byId(id_boxright);
	frontbox = Body::byId(id_boxfront);
	backbox = Body::byId(id_boxback);
	topbox = Body::byId(id_topbox);
	boxbas = Body::byId(id_boxbas);

	if(firstIt)
	{
		it_begin=scene->iter;
		H0=topbox->state->pos.y();
		X0=topbox->state->pos.x();
		Vector3r F_sup=scene->forces.getForce(id_topbox);
		Fn0=F_sup.y();
		Ft0=F_sup.x();

		Real	OnlySsInt=0	// the half number of real sphere-sphere (only) interactions, at the beginning of the perturbation
			,TotInt=0	// the half number of all the real interactions, at the beginning of the perturbation
			;

		InteractionContainer::iterator ii    = scene->interactions->begin();
		InteractionContainer::iterator iiEnd = scene->interactions->end();
        	for(  ; ii!=iiEnd ; ++ii )
        	{
        		if ((*ii)->isReal())
                	{
				TotInt++;
				const shared_ptr<Body>& b1 = Body::byId( (*ii)->getId1() );
				const shared_ptr<Body>& b2 = Body::byId( (*ii)->getId2() );
				if ( (b1->isDynamic()) && (b2->isDynamic()) )
					OnlySsInt++;
			}
		}

		coordSs0 = OnlySsInt/8590;	// 8590 is the number of spheres in the CURRENT case
		coordTot0 = TotInt / 8596;	// 8596 is the number of bodies in the CURRENT case

		firstIt=false;
	}


	if ( (scene->iter-it_begin) < nbre_iter)
	{	letDisturb();
	}
	else if ( (scene->iter-it_begin) == nbre_iter)
	{
		stopMovement();
		string fileName=Key + "DR"+boost::lexical_cast<string> (nbre_iter)+"ItAtV_"+boost::lexical_cast<string> (v)+"done.xml";
// 		Omega::instance().saveSimulation ( fileName );
		saveData();
	}
}

void Disp2DPropLoadEngine::letDisturb()
{

	const Real& dt = scene->dt;
	dgamma=cos(theta*Mathr::PI/180.0)*v*dt;
	dh=sin(theta*Mathr::PI/180.0)*v*dt;

	Real Ysup = topbox->state->pos.y();
	Real Ylat = leftbox->state->pos.y();

// 	Changes in vertical and horizontal position :
	topbox->state->pos += Vector3r(dgamma,dh,0);

	leftbox->state->pos += Vector3r(dgamma/2.0,dh/2.0,0);
	rightbox->state->pos += Vector3r(dgamma/2.0,dh/2.0,0);

	Real Ysup_mod = topbox->state->pos.y();
	Real Ylat_mod = leftbox->state->pos.y();

//	with the corresponding velocities :
	topbox->state->vel=Vector3r(dgamma/dt,dh/dt,0);

	leftbox->state->vel = Vector3r((dgamma/dt)/2.0,dh/(2.0*dt),0);

	rightbox->state->vel = Vector3r((dgamma/dt)/2.0,dh/(2.0*dt),0);

//	Then computation of the angle of the rotation to be done :
	computeAlpha();
	if (alpha == Mathr::PI/2.0)	// Case of the very beginning
	{
		dalpha = - atan( dgamma / (Ysup_mod -Ylat_mod) );
	}
	else
	{
		Real A = (Ysup_mod - Ylat_mod) * 2.0*tan(alpha) / (2.0*(Ysup - Ylat) + dgamma*tan(alpha) );
		dalpha = atan( (A - tan(alpha))/(1.0 + A * tan(alpha)));
	}

	Quaternionr qcorr(AngleAxisr(dalpha,Vector3r::UnitZ()));
	if(LOG)
		cout << "Quaternion associe a la rotation incrementale : " << qcorr.w() << " " << qcorr.x() << " " << qcorr.y() << " " << qcorr.z() << endl;

// On applique la rotation en changeant l'orientation des plaques, leurs vang et en affectant donc alpha
	leftbox->state->ori = qcorr*leftbox->state->ori;
	leftbox->state->angVel = Vector3r(0,0,1)*dalpha/dt;

	rightbox->state->ori = qcorr*leftbox->state->ori;
	rightbox->state->angVel = Vector3r(0,0,1)*dalpha/dt;

}



void Disp2DPropLoadEngine::computeAlpha()
{
	Quaternionr orientationLeftBox,orientationRightBox;
	orientationLeftBox = leftbox->state->ori;
	orientationRightBox = rightbox->state->ori;
	if(orientationLeftBox!=orientationRightBox)
	{
		cout << "WARNING !!! your lateral boxes have not the same orientation, you're not in the case of a box imagined for creating these engines" << endl;
	}
	AngleAxisr aa(orientationLeftBox);
	alpha=Mathr::PI/2.0-aa.angle();		// right if the initial orientation of the body (on the beginning of the simulation) is q =(1,0,0,0) = FromAxisAngle((0,0,1),0)
}


void Disp2DPropLoadEngine::stopMovement()
{
	// annulation de la vitesse de la plaque du haut
	topbox->state->vel	=  Vector3r(0,0,0);

	// de la plaque gauche
	leftbox->state->vel	=  Vector3r(0,0,0);
	leftbox->state->angVel	=  Vector3r(0,0,0);

	// de la plaque droite
	rightbox->state->vel	=  Vector3r(0,0,0);
	rightbox->state->angVel	=  Vector3r(0,0,0);
}


void Disp2DPropLoadEngine::saveData()
{
	Real Xleft = leftbox->state->pos.x() + (SUDODEM_CAST<Box*>(leftbox->shape.get()))->extents.x();

	Real Xright = rightbox->state->pos.x() - (SUDODEM_CAST<Box*>(rightbox->shape.get()))->extents.x();

	Real Zfront = frontbox->state->pos.z() - SUDODEM_CAST<Box*>(frontbox->shape.get())->extents.z();
	Real Zback = backbox->state->pos.z() + (SUDODEM_CAST<Box*>(backbox->shape.get()))->extents.z();

	Real Scontact = (Xright-Xleft)*(Zfront-Zback);	// that's so the value of section at the middle of the height of the box

	InteractionContainer::iterator ii    = scene->interactions->begin();
        InteractionContainer::iterator iiEnd = scene->interactions->end();

	Real	OnlySsInt=0	// the half number of real sphere-sphere (only) interactions, at the end of the perturbation
		,TotInt=0	// the half number of all the real interactions, at the end of the perturbation
		;
        for(  ; ii!=iiEnd ; ++ii )
        {
        	if ((*ii)->isReal())
                {
			TotInt++;
			const shared_ptr<Body>& b1 = Body::byId( (*ii)->getId1() );
			const shared_ptr<Body>& b2 = Body::byId( (*ii)->getId2() );
			if ( (b1->isDynamic()) && (b2->isDynamic()) )
				OnlySsInt++;
		}
	}

	Real	coordSs = OnlySsInt/8590,	// 8590 is the number of spheres in the CURRENT case
		coordTot = TotInt / 8596;	// 8596 is the number of bodies in the CURRENT case

	Vector3r F_sup = scene->forces.getForce(id_topbox);

	Real	dFn=F_sup.y()-Fn0	// OK pour le signe
		,dFt=(F_sup.x()-Ft0)
		,du=-( topbox->state->pos.y() - H0 )	// OK pour le signe (>0 en compression)
		,dgamma=topbox->state->pos.x()-X0
		,SigN0 = (Fn0/Scontact)/1000	// en kPa, pour comparer à Fortran
		,Tau0 = -(Ft0/Scontact)/1000	// en kPa, pour comparer à Fortran, Ok pour le signe, cf p. SudoDEM29
		,dSigN= (dFn/Scontact)/1000	// Ok pour le signe
		,dTau = -(dFt/Scontact)/1000	// Ok pour le signe, idem que Tau0
		,d2W = dSigN * du + dTau * dgamma
		;

	ofile << boost::lexical_cast<string> (theta) << " "<< boost::lexical_cast<string> (dTau) << " " << boost::lexical_cast<string> (dSigN) << " "
		<< boost::lexical_cast<string> (dgamma)<<" " << boost::lexical_cast<string> (du) << " " << boost::lexical_cast<string> (Tau0) << " "
		<< boost::lexical_cast<string> (SigN0) << " " << boost::lexical_cast<string> (d2W) << " "
		<< boost::lexical_cast<string> (coordSs0) << " " << boost::lexical_cast<string> (coordTot0) << " "
		<< boost::lexical_cast<string> (coordSs) << " " << boost::lexical_cast<string> (coordTot) <<endl;
}
