#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/common/Aabb.hpp>
#include"PolySuperellipsoid.hpp"

#include<time.h>

#define _USE_MATH_DEFINES
#define TIGHT_BOX

SUDODEM_PLUGIN(/* self-contained in hpp: */ (PolySuperellipsoid) (PolySuperellipsoidGeom)
							(Bo1_PolySuperellipsoid_Aabb) (PolySuperellipsoidPhys)
							(PolySuperellipsoidMat) (Ip2_PolySuperellipsoidMat_PolySuperellipsoidMat_PolySuperellipsoidPhys)
							(PolySuperellipsoidLaw)
	/* some code in cpp (this file): */
	#ifdef SUDODEM_OPENGL
		(Gl1_PolySuperellipsoid) (Gl1_PolySuperellipsoidGeom) /*(Gl1_PolySuperellipsoidPhys)*/
	#endif
	);

//***************************************************************************************
/*Constractor*/

PolySuperellipsoid::PolySuperellipsoid(Vector2r eps_in, Vector6r halflen)
{
	createIndex();

	eps1=eps_in[0];
	eps2=eps_in[1];
  rxyz = halflen;
	rxyz_ref = halflen;
  eps = eps_in;
//cout<<"constructor"<<endl;
  init = false;//this is very important
	//initialize to calculate some basic geometric parameters.
	Initial();
	//get surface
	//Surface = getSurface();
}


void PolySuperellipsoid::Initial()
{       //cout<<"watch init"<<init<<endl;
        if (init) return;
        //
        //rx = rxyz(0);
        //ry = rxyz(1);
        //rz = rxyz(2);
        eps1 = eps(0);
        eps2 = eps(1);//
        //
        //cout<<"rxyz"<<rxyz(0)<<" "<<rxyz(1)<<" "<<rxyz(2)<<endl;
        //cout<<"eps"<<eps(0)<<" "<<eps(1)<<endl;

	//The volume and moments of inetia can be expressed in terms of Beta fuction.
	//The detailed derivations are shown in the reference, i.e., A. Jaklič, A. Leonardis,
	//F. Solina, Segmentation and Recovery of PolySuperellipsoid, Springer Netherlands, Dordrecht, 2000.
	//compute volume. moment of inertia here

		//get the mass center
		std::vector<Vector3r> centers,inertias;
		//std::vector<double> vols;
		massCenter = Vector3r(0,0,0);
		Inertia = Vector3r::Zero();
		Volume = 0.0;
		for(int k=0;k<2;k++){
			for(int j=0;j<2;j++){
				for(int i=0;i<2;i++){
					Vector3r center;
					double beta1,beta2,rx,ry,rz;
					//set rx,ry,rz
					rx = rxyz[i];
					ry = rxyz[j+2];
					rz = rxyz[k+4];
					beta1 = Mathr::beta(0.5*eps1,eps1)*Mathr::beta(0.5*eps2,1.5*eps2);
					beta2 = Mathr::beta(0.5*eps1,0.5*eps1)*Mathr::beta(0.5*eps2,eps2);
					center[0] = (1-2*i)*rx*3.0/4.0*beta1/beta2;
					center[1] = (1-2*j)*ry*3.0/4.0*beta1/beta2;
					center[2] = (1-2*k)*rz*3.0/4.0*Mathr::beta(eps2,eps2)/Mathr::beta(eps2,0.5*eps2);
					centers.push_back(center);
					double a,v;//some coefficients
					a = 0.125*rx*ry*rz*eps1*eps2;

					beta1 = Mathr::beta(1.5*eps1, 0.5*eps1)*Mathr::beta(0.5*eps2,2.*eps2+1);
					beta2 = Mathr::beta(0.5*eps1, 0.5*eps1+1)*Mathr::beta(1.5*eps2, eps2+1);
					v = 2.0*a*Mathr::beta(eps1*0.5,eps1*0.5)*Mathr::beta(eps2*0.5+1,eps2);     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
					//vols.push_back(v);
					Volume += v;
					massCenter += center*v;
					Inertia[0] += 0.5*a*(pow(ry, 2.)*beta1 + 4.*pow(rz, 2.)*beta2);	//I_xx
					Inertia[1] += 0.5*a*(pow(rx, 2.)*beta1 + 4.*pow(rz, 2.)*beta2);	//I_yy
					Inertia[2] += 0.5*a*(pow(rx, 2.) + pow(ry, 2.))*beta1;			//I_zz

				}
			}
		}
		massCenter /= Volume;
		massCenter_ref = massCenter;
		//applying the parallel theorem:I_xx = I_{xx}^c + M*l_x^2 = I_{xx}^c + M*(y^2+z^2)
		//The formular is wrong in the reference(John F. Peters, Mark A. Hopkins, Raju Kala, Ronald E. Wahl, (2009) "A poly‐ellipsoid particle for non‐spherical discrete element method", Engineering Computations, Vol. 26 Issue: 6, pp.645-657)
	  //also in Boning Zhang, Richard Regueiro, Andrew Druckrey, Khalid Alshibli, (2018) "Construction of polyellipsoidal grain shapes from SMT imaging on sand, and the development of a new DEM contact detection algorithm", Engineering Computations, Vol. 35 Issue: 2, pp.733-771
		//we need the moment of inertia with respect to the mass center from that with respect to the geometric center.
		//so the the right one is below
		//cout<<"inertia1="<<Inertia<<endl;
		Inertia[0] = Inertia[0] - (pow(massCenter[1],2)+pow(massCenter[2],2))*Volume;
		Inertia[1] = Inertia[1] - (pow(massCenter[0],2)+pow(massCenter[2],2))*Volume;
		Inertia[2] = Inertia[2] - (pow(massCenter[0],2)+pow(massCenter[1],2))*Volume;
		//cout<<"inertia2="<<Inertia<<endl;
	double rx,ry,rz,rx1,ry1,rz1,rx2,ry2,rz2;
	rx1 = rxyz[0]-massCenter[0];
	rx2 = rxyz[1]+massCenter[0];
	ry1 = rxyz[2]-massCenter[1];
	ry2 = rxyz[3]+massCenter[1];
	rz1 = rxyz[4]-massCenter[2];
	rz2 = rxyz[5]+massCenter[2];
	rx = max(rx1,rx2);
	ry = max(ry1,ry2);
	rz = max(rz1,rz2);
	r_max = (rx > ry) ? rx : ry;
	r_max = r_max > rz ? r_max : rz;
	if ((eps1 < 1.0) || (eps2 < 1.0)){
					r_max *= pow(2.0,0.5);
	}
	r_max_ref = r_max;
	Orientation = Quaternionr::Identity();
	//rotation matrices
  rot_mat2local = (Orientation).conjugate().toRotationMatrix();//to particle's system
	rot_mat2global = (Orientation).toRotationMatrix();//to global system

	//
	Position = Vector3r::Zero();
	//test
	//Quaternionr Rot(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	//Rot.normalize();
	//Orientation =Rot;
	//test_end
	ContactNorm = Vector3r(0.,0.,0.);
	CurrentP = Vector3r(0.,0.,0.);

        //initialization done
	init = true;
}


//****************************************************************************************
/* Destructor */

PolySuperellipsoid::~PolySuperellipsoid(){}
PolySuperellipsoidGeom::~PolySuperellipsoidGeom(){}

//****************************************************************************************
bool PolySuperellipsoid::isInside(Vector3r p)//check point p is inside the surface using the inside-outside function
{
    //double alpha1,alpha2;
    p = rot_mat2local*p + massCenter;
    return 1 > pow( pow(fabs(p(0)/rxyz[(p[0]>0?0:1)]),2.0/eps1) + pow(fabs(p(1)/rxyz[(p[1]>0?2:3)]),2.0/eps1), eps1/eps2) + pow(fabs(p(2)/rxyz[(p[2]>0?4:5)]),2.0/eps2);
}
//****************************************************************************************
Vector3r PolySuperellipsoid::getSurface(Vector2r phi) const//in local coordinates system
{
		assert(phi[0]>=0);
		double x,y,z;
		x = Mathr::Sign(cos(phi(0)))*pow(fabs(cos(phi(0))), eps[0])*pow(fabs(cos(phi(1))), eps[1]);
		y = Mathr::Sign(sin(phi(0)))*pow(fabs(sin(phi(0))), eps[0])*pow(fabs(cos(phi(1))), eps[1]);
		z = Mathr::Sign(sin(phi(1)))*pow(fabs(sin(phi(1))), eps[1]);
		x *= rxyz[(x>0?0:1)];
		y *= rxyz[(y>0?2:3)];
		z *= rxyz[(z>0?4:5)];
		return Vector3r(x,y,z);
	//std::cout<<"phi="<<phi<<"surf"<<Surf<<std::endl;
	//Matrix3r A = Orientation.toRotationMatrix();
	//Surf = A*Surf + Position;
}
Vector3r PolySuperellipsoid::getSurfaceMC(Vector2r phi) const//in local coordinates system with the origin at the mass center
{
		return getSurface(phi) - massCenter;
}

Vector3r PolySuperellipsoid::getNormal(Vector2r phi)
{
	Vector3r n;
	n(0) = Mathr::Sign(cos(phi(0))) *pow(fabs(cos(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(1) = Mathr::Sign(sin(phi(0))) *pow(fabs(sin(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(2) = Mathr::Sign(sin(phi(1))) *pow(fabs(sin(phi(1))), 2.-eps2);
	n(0) /= rxyz[(n(0)>0?0:1)];
	n(1) /= rxyz[(n(1)>0?2:3)];
	n(2) /= rxyz[(n(2)>0?4:5)];
	return n;
}

//*************************************************************************************
Vector2r PolySuperellipsoid::Normal2Phi(Vector3r n)
{       //n.normalize();
    //cout<<"contact normal: "<<n<<endl;
	Vector2r phi;
	/*
	phi(0) = atan(Mathr::Sign(n(1))/Mathr::Sign(n(0))*pow((fabs(ry*n(1)/rx/n(0))),1/(2.-eps1)));
	if (fabs(rx*n(0))>fabs(ry*n(1)))
	{
		phi(1) = atan(Mathr::Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(cos(phi(0))),2.-eps1) / fabs(rx*n(0))), 1/(2.-eps2) ) );
	}
	else{
		phi(1) = atan(Mathr::Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(sin(phi(0))),2.-eps1) / fabs(ry*n(1))), 1/(2.-eps2) ) );
	}
	*/
	//use atan2(y,x), which returns arctan(y/x).
	double rx,ry,rz;
	rx = rxyz[(n(0)>0?0:1)];
	ry = rxyz[(n(1)>0?2:3)];
	rz = rxyz[(n(2)>0?4:5)];
	phi(0) = atan2(Mathr::Sign(n(1))*pow((fabs(ry*n(1))),1/(2.-eps1)), Mathr::Sign(n(0))*pow((fabs(rx*n(0))),1/(2.-eps1))); //atan2 returns in [-pi,pi]//atan(x) returns in [-pi/2, pi/2]

    if (1)//(fabs(rx*n(0))>fabs(ry*n(1)))
	{
		phi(1) = atan2(Mathr::Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(cos(phi(0))),2.-eps1)), 1/(2.-eps2) ), pow(fabs(rx*n(0)), 1/(2.-eps2) ) );
		//std::cout<<"phi1==="<<phi(1)<<std::endl;
	}
	else{
		phi(1) = atan2(Mathr::Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(sin(phi(0))),2.-eps1)), 1/(2.-eps2) ), pow(fabs(ry*n(1)), 1/(2.-eps2) ) );
		//std::cout<<"phi1=111=="<<phi(1)<<std::endl;
	}

	/*if (1)//(fabs(rx*n(0))>fabs(ry*n(1)))
	{
		phi(1) = (Mathr::Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(cos(phi(0))),2.-eps1)), 1/(2.-eps2) ))/ pow(fabs(rx*n(0)), 1/(2.-eps2) ) ;
        phi(1) = atan(phi(1));
		//std::cout<<"phi1==="<<phi(1)<<std::endl;
	}
	else{
		phi(1) = atan2(Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(sin(phi(0))),2.-eps1)), 1/(2.-eps2) ), pow(fabs(ry*n(1)), 1/(2.-eps2) ) );
		//std::cout<<"phi1=111=="<<phi(1)<<std::endl;
	}*/
	return phi;
}

Vector3r PolySuperellipsoid::Normal2SurfaceMC(Vector3r n) const
{       //n.normalize();
    //cout<<"contact normal: "<<n<<endl;
	//Vector2r phi;
	//use atan2(y,x), which returns arctan(y/x).
	double rx,ry,rz;
	rx = rxyz[(n(0)>0?0:1)];
	ry = rxyz[(n(1)>0?2:3)];
	rz = rxyz[(n(2)>0?4:5)];
	double x,y,z,cos_phi0,sin_phi0,cos_phi1,sin_phi1;
	//phi(0) = atan2(Mathr::Sign(n(1))*pow((fabs(ry*n(1))),1/(2.-eps1)), Mathr::Sign(n(0))*pow((fabs(rx*n(0))),1/(2.-eps1))); //atan2 returns in [-pi,pi]//atan(x) returns in [-pi/2, pi/2]
	if(n(0)==0){
		cos_phi0 = 0.0;
		sin_phi0 = 1.0;
		if(n(1)==0){
			cos_phi1 = 0.0;
			sin_phi1 = 1.0;
		}else{
			double phi1_cos2 = 1.0/(1.0 + pow(fabs(rz*n(2)* pow(sin_phi0,2.-eps1)/ry/n(1)),2/(2.-eps2)));
			cos_phi1 = pow(phi1_cos2,0.5);
			sin_phi1 = pow(1.0 - phi1_cos2, 0.5);
		}
	}else{
		double phi0_cos2 = 1.0/(1.0 + pow((fabs(ry*n(1)/rx/n(0))),2/(2.-eps1)));
		cos_phi0 = pow(phi0_cos2,0.5);
		sin_phi0 = pow(1.0 - phi0_cos2, 0.5);
		double phi1_cos2 = 1.0/(1.0 + pow(fabs(rz*n(2)* pow(cos_phi0,2.-eps1)/rx/n(0)),2/(2.-eps2)));
		cos_phi1 = pow(phi1_cos2,0.5);
		sin_phi1 = pow(1.0 - phi1_cos2, 0.5);

	}

	x = Mathr::Sign(n(0))*rx*pow(cos_phi0, eps[0])*pow(cos_phi1, eps[1]);
	y = Mathr::Sign(n(1))*ry*pow(sin_phi0, eps[0])*pow(cos_phi1, eps[1]);
	z = Mathr::Sign(n(2))*rz*pow(sin_phi1, eps[1]);

	return Vector3r(x,y,z) - massCenter;
}
Vector3r PolySuperellipsoid::Normal2SurfaceMCgl(Vector3r gn) const
{       //n.normalize();
  Vector3r n = rot_mat2local*gn;//normal at global sys to local sys
	double rx,ry,rz;
	rx = rxyz[(n(0)>0?0:1)];
	ry = rxyz[(n(1)>0?2:3)];
	rz = rxyz[(n(2)>0?4:5)];
	double x,y,z,cos_phi0,sin_phi0,cos_phi1,sin_phi1;
	//phi(0) = atan2(Mathr::Sign(n(1))*pow((fabs(ry*n(1))),1/(2.-eps1)), Mathr::Sign(n(0))*pow((fabs(rx*n(0))),1/(2.-eps1))); //atan2 returns in [-pi,pi]//atan(x) returns in [-pi/2, pi/2]
	if(n(0)==0){
		cos_phi0 = 0.0;
		sin_phi0 = 1.0;
		if(n(1)==0){
			cos_phi1 = 0.0;
			sin_phi1 = 1.0;
		}else{
			double phi1_cos2 = 1.0/(1.0 + pow(fabs(rz*n(2)* pow(sin_phi0,2.-eps1)/ry/n(1)),2/(2.-eps2)));
			cos_phi1 = pow(phi1_cos2,0.5);
			sin_phi1 = pow(1.0 - phi1_cos2, 0.5);
		}
	}else{
		double phi0_cos2 = 1.0/(1.0 + pow((fabs(ry*n(1)/rx/n(0))),2/(2.-eps1)));
		cos_phi0 = pow(phi0_cos2,0.5);
		sin_phi0 = pow(1.0 - phi0_cos2, 0.5);
		double phi1_cos2 = 1.0/(1.0 + pow(fabs(rz*n(2)* pow(cos_phi0,2.-eps1)/rx/n(0)),2/(2.-eps2)));
		cos_phi1 = pow(phi1_cos2,0.5);
		sin_phi1 = pow(1.0 - phi1_cos2, 0.5);

	}

	x = Mathr::Sign(n(0))*rx*pow(cos_phi0, eps[0])*pow(cos_phi1, eps[1]);
	y = Mathr::Sign(n(1))*ry*pow(sin_phi0, eps[0])*pow(cos_phi1, eps[1]);
	z = Mathr::Sign(n(2))*rz*pow(sin_phi1, eps[1]);
	return rot_mat2global*(Vector3r(x,y,z) - massCenter);
}

Vector3r PolySuperellipsoid::support(Vector3r gn) const
{       //n.normalize();
  Vector3r n = rot_mat2local*gn;//normal at global sys to local sys
	double rx,ry,rz;
	rx = rxyz[(n(0)>0?0:1)];
	ry = rxyz[(n(1)>0?2:3)];
	rz = rxyz[(n(2)>0?4:5)];
	double x,y,z,cos_phi0,sin_phi0,cos_phi1,sin_phi1;
	//phi(0) = atan2(Mathr::Sign(n(1))*pow((fabs(ry*n(1))),1/(2.-eps1)), Mathr::Sign(n(0))*pow((fabs(rx*n(0))),1/(2.-eps1))); //atan2 returns in [-pi,pi]//atan(x) returns in [-pi/2, pi/2]
	if(n(0)==0){
		cos_phi0 = 0.0;
		sin_phi0 = 1.0;
		if(n(1)==0){
			cos_phi1 = 0.0;
			sin_phi1 = 1.0;
		}else{
			double phi1_cos2 = 1.0/(1.0 + pow(fabs(rz*n(2)* pow(sin_phi0,2.-eps1)/ry/n(1)),2/(2.-eps2)));
			cos_phi1 = pow(phi1_cos2,0.5);
			sin_phi1 = pow(1.0 - phi1_cos2, 0.5);
		}
	}else{
		double phi0_cos2 = 1.0/(1.0 + pow((fabs(ry*n(1)/rx/n(0))),2/(2.-eps1)));
		cos_phi0 = pow(phi0_cos2,0.5);
		sin_phi0 = pow(1.0 - phi0_cos2, 0.5);
		double phi1_cos2 = 1.0/(1.0 + pow(fabs(rz*n(2)* pow(cos_phi0,2.-eps1)/rx/n(0)),2/(2.-eps2)));
		cos_phi1 = pow(phi1_cos2,0.5);
		sin_phi1 = pow(1.0 - phi1_cos2, 0.5);

	}

	x = Mathr::Sign(n(0))*rx*pow(cos_phi0, eps[0])*pow(cos_phi1, eps[1]);
	y = Mathr::Sign(n(1))*ry*pow(sin_phi0, eps[0])*pow(cos_phi1, eps[1]);
	z = Mathr::Sign(n(2))*rz*pow(sin_phi1, eps[1]);
	return rot_mat2global*(Vector3r(x,y,z) - massCenter);
}
//****************************************************************************************
/* AaBb overlap checker  */

void Bo1_PolySuperellipsoid_Aabb::go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*){

	PolySuperellipsoid* t=static_cast<PolySuperellipsoid*>(ig.get());
	//t->rot_mat2local = (se3.orientation).conjugate().toRotationMatrix();//to particle's system
	//t->rot_mat2global = (se3.orientation).toRotationMatrix(); //to global system


	//find the maximal axis
	//double r_max = t->getrx() > t->getry() ? t->getrx() : t->getry();
	//r_max = r_max > t->getrz() ? r_max : t->getrz();
        //cout<<"pos="<<t_pos(0)<<" "<<t_pos(1)<<" "<<t_pos(2)<<endl;
        //cout<<"px="<<p_x(0)<<" py="<<p_y(1)<<" pz="<<p_z(2)<<endl;
	//cout<<"AABB"<<endl;
	if (!t->IsInitialized()){
	        t->Initial();
	        t->rot_mat2local = (se3.orientation).conjugate().toRotationMatrix();//to particle's system
	        t->rot_mat2global = (se3.orientation).toRotationMatrix(); //to global system
	}
	double r_max = t->getr_max();
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
        //Vector3r mincoords(-p_x(0),-p_y(1),-p_z(2)), maxcoords(p_x(0),p_y(1),p_z(2));
	#ifdef TIGHT_BOX
	//use tight box
	Matrix3r rot_mat1 = t->rot_mat2local;//(se3.orientation).conjugate().toRotationMatrix();//to particle's system
	Matrix3r rot_mat2 = t->rot_mat2global;//(se3.orientation).toRotationMatrix();//to global system
	//find the furthest point in a direction
	Vector3r n_x_min = Vector3r(-1,0,0);
	Vector3r n_x_max = Vector3r(1,0,0);
	Vector3r n_y_min = Vector3r(0,-1,0);
	Vector3r n_y_max = Vector3r(0,1,0);
	Vector3r n_z_min = Vector3r(0,0,-1);
	Vector3r n_z_max = Vector3r(0,0,1);
	//x axis
	Vector2r phi = t->Normal2Phi(rot_mat1*n_x_min);//phi in particle's local system
	Vector3r p_x_min = rot_mat2*( t->getSurfaceMC(phi));	//
	phi = t->Normal2Phi(rot_mat1*n_x_max);//phi in particle's local system
	Vector3r p_x_max = rot_mat2*( t->getSurfaceMC(phi));	//
	//y axis
	phi = t->Normal2Phi(rot_mat1*n_y_min);//phi in particle's local system
	Vector3r p_y_min = rot_mat2*( t->getSurfaceMC(phi));	//
	phi = t->Normal2Phi(rot_mat1*n_y_max);//phi in particle's local system
	Vector3r p_y_max = rot_mat2*( t->getSurfaceMC(phi));	//
	//z axis
	phi = t->Normal2Phi(rot_mat1*n_z_min);//phi in particle's local system
	Vector3r p_z_min = rot_mat2*( t->getSurfaceMC(phi));	//
	phi = t->Normal2Phi(rot_mat1*n_z_max);//phi in particle's local system
	Vector3r p_z_max = rot_mat2*( t->getSurfaceMC(phi));	//
  Vector3r mincoords(p_x_min[0],p_y_min[1],p_z_min[2]), maxcoords(p_x_max[0],p_y_max[1],p_z_max[2]);
	//using tight box ends
	#else
	Vector3r mincoords(-r_max,-r_max,-r_max), maxcoords(r_max,r_max,r_max);
	#endif
	//the aabb should be optimized in future.

	aabb->min=se3.position+mincoords;
	aabb->max=se3.position+maxcoords;
}

//**********************************************************************************
/* Plotting */

#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	bool Gl1_PolySuperellipsoid::wire;
	int Gl1_PolySuperellipsoid::slices;
	int Gl1_PolySuperellipsoid::pre_slices = -1;
    int Gl1_PolySuperellipsoid::glSolidList=-1;
    int Gl1_PolySuperellipsoid::glWireList=-1;
	//vector<Vector3r> Gl1_PolySuperellipsoid::vertices;

void Gl1_PolySuperellipsoid::initSolidGlList(PolySuperellipsoid* t)
{
	//Generate the list. Only once for each qtView, or more if quality is modified.
	glDeleteLists(t->m_glSolidListNum,1);
	t->m_glSolidListNum = glGenLists(1);
	glNewList(t->m_glSolidListNum,GL_COMPILE);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);//glShadeModel(GL_FLAT);

	int i,j,w=2*slices,h=slices;
	double a=0.0,b=0.0;
	double hStep=M_PI/(h-1);
	double wStep=2*M_PI/w;
	//Vector3r p = t->getPosition();
	//std::cerr << "Gl1_PolySuperellipsoid:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";

	vector<Vector3r> vertices;
    Vector3r vert;
	vertices.clear();
	for(a=0.0,i=0;i<h;i++,a+=hStep)
	  for(b=0.0,j=0;j<w;j++,b+=wStep)
	  {
		vertices.push_back(t->getSurfaceMC(Vector2r(b,a-M_PI_2l)));
	  }

   /* glBegin(GL_QUADS);
		for(i=0;i<h-1;i++)
		{
		    for(j=0;j<w-1;j++)
		    {
			//drawSlice(vertices[i*w+j],vertices[i*w+j+1],vertices[(i+1)*w+j+1],vertices[(i+1)*w+j]);
            vert = vertices[i*w+j];
            glVertex3v(vertices[i*w+j]);
            vert = vertices[i*w+j+1];
            glVertex3v(vertices[i*w+j+1]);
            vert = vertices[(i+1)*w+j+1];
            glVertex3v(vertices[(i+1)*w+j+1]);
            vert = vertices[(i+1)*w+j];
            glVertex3v(vertices[(i+1)*w+j]);
		    }
		     //drawSlice(vertices[i*w+j],vertices[i*w],vertices[(i+1)*w],vertices[(i+1)*w+j]);
            vert = vertices[i*w+j];
            glVertex3v(vertices[i*w+j]);
            vert = vertices[i*w];
            glVertex3v(vertices[i*w]);
            vert = vertices[(i+1)*w];
            glVertex3v(vertices[(i+1)*w]);
            vert = vertices[(i+1)*w+j];
            glVertex3v(vertices[(i+1)*w+j]);
		}

    glEnd();*/
    //////////////////////

    //FIXME:the two polar points should be just single points.
    //polygons of surface slices

    glBegin(GL_TRIANGLE_STRIP);
    for(i=0;i<h-1;i++)
    {   //Use TRIANGLE_STRIP for faster display of adjacent facets
        //glBegin(GL_TRIANGLE_STRIP);
        for(j=0;j<w-1;j++)
        {
            vert = vertices[(i+1)*w+j];vert.normalize();
            glNormal3v(vert);glVertex3v(vertices[(i+1)*w+j]);
            vert = vertices[i*w+j];vert.normalize();
            glNormal3v(vert);glVertex3v(vertices[i*w+j]);
        }
            vert = vertices[(i+1)*w];vert.normalize();
            glNormal3v(vert);glVertex3v(vertices[(i+1)*w]);
            vert = vertices[i*w];vert.normalize();
            glNormal3v(vert);glVertex3v(vertices[i*w]);

        //glEnd();
    }
    glEnd();
	glEndList();
}

void Gl1_PolySuperellipsoid::initWireGlList(PolySuperellipsoid* t)
{
	//Generate the "no-stripes" display list, each time quality is modified
	glDeleteLists(t->m_glWireListNum,1);
	t->m_glWireListNum = glGenLists(1);
	glNewList(t->m_glWireListNum,GL_COMPILE);
    glDisable(GL_LIGHTING);
	int i,j,w=2*slices,h=slices;
	double a=0.0,b=0.0;
	double hStep=M_PI/(h-1);
	double wStep=2*M_PI/w;
	//Vector3r p = t->getPosition();
	//std::cerr << "Gl1_PolySuperellipsoid:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";

	vector<Vector3r> vertices;
    Vector3r vert;
	vertices.clear();
	for(a=0.0,i=0;i<h;i++,a+=hStep)
	  for(b=0.0,j=0;j<w;j++,b+=wStep)
	  {
		vertices.push_back(t->getSurfaceMC(Vector2r(b,a-M_PI_2l)));
	  }
   /* glBegin(GL_LINE_LOOP);
		for(i=0;i<h-1;i++)
		{
		    for(j=0;j<w-1;j++)
		    {
			//drawSlice(vertices[i*w+j],vertices[i*w+j+1],vertices[(i+1)*w+j+1],vertices[(i+1)*w+j]);
            vert = vertices[i*w+j];
            glVertex3v(vertices[i*w+j]);
            vert = vertices[i*w+j+1];
            glVertex3v(vertices[i*w+j+1]);
            vert = vertices[(i+1)*w+j+1];
            glVertex3v(vertices[(i+1)*w+j+1]);
            vert = vertices[(i+1)*w+j];
            glVertex3v(vertices[(i+1)*w+j]);
		    }
		     //drawSlice(vertices[i*w+j],vertices[i*w],vertices[(i+1)*w],vertices[(i+1)*w+j]);
            vert = vertices[i*w+j];
            glVertex3v(vertices[i*w+j]);
            vert = vertices[i*w];
            glVertex3v(vertices[i*w]);
            vert = vertices[(i+1)*w];
            glVertex3v(vertices[(i+1)*w]);
            vert = vertices[(i+1)*w+j];
            glVertex3v(vertices[(i+1)*w+j]);
		}

    glEnd();*/

    glBegin(GL_LINE_STRIP);
    //draw the latitudes first
    for(i=1;i<h-1;i++)
    {
        //glBegin(GL_LINE_LOOP);
        for(j=0;j<w;j++)
        {	glVertex3v(vertices[i*w+j]);

        }
        glVertex3v(vertices[i*w]);
        //glEnd();
    }
    glVertex3v(vertices[(h-1)*w]);//the top pole

    //downwards the bottom pole at the opposite side
    for(i=h-2;i>=0;i--)
    {
        glVertex3v(vertices[i*w+slices]);////FIXME: it is assumed that w = 2*slices; the opposite column is the "slices".
    }

    //draw the rest longitudes
    for(j=1;j<slices;j++)
    {
        for(i=1;i<h;i++)
        {
            glVertex3v(vertices[i*w+j]);
        }
        for(i=h-1;i>=0;i--)
        {
            glVertex3v(vertices[i*w+j+slices]);
        }

    }
    glEnd();
	glEndList();
}


	void Gl1_PolySuperellipsoid::drawSlice(Vector3r &p1,Vector3r &p2,Vector3r &p3,Vector3r &p4)
	{       glDisable(GL_LIGHTING);
		if(!wire)
		{
		  //glEnable(GL_LIGHTING);
		  glBegin(GL_QUADS);
		 }
		else{
		  //glDisable(GL_LIGHTING);
		  glBegin(GL_LINE_LOOP);
		 }
		//  break;
		//}
		//glBegin(GL_LINE_LOOP);
		//glColor3f(0,1,0);
                //glColor3f(color(0),color(1),color(2));
		glVertex3v(p1);glVertex3v(p2);glVertex3v(p3);glVertex3v(p4);
		glEnd();
	}

	void Gl1_PolySuperellipsoid::go(const shared_ptr<Shape>& cm, const shared_ptr<State>&state,bool wire2,const GLViewInfo&)
	{
		glMaterialv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,Vector3r(cm->color[0],cm->color[1],cm->color[2]));
		glColor3v(cm->color);
		PolySuperellipsoid* t=static_cast<PolySuperellipsoid*>(cm.get());
		if (!t->IsInitialized()){
	                t->Initial();
	                t->rot_mat2local = state->ori.conjugate().toRotationMatrix();//to particle's system
	                t->rot_mat2global = state->ori.toRotationMatrix(); //to global system
	        }


        /*bool somethingChanged = (pre_slices!=slices || glIsList(glSolidList)!=GL_TRUE);
		if (somethingChanged) {
            //std::cout<<"haha preslices"<<pre_slices<<"slices"<<slices<<std::endl;
            initSolidGlList(t); initWireGlList(t);pre_slices=slices;}
        if(!wire)
        {
            glCallList(glSolidList);
        }else{
            glCallList(glWireList);
        }*/
        bool somethingChanged = (t->m_GL_slices!=slices|| glIsList(t->m_glSolidListNum)!=GL_TRUE);
		if (somethingChanged) {
            //std::cout<<"haha preslices"<<t->m_GL_slices<<"slices"<<slices<<std::endl;
            initSolidGlList(t); initWireGlList(t);t->m_GL_slices=slices;}
        if(wire||wire2)
        {
            glCallList(t->m_glWireListNum);//glCallList(t->m_glListNum);
        }else{
            glCallList(t->m_glSolidListNum);
        }




	}
	//
	void Gl1_PolySuperellipsoidGeom::go(const shared_ptr<IGeom>& ig, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool){
	        PolySuperellipsoidGeom* pg=static_cast<PolySuperellipsoidGeom*>(ig.get());
	        Vector3r point1,point2;

					point1 = (!scene->isPeriodic ? pg->point1 : scene->cell->wrapShearedPt(pg->point1));
	        point2 = (!scene->isPeriodic ? pg->point2 : scene->cell->wrapShearedPt(pg->point2));
	        //testing

	        Vector3r point11,point22;
	        point11 = (!scene->isPeriodic ? pg->point11 : scene->cell->wrapShearedPt(pg->point11));
	        point22 = (!scene->isPeriodic ? pg->point22 : scene->cell->wrapShearedPt(pg->point22));

					/*
	        point1 = pg->point1;
	        point2 = pg->point2;
	        //testing

	        Vector3r point11,point22;
	        point11 = pg->point11;
	        point22 = pg->point22;
					*/
	        //draw
	        glColor3f(0.5, 1., 1.0);
          glPointSize(10.0f);
          //the two closest points
          glBegin(GL_POINTS);
          glVertex3f(point1(0),point1(1),point1(2));
          glVertex3f(point2(0),point2(1),point2(2));
          glEnd();

          //Line connecting the two closest points
          glLineStipple (1, 0x0F0F);
          glBegin(GL_LINES);glLineWidth (3.0);

          glVertex3f(point1(0),point1(1),point1(2)); glVertex3f(point2(0),point2(1),point2(2));
          //
          glVertex3f(point1(0),point1(1),point1(2)); glVertex3f(point11(0),point11(1),point11(2));
          glVertex3f(point2(0),point2(1),point2(2)); glVertex3f(point22(0),point22(1),point22(2));

          glEnd();

          //


	}

#endif
//**********************************************************************************
//!Precompute data needed for rotating tangent vectors attached to the interaction

void PolySuperellipsoidGeom::precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r&
currentNormal, bool isNew, const Vector3r& shift2){

	if(!isNew) {
		orthonormal_axis = normal.cross(currentNormal);
		Real angle = scene->dt*0.5*currentNormal.dot(rbp1.angVel + rbp2.angVel);
		twist_axis = angle*currentNormal;}
        //if(isnan(orthonormal_axis[0])){std::cerr << "NAN000, normal "<<normal<<"currentN "<<currentNormal<< "\n";}
	else twist_axis=orthonormal_axis=Vector3r::Zero();
	//Update contact normal
	normal=currentNormal;
	//Precompute shear increment
	Vector3r c1x = (contactPoint - rbp1.pos);
	Vector3r c2x = (contactPoint - rbp2.pos - shift2);
	//std::cout<<"precompute---"<<c1x<<" "<<c2x<<std::endl;
	Vector3r relativeVelocity = (rbp2.vel+rbp2.angVel.cross(c2x)) - (rbp1.vel+rbp1.angVel.cross(c1x));
	//keep the shear part only
	relativeVn = normal.dot(relativeVelocity)*normal;
	relativeVs = relativeVelocity-relativeVn;
	shearInc = relativeVs*scene->dt;
}
//**********************************************************************************
Vector3r& PolySuperellipsoidGeom::rotate(Vector3r& shearForce) const {
	// approximated rotations
    //if(isnan(shearForce[0])){std::cerr << "NAN"<< "\n";}
	//std::cerr << "rto:"<<gettid4()<<orthonormal_axis[0]<<orthonormal_axis[1]<<orthonormal_axis[2]<< "\n";
	//std::cerr << "rtt:"<<gettid4()<<twist_axis[0]<<twist_axis[1]<<twist_axis[2]<< "\n";
	//std::cerr << "rtn:"<<gettid4()<<normal[0]<<normal[1]<<normal[2]<< "\n";
	//char filename[20];
	//sprintf(filename,"%d",gettid4());
	//FILE * fin = fopen(filename,"a");
    //if(shearForce.squaredNorm()<1E-18){return shearForce;}
	//fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
	//fprintf(fin,"shear_force1\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	//fprintf(fin,"orthonormal_axis\t%e\t%e\t%e\n",orthonormal_axis[0],orthonormal_axis[1],orthonormal_axis[2]);
	//fprintf(fin,"twist_axis\t%e\t%e\t%e\n",twist_axis[0],twist_axis[1],twist_axis[2]);
	//fprintf(fin,"normal\t%e\t%e\t%e\n",normal[0],normal[1],normal[2]);
	shearForce -= shearForce.cross(orthonormal_axis);
    //if(isnan(shearForce[0])){std::cerr << "NAN1:orthonormal_axis"<<orthonormal_axis<< "\n";}
	//fprintf(fin,"shear_force2\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	shearForce -= shearForce.cross(twist_axis);
    //if(isnan(shearForce[0])){std::cerr << "NAN2:"<<shearForce<< "\n";}
	//fprintf(fin,"shear_force3\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	//NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	shearForce -= normal.dot(shearForce)*normal;
    //if(isnan(shearForce[0])){std::cerr << "NAN3:"<<shearForce<< "\n";}
	//fprintf(fin,"shear_force4\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	//fclose(fin);
	return shearForce;
}


//**********************************************************************************
/* Material law, physics */

void Ip2_PolySuperellipsoidMat_PolySuperellipsoidMat_PolySuperellipsoidPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<PolySuperellipsoidMat>& mat1 = SUDODEM_PTR_CAST<PolySuperellipsoidMat>(b1);
	const shared_ptr<PolySuperellipsoidMat>& mat2 = SUDODEM_PTR_CAST<PolySuperellipsoidMat>(b2);
	interaction->phys = shared_ptr<PolySuperellipsoidPhys>(new PolySuperellipsoidPhys());
	const shared_ptr<PolySuperellipsoidPhys>& contactPhysics = SUDODEM_PTR_CAST<PolySuperellipsoidPhys>(interaction->phys);
    const Body::id_t id= interaction->id1;
	//Real Kna 	= mat1->Kn;
	//Real Knb 	= mat2->Kn;
	//Real Ksa 	= mat1->Ks;
	//Real Ksb 	= mat2->Ks;
    if (id>5){
        Real frictionAngle = std::min(mat1->frictionAngle,mat2->frictionAngle);
        contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);
        contactPhysics->kn = 2.0*mat1->Kn*mat2->Kn/(mat1->Kn+mat2->Kn);
	    contactPhysics->ks = 2.0*mat1->Ks*mat2->Ks/(mat1->Ks+mat2->Ks);

    }else{//walls, and the wall number is less than 6 by default.
        contactPhysics->tangensOfFrictionAngle = std::tan(mat1->frictionAngle);
        contactPhysics->kn = mat1->Kn;
	    contactPhysics->ks = mat1->Ks;
    }

	contactPhysics->betan = std::max(mat1->betan,mat2->betan);//
	contactPhysics->betas = std::max(mat1->betas,mat2->betas);


};

//**************************************************************************************
#if 1
Real PolySuperellipsoidLaw::getPlasticDissipation() {return (Real) plasticDissipation;}
void PolySuperellipsoidLaw::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}
Real PolySuperellipsoidLaw::elasticEnergy()
{
	Real energy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		FrictPhys* phys = dynamic_cast<FrictPhys*>(I->phys.get());
		if(phys) {
			energy += 0.5*(phys->normalForce.squaredNorm()/phys->kn + phys->shearForce.squaredNorm()/phys->ks);}
	}
	return energy;
}
#endif


//**************************************************************************************
// Apply forces on PolySuperellipsoid in collision based on geometric configuration
bool PolySuperellipsoidLaw::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){

		if (!I->geom) {return true;}
		const shared_ptr<PolySuperellipsoidGeom>& contactGeom(SUDODEM_PTR_DYN_CAST<PolySuperellipsoidGeom>(I->geom));
		if(!contactGeom) {return true;}
		const Body::id_t idA=I->getId1(), idB=I->getId2();
		const shared_ptr<Body>& A=Body::byId(idA), B=Body::byId(idB);

    State* de1 = Body::byId(idA,scene)->state.get();
	  State* de2 = Body::byId(idB,scene)->state.get();

		PolySuperellipsoidPhys* phys = dynamic_cast<PolySuperellipsoidPhys*>(I->phys.get());

		Matrix3r cellHsize; Vector3r ab_min,ab_max;
    ab_min = B->bound->min;ab_max = B->bound->max;
    if(scene->isPeriodic){
      cellHsize=scene->cell->hSize;
      Vector3r shift2=cellHsize*I->cellDist.cast<Real>();
      ab_min += shift2;
      ab_max += shift2;
    }
		//erase the interaction when aAbB shows separation, otherwise keep it to be able to store previous separating plane for fast detection of separation
		if (A->bound->min[0] >= ab_max[0] || ab_min[0] >= A->bound->max[0] || A->bound->min[1] >= ab_max[1] || ab_min[1] >= A->bound->max[1]  || A->bound->min[2] >= ab_max[2] || ab_min[2] >= A->bound->max[2])  {
		        phys->normalForce = Vector3r(0.,0.,0.); phys->shearForce = Vector3r(0.,0.,0.);
			scene->interactions->requestErase(I);
			return false;
		}
/*
		//erase the interaction when aAbB shows separation, otherwise keep it to be able to store previous separating plane for fast detection of separation
		if (A->bound->min[0] >= B->bound->max[0] || B->bound->min[0] >= A->bound->max[0] || A->bound->min[1] >= B->bound->max[1] || B->bound->min[1] >= A->bound->max[1] || A->bound->min[2] >= B->bound->max[2] || B->bound->min[2] >= A->bound->max[2])  {
		        phys->normalForce = Vector3r(0.,0.,0.); phys->shearForce = Vector3r(0.,0.,0.);
			scene->interactions->requestErase(I);
			return false;
		}
*/
		//zero penetration depth means no interaction force
		if(!(contactGeom->PenetrationDepth > 1E-18) ) {
            phys->normalForce = Vector3r(0.,0.,0.); phys->shearForce = Vector3r(0.,0.,0.);
            return true;
        }
		Vector3r normalForce=contactGeom->normal*contactGeom->PenetrationDepth*phys->kn;
		Vector3r shearForce1=Vector3r::Zero();//zhswee deprecated PolySuperellipsoidLaw.shearForce
		//shear force: in case the polyhdras are separated and come to contact again, one should not use the previous shear force
		//std::cerr << "p1:"<<gettid4()<<shearForce[0]<< "\n";
		///////////////////////
		bool useDamping=(phys->betan > 1e-5 || phys->betas > 1e-5);//using viscous damping?
		// tangential and normal stiffness coefficients, recomputed from betan,betas at every step
        Real cn=0, cs=0;
        Vector3r normalViscous = Vector3r::Zero();
        Vector3r shearViscous = Vector3r::Zero();
        if (useDamping){//

            Real mbar = (!A->isDynamic() && B->isDynamic()) ? de2->mass : ((!B->isDynamic() && A->isDynamic()) ? de1->mass : (de1->mass*de2->mass / (de1->mass + de2->mass))); // get equivalent mass if both bodies are dynamic, if not set it equal to the one of the dynamic body
	        //Real mbar = de1->mass*de2->mass / (de1->mass + de2->mass); // equivalent mass
	        Real Cn_crit = 2.*sqrt(mbar*phys->kn); // Critical damping coefficient (normal direction)
	        Real Cs_crit = 2.*sqrt(mbar*phys->ks); // Critical damping coefficient (shear direction)
	        // Note: to compare with the analytical solution you provide cn and cs directly (since here we used a different method to define c_crit)
	        cn = Cn_crit*phys->betan; // Damping normal coefficient
	        cs = Cs_crit*phys->betas; // Damping tangential coefficient
	        //get the relative velocity of two bodies at the contact
	        //it is more complecated when the particle is non-spherical.
	        //this has been done when calculating the shear increamental displacement (please see the function 'precompute')
	        normalViscous = cn*contactGeom->relativeVn;
	        //std::cerr << "normalViscous:"<<normalViscous(0)<<" "<<normalViscous(1)<<" "<<normalViscous(2)<< "\n";
	        Vector3r normTemp = normalForce - normalViscous; // temporary normal force
	        // viscous force should not exceed the value of current normal force, i.e. no attraction force should be permitted if particles are non-adhesive
	        //std::cerr << "nF1:"<<normalForce.norm()<< "\n";
	        if (normTemp.dot(contactGeom->normal)<0.0){

				normalViscous = normalForce; // normal viscous force is such that the total applied force is null - it is necessary to compute energy correctly!
				normalForce = Vector3r::Zero();
			}
			else{normalForce -= normalViscous;}
	                //std::cerr << "nF2:"<<normalForce.norm()<< "\n";
        }
        //if(isnan(phys->shearForce[0])){std::cerr << "NAN, phys->shearForce:"<<shearForce1<< "\n";}
		if (contactGeom->isShearNew){
			//shearForce = Vector3r::Zero();
			shearForce1 = Vector3r(0.,0.,0.);
		}else{
			//shearForce = contactGeom->rotate(shearForce);

			shearForce1 = contactGeom->rotate(phys->shearForce);
			//std::cerr << "p2:"<<gettid4() <<shearForce[0]<< "\n";
        }
		const Vector3r& shearDisp = contactGeom->shearInc;
		shearForce1 -= phys->ks*shearDisp;
        //if(isnan(shearForce1[0])){std::cerr << "NAN, shearDisp:"<<shearDisp<< "\n";}
                //std::cerr << "p3:"<<gettid4() <<shearDisp[0]<< "\n";
		Real maxFs = normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);
		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if( shearForce1.squaredNorm() > maxFs ){
			Real ratio = sqrt(maxFs) / shearForce1.norm();
			shearForce1 *= ratio;}
		else{
            if (useDamping){ // add current contact damping if we do not slide and if damping is requested
		            shearViscous = cs*contactGeom->relativeVs; // get shear viscous component
		            //std::cerr << "sF1:"<<shearForce1(0)<<" "<<shearForce1(1)<<" "<<shearForce1(2)<< "\n";
		            //std::cerr << "sF1:"<<shearForce1.norm()<< "\n";
		            //shearForce1 -= shearViscous;
		            //std::cerr << "sF2:"<<shearForce1.norm()<< "\n";
		    }
        }
        //if(isnan(shearForce1[0])){std::cerr << "NAN, shearViscous:"<<shearViscous<< "\n";}
		//std::cerr << "using viscous damping:"<<useDamping<< "\n";
		/*
		if (likely(!scene->trackEnergy  && !traceEnergy)){//Update force but don't compute energy terms (see below))
			// PFC3d SlipModel, is using friction angle. CoulombCriterion
			if( shearForce.squaredNorm() > maxFs ){
				Real ratio = sqrt(maxFs) / shearForce.norm();
				shearForce *= ratio;}
		} else {
			//almost the same with additional Vector3r instatinated for energy tracing,
			//duplicated block to make sure there is no cost for the instanciation of the vector when traceEnergy==false
			if(shearForce.squaredNorm() > maxFs){
				Real ratio = sqrt(maxFs) / shearForce.norm();
				Vector3r trialForce=shearForce;//store prev force for definition of plastic slip
				//define the plastic work input and increment the total plastic energy dissipated
				shearForce *= ratio;
				Real dissip=((1/phys->ks)*(trialForce-shearForce)).dot(shearForce);
				if (traceEnergy) plasticDissipation += dissip;
				else if(dissip>0) scene->energy->add(dissip,"plastDissip",plastDissipIx,false);
			}
			// compute elastic energy as well
			scene->energy->add(0.5*(normalForce.squaredNorm()/phys->kn+shearForce.squaredNorm()/phys->ks),"elastPotential",elastPotentialIx,true);
		}*/



		/*
		FILE * fin = fopen("Forces.dat","a");
		fprintf(fin,"************** IDS %d %d **************\n",idA, idB);
		Vector3r T = (B->state->pos-contactGeom->contactPoint).cross(F);
		fprintf(fin,"volume\t%e\n",contactGeom->penetrationVolume);
		fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
		fprintf(fin,"shear_force\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
		fprintf(fin,"total_force\t%e\t%e\t%e\n",F[0],F[1],F[2]);
		fprintf(fin,"torsion\t%e\t%e\t%e\n",T[0],T[1],T[2]);
		fprintf(fin,"A\t%e\t%e\t%e\n",A->state->pos[0],A->state->pos[1],A->state->pos[2]);
		fprintf(fin,"B\t%e\t%e\t%e\n",B->state->pos[0],B->state->pos[1],B->state->pos[2]);
		fprintf(fin,"centroid\t%e\t%e\t%e\n",contactGeom->contactPoint[0],contactGeom->contactPoint[1],contactGeom->contactPoint[2]);
		fclose(fin);
		*/
		if(isnan(shearForce1[0])||isnan(normalForce[0])){
		  //char filename[20];
          time_t tmNow = time(NULL);
          char tmp[64];
          strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&tmNow) );
          const char* filename = "PolySuperellipsoidLaw_err.log";
		  //sprintf(filename,"%d",gettid4());
		  FILE * fin = fopen(filename,"a");
          fprintf(fin,"\n*******%s*******\n",tmp);
          Vector3r normal = contactGeom->normal;
          Real depth = contactGeom->PenetrationDepth;
		  //fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
		  fprintf(fin,"shear_force\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
		  fprintf(fin,"shear_force1\t%e\t%e\t%e\n",shearForce1[0],shearForce1[1],shearForce1[2]);
          //fprintf(fin,"F\t%e\t%e\t%e\n",F[0],F[1],F[2]);
		  fprintf(fin,"ks\t%emaxFs\t%esfn\t%e\n",phys->ks,maxFs,shearForce.squaredNorm());
		  fprintf(fin,"shearDisp\t%e\t%e\t%e\n",shearDisp[0],shearDisp[1],shearDisp[2]);
		  fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
          fprintf(fin,"contact normal \t%e\t%e\t%e\n",normal[0],normal[1],normal[2]);
          fprintf(fin,"depth\t%e\n",depth);
		  fclose(fin);
		  //throw runtime_error("P4 #");
          LOG_WARN("Error in computation of contact forces");
          //scene->stopAtIter = scene->iter + 1;//Not fatal error for only one iteration, so just pass it but with a warning.
          normalForce = Vector3r::Zero();
          shearForce1 = Vector3r::Zero();
		}
		Vector3r F = -normalForce-shearForce1+shearViscous;	//the normal force acting on particle A (or 1) is equal to the normal force.
		if (contactGeom->PenetrationDepth != contactGeom->PenetrationDepth) exit(1);
		scene->forces.addForce (idA,F);
		scene->forces.addForce (idB, -F);
		scene->forces.addTorque(idA, -(A->state->pos-contactGeom->contactPoint).cross(F));
		scene->forces.addTorque(idB, (B->state->pos-contactGeom->contactPoint).cross(F));
		//needed to be able to acces interaction forces in other parts of sudodem
		shearForce=shearForce1;
		phys->normalForce = normalForce;
		phys->shearForce = shearForce1;
		//may need to add phys->shearViscous and phys->normalViscous
		return true;
}
