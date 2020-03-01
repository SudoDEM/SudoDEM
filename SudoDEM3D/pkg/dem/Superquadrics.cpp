#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/common/Aabb.hpp>
#include"Superquadrics.hpp"

#include<time.h>

#define _USE_MATH_DEFINES


SUDODEM_PLUGIN(/* self-contained in hpp: */ (Superquadrics) (SuperquadricsGeom) (SuperquadricsGeom2) (Bo1_Superquadrics_Aabb) (SuperquadricsPhys) (SuperquadricsHertzMindlinPhys) (SuperquadricsMat) (SuperquadricsMat2) (Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys) (Ip2_SuperquadricsMat2_SuperquadricsMat2_HertzMindlinPhys) (SuperquadricsLaw) (SuperquadricsLaw2)
	/* some code in cpp (this file): */
	#ifdef SUDODEM_OPENGL
		(Gl1_Superquadrics) (Gl1_SuperquadricsGeom) /*(Gl1_SuperquadricsPhys)*/
	#endif
	);


////////////////////////////////some auxiliary functions/////////////////////////////////////////////////////////////////
double cot(double x){return tan(M_PI_2l-x);}


//***************************************************************************************
/*Constractor*/

Superquadrics::Superquadrics(double x, double y, double z, double ep1, double ep2)
{       createIndex();
	rx=x;
	ry=y;
	rz=z;
	eps1=ep1;
	eps2=ep2;
        rxyz = Vector3r(x,y,z);
        eps = Vector2r(ep1,ep2);
	//cout<<"constructor"<<endl;
        init = false;//this is very important
	//initialize to calculate some basic geometric parameters.
	Initial();
	//get surface
	//Surface = getSurface();
}
//sssssssssssssssssssssss
/*sssssssssssssss*/
void Superquadrics::Initial()
{       //cout<<"watch init"<<init<<endl;
        if (init) return;
        //
        rx = rxyz(0);
        ry = rxyz(1);
        rz = rxyz(2);
        eps1 = eps(0);
        eps2 = eps(1);//
        //
        //cout<<"rxyz"<<rxyz(0)<<" "<<rxyz(1)<<" "<<rxyz(2)<<endl;
        //cout<<"eps"<<eps(0)<<" "<<eps(1)<<endl;
        r_max = (rx > ry) ? rx : ry;
        r_max = r_max > rz ? r_max : rz;
        if ((eps1 < 1.0) || (eps2 < 1.0)){
                r_max *= pow(2.0,0.5);
        }
	//The volume and moments of inetia can be expressed in terms of Beta fuction.
	//The detailed derivations are shown in the reference, i.e., A. Jaklič, A. Leonardis,
	//F. Solina, Segmentation and Recovery of Superquadrics, Springer Netherlands, Dordrecht, 2000.
	double a,beta1,beta2;//some coefficients
	a = rx*ry*rz*eps1*eps2;
	beta1 = Mathr::beta(1.5*eps1, 0.5*eps1)*Mathr::beta(0.5*eps2,2.*eps2+1);
	beta2 = Mathr::beta(0.5*eps1, 0.5*eps1+1)*Mathr::beta(1.5*eps2, eps2+1);
	Volume = 2.*a*Mathr::beta(eps1*0.5,eps1*0.5)*Mathr::beta(eps2*0.5+1,eps2);     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
	Inertia(0) = 0.5*a*(pow(ry, 2.)*beta1 + 4.*pow(rz, 2.)*beta2);	//I_xx
	Inertia(1) = 0.5*a*(pow(rx, 2.)*beta1 + 4.*pow(rz, 2.)*beta2);	//I_yy
	Inertia(2) = 0.5*a*(pow(rx, 2.) + pow(ry, 2.))*beta1;			//I_zz

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

Superquadrics::~Superquadrics(){}
SuperquadricsGeom::~SuperquadricsGeom(){}
SuperquadricsGeom2::~SuperquadricsGeom2(){}
//****************************************************************************************
bool Superquadrics::isInside(Vector3r p)//check point p is inside the surface using the inside-outside function
{
    //double alpha1,alpha2;
    p = rot_mat2local*p;
    return 1 > pow( pow(fabs(p(0)/rx),2.0/eps1) + pow(fabs(p(1)/ry),2.0/eps1), eps1/eps2) + pow(fabs(p(2)/rz),2.0/eps2);
}
//****************************************************************************************
Vector3r Superquadrics::getSurface(Vector2r phi) const//in local coordinates system
{
	Vector3r Surf;
	Surf(0) = Mathr::Sign(cos(phi(0)))*rx*pow(fabs(cos(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	Surf(1) = Mathr::Sign(sin(phi(0)))*ry*pow(fabs(sin(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	Surf(2) = Mathr::Sign(sin(phi(1)))*rz*pow(fabs(sin(phi(1))), eps2);
	//std::cout<<"phi="<<phi<<"surf"<<Surf<<std::endl;
	//Matrix3r A = Orientation.toRotationMatrix();
	//Surf = A*Surf + Position;
	return Surf;
}

void Superquadrics::getCurrentP(Vector2r phi)  //get the current p on the surface
{
	CurrentP(0) = Mathr::Sign(cos(phi(0)))*rx*pow(fabs(cos(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	CurrentP(1) = Mathr::Sign(sin(phi(0)))*ry*pow(fabs(sin(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	CurrentP(2) = Mathr::Sign(sin(phi(1)))*rz*pow(fabs(sin(phi(1))), eps2);

}

//Get principal curvatures at a given point
Vector2r Superquadrics::Curvatures(Vector2r phi,Vector3r &CurvatureDir_max,Vector3r &CurvatureDir_min,bool &curvatureMaxflag)
{   curvatureMaxflag = true;//using CurvatureDir_max first by default.
    //Mat32r DerivateR;// = Derivate_P2Phi(phi);
	//std::cout<<"phi"<<phi(0)<<'\t'<<phi(1)<<std::endl;
	//std::cout<<"a"<<rx<<"b"<<ry<<"c"<<rz<<"eps1"<<eps1<<"eps2"<<eps2<<std::endl;
	//check phi1 and phi2 to avoid numerical errors "NAN", i.e., sin(phi(0)) = 0 and sin(phi(1))=0 would cast a NAN in the following computation.
	if (fabs(phi(0)) < 1E-20){phi(0)=1E-20;}
	if (fabs(phi(1)) < 1E-20){phi(1)=1E-20;}
    Vector3r R1,R2,R12,R11,R22;
    double  a,b,c,d;//,a1,b1,c1,d1;
	double s1,s2,c1,c2;
	s1 = sin(phi(0));s2 = sin(phi(1));
	c1 = cos(phi(0));c2 = cos(phi(1));
	double s_s1,s_s2,s_c1,s_c2;
	s_s1 = Mathr::Sign(s1); s_s2 = Mathr::Sign(s2);
	s_c1 = Mathr::Sign(c1); s_c2 = Mathr::Sign(c2);
	double abs_s1,abs_s2,abs_c1,abs_c2;
	abs_s1 = fabs(s1); abs_s2 = fabs(s2);
	abs_c1 = fabs(c1); abs_c2 = fabs(c2);
	double s10,s20,c10,c20;//f^{eps-0}
	s10 = s_s1*pow(abs_s1,eps1);
	s20 = s_s2*pow(abs_s2,eps2);
	c10 = s_c1*pow(abs_c1,eps1);
	c20 = s_c2*pow(abs_c2,eps2);
	double s11,s21,c11,c21;//f^{eps-1}//not used cause eps1=1 would prompt a bug
	s11 = s_s1*pow(abs_s1,eps1-1.0);
	s21 = s_s2*pow(abs_s2,eps2-1.0);
	c11 = s_c1*pow(abs_c1,eps1-1.0);
	c21 = s_c2*pow(abs_c2,eps2-1.0);
	double s12,s22,c12,c22;//f^(eps-2)
	s12 = s_s1*pow(abs_s1,eps1-2.0);
	s22 = s_s2*pow(abs_s2,eps2-2.0);
	c12 = s_c1*pow(abs_c1,eps1-2.0);
	c22 = s_c2*pow(abs_c2,eps2-2.0);

	R1 = Vector3r(-rx*eps1*c12*c1*c20*s1,ry*eps1*c1*c20*s12*s1,0);
	R2 = Vector3r(-rx*eps2*c10*c22*c2*s2, -ry*eps2*c22*c2*s10*s2, rz*eps2*c2*s22*s2);
	R12 = Vector3r(a*eps1*eps2*c12*c1*c22*c2*s1*s2,-ry*eps1*eps2*c1*c22*c2*s12*s1*s2,0);
	R11 = Vector3r(rx*eps1*c12*s20*s1*s1*(eps1-1.0)-rx*eps1*c10*c20,
				  ry*eps1*c1*c1*c20*s12*(eps1-1.0)-ry*eps1*c20*s10,
				  0);
	R22 = Vector3r(rx*eps2*c10*c22*s2*s2*(eps2-1.0)-rx*eps2*c10*c20,
				   ry*eps2*s10*s1*s1*(eps2-1.0)-ry*eps2*c20*s10,
				   rz*eps2*c2*c2*s22*(eps2-1.0)-rz*eps2*s20);


    //std::cout<<"Der_Rphi1"<<R1<<std::endl;
    //std::cout<<"Der_Rphi2"<<R2<<std::endl;
	//std::cout<<"R11"<<R11<<std::endl;
    //std::cout<<"R22"<<R22<<std::endl;
	//std::cout<<"R12"<<R12<<std::endl;

    Vector3r RR = R1.cross(R2);
	Vector3r n = RR/RR.norm();
	//std::cout<<"n"<<n<<std::endl;
	double L,N,E,F,G,K,H,k1,k2;
	L = R11.dot(n);
	N = R22.dot(n);
	E = R1.dot(R1);
	F = R1.dot(R2);
	G = R2.dot(R2);
    K = L*N/(E*G-F*F);
	H = 0.5*(E*N+G*L)/(E*G-F*F);//Caution: calculated curvatures are negative by defaut.
	double tmp0 = sqrt(fabs(H*H-K));//pow(H*H-K,0.5);//H2>=K; this condition may not be fullfilled for numerical errors when H^2 is very close to K;
	k1 = H+tmp0;
	k2 = H- tmp0;//should fix the bug, sqrt(a), a>=0
    /*if(isnan(k1)||isnan(k2)){//for debugging
		std::cout<<std::setprecision(12)<<"phi"<<phi(0)<<'\t'<<phi(1)<<std::endl;
		std::cout<<"a"<<rx<<"b"<<ry<<"c"<<rz<<"eps1"<<eps1<<"eps2"<<eps2<<std::endl;
		std::cout<<"Der_Rphi1"<<R1<<std::endl;
		std::cout<<"Der_Rphi2"<<R2<<std::endl;
		std::cout<<"R11"<<R11<<std::endl;
		std::cout<<"R22"<<R22<<std::endl;
		std::cout<<"R12"<<R12<<std::endl;
		std::cout<<"NNNN  L"<<L<<" N"<<N<<" E"<<E<<" F"<<F<<" G"<<G<<" K"<<K<<" H"<<H<<'\t'<<std::endl;
		std::cout<<"k1 "<<k1<<" k2"<<k2<<std::endl;
		std::cout<<"tmp0 "<<tmp0<<" H*H-K="<<H*H-K<<"tmp1"<<sqrt(H*H-K)<<std::endl;

	}*/
	//std::cout<<"NNNN  L"<<L<<" N"<<N<<" E"<<E<<" F"<<F<<" G"<<G<<" K"<<K<<" H"<<H<<'\t'<<std::endl;
	//std::cout<<"k1 "<<k1<<" k2"<<k2<<std::endl;
	//BUGS:when rx equals ry, I cann't calculate the direction of principal curvatures. In other words,
	//double lambda_max =  k1*F/(N-k1*G);//direction of maximum principal curvature
    //cout<<"lambda_max="<<lambda_max<<endl;
	//Vector3r dir1 = R2*lambda_max+R1;//transfer the direction in the parametric space to R3xyz.
    Vector3r dir1 = R1*(N-k1*G)+R2*k1*F;
    if (dir1.norm()==0){//this is very tricky to check a variable equal to zero.
        //if yes, dir1 isn't able to provide any info of direction
        //we need to compute dir2, i.e., computing direction according to the minimum principal curvature
        curvatureMaxflag = false;//using dir2

    }else{
        dir1.normalize();
        CurvatureDir_max = rot_mat2global*dir1; //restore the direction of maximum principal curvature.
        //cout<<"CurvatureDir_max"<<CurvatureDir_max<<endl;
    }
	Vector3r dir2 = R1*(N-k2*G)+R2*k2*F;
    dir2.normalize();
    CurvatureDir_min = rot_mat2global*dir2;
    //cout<<"CurvatureDir_min"<<CurvatureDir_min<<endl;
	/*double lambda_min =  k2*F/(N-k2*G);
	Vector3r dir2 = R2*lambda_min+R1;
	dir2.normalize();*/
	//Matrix3r Rot = Orientation.toRotationMatrix();//to global system
	//CurvatureDir_max1 = rot_mat2global*dir1; //restore the direction of maximum principal curvature.
    //cout<<"CurvatureDir_max"<<CurvatureDir_max<<endl;
	//CurvatureDir_min = rot_mat2global*dir2;
	//std::cout<<"dir"<<dir1<<'\t'<<"dir2"<<dir2<<'\t'<<"dir1*dir2="<<dir1.dot(dir2)<<std::endl;
    //std::cout<<"direction=="<<lambda_max<<std::endl;
    //test an ellipsoid
    /*
    getCurrentP(phi);
    Vector3d m=CurrentP;
	double tmp = (m(0)*m(0)/pow(rx,4)+m(1)*m(1)/pow(ry,4)+m(2)*m(2)/pow(rz,4));
    double k=pow(rx*ry*rz*tmp,-2);
	double h=-(m(0)*m(0)+m(1)*m(1)+m(2)*m(2)-rx*rx-ry*ry-rz*rz)/(2*pow(rx*ry*rz,2)*pow(tmp,1.5));

    double x1=(h+sqrt(pow(h,2)-k));
    double x2=(h-sqrt(pow(h,2)-k));
    std::cout<<"ellipsoid curvatures="<<x1<<'\t'<<x2<<k1<<'\t'<<k2<<std::endl;
    */
    //std::cout<<"H K h k    "<<H<<'\t'<<K<<'\t'<<h<<'\t'<<k<<std::endl;
    //std::cout<<"aaaa"<<h1+sqrt(h1*h1-k1)<<'\t'<<h1-sqrt(h1*h1-k)<<std::endl;
    return Vector2r(k1,k2);
}


Vector3r Superquadrics::getNormal(Vector2r phi)
{
	Vector3r n;
	n(0) = Mathr::Sign(cos(phi(0))) /rx *pow(fabs(cos(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(1) = Mathr::Sign(sin(phi(0))) /ry *pow(fabs(sin(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(2) = Mathr::Sign(sin(phi(1))) /rz *pow(fabs(sin(phi(1))), 2.-eps2);

	return n;
}


Vector2r Superquadrics::Normal2Phi(Vector3r n)
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
		phi(1) = atan2(Mathr::Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(sin(phi(0))),2.-eps1)), 1/(2.-eps2) ), pow(fabs(ry*n(1)), 1/(2.-eps2) ) );
		//std::cout<<"phi1=111=="<<phi(1)<<std::endl;
	}*/
	return phi;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
////the following fuctions aim to calculating the derivates of Points to Phi.
////////////////////////////////////////
//
/*Mat32r Superquadrics::Derivate_P2Phi(Vector2r phi)  //3*2 matrix
{	Vector3r P;
	P(0) = Mathr::Sign(cos(phi(0)))*rx*pow(fabs(cos(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	P(1) = Mathr::Sign(sin(phi(0)))*ry*pow(fabs(sin(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	P(2) = Mathr::Sign(sin(phi(1)))*rz*pow(fabs(sin(phi(1))), eps2);
	Mat32r p2phi;
	p2phi <<
			P(0)*eps1*fabs(tan(phi(0))), P(0)*eps2*fabs(tan(phi(1))),
			P(1)*eps1*fabs(cot(phi(0))), P(1)*eps2*fabs(tan(phi(1))),
			P(2)*0.					  , P(2)*eps2*fabs(cot(phi(1)));

	return p2phi;
}*/

Mat32r Superquadrics::Derivate_P2Phi(Vector2r phi)  //3*2 matrix
{	Vector3r P;
    /*
	P(0) = Mathr::Sign(cos(phi(0)))*rx*pow(fabs(cos(phi(0))), eps1-1)*pow(fabs(cos(phi(1))), eps2-1);
	P(1) = Mathr::Sign(sin(phi(0)))*ry*pow(fabs(sin(phi(0))), eps1-1)*pow(fabs(cos(phi(1))), eps2-1);
	P(2) = Mathr::Sign(sin(phi(1)))*rz*pow(fabs(sin(phi(1))), eps2-1);
	Mat32r p2phi;
	p2phi <<
			-P(0)*eps1*sin(phi(0))*fabs(cos(phi(1))), -P(0)*eps2*fabs(cos(phi(0)))*sin(phi(1)),
			P(1)*eps1*cos(phi(0))*fabs(cos(phi(1))), -P(1)*eps2*fabs(sin(phi(0)))*sin(phi(1)),
			0.					  , P(2)*eps2*cos(phi(1));
    */

    P(0) = Mathr::Sign(cos(phi(0)))*rx*pow(fabs(cos(phi(0))), eps1-1)*pow(fabs(cos(phi(1))), eps2-1);
	P(1) = Mathr::Sign(sin(phi(0)))*ry*pow(fabs(sin(phi(0))), eps1-1)*pow(fabs(cos(phi(1))), eps2-1);
	P(2) = Mathr::Sign(sin(phi(1)))*rz*pow(fabs(sin(phi(1))), eps2-1);
	Mat32r p2phi;
	p2phi <<
			-Mathr::Sign(phi(0))*P(0)*eps1*sin(phi(0))*fabs(cos(phi(1))), -P(0)*eps2*fabs(cos(phi(0)))*sin(phi(1)),
			Mathr::Sign(phi(0))*P(1)*eps1*cos(phi(0))*fabs(cos(phi(1))), - P(1)*eps2*fabs(sin(phi(0)))*sin(phi(1)),
			0.					  , Mathr::Sign(phi(1))*P(2)*eps2*cos(phi(1));

	return p2phi;
}

Vector3r Superquadrics::Derivate_P2Phi_12(Vector2r phi)
{
	Vector3r Surf;
	Surf = CurrentP;//getSurface(phi);  //Caution: the current p should be updated before invoking this function.
	Surf(0) *= eps1*eps2*fabs(tan(phi(0))*tan(phi(1)));
	Surf(1) *= eps1*eps2*fabs(cot(phi(0))*tan(phi(1)));
	Surf(2) *= 0.;
	return Surf;
}

Vector3r Superquadrics::Derivate_P2Phi_11(Vector2r phi)
{
	Vector3r Surf;
	Surf = CurrentP;//getSurface(phi);  //Caution: the current p should be updated before invoking this function.
	Surf(0) *= eps1*(1.+(eps1-1.)*pow(fabs(tan(phi(0))),2.));
	Surf(1) *= eps1*(1.+(eps1-1.)*pow(fabs(cot(phi(0))),2.));
	Surf(2) *= 0.;
	return Surf;
}

Vector3r Superquadrics::Derivate_P2Phi_22(Vector2r phi)
{
	Vector3r Surf;
	Surf = CurrentP;//getSurface(phi);  //Caution: the current p should be updated before invoking this function.
	Surf(0) *= eps2*(1.+(eps2-1.)*pow(fabs(tan(phi(1))),2.));
	Surf(1) *= eps2*(1.+(eps2-1.)*pow(fabs(tan(phi(0))),2.));
	Surf(2) *= eps2*(1.+(eps2-1.)*pow(fabs(cot(phi(0))),2.));
	return Surf;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
////the following fuctions aim to calculating the derivates of Phi to N.
////////////////////////////////////////

Mat23r Superquadrics::Derivate_Phi2N(Vector2r phi, Vector3r n)
{
	Mat23r phi2n;		//2*3 matrix
	double a = Mathr::Sign(n(1))*pow(fabs(ry*n(1)/rx/n(0)), 1/(2.-eps1));
	double b = a/(1.+ pow(a,2.));


	phi2n(0,0) = -b/(2.-eps1)/fabs(n(0));		//phi1 to n1
	phi2n(0,1) = b/(2.-eps1)/fabs(n(1));		//phi1 to n2
	phi2n(0,2) = 0.;							//phi1 to n3

	//phi2 to n
	if(1)// (fabs(rx*n(0))>fabs(ry*n(1)))		//|r1n1|>|r2n2|
	{
		a = Mathr::Sign(n(2)) * pow(rz*fabs(n(2))* pow(fabs(cos(phi(0))),2.-eps1) / fabs(rx*n(0)), 1/(2.-eps2) ) ;
		b = a/(1.+pow(a,2.));
		//phi2n(1,0) = -b*((2.-eps1)/(2.-eps2)*tan(phi(0))*phi2n(0,0) + 1/(2.-eps2)/n(0));
        phi2n(1,0) = 1/(2.-eps2)*b*((eps1-2.0)*tan(phi(0))*phi2n(0,0) - 1/n(0));

		//phi2n(1,1) = -b*((2.-eps1)/(2.-eps2)*tan(phi(0))*phi2n(0,1) );
        phi2n(1,1) = -b*((2.-eps1)/(2.-eps2)*tan(phi(0))*phi2n(0,1) );
		phi2n(1,2) = Mathr::Sign(phi(1))*b*(1./(2.-eps2)/fabs(n(2)));
	}
	else{
		a = Mathr::Sign(n(2)) * pow(rz*n(2)* pow(fabs(sin(phi(0))),2.-eps1) / fabs(ry*n(1)), 1/(2.-eps2) ) ;
		b = a/(1.+pow(a,2.));
		phi2n(1,0) = b*((2.-eps1)/(2.-eps2)*cot(phi(0))*phi2n(0,0) );
		phi2n(1,1) = b*((2.-eps1)/(2.-eps2)*cot(phi(0))*phi2n(0,1) - 1/(2.-eps2)/fabs(n(1)));
		phi2n(1,2) = b*(1./(2.-eps2)/fabs(n(2)));
	}
    /*
    std::cout<<"phi2n==="<<phi2n<<std::endl;
    std::cout<<"n==="<<n<<std::endl;
    //validate for ellipsoids
    //test1

    a = rz*n(2)* cos(phi(0)) /rx/n(0);
	b = a/(1.+pow(a,2.));
	//phi2n(1,0) = -b*((2.-eps1)/(2.-eps2)*tan(phi(0))*phi2n(0,0) + 1/(2.-eps2)/n(0));
    double p0 = 1/(2.-eps2)*b*((eps1-2.0)*tan(phi(0))*phi2n(0,0) - 1/n(0));

	//phi2n(1,1) = -b*((2.-eps1)/(2.-eps2)*tan(phi(0))*phi2n(0,1) );
    double p1 = -b*(tan(phi(0))*phi2n(0,1));
	double p2 = b*(1./fabs(n(2)));
    std::cout<<"p0="<<p0<<"p1="<<p1<<"p2="<<p2<<std::endl;
    //test2
    double a1 = rx*n(0)/ry/n(1)+ry*n(1)/rx/n(0);
    std::cout<<"phi00="<<-1.0/a1/n(0)<<"phi01="<<1/a1/n(1)<<std::endl;
    double b1 = rz*n(2)*cos(phi(0))/rx/n(0);
    double b2 = b1/(1.0+b1*b1)*(-1.0/n(0)+tan(phi(0))/a1/n(0));
    double c1 = b1/(1.0+b1*b1)*(-tan(phi(0))/a1/n(1));
    double d1 = b1/(1.0+b1*b1)/n(2);
    std::cout<<"phi10="<<b2<<"phi11="<<c1<<"phi13="<<d1<<std::endl;
    */
	return phi2n;
}
//////////////////derivates of contact normal to alpha_i, i=1,2
Mat32r Superquadrics::Derivate_C2Alpha(Vector2r alpha, Matrix3r RotationMat)
{
	Mat32r c2alpha;
	c2alpha<<
			-sin(alpha(0))*cos(alpha(1)), -cos(alpha(0))*sin(alpha(1)),
			cos(alpha(0))*cos(alpha(1)), -sin(alpha(0))*sin(alpha(1)),
			0., cos(alpha(1));

	c2alpha = RotationMat * c2alpha;

	return c2alpha;
}
Vector3r Superquadrics::P_alpha12(Vector2r para, Vector3r &n,Matrix3r globalContact,Matrix3r RotationMat, int sign)
{
	//para: [alpha1, alpha2]
	//point:[x1,x2,x3], the global Cartesian coordinates of the current point
	//m:number of parameters (alpha1, alpha2), and it is constantly equal to 2 here.
	//n:number of non-linear functions. These fuctions are p(alpha1,alpha2)_i, i=1,2,3. Thus, it is 3 constantly.
	//RotationMat: rotation matrix to particle's local system, equal to Orientation.conjugate().toRotationMatrix();
	 //globalContact: Rotation matrix of the unit base vectors e_i (i=1,2,3) where the contact vector is defined.
	//double alpha1=para(0), alpha2=para(1);
	//int sign = 1;

	Vector3r contact(0.,0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.
	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	contact = globalContact*contact;	   //to global
	n = Mathr::Sign(sign)*contact;			//return the normal vector in the global system
	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact = RotationMat*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	Vector2r phi= Normal2Phi( Mathr::Sign(sign)*contact);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	return getSurface( phi ); //in particle's local
}

Vector3r Superquadrics::P_alpha12(Vector2r para,Matrix3r globalContact,int sign)
{
	//para: [alpha1, alpha2]
	//point:[x1,x2,x3], the global Cartesian coordinates of the current point
	//m:number of parameters (alpha1, alpha2), and it is constantly equal to 2 here.
	//n:number of non-linear functions. These fuctions are p(alpha1,alpha2)_i, i=1,2,3. Thus, it is 3 constantly.
	//globalContact: Rotation matrix of the unit base vectors e_i (i=1,2,3) where the contact vector is defined.
	double alpha1=para(0), alpha2=para(1);
	//int sign = 1;

	Vector3r contact(0.,0.,0.);

	contact(0) = cos(alpha1)*cos(alpha2);  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(alpha1)*cos(alpha2);  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(alpha2);			   //e_i is pointing in the direction from the first particle center to the second particle center.
	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	contact = globalContact*contact;	   //to global
	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact = Orientation.conjugate().toRotationMatrix()*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	Vector2r phi= Normal2Phi( Mathr::Sign(sign)*contact);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	return getSurface( phi );
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	//std::cout<<"sign= "<<sign<<"\t point= "<<p<<std::endl;

}

Vector3r Superquadrics::P_alpha12(Vector2r para,Matrix3r globalContact,Matrix3r RotationMat,int sign)
{
	//para: [alpha1, alpha2]
	//point:[x1,x2,x3], the global Cartesian coordinates of the current point
	//m:number of parameters (alpha1, alpha2), and it is constantly equal to 2 here.
	//n:number of non-linear functions. These fuctions are p(alpha1,alpha2)_i, i=1,2,3. Thus, it is 3 constantly.
	//RotationMat: rotation matrix to particle's local system, equal to Orientation.conjugate().toRotationMatrix();
	//globalContact: Rotation matrix of the unit base vectors e_i (i=1,2,3) where the contact vector is defined.
	//double alpha1=para(0), alpha2=para(1);
	//int sign = 1;

	Vector3r contact(0.,0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.
	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	contact = globalContact*contact;	   //to global
	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact = RotationMat*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	Vector2r phi= Normal2Phi( Mathr::Sign(sign)*contact);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	return getSurface( phi );

}


Mat32r Superquadrics::JacFunc(Vector2r para, Matrix3r globalContact,int sign)
{
	//para: [alpha1, alpha2]
	//jac:Jacobian of those non-linear functions, p(alpha1, alpha2)_i, i=1,2,3
	//m:number of parameters (alpha1, alpha2), and it is constantly equal to 2 here.
	//n:number of non-linear functions. These fuctions are p(alpha1,alpha2)_i, i=1,2,3. Thus, it is 3 constantly.
	//double alpha1=para(0), alpha2=para(1);
	//int sign = 1;

	Vector3r p(0.,0.,0.), contact(0.,0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	contact = globalContact*contact;	   //to global
	Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	contact = Rot*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	Vector2r phi = Normal2Phi( Mathr::Sign(sign)*contact);

	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	Mat32r d_p2phi = Rot.inverse()*Derivate_P2Phi(phi);//to global system

	Mat23r d_phi2n = Derivate_Phi2N(phi, Mathr::Sign(sign)*contact);

	Mat32r c2alpha = Mathr::Sign(sign)*Rot*Derivate_C2Alpha(para, globalContact);

	return d_p2phi*d_phi2n*c2alpha;
	//std::cout<<"JacMat"<<JacMat<<std::endl;
	/*
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<2; j++)
			jac[i*2+j] = JacMat(i,j);
	}
	*/
}

//**********************
Mat32r Superquadrics::JacFunc(Vector2r para, Matrix3r globalContact, Matrix3r RotationMat, int sign)
{
	//para: [alpha1, alpha2]
	//jac:Jacobian of those non-linear functions, p(alpha1, alpha2)_i, i=1,2,3
	//m:number of parameters (alpha1, alpha2), and it is constantly equal to 2 here.
	//n:number of non-linear functions. These fuctions are p(alpha1,alpha2)_i, i=1,2,3. Thus, it is 3 constantly.
	//RotationMat: rotation matrix to particle's local system, equal to Orientation.conjugate().toRotationMatrix();
	//double alpha1=para(0), alpha2=para(1);
	//int sign = 1;

	Vector3r p(0.,0.,0.), contact(0.,0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	contact = globalContact*contact;	   //to global
	//Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	contact = RotationMat*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	Vector2r phi = Normal2Phi( Mathr::Sign(sign)*contact);

	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	Mat32r d_p2phi = RotationMat.inverse()*Derivate_P2Phi(phi);//to global system

	Mat23r d_phi2n = Derivate_Phi2N(phi, Mathr::Sign(sign)*contact);

	Mat32r c2alpha = Mathr::Sign(sign)*RotationMat*Derivate_C2Alpha(para, globalContact);

	return d_p2phi*d_phi2n*c2alpha;
}
///this function is for testing...
Vector3r Superquadrics::N_alpha12(Vector2r para, Matrix3r globalContact,int sign)
{	double alpha1=para(0), alpha2=para(1);
	Vector3r p(0.,0.,0.), contact(0.,0.,0.);

	contact(0) = cos(alpha1)*cos(alpha2);  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(alpha1)*cos(alpha2);  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(alpha2);			   //e_i is pointing in the direction from the first particle center to the second particle center.

	contact = globalContact*contact;	   //to global
	Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	contact = Rot*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	Vector2r phi = Normal2Phi( Mathr::Sign(sign)*contact);
	Vector3r N = Rot.inverse()*getNormal(phi);//to global system
	N.normalize();
	//Vector3r p1 = getSurface( phi );
	return N;
	//return p1+0.5*N;
}

//*************************************************************************************
//****************************************************************************************
/* AaBb overlap checker  */

void Bo1_Superquadrics_Aabb::go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*){

	Superquadrics* t=static_cast<Superquadrics*>(ig.get());
	//t->rot_mat2local = (se3.orientation).conjugate().toRotationMatrix();//to particle's system
	//t->rot_mat2global = (se3.orientation).toRotationMatrix(); //to global system
	/*
	Matrix3r rot_mat1 = (se3.orientation).conjugate().toRotationMatrix();//to particle's system
	Matrix3r rot_mat2 = (se3.orientation).toRotationMatrix();//to global system
	Vector3r t_pos = se3.position;
	//find the furthest point in a direction
	//Vector3r n_x_min = Vector3r(-1,0,0);
	Vector3r n_x_max = Vector3r(1,0,0);
	//Vector3r n_y_min = Vector3r(0,-1,0);
	Vector3r n_y_max = Vector3r(0,1,0);
	//Vector3r n_z_min = Vector3r(0,0,-1);
	Vector3r n_z_max = Vector3r(0,0,1);
	//using the symmetricity of the shape, thus x_max -pos_x = pos_x - x_min
	//x axis
	Vector2r phi = t->Normal2Phi(rot_mat1*n_x_max);//phi in particle's local system
	Vector3r p_x = rot_mat2*( t->getSurface(phi));	//
	//y axis
	phi = t->Normal2Phi(rot_mat1*n_y_max);//phi in particle's local system
	Vector3r p_y = rot_mat2*( t->getSurface(phi));	//
	//z axis
	phi = t->Normal2Phi(rot_mat1*n_z_max);//phi in particle's local system
	Vector3r p_z = rot_mat2*( t->getSurface(phi));	//
	*/
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
	Vector3r mincoords(-r_max,-r_max,-r_max), maxcoords(r_max,r_max,r_max);
	//the aabb should be optimized in future.

	aabb->min=se3.position+mincoords;
	aabb->max=se3.position+maxcoords;
}

//**********************************************************************************
/* Plotting */

#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	bool Gl1_Superquadrics::wire;
	int Gl1_Superquadrics::slices;
	int Gl1_Superquadrics::pre_slices = -1;
    int Gl1_Superquadrics::glSolidList=-1;
    int Gl1_Superquadrics::glWireList=-1;
	//vector<Vector3r> Gl1_Superquadrics::vertices;

void Gl1_Superquadrics::initSolidGlList(Superquadrics* t)
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
	//std::cerr << "Gl1_Superquadrics:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";

	vector<Vector3r> vertices;
    Vector3r vert;
	vertices.clear();
	for(a=0.0,i=0;i<h;i++,a+=hStep)
	  for(b=0.0,j=0;j<w;j++,b+=wStep)
	  {
		vertices.push_back(t->getSurface(Vector2r(b,a-M_PI_2l)));
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

void Gl1_Superquadrics::initWireGlList(Superquadrics* t)
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
	//std::cerr << "Gl1_Superquadrics:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";

	vector<Vector3r> vertices;
    Vector3r vert;
	vertices.clear();
	for(a=0.0,i=0;i<h;i++,a+=hStep)
	  for(b=0.0,j=0;j<w;j++,b+=wStep)
	  {
		vertices.push_back(t->getSurface(Vector2r(b,a-M_PI_2l)));
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


	void Gl1_Superquadrics::drawSlice(Vector3r &p1,Vector3r &p2,Vector3r &p3,Vector3r &p4)
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

	void Gl1_Superquadrics::go(const shared_ptr<Shape>& cm, const shared_ptr<State>&state,bool wire2,const GLViewInfo&)
	{
		glMaterialv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,Vector3r(cm->color[0],cm->color[1],cm->color[2]));
		glColor3v(cm->color);
		Superquadrics* t=static_cast<Superquadrics*>(cm.get());
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
	void Gl1_SuperquadricsGeom::go(const shared_ptr<IGeom>& ig, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool){
	        SuperquadricsGeom* pg=static_cast<SuperquadricsGeom*>(ig.get());
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
///function for calculation of contact force for Hertz-Mindlin model
void SuperquadricsGeom2::HertzMindlin(Vector2r curvature1, Vector2r curvature2, double theta){
    //double rho11, double rho12, double rho21, double rho22,
	double rho11(0), rho12(0), rho21(0), rho22(0);//curvatures from two surfaces
    rho11 = fabs(curvature1(0));
    rho12 = fabs(curvature1(1));
    rho21 = fabs(curvature2(0));
    rho22 = fabs(curvature2(1));
    if (rho11*rho12 <1E-18){//K=0 Fix me!
        rho11=rho12=0.0;
        cout<<"rho1112===00000000000"<<endl;
    }
    if(rho21*rho22 <1E-18){//K=0 Fix me!
        rho21=rho22=0.0;
        cout<<"rho2122===00000000000"<<endl;
    }
	//double theta(0);//angle between two maximum principal curvatures
    //cout<<"rho==="<<rho11<<'\t'<<rho12<<'\t'<<rho21<<'\t'<<rho22<<" angle "<<theta<<endl;
	double sum_rho = rho11+rho12+rho21+rho22;
    //cout<<"ffff"<<pow(rho11-rho12,2)+pow(rho21-rho22,2)+2*(fabs(rho11-rho12))*fabs(rho21-rho22)*cos(2*theta)<<endl;
	double F_rho = -sqrt(pow(rho11-rho12,2)+pow(rho21-rho22,2)+2*(fabs(rho11-rho12))*fabs(rho21-rho22)*cos(2*theta));
    //cout<<"sum_rho,F_rho==="<<sum_rho<<'\t'<<F_rho<<endl;
	double A = (sum_rho+F_rho)/4.;
	double B = (sum_rho-F_rho)/4.;
	/*if(isnan(A)||isnan(B)){//for debugging
		cout<<"rho==="<<rho11<<'\t'<<rho12<<'\t'<<rho21<<'\t'<<rho22<<" angle "<<theta<<endl;
		cout<<"sum_rho,F_rho==="<<sum_rho<<'\t'<<F_rho<<endl;
		cout<<"B,A==="<<B<<'\t'<<A<<endl;
	}*/
    curvatures_sum = A+B;//save info
    //cout<<"B,A==="<<B<<'\t'<<A<<endl;
	//k=(B/A)^gamma'
	//approximate numerical method to calculate k
	//see the reference by J. F. Antoine et al. "Approximate Analystical Model for HERIZIAN Elliptical Contact Problems"

	double x_squ = pow(log10(B/A),2);//log10
	double gamma1 = 2*(1+x_squ*(0.40227436+x_squ*(3.7491752e-2+x_squ*(7.4855761e-4+x_squ*2.1667028e-6))));
	double gamma2 = 3*(1+x_squ*(0.42678878+x_squ*(4.2605401e-2+x_squ*(9.0786922e-4+x_squ*2.7868927e-6))));

	//double k = pow((B/A), gamma1/gamma2);
    alpha = pow((B/A), gamma1/gamma2);

	double m1 = 1/pow(alpha,2);
	//cout<<"m1"<<m1<<"B/A"<<B/A<<endl;
	K = (1.3862944+m1*(0.1119723+0.0725296*m1))-log(m1)*(0.5+m1*(0.1213478+0.0288729*m1));
	E = (1+m1*(0.4630151+0.1077812*m1))-log(m1)*(0.2452727+0.0412496*m1)*m1;
    //cout<<"E_K F_K   k"<<E<<'\t'<<K<<'\t'<<alpha<<endl;
	//double Fn = 4./3.*M_PIl*k*sqrt(E_k/pow(F_k,3)/sum_rho);

	//cout<<"Fn============="<<Fn<<endl;

	//Mindlin tangential stiffness coefficient
	//angle between rho1 and semimajor of the contact ellipse
}
//**********************************************************************************
//!Precompute data needed for rotating tangent vectors attached to the interaction

void SuperquadricsGeom::precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r&
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
Vector3r& SuperquadricsGeom::rotate(Vector3r& shearForce) const {
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

void Ip2_SuperquadricsMat_SuperquadricsMat_SuperquadricsPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<SuperquadricsMat>& mat1 = SUDODEM_PTR_CAST<SuperquadricsMat>(b1);
	const shared_ptr<SuperquadricsMat>& mat2 = SUDODEM_PTR_CAST<SuperquadricsMat>(b2);
	interaction->phys = shared_ptr<SuperquadricsPhys>(new SuperquadricsPhys());
	const shared_ptr<SuperquadricsPhys>& contactPhysics = SUDODEM_PTR_CAST<SuperquadricsPhys>(interaction->phys);
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
/* Hertz-Mindlin Material law, physics */
void Ip2_SuperquadricsMat2_SuperquadricsMat2_HertzMindlinPhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;

	const shared_ptr<SuperquadricsGeom2>& contactGeom(SUDODEM_PTR_DYN_CAST<SuperquadricsGeom2>(interaction->geom));
	//if(!contactGeom) {return;}

	const shared_ptr<SuperquadricsMat2>& mat1 = SUDODEM_PTR_CAST<SuperquadricsMat2>(b1);
	const shared_ptr<SuperquadricsMat2>& mat2 = SUDODEM_PTR_CAST<SuperquadricsMat2>(b2);
	interaction->phys = shared_ptr<SuperquadricsHertzMindlinPhys>(new SuperquadricsHertzMindlinPhys());
	const shared_ptr<SuperquadricsHertzMindlinPhys>& contactPhysics = SUDODEM_PTR_CAST<SuperquadricsHertzMindlinPhys>(interaction->phys);
    const Body::id_t id= interaction->id1;
    //zero penetration depth means no interaction force
	//if(!(contactGeom->PenetrationDepth > 1E-18) ) return;
    if (id>5){

		contactPhysics->nuab 	= 0.5*(mat1->nu+mat2->nu);//average Possion's ratio
		contactPhysics->Gab 	= 1.0/((1-mat1->nu)/mat1->G+(1-mat2->nu)/mat2->G);//average shear modulus
		contactPhysics->mu1 = 1./mat1->G+1./mat2->G;
		contactPhysics->mu2 = mat1->nu/mat1->G+mat2->nu/mat2->G;
		Real frictionAngle = std::min(mat1->frictionAngle,mat2->frictionAngle);
    	contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle);


	    //cout<<"ksx"<<ksx<<"ksy"<<ksy<<"ks"<<ks<<"Kn"<<Kn<<"Kn/ks"<<Kn/ks<<endl;
	    //compute Hertz normal contact stiffness

        //contactPhysics->kn = 3e4;//Kn;//2.0*mat1->Kn*mat2->Kn/(mat1->Kn+mat2->Kn);
	    //contactPhysics->ks = 3e4;//ks;//2.0*mat1->Ks*mat2->Ks/(mat1->Ks+mat2->Ks);

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
Real SuperquadricsLaw::getPlasticDissipation() {return (Real) plasticDissipation;}
void SuperquadricsLaw::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}
Real SuperquadricsLaw::elasticEnergy()
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
// Apply forces on superquadrics in collision based on geometric configuration
bool SuperquadricsLaw::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){

		if (!I->geom) {return true;}
		const shared_ptr<SuperquadricsGeom>& contactGeom(SUDODEM_PTR_DYN_CAST<SuperquadricsGeom>(I->geom));
		if(!contactGeom) {return true;}
		const Body::id_t idA=I->getId1(), idB=I->getId2();
		const shared_ptr<Body>& A=Body::byId(idA), B=Body::byId(idB);

                State* de1 = Body::byId(idA,scene)->state.get();
	        State* de2 = Body::byId(idB,scene)->state.get();

		SuperquadricsPhys* phys = dynamic_cast<SuperquadricsPhys*>(I->phys.get());

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
		Vector3r shearForce1=Vector3r::Zero();//zhswee deprecated SuperquadricsLaw.shearForce
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
          const char* filename = "SuperquadricsLaw_err.log";
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

//**************************************************************************************
// Apply forces on superquadrics in collision based on geometric configuration
bool SuperquadricsLaw2::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){

		if (!I->geom) {return true;}
		const shared_ptr<SuperquadricsGeom2>& contactGeom(SUDODEM_PTR_DYN_CAST<SuperquadricsGeom2>(I->geom));
		if(!contactGeom) {return true;}
		const Body::id_t idA=I->getId1(), idB=I->getId2();
		const shared_ptr<Body>& A=Body::byId(idA), B=Body::byId(idB);

        State* de1 = Body::byId(idA,scene)->state.get();
	    State* de2 = Body::byId(idB,scene)->state.get();

		SuperquadricsHertzMindlinPhys* phys = dynamic_cast<SuperquadricsHertzMindlinPhys*>(I->phys.get());

		//erase the interaction when aAbB shows separation, otherwise keep it to be able to store previous separating plane for fast detection of separation
		if (A->bound->min[0] >= B->bound->max[0] || B->bound->min[0] >= A->bound->max[0] || A->bound->min[1] >= B->bound->max[1] || B->bound->min[1] >= A->bound->max[1] || A->bound->min[2] >= B->bound->max[2] || B->bound->min[2] >= A->bound->max[2])  {
		        phys->normalForce = Vector3r(0.,0.,0.); phys->shearForce = Vector3r(0.,0.,0.);
			scene->interactions->requestErase(I);
			return false;
		}

		//zero penetration depth means no interaction force
		if(!(contactGeom->PenetrationDepth > 1E-18) ) {
            phys->normalForce = Vector3r(0.,0.,0.); phys->shearForce = Vector3r(0.,0.,0.);
            return true;
        }
		////////

    	if (idA>5){
			Real depth = contactGeom->PenetrationDepth;
			Real curvatures_sum = contactGeom->curvatures_sum;
			if (contactGeom->isSphere){

				Real tmp = sqrt(depth/curvatures_sum);
				phys->kn = 8./3.*phys->Gab*tmp;
				phys->ks = 8./(2.0*phys->mu1-phys->mu2)*tmp;
				//cout<<"geom-Kn="<<Kn<<"Ks"<<endl;
			}else{

				//contact geometry info


				Real E = contactGeom->E;
				Real K = contactGeom->K;
				Real alpha = contactGeom->alpha;
				//cout<<"K"<<K<<" E"<<E<<" alpha"<<alpha<<"depth"<<depth<<endl;
				Real a = sqrt(depth*E/curvatures_sum/K);//semi-length of contact ellipse
				phys->kn = 4./3.*M_PIl*a*phys->Gab*alpha/K;
				Real e2 = 1-1/pow(alpha,2.0);
				Real ksx = phys->mu1*K-phys->mu2*(K-E)/e2;
				Real ksy = phys->mu1*K+phys->mu2*((1-e2)*K-E)/e2;
				phys->ks = M_PIl*alpha*a*(1./ksx+1./ksy);
			}
		}
		Vector3r normalForce=contactGeom->normal*contactGeom->PenetrationDepth*phys->kn;
		//std::cerr << "Kn"<<phys->kn<<" ks"<<phys->ks<< "\n";
		Vector3r shearForce1=Vector3r::Zero();//zhswee deprecated SuperquadricsLaw.shearForce
		//shear force: in case the polyhdras are separated and come to contact again, one should not use the previous shear force
		//std::cerr << "p1:"<<gettid4()<<shearForce[0]<< "\n";
		///////////////////////
		bool useDamping=(phys->betan >0.00001 || phys->betas > 0.00001);//using viscous damping?
		// tangential and normal stiffness coefficients, recomputed from betan,betas at every step
	        Real cn=0, cs=0;
	        Vector3r normalViscous = Vector3r::Zero();
	        Vector3r shearViscous = Vector3r::Zero();
        if (useDamping){//
            //std::cerr << "viscous damping"<< "\n";
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
			}else{normalForce -= normalViscous;}
                //std::cerr << "nF2:"<<normalForce.norm()<< "\n";
        }
		if (contactGeom->isShearNew){
			//shearForce = Vector3r::Zero();
			shearForce1 = Vector3r(0.,0.,0.);
		}else{
			//shearForce = contactGeom->rotate(shearForce);
			shearForce1 = contactGeom->rotate(phys->shearForce);
		}
			//std::cerr << "p2:"<<gettid4() <<shearForce[0]<< "\n";
		const Vector3r& shearDisp = contactGeom->shearInc;
		//std::cerr << "shearDisplacement="<<shearDisp.norm()<< "\n";
		shearForce1 -= phys->ks*shearDisp;
                //std::cerr << "p3:"<<gettid4() <<shearDisp[0]<< "\n";
		Real maxFs = normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);
		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if( shearForce1.squaredNorm() > maxFs ){
			Real ratio = sqrt(maxFs) / shearForce1.norm();
			shearForce1 *= ratio;}
		else if (useDamping){ // add current contact damping if we do not slide and if damping is requested
		        shearViscous = cs*contactGeom->relativeVs; // get shear viscous component
		        //std::cerr << "sF1:"<<shearForce1(0)<<" "<<shearForce1(1)<<" "<<shearForce1(2)<< "\n";
		        //std::cerr << "sF1:"<<shearForce1.norm()<< "\n";
		        shearForce1 -= shearViscous;
		        //std::cerr << "sF2:"<<shearForce1.norm()<< "\n";

		}
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
		  //char filename[50];
          const char* filename = "SuperquadricLaw2_err.log";
		  //sprintf(filename,"%d",gettid4());
		  FILE * fin = fopen(filename,"a");
          Vector3r normal = contactGeom->normal;
          Real depth = contactGeom->PenetrationDepth;
		  //fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
		  fprintf(fin,"shear_force\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
		  fprintf(fin,"shear_force1\t%e\t%e\t%e\n",shearForce1[0],shearForce1[1],shearForce1[2]);
		  //fprintf(fin,"F\t%e\t%e\t%e\n",F[0],F[1],F[2]);
			//Hertz-Mindlin parameters
		  fprintf(fin,"E\t%eK\t%ealpha\t%ecurvatures_sum\t%e\n",contactGeom->E,contactGeom->K,contactGeom->alpha,contactGeom->curvatures_sum);
		  //
		  fprintf(fin,"ks\t%emaxFs\t%esfn\t%e\n",phys->ks,maxFs,shearForce.squaredNorm());
		  fprintf(fin,"shearDisp\t%e\t%e\t%e\n",shearDisp[0],shearDisp[1],shearDisp[2]);
		  fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
          fprintf(fin,"contact normal \t%e\t%e\t%e\n",normal[0],normal[1],normal[2]);
          fprintf(fin,"depth\t%e\n",depth);
		  fclose(fin);
		  //throw runtime_error("P4 #");
          LOG_WARN("Error in computation of contact forces");
          scene->stopAtIter = scene->iter + 1;
          normalForce = Vector3r::Zero();
          shearForce1 = Vector3r::Zero();
		}
		Vector3r F = -normalForce-shearForce1;	//the normal force acting on particle A (or 1) is equal to the normal force.
		if (contactGeom->PenetrationDepth != contactGeom->PenetrationDepth) exit(1);
		scene->forces.addForce (idA,F);
		scene->forces.addForce (idB, -F);
		scene->forces.addTorque(idA, -(A->state->pos-contactGeom->contactPoint).cross(F));
		scene->forces.addTorque(idB, (B->state->pos-contactGeom->contactPoint).cross(F));//FIXME:shoud add shift2 to B's pos when considering periodic boundary
		//needed to be able to acces interaction forces in other parts of sudodem
		shearForce=shearForce1;
		phys->normalForce = normalForce;
		phys->shearForce = shearForce1;
		//may need to add phys->shearViscous and phys->normalViscous
		return true;
}
