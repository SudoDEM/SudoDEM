#include "Superquadrics.hpp"
#include "Superquadrics_Ig2.hpp"
#include <cmath>
//#include <sys/syscall.h>
#define _USE_MATH_DEFINES

#define _USE_LM_OPTIMIZATION
#define _DEBUG_LM_NM1
//LM: Levenberg-Marquardt method
//NM: Nerld-Mead simplex method

#include<time.h>

SUDODEM_PLUGIN(/* self-contained in hpp: */ (Ig2_Superquadrics_Superquadrics_SuperquadricsGeom)
        (Ig2_Superquadrics_Superquadrics_SuperquadricsGeom2)
		(Ig2_Wall_Superquadrics_SuperquadricsGeom)
        (Ig2_Wall_Superquadrics_SuperquadricsGeom2)
		(Ig2_Facet_Superquadrics_SuperquadricsGeom)
	   );

//**********************************************************************************
/*some auxiliary functions*/

//////////////////////////////////////////////
#define LM_INFO_SZ    	 10
#define LM_ERROR         -1
#define LM_INIT_MU    	 1E-03
#define LM_STOP_THRESH	 1E-17

#define EPSILON       1E-12
#define ONE_THIRD     0.3333333333333334 /* 1.0/3.0 */
#define LM_REAL_MIN -DBL_MAX
#define LM_CNST(x) (x)
#define LM_REAL_MAX DBL_MAX
#define LM_FINITE finite
#define LM_ISINF(x) isinf(x)

//using namespace std;

template <typename T>
Vector3r Disfunction(Vector2r para,Vector3r &center, T *particle1, T *particle2,
		     Matrix3r RotationMat1,Matrix3r RotationMat2,
		     Vector3r position1, Vector3r position2,
		     Matrix3r globalContact, bool &stop_flag)
{	Vector3r d,p1,p2,n1,n2;
	stop_flag = false;
	p2 = particle2->P_alpha12(para, n2, globalContact, RotationMat2, -1);
	p1 = particle1->P_alpha12(para, n1, globalContact, RotationMat1, 1);
	//to global system
	p2 = RotationMat2.inverse()*p2 + position2;
	p1 = RotationMat1.inverse()*p1 + position1;
	d = p2 -p1;
	center = 0.5*(p1 + p2);
	//we should to check whether the angle between n1 and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (n1.dot(d) > 0)  //the angle < pi/2.
	{
		stop_flag = true;
	}
	return d;//CAUTION:d is the vector from point 1 on the particle1 to point 2 on the particle2.

}
template <typename T>
void Disfunction(Vector2r para,Vector3r &center, T *particle1, T *particle2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2,
		     Matrix3r globalContact, bool &stop_flag)
{
	Vector3r d;
	stop_flag = false;
	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
        Vector2r phi(0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	contact = globalContact*contact;	   //to global
	//we should to check whether the angle between n1 (contact) and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (contact.dot(d) > 0)  //the angle < pi/2.
	{
		stop_flag = true;
		return ;//FIXME:May not be zero
	}
	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi( contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1 + position1;
	//point 2 on Particle 2
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	phi = particle2->Normal2Phi( (-1)*contact);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
        //get the coordinates of the current point in the global system
	p2 = particle2->getSurface( phi );
	p2 = particle2->rot_mat2global*p2 + position2;

	d = p2 -p1;
	center = 0.5*(p1 + p2);

	return ;//CAUTION:d is the vector from point 1 on the particle1 to point 2 on the particle2.

}
/////

template <typename T>
double Disfunction(Vector2r para,Vector3r &center, T *particle1, T *particle2,
		     Matrix3r RotationMat1,Matrix3r RotationMat2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2,
		     Matrix3r globalContact, bool &stop_flag)
{	Vector3r d,n1,n2;
	stop_flag = false;
	p2 = particle2->P_alpha12(para, n2, globalContact, RotationMat2, -1);
	p1 = particle1->P_alpha12(para, n1, globalContact, RotationMat1, 1);
	//to global system
	//p2 = RotationMat2.inverse()*p2 + position2;
	//p1 = RotationMat1.inverse()*p1 + position1;
        p2 = particle2->rot_mat2global*p2 + position2;
	p1 = particle1->rot_mat2global*p1 + position1;
	d = p2 -p1;
	center = 0.5*(p1 + p2);
	//we should to check whether the angle between n1 and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (n1.dot(d) > 0)  //the angle < pi/2.
	{
		stop_flag = true;
	}
	return d.norm();//CAUTION:d is the vector from point 1 on the particle1 to point 2 on the particle2.

}

template <typename T>//optimal function to compute the distance
double Disfunction(Vector2r para,Vector3r &center, T *particle1, T *particle2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2,
		     Matrix3r globalContact, bool &stop_flag)
{	Vector3r d;
	stop_flag = false;
	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
        Vector2r phi(0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	contact = globalContact*contact;	   //to global

	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi( contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1 + position1;
	//point 2 on Particle 2
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	phi = particle2->Normal2Phi( (-1)*contact);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
        //get the coordinates of the current point in the global system
	p2 = particle2->getSurface( phi );
	p2 = particle2->rot_mat2global*p2 + position2;

	d = p2 -p1;
	center = 0.5*(p1 + p2);
	//we should to check whether the angle between n1 (contact) and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (contact.dot(d) > 0)  //the angle < pi/2.
	{
		stop_flag = true;
		return 0;//FIXME:May not be zero
	}
	return d.norm();//CAUTION:d is the vector from point 1 on the particle1 to point 2 on the particle2.

}


template <typename T>
Mat32r Jacfunction(Vector2r para, T *particle1, T *particle2,
		   Matrix3r RotationMat1,Matrix3r RotationMat2,
		   Matrix3r globalContact)
{
	Mat32r JacMat1, JacMat2;
	JacMat1 = particle1->JacFunc(para, globalContact, RotationMat1, 1);
	JacMat2 = particle2->JacFunc(para, globalContact, RotationMat2, -1);
	//std::cout<<"JacMat1"<<JacMat1<<std::endl;

	return JacMat2 - JacMat1;
}

template <typename T>
Mat32r Jacfunction(Vector2r para, T *particle1, T *particle2,
		   Matrix3r globalContact)
{
	Mat32r JacMat1, JacMat2;

	Vector3r p(0.,0.,0.), contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
        Vector2r phi;
	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	contact = globalContact*contact;	   //to global
	//Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	phi = particle1->Normal2Phi( contact1);

	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	Mat32r d_p2phi = particle1->rot_mat2global*particle1->Derivate_P2Phi(phi);//to global system
    //cout<<"d_p2phi="<<d_p2phi<<endl;
	Mat23r d_phi2n = particle1->Derivate_Phi2N(phi, contact1);
    //cout<<"d_phi2n="<<d_phi2n<<endl;
	Mat32r c2alpha = particle1->rot_mat2local*particle1->Derivate_C2Alpha(para, globalContact);
    //cout<<"c2alpha ="<<c2alpha<<endl;
	JacMat1 = d_p2phi*d_phi2n*c2alpha;
    //cout<<"Jacobian mat1="<<JacMat1<<endl;
	//
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	phi = particle2->Normal2Phi( (-1)*contact2);

	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	d_p2phi = particle2->rot_mat2global*particle2->Derivate_P2Phi(phi);//to global system

	d_phi2n = particle2->Derivate_Phi2N(phi, (-1)*contact2);

	c2alpha = (-1)*particle2->rot_mat2local*particle2->Derivate_C2Alpha(para, globalContact);

	JacMat2 = d_p2phi*d_phi2n*c2alpha;
    //cout<<"Jacobian mat2="<<JacMat2<<endl;
	return JacMat2 - JacMat1;
}

template <typename T>
void CheckDC(Vector2r &para, T *particle1, T *particle2, Matrix3r globalContact)
{
	Vector3r p(0.,0.,0.), contact(0.,0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	contact = globalContact*contact;	   //to global


	Vector3r d;


	d = particle2->P_alpha12(para, globalContact,-1) - particle1->P_alpha12(para, globalContact,1);

	double phi = acos(d.dot(contact)/(d.norm()*contact.norm()));//phi in [0, pi]
	cout<<"check whether d is along c:--phi="<<phi<<endl;
}

template <typename T>
double Check_DC(Vector2r &para, T *particle1, T *particle2, Matrix3r globalContact)
{
	Vector3r p(0.,0.,0.), contact(0.,0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	contact = globalContact*contact;	   //to global


	Vector3r d;


	d = particle2->P_alpha12(para, globalContact,-1) - particle1->P_alpha12(para, globalContact,1);

	double phi = acos(d.dot(contact)/(d.norm()*contact.norm()));//phi in [0, pi]
	//cout<<"check whether d is along c:--phi="<<phi<<endl;
	return phi;
}
Matrix3r CRotationMat(Vector3r position1, Vector3r position2)//positions of particle 1 and particle 2.
{
	//Matrix3r globalContact;
	Vector3r e1,n1,x1(1.,0,0);
	Quaternionr p;
	e1 = position2 - position1;		//vector with direction from the center of particle 1 to that of particle 2.
	//calculat the angle between e1 and x1
	//double phi = acos(x1.dot(e1)/(x1.norm()*e1.norm()));//phi in [0, pi]
	double cos_phi = x1.dot(e1)/(x1.norm()*e1.norm());
	double cos_halfphi = sqrt(0.5*(1+cos_phi));
	double sin_halfphi = sqrt(0.5*(1-cos_phi));
	if(abs(cos_phi)>1.0-1e-10){//FIX BUG: e1 is very close to x

		p = Quaternionr(cos_halfphi,0,0,sin_halfphi);
	}else{
		n1 = x1.cross(e1);
		n1.normalize();//caution:n1 must be normalized as the normal unit vector
		n1 *=sin_halfphi;

		p = Quaternionr(cos_halfphi,n1(0),n1(1),n1(2));
	}

	//p.normalize();
	return p.toRotationMatrix();
}

////////////////optimazation
//Vector3r (*func)(Vector2r p, Vector3r &center, T *particle1, T *particle2,
//		   Matrix3r RotationMat1,Matrix3r RotationMat2,
//		   Vector3r position1, Vector3r position2,
//		   Matrix3r globalContact, bool &stop_flag)
template <typename T>	//return distance
void lmder_z(
  void (*func)(Vector2r p, Vector3r &center, T *particle1, T *particle2,
		   Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2,
		   Matrix3r globalContact, bool &stop_flag), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  Mat32r (*jacf)(Vector2r p, T *particle1, T *particle2, Matrix3r globalContact),  /* function to evaluate the Jacobian \part x / \part p */
  Vector2r &p,         /* I/O: initial parameter estimates. On output has the estimated solution */
  int itmax,          /* I: maximum number of iterations */
  double opts[4],    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
                       * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used
                       */
  double info[LM_INFO_SZ],//,
					           /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]= ||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small Dp
                      *                                 3 - stopped by itmax
                      *                                 4 - singular matrix. Restart from current p with increased mu
                      *                                 5 - no further error reduction is possible. Restart with increased mu
                      *                                 6 - stopped by small ||e||_2
                      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                      * info[7]= # function evaluations
                      * info[8]= # Jacobian evaluations
                      * info[9]= # linear systems solved, i.e. # attempts for reducing error
                      */

T *particle1, T *particle2,
Vector3r position1, Vector3r position2,
Vector3r &p1,Vector3r &p2,
//Matrix3r RotationMat1,Matrix3r RotationMat2,
Matrix3r globalContact, Vector3r &center, bool &touching
)
{
	register int k;
	touching = true;
	Vector3r e, p_f(0.,0.,0.), pDp_f(0.,0.,0.),deltaF,hx(0.,0.,0.);
	Vector2r jacTf, diag_jacTjac(0,0);
	Mat32r jac;
	Matrix2d jacTjac;
	Vector2r Dp;
	Vector2r pDp;
	bool stop_flag(false);

	register double mu(0.),  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
	double jacTf_inf(0.); /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
	double p_L2, Dp_L2=LM_REAL_MAX, dF, dL;
	double p_fL2(0.), pDp_fL2(0.),pDp_feL2(0.),pDp_L2=1.0; //zhswee 2-norm of f(x)
	double tau, eps1, eps2, eps2_sq, eps3, init_p_fL2;
	int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;



////////////////////check input parameters//////////////
	if(opts){
		tau=opts[0];
		eps1=opts[1];
	  	eps2=opts[2];
	  	eps2_sq=opts[2]*opts[2];
      	eps3=opts[3];
  	}
  	else{ // use default values
	  	tau=LM_CNST(LM_INIT_MU);
	  	eps1=LM_CNST(LM_STOP_THRESH);
	  	eps2=LM_CNST(LM_STOP_THRESH);
	  	eps2_sq=LM_CNST(LM_STOP_THRESH)*LM_CNST(LM_STOP_THRESH);
      	eps3=LM_CNST(LM_STOP_THRESH);
  	}


////////////////////////////////////////////////////////////
  /* compute f(p) and its L2 norm */
  	(*func)(p, center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev=1;
        p_f = p2 - p1;//f(p) = p2 - p1, and
	//cout<<"pf="<<p_f(0)<<"\t"<<p_f(1)<<"\t"<<p_f(2)<<endl;
	p_fL2=p_f.norm();


  	init_p_fL2=p_fL2;
  	if(!LM_FINITE(p_fL2)) stop=7;
    //compute an approximation of Jacobian at p'
    //Vector2r p_1 = p + Vector2r(0.001,0.001);
    /*Vector3r p_f1,p_f2;
    (*func)(Vector2r(p(0)+0.001,p(1)), center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev++;
    p_f1 = p2 -p1;
    p_f1 = p_f1 - p_f;
    (*func)(Vector2r(p(0),p(1)+0.001), center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev++;
    p_f2 = p2 -p1;
    p_f2 = p_f2 - p_f;
    //Compute the Jacobian J at p,  J^T J,  J^T f,  ||J^T f||_inf and ||p||^2.
    //jac = (*jacf)(p, particle1, particle2, globalContact); ++njev;
    for(int i1=0; i1<3;i1++){
        jac(i1,0) = p_f1(i1)/0.001;
        jac(i1,1) = p_f2(i1)/0.001;
    }
    cout<<"p_f="<<p_f<<"p_f1="<<p_f1<<"jac="<<jac<<endl;
    */
  	for(k=0; k<itmax && !stop; ++k){
    /* Note that p and p_f have been updated at a previous iteration */

    	if(p_fL2<=eps3){ /* the distance is pretty small */
      		stop=6;
      		break;
    	}

		if (stop_flag){//stop the rutine
			stop=9;
			touching = false;
			break;
		}

    	//Compute the Jacobian J at p,  J^T J,  J^T f,  ||J^T f||_inf and ||p||^2.
    	jac = (*jacf)(p, particle1, particle2, globalContact); ++njev;
		//cout<<"Jacobian mat="<<jac<<endl;
    	/* J^T J, J^T f */
	 	jacTjac = jac.transpose()*jac;
	 	jacTf = jac.transpose()*(-1.0*p_f);//jac.transpose()*p_f;
 		//cout<<"jacTf="<<jacTf<<endl;
		//cout<<"p2"<<p<<endl;
	  	/* Compute ||J^T f||_inf and ||p||^2 */
		jacTf_inf = (fabs(jacTf(0)) > fabs(jacTf(1)) ? fabs(jacTf(0)) : fabs(jacTf(1)));
		p_L2 = p.norm();
		//cout<<"jacTf_inf"<<jacTf_inf<<endl;
    	/* check for convergence */
    	if((jacTf_inf <= eps1)){
      		Dp_L2=0.0; /* no increment for p in this case */
      		stop=1;
      		break;
    	}
		deltaF = pDp_f - p_f;
		pDp_feL2 = fabs(pDp_f.norm() - p_f.norm());
		/*if(pDp_feL2 <=eps2*p_fL2)
		{
			stop=8;
			break;
		}*/
   		/* compute initial damping factor */
    	if(k==0){
			tmp = (jacTjac(0,0) > jacTjac(1,1) ? jacTjac(0,0) : jacTjac(1,1));//find max diagonal element
      		mu=tau*tmp;
    	}
		diag_jacTjac(0) = jacTjac(0,0);diag_jacTjac(1) = jacTjac(1,1);
    	/* determine increment using adaptive damping */
   	 	while(1){
      		/* augment normal equations */
	  		//jacTjac += mu*Matrix2d::Identity();
			jacTjac(0,0) +=mu;
			jacTjac(1,1) +=mu;
      		/* solve augmented equations */
			//cout<<"jacTjac= "<<jacTjac<<endl;
	  		//Dp = jacTjac.inverse()*jacTf;//jacTjac.inverse()*(-jacTf);//caution:sigular
            double n_a = jacTjac(0,0)*jacTjac(1,1) - jacTjac(0,1)*jacTjac(0,1);
            Dp(0) = (jacTf(0)*jacTjac(1,1) - jacTf(1)*jacTjac(0,1))/n_a;
            Dp(1) = (jacTf(1)*jacTjac(0,0) - jacTf(0)*jacTjac(0,1))/n_a;
 			//cout<<"Dp="<<Dp(0)<<"\t"<<Dp(1)<<endl;

			Dp_L2 = Dp.norm();

        	//cout<<"pDp="<<pDp(0)<<"\t"<<pDp(1)<<endl;

			//std::cout<<"solved --- pDp-"<<pDp[0]<<"\t"<<pDp[1]<<std::endl;
        	//if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        	/*if(Dp_L2<=eps2*p_L2){ //relative change in p is small, stop
          		stop=2;
		   	//	printf("Dp_L2= %.9g \n", Dp_L2);
          		break;
        	}*/

            //test starts
            /*contact = globalContact*contact;	   //to global
	        d = pt2 - pt1;
	        check_dc = acos(d.dot(contact)/(d.norm()*contact.norm()));//phi in [0, pi]
		    if (  check_dc > 2.8 ){//
		            info[0] = nfunc;	//number of function invoking
			    info[1] = y_tmp;		//the shortest distance
		            return p_tmp;

		    }*/
            //test ends
            pDp = p + Dp;
            cout<<"alpha="<<pDp<<endl;
            pDp_L2 = pDp.norm();
       		if(Dp_L2>=(p_L2+eps2)/(LM_CNST(EPSILON)*LM_CNST(EPSILON))){ /* almost singular */
       			//if(Dp_L2>=(p_L2+eps2)/LM_CNST(EPSILON)){ /* almost singular */
         		stop=4;
         		break;
       		}

        	 (*func)(pDp, center, particle1, particle2,position1, position2, p1, p2,globalContact, stop_flag); ++nfev; /* evaluate function at p + Dp */
                 pDp_f = p2 - p1;
			//cout<<"pDp_f="<<pDp_f(0)<<"\t"<<pDp_f(1)<<"\t"<<pDp_f(2)<<endl;

			pDp_fL2 = pDp_f.norm();
			//cout<<"dis p="<<p_fL2<<endl;
			//cout<<"dis pDp="<<pDp_fL2<<endl;

        	if(!LM_FINITE(pDp_fL2)){ /* sum of squares is not finite, most probably due to a user error.
                                  * This check makes sure that the inner loop does not run indefinitely.
                                  * Thanks to Steve Danauskas for reporting such cases
                                  */
          		stop=7;
          		break;
        	}


			dL = Dp.transpose()*(mu*Dp + jacTf);//0.5*Dp.transpose()*(mu*Dp - jacTf);
        	dF=p_fL2 - pDp_fL2;
			//cout<<"dF/dL:"<<dF/dL<<endl;
        	/*if(dL>0.0 && dF>0.0){ // reduction in error, increment is accepted
          		tmp=(LM_CNST(2.0)*dF/dL-LM_CNST(1.0));
          		tmp=LM_CNST(1.0)-tmp*tmp*tmp;
		  		std::cout<<"tmp:"<<tmp<<std::endl;
          		mu=mu*( (tmp>=LM_CNST(ONE_THIRD))? tmp : LM_CNST(ONE_THIRD) );
          		nu=2;
				p = pDp;
				//e = hx;
		  		p_f = pDp_f;

		  		p_fL2=pDp_fL2;
          		break;
        	}*/
    		tmp = dF/dL;
            cout<<"damping factor = "<< mu<<" tho= "<<tmp<<" ht"<<Dp.transpose()*Dp<<"hh"<<Dp.transpose()*jacTf<<endl;
			if(tmp>0.0){//dL is greater than zero.
				tmp = 1.0 - pow(2*tmp-1., 3.);
				//mu/=4;
				mu=mu*( (tmp>=LM_CNST(ONE_THIRD))? tmp : LM_CNST(ONE_THIRD) );
				nu=2;

				p = pDp;
				//e = hx;
		  		p_f = pDp_f;

		  		p_fL2=pDp_fL2;

            //Jacobian
            //(*func)(p, center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev=1;
            /*(*func)(Vector2r(p(0)+Dp(0),p(1)), center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev++;
        p_f1 = p2 -p1;
        p_f1 = p_f1 - p_f;
        (*func)(Vector2r(p(0),p(1)+Dp(1)), center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev++;
        p_f2 = p2 -p1;
        p_f2 = p_f2 - p_f;
        //Compute the Jacobian J at p,  J^T J,  J^T f,  ||J^T f||_inf and ||p||^2.
        //jac = (*jacf)(p, particle1, particle2, globalContact); ++njev;
        for(int i1=0; i1<3;i1++){
            jac(i1,0) = p_f1(i1)/Dp(0);
            jac(i1,1) = p_f2(i1)/Dp(1);
        }*/



          		break;
			}
      /* if this point is reached, either the linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      		mu*=nu;
      		nu2=nu<<1; // 2*nu;
      		if(nu2<=nu){ // nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case
        		stop=5;
        		break;
      		}
      		nu=nu2;
			//mu *=10.;
			jacTjac(0,0) = diag_jacTjac(0);jacTjac(1,1) = diag_jacTjac(1);//restore diagonal J^T J entries

    	} /* inner loop */
		//cout<<"mu:"<<mu<<endl;

  }
  if(k>=itmax) stop=3;

  if(info){
    info[0]=init_p_fL2;
    info[1]=p_fL2;
    info[2]=jacTf_inf;
    info[3]=Dp_L2;
    //for(i=0, tmp=LM_REAL_MIN; i<m; ++i)
      //if(tmp<jacTjac[i*m+i]) tmp=jacTjac[i*m+i];
    info[4]=p_fL2;
    info[5]=(double)k;
    info[6]=(double)stop;
    info[7]=(double)nfev;
    info[8]=(double)njev;
    info[9]=(double)nlss;
  }

  return;//(stop!=4 && stop!=7)?  k : LM_ERROR;
}

/*****************************************************/
template <typename T>
void Disfunction(Vector2r para,Vector3r &center, T *particle1, T *particle2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag)
{
	Vector3r d;
	//stop_flag = false;
	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
    Vector2r phi(0.,0.);
    //contact direction c is defined at the global coordinate system
	//contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	//contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	//contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.
    contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));
	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	//contact = globalContact*contact;	   //to global

	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi(contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1 + position1;

	//point 2 on Particle 2
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	phi = particle2->Normal2Phi( (-1)*contact2);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
    //cout<<"phi of particle 2 = "<<phi<<endl;
        //get the coordinates of the current point in the global system
	p2 = particle2->getSurface( phi );
    //cout<<"point of par 2= "<<p2<<endl;
	p2 = particle2->rot_mat2global*p2 + position2;
    /*if(isnan(p1(0)) || isnan(p2(0))){
        //cout<<"NAN ,p1"<<p1<<"p2="<<p2<<endl;//scene->stopAtIter = scene->iter + 1;
        FILE * fin = fopen("direction.dat","a");
		fprintf(fin,"p1\t%e\t%e\t%e\n",p1[0],p1[1],p1[2]);
		fprintf(fin,"p2\t%e\t%e\t%e\n",p2[0],p2[1],p2[2]);
        fprintf(fin,"para\t%e\t%e\n",para[0],para[1]);
        fprintf(fin,"contact\t%e\t%e\t%e\n",contact[0],contact[1],contact[2]);
		fclose(fin);
    }*/
	d = p2 -p1;
    #ifdef _DEBUG_LM_NM
    FILE * fin1 = fopen("para.dat","a");

    fprintf(fin1,"\t%e\t%e\n",para[0],para[1]);
	fclose(fin1);
    #endif
	center = 0.5*(p1 + p2);
	//we should to check whether the angle between n1 (contact) and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (contact.dot(d) > 0)  //the angle < pi/2.
	{
		stop_flag = true;
		return ;//FIXME:May not be zero
	}

	return ;//CAUTION:d is the vector from point 1 on the particle1 to point 2 on the particle2.

}
template <typename T>
Mat32r Jacfunction(Vector2r para, T *particle1, T *particle2)
{



	Mat32r JacMat1, JacMat2;

	Vector3r p(0.,0.,0.), contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
        Vector2r phi;
	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));

	//contact = globalContact*contact;	   //to global
	//Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	phi = particle1->Normal2Phi( contact1);
    //cout<<"phi1= "<<phi<<endl;
    //cout<<"local rot= "<<particle1->rot_mat2local<<endl;
    //cout<<"global rot= "<<particle1->rot_mat2global<<endl;
    //cout<<"local*global= "<<particle1->rot_mat2local*particle1->rot_mat2global<<endl;
	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	Mat32r d_p2phi = particle1->rot_mat2global*particle1->Derivate_P2Phi(phi);//to global system
    //test P2Phi
    /*Vector3r p10 = particle1->getSurface(Vector2r(phi(0),phi(1)));
    Vector3r p11 = particle1->getSurface(Vector2r(phi(0)+1e-5,phi(1)));
    Vector3r p12 = particle1->getSurface(Vector2r(phi(0),phi(1)+1e-5));
    Mat32r P12Phi_d;
    for(int i=0;i<3;i++){
        P12Phi_d(i,0) = (p11(i)-p10(i))/1e-5;
    }
    for(int i=0;i<3;i++){
        P12Phi_d(i,1) = (p12(i)-p10(i))/1e-5;
    }
    cout<<"Approx. d_p12phi="<<particle1->rot_mat2global*P12Phi_d<<endl;
    cout<<"d_p12phi="<<d_p2phi<<endl;
    */
	Mat23r d_phi2n = particle1->Derivate_Phi2N(phi, contact1);
    //cout<<"d_phi2n="<<d_phi2n<<endl;
    Mat32r c2alpha;
	c2alpha<<
			-sin(para(0))*cos(para(1)), -cos(para(0))*sin(para(1)),
			cos(para(0))*cos(para(1)), -sin(para(0))*sin(para(1)),
			0., cos(para(1));
	Mat32r c2alpha1 = particle1->rot_mat2local*c2alpha;
    //cout<<"c2alpha ="<<c2alpha<<endl;
	JacMat1 = d_p2phi*d_phi2n*c2alpha1;
    //test JacMat1
    /*
    Vector3r pa0 = SurfaceP(Vector2r(para(0),para(1)),particle1);
    Vector3r pa1 = SurfaceP(Vector2r(para(0)+1e-7,para(1)),particle1);
    Vector3r pa2 = SurfaceP(Vector2r(para(0),para(1)+1e-7),particle1);
    Mat32r JacMat1_d;
    for(int i=0;i<3;i++){
        JacMat1_d(i,0) = (pa1(i)-pa0(i))/1e-7;
    }
    for(int i=0;i<3;i++){
        JacMat1_d(i,1) = (pa2(i)-pa0(i))/1e-7;
    }
    cout<<"Jacobian mat1="<<JacMat1<<endl;
    cout<<"Approx. Jacobian mat1="<<JacMat1_d<<endl;


    */
	//
	contact2 = (-1.0)*particle2->rot_mat2local*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	phi = particle2->Normal2Phi(contact2);
    //cout<<"phi2= "<<phi<<endl;
	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	d_p2phi = particle2->rot_mat2global*particle2->Derivate_P2Phi(phi);//to global system
    //test P2Phi
    /*Vector3r p0 = particle2->getSurface(Vector2r(phi(0),phi(1)));
    Vector3r p1 = particle2->getSurface(Vector2r(phi(0)+1e-5,phi(1)));
    Vector3r p2 = particle2->getSurface(Vector2r(phi(0),phi(1)+1e-5));
    Mat32r P2Phi_d;
    for(int i=0;i<3;i++){
        P2Phi_d(i,0) = (p1(i)-p0(i))/1e-5;
    }
    for(int i=0;i<3;i++){
        P2Phi_d(i,1) = (p2(i)-p0(i))/1e-5;
    }
    cout<<"Approx. d_p2phi="<<particle2->rot_mat2global*P2Phi_d<<endl;
    cout<<"d_p2phi="<<d_p2phi<<endl;
    */
	d_phi2n = particle2->Derivate_Phi2N(phi, contact2);
    //test Phi2N
    /*
    Vector2r pn0 = particle2->Normal2Phi(Vector3r(contact2(0),contact2(1),contact2(2)));
    Vector2r pn1 = particle2->Normal2Phi(Vector3r(contact2(0)+1e-5,contact2(1),contact2(2)));
    Vector2r pn2 = particle2->Normal2Phi(Vector3r(contact2(0),contact2(1)+1e-5,contact2(2)));
    Vector2r pn3 = particle2->Normal2Phi(Vector3r(contact2(0),contact2(1),contact2(2)+1e-5));
    Mat23r Phi2N_d;
    for(int i=0;i<2;i++){
        Phi2N_d(i,0) = (pn1(i)-pn0(i))/1e-5;
    }
    for(int i=0;i<2;i++){
        Phi2N_d(i,1) = (pn2(i)-pn0(i))/1e-5;
    }
    for(int i=0;i<2;i++){
        Phi2N_d(i,2) = (pn3(i)-pn0(i))/1e-5;
    }
    cout<<"Approx. d_phi2n="<<Phi2N_d<<endl;
    cout<<"d_phi2n="<<d_phi2n<<endl;
    */
	c2alpha1 = (-1)*particle2->rot_mat2local*c2alpha;

	JacMat2 = d_p2phi*d_phi2n*c2alpha1;
    //test JacMat2
    /*
    Vector3r pa20 = SurfaceP2(Vector2r(para(0),para(1)),particle2);
    Vector3r pa21 = SurfaceP2(Vector2r(para(0)+1e-7,para(1)),particle2);
    Vector3r pa22 = SurfaceP2(Vector2r(para(0),para(1)+1e-7),particle2);
    Mat32r JacMat2_d;
    for(int i=0;i<3;i++){
        JacMat2_d(i,0) = (pa21(i)-pa20(i))/1e-7;
    }
    for(int i=0;i<3;i++){
        JacMat2_d(i,1) = (pa22(i)-pa20(i))/1e-7;
    }
    cout<<"Jacobian mat2="<<JacMat2<<endl;
    cout<<"Approx. Jacobian mat2="<<JacMat2_d<<endl;
    */
	return JacMat2 - JacMat1;
}

template <typename T>
Vector3r SurfaceP(Vector2r para,T *particle1)
{
	Vector3r d,p1;

	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
    Vector2r phi(0.,0.);
    //contact direction c is defined at the global coordinate system
    contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));

	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	//contact = globalContact*contact;	   //to global

	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi(contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1;

	return p1;
}

template <typename T>
Vector3r SurfaceP2(Vector2r para,T *particle1)
{
	Vector3r d,p1;

	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
    Vector2r phi(0.,0.);
    //contact direction c is defined at the global coordinate system
    contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));

	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	//contact = globalContact*contact;	   //to global

	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 =(-1.0)* particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi(contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1;

	return p1;
}


template <typename T>//simplized Jacbian matrix calculation
Mat32r JacfunctionSim(Vector2r para, T *particle1, T *particle2)
{



	Mat32r JacMat1, JacMat2;

	Vector3r p(0.,0.,0.), contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
        Vector2r phi;
    contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));

	//contact = globalContact*contact;	   //to global
	//Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	contact1 = contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	phi = particle1->Normal2Phi( contact1);
    //cout<<"phi1"<<phi<<endl;
	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	Mat32r d_p2phi = particle1->Derivate_P2Phi(phi);//to global system
    //test P2Phi
    /*Vector3r p10 = particle1->getSurface(Vector2r(phi(0),phi(1)));
    Vector3r p11 = particle1->getSurface(Vector2r(phi(0)+1e-5,phi(1)));
    Vector3r p12 = particle1->getSurface(Vector2r(phi(0),phi(1)+1e-5));
    Mat32r P12Phi_d;
    for(int i=0;i<3;i++){
        P12Phi_d(i,0) = (p11(i)-p10(i))/1e-5;
    }
    for(int i=0;i<3;i++){
        P12Phi_d(i,1) = (p12(i)-p10(i))/1e-5;
    }
    cout<<"Approx. d_p12phi="<<P12Phi_d<<endl;
    cout<<"d_p12phi="<<d_p2phi<<endl;
    */
    Mat23r d_phi2n = particle1->Derivate_Phi2N(phi, contact1);
    //test Phi2N
    /*Vector2r pn0 = particle1->Normal2Phi(Vector3r(contact1(0),contact1(1),contact1(2)));
    Vector2r pn1 = particle1->Normal2Phi(Vector3r(contact1(0)+1e-5,contact1(1),contact1(2)));
    Vector2r pn2 = particle1->Normal2Phi(Vector3r(contact1(0),contact1(1)+1e-5,contact1(2)));
    Vector2r pn3 = particle1->Normal2Phi(Vector3r(contact1(0),contact1(1),contact1(2)+1e-5));
    Mat23r Phi2N_d;
    for(int i=0;i<2;i++){
        Phi2N_d(i,0) = (pn1(i)-pn0(i))/1e-5;
    }
    for(int i=0;i<2;i++){
        Phi2N_d(i,1) = (pn2(i)-pn0(i))/1e-5;
    }
    for(int i=0;i<2;i++){
        Phi2N_d(i,2) = (pn3(i)-pn0(i))/1e-5;
    }
    cout<<"Approx. d_phi2n="<<Phi2N_d<<endl;
    cout<<"d_phi2n="<<d_phi2n<<endl;
    */

    Mat32r c2alpha;
	c2alpha<<
			-sin(para(0))*cos(para(1)), -cos(para(0))*sin(para(1)),
			cos(para(0))*cos(para(1)), -sin(para(0))*sin(para(1)),
			0., cos(para(1));
	Mat32r c2alpha1 = c2alpha;
    //cout<<"c2alpha ="<<c2alpha<<endl;
	JacMat1 = d_p2phi*d_phi2n*c2alpha1;
    //test JacMat1
    /*
    Vector3r pa0 = SurfaceP(Vector2r(para(0),para(1)),particle1);
    Vector3r pa1 = SurfaceP(Vector2r(para(0)+1e-7,para(1)),particle1);
    Vector3r pa2 = SurfaceP(Vector2r(para(0),para(1)+1e-7),particle1);
    Mat32r JacMat1_d;
    for(int i=0;i<3;i++){
        JacMat1_d(i,0) = (pa1(i)-pa0(i))/1e-7;
    }
    for(int i=0;i<3;i++){
        JacMat1_d(i,1) = (pa2(i)-pa0(i))/1e-7;
    }
    cout<<"Jacobian mat1="<<JacMat1<<endl;
    cout<<"Approx. Jacobian mat1="<<JacMat1_d<<endl;
    */
	//
	contact2 = -1.0*contact; //to particle's local system
	//contact.normalize();
	//std::cout<<"contact normal= "<<contact<<std::endl;
	phi = particle2->Normal2Phi( contact2);
    //if(phi(0)<0){phi(0)+=Mathr::TWO_PI;}
    //if(phi(1)<0){phi(1)+=Mathr::TWO_PI;}
    cout<<"phi2"<<phi<<endl;
	//get the coordinates of the current point in the global system
	//p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	d_p2phi = particle2->Derivate_P2Phi(phi);//to global system
    //test P2Phi
    /*Vector3r p0 = particle2->getSurface(Vector2r(phi(0),phi(1)));
    Vector3r p1 = particle2->getSurface(Vector2r(phi(0)+1e-5,phi(1)));
    Vector3r p2 = particle2->getSurface(Vector2r(phi(0),phi(1)+1e-5));
    Mat32r P2Phi_d;
    for(int i=0;i<3;i++){
        P2Phi_d(i,0) = (p1(i)-p0(i))/1e-5;
    }
    for(int i=0;i<3;i++){
        P2Phi_d(i,1) = (p2(i)-p0(i))/1e-5;
    }
    cout<<"Approx. d_p2phi="<<P2Phi_d<<endl;
    cout<<"d_p2phi="<<d_p2phi<<endl;
    */
	d_phi2n = particle2->Derivate_Phi2N(phi, contact2);
    //test Phi2N
    /*
    Vector2r pn0 = particle2->Normal2Phi(Vector3r(contact2(0),contact2(1),contact2(2)));
    Vector2r pn1 = particle2->Normal2Phi(Vector3r(contact2(0)+1e-5,contact2(1),contact2(2)));
    Vector2r pn2 = particle2->Normal2Phi(Vector3r(contact2(0),contact2(1)+1e-5,contact2(2)));
    Vector2r pn3 = particle2->Normal2Phi(Vector3r(contact2(0),contact2(1),contact2(2)+1e-5));
    Mat23r Phi2N_d;
    for(int i=0;i<2;i++){
        Phi2N_d(i,0) = (pn1(i)-pn0(i))/1e-5;
    }
    for(int i=0;i<2;i++){
        Phi2N_d(i,1) = (pn2(i)-pn0(i))/1e-5;
    }
    for(int i=0;i<2;i++){
        Phi2N_d(i,2) = (pn3(i)-pn0(i))/1e-5;
    }
    cout<<"Approx. d_phi2n="<<Phi2N_d<<endl;
    cout<<"d_phi2n="<<d_phi2n<<endl;
    */
	c2alpha1 = (-1)*c2alpha;

	JacMat2 = d_p2phi*d_phi2n*c2alpha1;
    //test JacMat2
    /*
    Vector3r pa0 = SurfaceP2(Vector2r(para(0),para(1)),particle2);
    Vector3r pa1 = SurfaceP2(Vector2r(para(0)+1e-7,para(1)),particle2);
    Vector3r pa2 = SurfaceP2(Vector2r(para(0),para(1)+1e-7),particle2);
    Mat32r JacMat2_d;
    for(int i=0;i<3;i++){
        JacMat2_d(i,0) = (pa1(i)-pa0(i))/1e-7;
    }
    for(int i=0;i<3;i++){
        JacMat2_d(i,1) = (pa2(i)-pa0(i))/1e-7;
    }
    cout<<"Jacobian mat2="<<JacMat2<<endl;
    cout<<"Approx. Jacobian mat2="<<JacMat2_d<<endl;
    */
	return JacMat2 - JacMat1;
}
//secent Jacobian matrix with approximation
template <typename T>
Mat32r JacfunctionApp(Vector2r para, Vector2r paraDelta, Vector3r &center, T *particle1, T *particle2,Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag)
{
    Vector3r p_f0,p_f1,p_f2;
    Vector2r para0;
    Disfunction(para,center,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f0 = p2 - p1;//f(p) = p2 - p1, and
    para0 = Vector2r(para(0)+paraDelta(0),para(1));
    Disfunction(para0,center,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f1 = p2 - p1;
    para0 = Vector2r(para(0),para(1)+paraDelta(1));
    Disfunction(para0,center,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f2 = p2 - p1;
    Mat32r JacMat;
    for(int i=0;i<3;i++){
        JacMat(i,0) = (p_f1(i)-p_f0(i))/paraDelta(0);
    }
    for(int i=0;i<3;i++){
        JacMat(i,1) = (p_f2(i)-p_f0(i))/paraDelta(1);
    }

	return JacMat;
}

template <typename T>	//return distance
void lmder2(
  void (*func)(Vector2r p, Vector3r &center, T *particle1, T *particle2,
		   Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  Mat32r (*jacf)(Vector2r p, T *particle1, T *particle2),  /* function to evaluate the Jacobian \part x / \part p */
  Vector2r &p,         /* I/O: initial parameter estimates. On output has the estimated solution */
  int itmax,          /* I: maximum number of iterations */
  double opts[4],    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
                       * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used
                       */
  double info[LM_INFO_SZ],//,
					           /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]= ||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small Dp
                      *                                 3 - stopped by itmax
                      *                                 4 - singular matrix. Restart from current p with increased mu
                      *                                 5 - no further error reduction is possible. Restart with increased mu
                      *                                 6 - stopped by small ||e||_2
                      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                      * info[7]= # function evaluations
                      * info[8]= # Jacobian evaluations
                      * info[9]= # linear systems solved, i.e. # attempts for reducing error
                      */

T *particle1, T *particle2,
Vector3r position1, Vector3r position2,
Vector3r &p1,Vector3r &p2, Vector3r &center, bool &touching, Vector3r &contact
)
{
	register int k;
	touching = true;
	Vector3r e, p_f(0.,0.,0.), pDp_f(0.,0.,0.),deltaF,hx(0.,0.,0.);
	Vector2r jacTf, diag_jacTjac(0,0);
	Mat32r jac;
	Matrix2d jacTjac;
	Vector2r Dp(1e-7,1e-7);
	Vector2r pDp;
	bool stop_flag(false);

	register double mu(0.),  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
	double jacTf_inf(0.); /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
	double p_L2, Dp_L2=LM_REAL_MAX, dF, dL;
	double p_fL2(0.), pDp_fL2(0.),pDp_feL2(0.),pDp_L2=1.0; //zhswee 2-norm of f(x)
	double tau, eps1, eps2, eps2_sq, eps3, init_p_fL2;
	int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;



////////////////////check input parameters//////////////
	if(opts){
		tau=opts[0];
		eps1=opts[1];
	  	eps2=opts[2];
	  	eps2_sq=opts[2]*opts[2];
      	eps3=opts[3];
  	}
  	else{ // use default values
	  	tau=LM_CNST(LM_INIT_MU);
	  	eps1=LM_CNST(LM_STOP_THRESH);
	  	eps2=LM_CNST(LM_STOP_THRESH);
	  	eps2_sq=LM_CNST(LM_STOP_THRESH)*LM_CNST(LM_STOP_THRESH);
      	eps3=LM_CNST(LM_STOP_THRESH);
  	}


////////////////////////////////////////////////////////////
  /* compute f(p) and its L2 norm */
  	(*func)(p, center, particle1, particle2,position1, position2,p1,p2, stop_flag); nfev=1;
        p_f = p2 - p1;//f(p) = p2 - p1, and
	//cout<<"pf="<<p_f(0)<<"\t"<<p_f(1)<<"\t"<<p_f(2)<<endl;
    //cout<<"para="<<p(0)<<"\t"<<p(1)<<"\t"<<p(2)<<endl;
	p_fL2=p_f.norm();


  	init_p_fL2=p_fL2;
  	if(!LM_FINITE(p_fL2)){ stop=7;}

  	for(k=0; k<itmax && !stop; ++k){
    /* Note that p and p_f have been updated at a previous iteration */

    	if(p_fL2<=eps3){ /* the distance is pretty small */
      		stop=6;
      		break;
    	}
        //break;
		if (stop_flag){//stop the rutine
			stop=9;
			touching = false;
			break;
		}
        {//check angle between contact direction and depth vector
        //Vector3r contact(0.,0.,0.);

        contact(0) = cos(p(0))*cos(p(1));  //the base vectors are e_i, i=1,2,3
        contact(1) = sin(p(0))*cos(p(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
        contact(2) = sin(p(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.
        //double angle = acos(p_f.dot(contact)/(p_f.norm()*contact.norm()));//phi in [0, pi]
        double angle = p_f.dot(contact)/(p_f.norm()*contact.norm());//phi in [0, pi]
        #ifdef _DEBUG_LM_NM
        cout<<"check whether d is along c:--angle="<<angle<<"para="<<p(0)<<" "<<p(1)<<endl;
        #endif
        if(fabs(angle) > 0.95)
        {stop = 8;
        //cout<<"angle2="<<DisContact(p, center, particle1, particle2,position1, position2,p1,p2, stop_flag)<<endl;
        break;}

        }
    	//Compute the Jacobian J at p,  J^T J,  J^T f,  ||J^T f||_inf and ||p||^2.
    	//jac = (*jacf)(p, particle1, particle2); ++njev;
        Mat32r jac0,jac1;
        jac0 = JacfunctionApp(p, Dp,center, particle1, particle2,position1, position2,p1,p2, stop_flag);nfev+=3;
        //jac0 = JacfunctionApp(p, Vector2r(1e-7,1e-7),center, particle1, particle2,position1, position2,p1,p2, stop_flag);nfev+=3;

        //jac1 = JacfunctionSim(p, particle1, particle2);
		//cout<<"Jacobian mat="<<jac<<endl;
        //cout<<"Approx. Jacobian mat="<<jac0<<endl;
        jac = jac0;
        //cout<<"Jacobian mat simplized ="<<jac1<<endl;
    	/* J^T J, J^T f */
	 	jacTjac = jac.transpose()*jac;
	 	jacTf = jac.transpose()*(-1.0*p_f);//jac.transpose()*p_f;
 		//cout<<"jacTf="<<jacTf<<endl;
		//cout<<"p2"<<p<<endl;
	  	/* Compute ||J^T f||_inf and ||p||^2 */
		jacTf_inf = (fabs(jacTf(0)) > fabs(jacTf(1)) ? fabs(jacTf(0)) : fabs(jacTf(1)));
		p_L2 = p.norm();
		//cout<<"jacTf_inf"<<jacTf_inf<<endl;
    	/* check for convergence */
    	if((jacTf_inf <= eps1)){
      		Dp_L2=0.0; /* no increment for p in this case */
      		stop=1;
      		break;
    	}
		deltaF = pDp_f - p_f;
		pDp_feL2 = fabs(pDp_f.norm() - p_f.norm());
		/*if(pDp_feL2 <=eps2*p_fL2)
		{
			stop=8;
			break;
		}*/
   		/* compute initial damping factor */
    	if(k==0){
			tmp = (jacTjac(0,0) > jacTjac(1,1) ? jacTjac(0,0) : jacTjac(1,1));//find max diagonal element
      		mu=tau*tmp;
    	}
		diag_jacTjac(0) = jacTjac(0,0);diag_jacTjac(1) = jacTjac(1,1);
    	/* determine increment using adaptive damping */
   	 	while(1){
      		/* augment normal equations */
	  		//jacTjac += mu*Matrix2d::Identity();
			jacTjac(0,0) +=mu;
			jacTjac(1,1) +=mu;
      		/* solve augmented equations */
			//cout<<"jacTjac= "<<jacTjac<<endl;
	  		//Dp = jacTjac.inverse()*jacTf;//jacTjac.inverse()*(-jacTf);//caution:sigular
            double n_a = jacTjac(0,0)*jacTjac(1,1) - jacTjac(0,1)*jacTjac(0,1);
            Dp(0) = (jacTf(0)*jacTjac(1,1) - jacTf(1)*jacTjac(0,1))/n_a;
            Dp(1) = (jacTf(1)*jacTjac(0,0) - jacTf(0)*jacTjac(0,1))/n_a;
 			//cout<<"Dp="<<Dp(0)<<"\t"<<Dp(1)<<endl;

			Dp_L2 = Dp.norm();

        	//cout<<"pDp="<<pDp(0)<<"\t"<<pDp(1)<<endl;

			//std::cout<<"solved --- pDp-"<<pDp[0]<<"\t"<<pDp[1]<<std::endl;
        	//if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        	if(Dp_L2 < 1e-5){//if(Dp_L2<=eps2*p_L2){ //relative change in p is small, stop
          		stop=2;
		   	//	printf("Dp_L2= %.9g \n", Dp_L2);
          		break;
        	}

            //test starts
            /*contact = globalContact*contact;	   //to global
	        d = pt2 - pt1;
	        check_dc = acos(d.dot(contact)/(d.norm()*contact.norm()));//phi in [0, pi]
		    if (  check_dc > 2.8 ){//
		            info[0] = nfunc;	//number of function invoking
			    info[1] = y_tmp;		//the shortest distance
		            return p_tmp;

		    }*/
            //test ends
            pDp = p + Dp;
            //cout<<"para="<<pDp<<endl;
            pDp_L2 = pDp.norm();
       		if(Dp_L2>=(p_L2+eps2)/(LM_CNST(EPSILON)*LM_CNST(EPSILON))){ /* almost singular */
       			//if(Dp_L2>=(p_L2+eps2)/LM_CNST(EPSILON)){ /* almost singular */
         		stop=4;
         		break;
       		}

        	 (*func)(pDp, center, particle1, particle2,position1, position2, p1, p2, stop_flag); ++nfev; /* evaluate function at p + Dp */
                 pDp_f = p2 - p1;
			//cout<<"pDp_f="<<pDp_f(0)<<"\t"<<pDp_f(1)<<"\t"<<pDp_f(2)<<endl;

			pDp_fL2 = pDp_f.norm();
			//cout<<"dis p="<<p_fL2<<endl;
			//cout<<"dis pDp="<<pDp_fL2<<endl;

        	if(!LM_FINITE(pDp_fL2)){ /* sum of squares is not finite, most probably due to a user error.
                                  * This check makes sure that the inner loop does not run indefinitely.
                                  * Thanks to Steve Danauskas for reporting such cases
                                  */
          		stop=7;
          		break;
        	}


			dL = Dp.transpose()*(mu*Dp + jacTf);//0.5*Dp.transpose()*(mu*Dp - jacTf);
        	dF=p_fL2 - pDp_fL2;
			//cout<<"dF/dL:"<<dF/dL<<endl;
        	/*if(dL>0.0 && dF>0.0){ // reduction in error, increment is accepted
          		tmp=(LM_CNST(2.0)*dF/dL-LM_CNST(1.0));
          		tmp=LM_CNST(1.0)-tmp*tmp*tmp;
		  		std::cout<<"tmp:"<<tmp<<std::endl;
          		mu=mu*( (tmp>=LM_CNST(ONE_THIRD))? tmp : LM_CNST(ONE_THIRD) );
          		nu=2;
				p = pDp;
				//e = hx;
		  		p_f = pDp_f;

		  		p_fL2=pDp_fL2;
          		break;
        	}*/
    		tmp = dF/dL;
            //cout<<"damping factor = "<< mu<<" tho= "<<tmp<<" ht"<<Dp.transpose()*Dp<<"hh"<<Dp.transpose()*jacTf<<endl;
			if(tmp>0.0){//dL is greater than zero.
				tmp = 1.0 - pow(2*tmp-1., 3.);
				//mu/=4;
				mu=mu*( (tmp>=LM_CNST(ONE_THIRD))? tmp : LM_CNST(ONE_THIRD) );
				nu=2;

				p = pDp;
				//e = hx;
		  		p_f = pDp_f;

		  		p_fL2=pDp_fL2;

            //Jacobian
            //(*func)(p, center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev=1;
            /*(*func)(Vector2r(p(0)+Dp(0),p(1)), center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev++;
        p_f1 = p2 -p1;
        p_f1 = p_f1 - p_f;
        (*func)(Vector2r(p(0),p(1)+Dp(1)), center, particle1, particle2,position1, position2,p1,p2,globalContact, stop_flag); nfev++;
        p_f2 = p2 -p1;
        p_f2 = p_f2 - p_f;
        //Compute the Jacobian J at p,  J^T J,  J^T f,  ||J^T f||_inf and ||p||^2.
        //jac = (*jacf)(p, particle1, particle2, globalContact); ++njev;
        for(int i1=0; i1<3;i1++){
            jac(i1,0) = p_f1(i1)/Dp(0);
            jac(i1,1) = p_f2(i1)/Dp(1);
        }*/



          		break;
			}
      /* if this point is reached, either the linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      		mu*=nu;
      		nu2=nu<<1; // 2*nu;
      		if(nu2<=nu){ // nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case
        		stop=5;
        		break;
      		}
      		nu=nu2;
			//mu *=10.;
			jacTjac(0,0) = diag_jacTjac(0);jacTjac(1,1) = diag_jacTjac(1);//restore diagonal J^T J entries

    	} /* inner loop */
		//cout<<"mu:"<<mu<<endl;

  }
  if(k>=itmax) stop=3;

  if(info){
    info[0]=init_p_fL2;
    info[1]=p_fL2;
    info[2]=jacTf_inf;
    info[3]=Dp_L2;
    //for(i=0, tmp=LM_REAL_MIN; i<m; ++i)
      //if(tmp<jacTjac[i*m+i]) tmp=jacTjac[i*m+i];
    info[4]=p_fL2;
    info[5]=(double)k;
    info[6]=(double)stop;
    info[7]=(double)nfev;
    info[8]=(double)njev;
    info[9]=(double)nlss;
  }

  return;//(stop!=4 && stop!=7)?  k : LM_ERROR;
}

template <typename T>
double DisContact(Vector2r para,Vector3r &center, T *particle1, T *particle2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag)
{
	Vector3r d;
	//stop_flag = false;
	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
    Vector2r phi(0.,0.);
    //contact direction c is defined at the global coordinate system
	//contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	//contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	//contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.
    contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));



	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi(contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1 + position1;

	//point 2 on Particle 2
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	phi = particle2->Normal2Phi( (-1)*contact2);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
    //cout<<"phi of particle 2 = "<<phi<<endl;
        //get the coordinates of the current point in the global system
	p2 = particle2->getSurface( phi );
    //cout<<"point of par 2= "<<p2<<endl;
	p2 = particle2->rot_mat2global*p2 + position2;
    if(isnan(p1(0)) || isnan(p2(0))){
        //cout<<"NAN ,p1"<<p1<<"p2="<<p2<<endl;//scene->stopAtIter = scene->iter + 1;
        FILE * fin = fopen("direction.dat","a");
		fprintf(fin,"p1\t%e\t%e\t%e\n",p1[0],p1[1],p1[2]);
		fprintf(fin,"p2\t%e\t%e\t%e\n",p2[0],p2[1],p2[2]);
        fprintf(fin,"para\t%e\t%e\n",para[0],para[1]);
        fprintf(fin,"contact\t%e\t%e\t%e\n",contact[0],contact[1],contact[2]);
		fclose(fin);
    }
	d = p2 -p1;
    FILE * fin1 = fopen("para.dat","a");

    fprintf(fin1,"\t%e\t%e\n",para[0],para[1]);
	fclose(fin1);
	center = 0.5*(p1 + p2);
	//we should to check whether the angle between n1 (contact) and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (contact.dot(d) > 0)  //the angle < pi/2.
	{
		stop_flag = true;
		return 1;//FIXME:May not be zero
	}
    double angle = d.dot(contact)/(d.norm()*contact.norm());//phi in [0, pi]
    //cout<<"check whether d is along c:--angle="<<angle<<endl;
	return angle;//minimia is -1.

}



///////////////////////nelder-mead simplex method////////////////////////////////
	//Extrapolates by a factor fac through the face of the simplex across from the
	//highest point, tries it, and replaces the high point if the new point is better.
	//double func(Vector2r, Vector3r &, T*, T*, Matrix3r,Matrix3r,Vector3r, Vector3r, Vector3r &,Vector3r &,Matrix3r,bool &),
	//      Vector3r &center, T *p1, T *p2,
	//      Matrix3r RotationMat1,Matrix3r RotationMat2,
	//      Vector3r position1, Vector3r position2,Vector3r &pt1,Vector3r &pt2,
	//      Matrix3r globalContact,bool &stop_flag)
template <class T>
double amotry(Vector2r (&p)[3], Vector3r &y, Vector2r &psum,
		const int ihi, const double fac,
	      double func(Vector2r, Vector3r &, T*, T*,Vector3r, Vector3r, Vector3r &,Vector3r &,Matrix3r,bool &),
	      Vector3r &center, T *p1, T *p2,
	      //Matrix3r RotationMat1,Matrix3r RotationMat2,//not used
	      Vector3r position1, Vector3r position2,Vector3r &pt1,Vector3r &pt2,
	      Matrix3r globalContact,bool &stop_flag)
{
	Vector2r ptry;
	double fac1=(1.0-fac)/2.0;
	double fac2=fac1-fac;
	ptry = psum*fac1 - p[ihi]*fac2;

	double ytry=func(ptry, center, p1, p2,position1, position2, pt1, pt2, globalContact,stop_flag); //Evaluate the function at the trial point.
	if (ytry < y(ihi)) {	//If it's better than the highest, then replace the highest.
		y(ihi)=ytry;
		psum += ptry - p[ihi];
		p[ihi] = ptry;
	}
	return ytry;
}
template <typename T>//The base class of T is Superquadrics (defined in Superquardrics.h)
Vector2r neldermead(
  	double (*func)(Vector2r p,Vector3r &center, T *particle1, T *particle2,
		       //Matrix3r RotationMat1,Matrix3r RotationMat2,//not used
		       Vector3r position1, Vector3r position2,Vector3r &pt1,Vector3r &pt2,
		       Matrix3r globalContact,bool &stop_flag), //function
	Vector2r point,
	Vector2r dels,
	double info[2], //
	T *p1,
	T *p2,
	//Matrix3r RotationMat1,Matrix3r RotationMat2,//not used
	Vector3r position1, Vector3r position2,Vector3r &pt1,Vector3r &pt2,
	Matrix3r globalContact, Vector3r &center,
	bool &touching
)
{
	const double ftol =1e-3; //the fractional convergence tolerance to be achieved in the function
	int nfunc;		//the number of function evaluations
    bool stop_flag(false);
	double fmin;	//function value at the minimum
	const int NMAX=5000;		//maximum allowed number of function evaluations
	const double TINY=1.0e-10;
	int ihi,ilo,inhi;
	Vector2r psum, pmin;
	Vector3r y;		//function values at the vertices of the simplex
	Vector2r p[3];	//current simplex

	p[0] = p[1] =p[2] =point;//each row is filled with the ndim dimentional vector (point)
	p[1](0) += dels(0);   //the first row is special. Add dels to the diagonal element of the sub matrix excluding the first row.
	p[2](1) += dels(1);
        Vector2r p_tmp;double y_tmp,check_dc(0.);
        p_tmp = p[0];
        y_tmp = 1000000.0;//
        Vector3r contact(0.,0.,0.),d;
	for (int i=0;i<3;i++){
		y(i)=func(p[i],center, p1, p2,position1, position2, pt1, pt2,globalContact,stop_flag);
		nfunc++;//function values
		if (y(i) < y_tmp){
		        p_tmp = p[i];
		        y_tmp = y(i);
		}
		if (stop_flag){//stop
		        touching = false;
		        info[0] = nfunc;	//number of function invoking
			info[1] = y(0);		//the shortest distance
		        return p[0];

		}
		check_dc = Check_DC(p[i], p1, p2, globalContact);
		//
	        contact(0) = cos(p[i](0))*cos(p[i](1));  //the base vectors are e_i, i=1,2,3
	        contact(1) = sin(p[i](0))*cos(p[i](1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	        contact(2) = sin(p[i](1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	        contact = globalContact*contact;	   //to global
	        d = pt2 - pt1;
	        check_dc = acos(d.dot(contact)/(d.norm()*contact.norm()));//phi in [0, pi]
		if (  check_dc > 2.8 ){//
		        info[0] = nfunc;	//number of function invoking
			info[1] = y_tmp;		//the shortest distance
		        return p_tmp;

		}
	}
	psum = p[0]+p[1]+p[2];
        //nfun:3
	for (;;) {
		ilo=0;
		if (stop_flag){//stop
		        touching = false;
		        info[0] = nfunc;	//number of function invoking
			info[1] = fmin;		//the shortest distance
		        return p[0];

		}
		//First, we must determine which point is the highest (worst), next-highest, and lowest (best),
		//by loping over the points in the simplex.
		ihi = y(0)>y(1) ? (inhi=1,0) : (inhi=0,1);
		for (int i=0;i<3;i++) {
			if (y(i) <= y(ilo)) ilo=i;
			if (y(i) > y(ihi)) {
				inhi=ihi;
				ihi=i;
			} else if (y(i) > y(inhi) && i != ihi) inhi=i;
		}
		//std::cout<<"y="<<y<<std::endl;
		//std::cout<<"ihi"<<ihi<<"inhi"<<inhi<<"ilo"<<ilo<<std::endl;
		//Compute the fractional range from highest to lowest and return if satisfactory
		double rtol=2.0*fabs(y(ihi)-y(ilo))/(fabs(y(ihi))+fabs(y(ilo))+TINY);
		//std::cout<<"rtol="<<rtol<<"ftol"<<ftol<<std::endl;
		if (rtol < ftol) {	//if returning, put the best point and value in slot 0.
			double tmp =y(0);//swap y(0) and y(il0)
			y(0) = y(ilo);
			y(ilo) = tmp;

			Vector2r tmp1 = p[ilo];//swap p[ilo] p[0]
			p[ilo] = p[0];
			p[0] = tmp1;
			//return data
			pmin=p[0];

			fmin=y(0);
			//return some info
			info[0] = nfunc;	//number of function invoking
			info[1] = fmin;		//the shortest distance
			return pmin;
		}
		if (nfunc >= NMAX) throw("NMAX exceeded");	//this is for test. It should be removed in the final rutine.
		//nfunc += 2;
		//Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
		//across from the high point, i.e., reflect the simplex from the high point.
		double ytry=amotry(p,y,psum,ihi,-1.0,func,center, p1, p2,position1, position2, pt1, pt2, globalContact,stop_flag);nfunc++;
		//nfunc:4
		if (stop_flag){//stop
		        touching = false;
		        info[0] = nfunc;	//number of function invoking
			info[1] = y(0);		//the shortest distance
		        return p[0];

		}
		if (ytry <= y(ilo))
			//Gives a result better than the best point, so try an additional extrapolation by a factor 2.
			{ytry=amotry(p,y,psum,ihi,2.0,func,center, p1, p2,position1, position2, pt1, pt2, globalContact,stop_flag);nfunc++;
			//nfunc:5
			if (stop_flag){//stop
		                touching = false;
		                info[0] = nfunc;	//number of function invoking
			        info[1] = y(0);		//the shortest distance
		                return p[0];

		        }

		        }
		else if (ytry >= y(inhi)) {
			//The reflected point is worse than the second-highest, so look for an intermediate lower point,
			//i.e., do a one-dimensional contraction.
			double ysave=y(ihi);
			ytry=amotry(p,y,psum,ihi,0.5,func,center, p1, p2,position1, position2, pt1, pt2, globalContact,stop_flag);nfunc++;
			//nfunc:6
			if (stop_flag){//stop
		                touching = false;
		                info[0] = nfunc;	//number of function invoking
			        info[1] = y(0);		//the shortest distance
		                return p[0];

		        }
			if (ytry >= ysave) {	//Can't seem to get rid of that high point.
				for (int i=0;i<3;i++) {	//Better contract around the lowest (best) point.
					if (i != ilo) {
						p[i] = psum = 0.5*(p[i]+p[ilo]);
						y(i)=func(psum,center, p1, p2,position1, position2, pt1, pt2, globalContact,stop_flag);nfunc++;
						if (stop_flag){//stop
		                                        touching = false;
		                                        info[0] = nfunc;	//number of function invoking
			                                info[1] = y(0);		//the shortest distance
		                                        return p[0];

		                                }
					}
				}
				//nfunc += ndim;	//Keep track of function evaluations.
				psum = p[0]+p[1]+p[2];//Recompute psum.
			}
		} //else --nfunc;	//Correct the evaluation count.
	}					//Go back for the test of doneness and the next iteration.
}
///////////////////////////
//optimal function to compute the distance
double Disfun(Vector2r para, Superquadrics *particle1, Superquadrics *particle2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2,
		     Matrix3r globalContact, bool &touching)
{	Vector3r d;
	//stop_flag = false;
	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
        Vector2r phi(0.,0.);

	contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	contact(2) = sin(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	//std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	contact = globalContact*contact;	   //to global
	//test
	//cout<<"rot mat0"<<particle1->rot_mat2local(0,0)<<" "<<particle1->rot_mat2local(0,1)<<" "<<particle1->rot_mat2local(0,2)<<endl;
	//cout<<"rot mat1"<<particle1->rot_mat2local(1,0)<<" "<<particle1->rot_mat2local(1,1)<<" "<<particle1->rot_mat2local(1,2)<<endl;
	//cout<<"rot mat2"<<particle1->rot_mat2local(2,0)<<" "<<particle1->rot_mat2local(2,1)<<" "<<particle1->rot_mat2local(2,2)<<endl;
	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
	contact1 = particle1->rot_mat2local*contact; //to particle's local system
	//std::cout<<"contact normal in local sys:  "<<contact<<std::endl;
	phi= particle1->Normal2Phi( contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	p1 = particle1->getSurface( phi );
	p1 = particle1->rot_mat2global*p1 + position1;
	//point 2 on Particle 2
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	phi = particle2->Normal2Phi( (-1)*contact2);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
        //get the coordinates of the current point in the global system

	p2 = particle2->getSurface( phi );
	p2 = particle2->rot_mat2global*p2 + position2;
    #ifdef _DEBUG_LM_NM
    FILE * fin1 = fopen("para_NM.dat","a");

    fprintf(fin1,"\t%e\t%e\n",para[0],para[1]);
	fclose(fin1);
    #endif
	d = p2 -p1;
	//we should to check whether the angle between n1 (contact) and d is less than pi/2.
	//If yes, then we should stop finding the shortest distance. In this case, the
	//two particles are not touching. Thus the further detection is not necessary.
	if (contact.dot(d) > 0)  //the angle < pi/2.
	{        //cout<<"t---"<<touching<<endl;
		touching = false;

	}

	return d.norm();//CAUTION:d is the vector from point 1 on the particle1 to point 2 on the particle2.

}


double amotry(Vector2r (&p)[3], Vector3r &y, Vector2r &psum,
		const int ihi, const double fac,
	       Superquadrics *p1, Superquadrics *p2,
	      Vector3r position1, Vector3r position2,Vector3r &pt1,Vector3r &pt2,
	      Matrix3r globalContact,bool &touching)
{
	Vector2r ptry;
	double fac1=(1.0-fac)/2.0;
	double fac2=fac1-fac;
	ptry = psum*fac1 - p[ihi]*fac2;

	double ytry=Disfun(ptry, p1, p2,position1, position2, pt1, pt2, globalContact,touching); //Evaluate the function at the trial point.
	if (ytry < y(ihi)) {	//If it's better than the highest, then replace the highest.
		y(ihi)=ytry;
		psum += ptry - p[ihi];
		p[ihi] = ptry;
	}
	return ytry;
}

Vector2r neldermead_simplex(
	Vector2r point,
	Vector2r dels,
	double info[2], //
	Superquadrics *p1,
	Superquadrics *p2,
	Vector3r position1, Vector3r position2,Vector3r &pt1,Vector3r &pt2,
	Matrix3r globalContact,
	bool &touching
)
{
	const double ftol = 1e-3;//1e-3 for packing; //the fractional convergence tolerance to be achieved in the function
	int nfunc;		//the number of function evaluations
    //bool stop_flag(false);
	double fmin;	//function value at the minimum
	const int NMAX=5000;		//maximum allowed number of function evaluations
	const double TINY=1.0e-10;
	int ihi,ilo,inhi;
	Vector2r psum, pmin;
	Vector3r y;		//function values at the vertices of the simplex
	Vector2r p[3];	//current simplex
    //std::cerr << "watch 4.01"<< "\n";
	p[0] = p[1] =p[2] =point;//each row is filled with the ndim dimentional vector (point)
	p[1](0) += dels(0);   //the first row is special. Add dels to the diagonal element of the sub matrix excluding the first row.
	p[2](1) += dels(1);
        Vector2r p_tmp;double y_tmp,check_dc(0.);
        p_tmp = p[0];
        y_tmp = 1000000.0;//
        Vector3r contact(0.,0.,0.),d;
	for (int i=0;i<3;i++){
		y(i)=Disfun(p[i],p1, p2,position1, position2, pt1, pt2,globalContact,touching);
		nfunc++;//function values
        //std::cerr << "watch 4.02"<< "\n";
		if (y(i) < y_tmp){
		        p_tmp = p[i];
		        y_tmp = y(i);
		}
		if (!touching){//stop
		        //touching = false;

		        info[0] = nfunc;	//number of function invoking
			info[1] = y(0);		//the shortest distance
		        return p[0];

		}
		/*
		//
	        contact(0) = cos(p[i](0))*cos(p[i](1));  //the base vectors are e_i, i=1,2,3
	        contact(1) = sin(p[i](0))*cos(p[i](1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	        contact(2) = sin(p[i](1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	        contact = globalContact*contact;	   //to global
	        d = pt2 - pt1;

	        check_dc = acos(d.dot(contact)/(d.norm()*contact.norm()));//phi in [0, pi]
		if (  check_dc > 2.8 ){//
		        info[0] = nfunc;	//number of function invoking
			info[1] = y_tmp;		//the shortest distance
		        return p_tmp;

		}
		*/
	}
	psum = p[0]+p[1]+p[2];
        //nfun:3
    //std::cerr << "watch 4.03"<< "\n";
	for (;;) {
		ilo=0;
		if (!touching){//stop
		        //touching = false;
		        info[0] = nfunc;	//number of function invoking
			info[1] = fmin;		//the shortest distance
		        return p[0];

		}
		//First, we must determine which point is the highest (worst), next-highest, and lowest (best),
		//by loping over the points in the simplex.
		ihi = y(0)>y(1) ? (inhi=1,0) : (inhi=0,1);
		for (int i=0;i<3;i++) {
			if (y(i) <= y(ilo)) ilo=i;
			if (y(i) > y(ihi)) {
				inhi=ihi;
				ihi=i;
			} else if (y(i) > y(inhi) && i != ihi) inhi=i;
		}
		//std::cout<<"y="<<y<<std::endl;
		//std::cout<<"ihi"<<ihi<<"inhi"<<inhi<<"ilo"<<ilo<<std::endl;
		//Compute the fractional range from highest to lowest and return if satisfactory
		double rtol=2.0*fabs(y(ihi)-y(ilo))/(fabs(y(ihi))+fabs(y(ilo))+TINY);
		//std::cout<<"rtol="<<rtol<<"ftol"<<ftol<<std::endl;
		if (rtol < ftol) {	//if returning, put the best point and value in slot 0.
			double tmp =y(0);//swap y(0) and y(il0)
			y(0) = y(ilo);
			y(ilo) = tmp;

			Vector2r tmp1 = p[ilo];//swap p[ilo] p[0]
			p[ilo] = p[0];
			p[0] = tmp1;
			//return data
			pmin=p[0];

			fmin=y(0);
			//return some info
			info[0] = nfunc;	//number of function invoking
			info[1] = fmin;		//the shortest distance
			return pmin;
		}
		if (nfunc >= NMAX) {
            //std::cerr << "watch 4.03"<< "\n";
            //LOG_WARN("NMAX ecceeded "<<NMAX<<" of ");
            info[0] = nfunc;	//number of function invoking
			info[1] = y(0);		//the shortest distance

			LOG_WARN("NMAX ecceeded");
             return p[0];
            //throw runtime_error("NMAX exceeded");	//this is for test. It should be removed in the final rutine.
        }
		//nfunc += 2;
		//Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
		//across from the high point, i.e., reflect the simplex from the high point.
		double ytry=amotry(p,y,psum,ihi,-1.0, p1, p2,position1, position2, pt1, pt2, globalContact,touching);nfunc++;
		//nfunc:4
		if (!touching){//stop
		        //touching = false;
		        info[0] = nfunc;	//number of function invoking
			info[1] = y(0);		//the shortest distance
		        return p[0];

		}
		if (ytry <= y(ilo))
			//Gives a result better than the best point, so try an additional extrapolation by a factor 2.
			{ytry=amotry(p,y,psum,ihi,2.0, p1, p2,position1, position2, pt1, pt2, globalContact,touching);nfunc++;
			//nfunc:5
			if (!touching){//stop
		                //touching = false;
		                info[0] = nfunc;	//number of function invoking
			        info[1] = y(0);		//the shortest distance
		                return p[0];

		        }

		        }
		else if (ytry >= y(inhi)) {
			//The reflected point is worse than the second-highest, so look for an intermediate lower point,
			//i.e., do a one-dimensional contraction.
			double ysave=y(ihi);
			ytry=amotry(p,y,psum,ihi,0.5, p1, p2,position1, position2, pt1, pt2, globalContact,touching);nfunc++;
			//nfunc:6
			if (!touching){//stop
		                //touching = false;
		                info[0] = nfunc;	//number of function invoking
			        info[1] = y(0);		//the shortest distance
		                return p[0];

		        }
			if (ytry >= ysave) {	//Can't seem to get rid of that high point.
				for (int i=0;i<3;i++) {	//Better contract around the lowest (best) point.
					if (i != ilo) {
						p[i] = psum = 0.5*(p[i]+p[ilo]);
						y(i)=Disfun(psum, p1, p2,position1, position2, pt1, pt2, globalContact,touching);nfunc++;
						if (!touching){//stop
		                                        //touching = false;
		                                        info[0] = nfunc;	//number of function invoking
			                                info[1] = y(0);		//the shortest distance
		                                        return p[0];

		                                }
					}
				}
				//nfunc += ndim;	//Keep track of function evaluations.
				psum = p[0]+p[1]+p[2];//Recompute psum.
			}
		} //else --nfunc;	//Correct the evaluation count.
	}					//Go back for the test of doneness and the next iteration.
}

//**********************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Superquadricss. */

bool Ig2_Superquadrics_Superquadrics_SuperquadricsGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){
    //std::cerr << "watch 1"<< "\n";
	//get polyhedras
	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
	Vector3r position1 = se31.position;
	Vector3r position2 = se32.position + shift2;
    //std::cerr << "watch 2"<< "\n";
	Superquadrics* A = static_cast<Superquadrics*>(cm1.get());
	Superquadrics* B = static_cast<Superquadrics*>(cm2.get());
        //cout<<"watch point"<<endl;
	bool isNew = !interaction->geom;
	//std::cerr << "pos1:" <<position1 << "\n";
	//std::cerr << "pos2:" <<position2 << "\n";

	//Vector3r p = A->getPosition();
	//std::cerr << "Gl1_Superquadrics:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";
        //A->rot_mat2local = (se31.orientation).conjugate().toRotationMatrix();//to particle's system
	//A->rot_mat2global = (se31.orientation).toRotationMatrix(); //to global system
	//B->rot_mat2local = (se32.orientation).conjugate().toRotationMatrix();//to particle's system
	//B->rot_mat2global = (se32.orientation).toRotationMatrix(); //to global system
	shared_ptr<SuperquadricsGeom> bang;
	Vector2r dels(0.01,0.01);
	Vector2r para(1e-5,1e-5);
    //std::cerr << "watch 3"<< "\n";

	if (isNew) {
		// new interaction
		//cout<<"new intersection"<<endl;
		bang=shared_ptr<SuperquadricsGeom>(new SuperquadricsGeom());
		//bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
		dels = Vector2r(0.1,0.1);
        para = bang->contactAngle;
        //cout<<"para0="<<para<<endl;
        #ifdef _USE_LM_OPTIMIZATION
            //cout<<"new geom!..."<<endl;
            Vector3r d0 = position2 - position1;
            //atan2(y,x): arc tangent of y/x
            para(0) = atan2(d0(1),d0(0));
            //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
            para(1) = atan(d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
            para(0) += 1e-5;para(1) += 1e-5;
            //cout<<"new para="<<para<<endl;
            //to do ...
        #endif
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperquadricsGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
        para = bang->contactAngle;
        //cout<<"para="<<para<<endl;
		//Vector3r v1 = state1.vel;
		//Vector3r v2 = state2.vel;//maybe a bug
		//cout<<"vvv1"<<v1(0)<<" "<<v1(1)<<" "<<v1(2)<<"  "<<v1.isZero(0)<<endl;
		//cout<<"vvv2"<<v2(0)<<" "<<v2(1)<<" "<<v2(2)<<"  "<<v2.isZero(0)<<endl;

		//if (v1.isZero(0)&&v2.isZero(0)){
		//        bang->PenetrationDepth=0;

               //         return true;
		//}
		dels = Vector2r(0.001,0.001);
	}
	Vector3r contactNormal(0,0,0), center(0,0,0);
  double pDepth(0.0);
	Vector3r p1,p2;
	//if the two particles are spherical
	if (A->isSphere && B->isSphere){
	        //particles are spherical
	        Vector3r dist = position2 -position1;
	        double r_sum = A->getr_max()+B->getr_max();
	        double depth = r_sum - dist.norm();
	        if (depth < 1E-18){//no contact
                bang->PenetrationDepth=0;
                bang->isShearNew = true;
                return true;
            }
            //compute contact geometric quantities
            dist.normalize();
            p1 = position1 + A->getr_max()*dist;
            p2 = position2 - B->getr_max()*dist;
            contactNormal = p1 - p2;
	        contactNormal.normalize();
	        bang->point1 = p1;
            bang->point2 = p2;

	        bang->contactPoint=0.5*(p1+p2);
	        bang->PenetrationDepth = depth;
	}else{
	//contact detection using bounding spheres
        double r_sum = A->getr_max()+B->getr_max();

        if ((position1 - position2).norm() > r_sum){//no contact
                bang->PenetrationDepth=0;
                bang->isShearNew = true;
                return true;
        }
        //Matrix3r globalContact = Matrix3r::Identity();
        //cout<<"watch point2"<<endl;
	    //globalContact = CRotationMat(position1, position2);
	    //cout<<"globalContact="<<globalContact<<endl;
	    //std::cerr << "watchpoint 1: " << "\n";
	    //para should be noted
	    //para = bang->contactAngle;

	    //nelder-mead simplex method
	    //Vector2r pmin;
	    bool touching(true);
	    double nm_info[2];

	    //std::cerr << "watch 4"<< "\n";

	    //CheckDC(para, A, B, globalContact);
	    //para = neldermead(Disfunction, para, dels,nm_info, A,B,rot_mat1,rot_mat2,position1,position2,p1,p2,globalContact,center,touching);//not used
	    //para =neldermead(Disfunction, para, dels,nm_info, A,B,position1,position2,p1,p2,globalContact,center,touching);//
        #ifdef _USE_LM_OPTIMIZATION

	    //Levenberg-Marquadt method
        //para(0) = 0.001;para(1)=0.001;
	    double opts[5], info[10];
        //cout<<"para="<<para<<endl;
	    /* optimization control parameters */
	    //opts[0]=0.01; opts[1]=1E-15; opts[2]=1E-5; opts[3]=1E-20;opts[4]=1e-6;
        opts[0]=1.0; opts[1]=1E-15; opts[2]=1E-3; opts[3]=1E-20;opts[4]=1e-6;


	    /* invoke the optimization function */
	    //lmder_z(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2, globalContact,center,touching);

        lmder2(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching,contactNormal);
        if(isnan(para[0]) || isnan(para[1]) || 7==info[6]){
            //FIXME: .
            if(touching){
                cout<<"Para or d is NAN. RESET para."<<endl;
                Vector3r d0 = position2 - position1;
                //atan2(y,x): arc tangent of y/x
                para(0) = atan2(d0(1),d0(0));
                //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
                para(1) = atan(d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
                para(0) += 1e-5;para(1) += 1e-5;
                //may be not safe
                lmder2(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching,contactNormal);
                //std::cerr << "id1="<<interaction->getId1()<<"id2="<<interaction->getId2()<<"\n";
                //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
                //scene->stopAtIter = scene->iter + 1;
            }
        }
        pDepth = (p1-p2).dot(contactNormal);
        if(touching){
            //check if the contact depth is not too small.
            if(10.0*(p2-p1).norm() > A->getr_max()){//the contact depth is tremendously large which might be at local minimia. The common contact depth should be 3 order of magnitude less than particle size.FIXME:the coefficient of 10.0 may be not the best try.
                if(!(A->isInside(p2-position1)) || !(B->isInside(p1-position2))){
                  //fixed:~(1>0.6)=-2 (random?);!(1>0.6)=0;so use ! for not operator.
                    //FIXME:the inside-outside function may not be sufficiently accurate for tiny contact depth during computation.

                    //cout<<"local minimia!"<<endl;
                    //cout<<"contact depth="<<(p2-p1).norm()<<endl;
                    //std::cerr << "id1="<<interaction->getId1()<<"id2="<<interaction->getId2()<<"\n";
                    //scene->stopAtIter = scene->iter + 1;
                    Vector3r d0 = position2 - position1;
                    //atan2(y,x): arc tangent of y/x
                    para(0) = atan2(d0(1),d0(0));
                    //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
                    para(1) = atan(d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
                    para(0) += 1e-5;para(1) += 1e-5;
                    //may be not safe
                    lmder2(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching,contactNormal);
                    pDepth = (p1-p2).dot(contactNormal);
                }
            }
        }
	    //contactNormal = -contactNormal;//make contact normal along the direction from particle 1 to particle 2.
        //#ifdef _DEBUG_LM_NM
	    //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
        //#endif
	    //printf("Best fit parameters: %.7g %.7g \n", para(0), para(1));


        //std::cerr << "touching: " <<touching<< "\n";
        //std::cerr << "watch 4.1"<< "\n";

        #else //use Nelder-Mead
        Matrix3r globalContact = Matrix3r::Identity();
        //cout<<"watch point2"<<endl;
	    globalContact = CRotationMat(position1, position2);
	    para =neldermead_simplex(para, dels,nm_info, A,B,position1,position2,p1,p2,globalContact,touching);//
        //CheckDC(para, A, B, globalContact);
        //para = pmin;
        //#ifdef _DEBUG_LM_NM
	    printf("Nelder-Mead method returned in %g iter, distance %g \n", nm_info[0], nm_info[1]);
        //#endif
	    //printf("Best fit parameters of nelder mead: %.7g %.7g \n", para(0), para(1));
	    //printf("Nelder-Mead: %g iter, distance %g, parameters %.7g %.7g\n", nm_info[0], nm_info[1], para(0), para(1));
        #endif

	    if (!touching){
	            bang->contactAngle = para;
	            bang->PenetrationDepth=0;
	            bang->isShearNew = true;
	            return true;}

	    //for testing the two closest points
	    //std::cerr << "watch 5"<< "\n";
	    //contactNormal = p1 - p2;
	    //contactNormal.normalize();
	    //cout<<"contact normal"<<contactNormal(0)<<" "<<contactNormal(1)<<" "<<contactNormal(2)<<endl;
	    //Real norm=contactNormal.norm();  // normal is unit vector now
	    //Vector3r contactNormal1 = contactNormal/norm;
	    //cout<<"contact normal111"<<contactNormal1(0)<<" "<<contactNormal1(1)<<" "<<contactNormal1(2)<<endl;
	    //std::cerr << "watchpoint 2: " << "\n";

	    //if((se32.position-centroid).dot(normal)<0) normal*=-1;
            bang->point1 = p1;
            bang->point2 = p2;

/*
	    Mat32r JacMat1, JacMat2;

	    Vector3r contact1(0.,0.,0.),contact2(0.,0.,0.);
        //Vector2r phi;
	    //contact(0) = cos(para(0))*sin(para(1));  //the base vectors are e_i, i=1,2,3
	    //contact(1) = sin(para(0))*sin(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	    //contact(2) = cos(para(1));

	    //contact = globalContact*contact;	   //to global
	    //Matrix3r Rot = Orientation.conjugate().toRotationMatrix();
	    contact1 = A->rot_mat2local*contact; //to particle's local system
	    //contact.normalize();
	    //std::cout<<"contact normal= "<<contact<<std::endl;
	    phi = A->Normal2Phi( contact1);
        //cout<<"phi1= "<<phi<<endl;
        //cout<<"local rot= "<<particle1->rot_mat2local<<endl;
        //cout<<"global rot= "<<particle1->rot_mat2global<<endl;
        //cout<<"local*global= "<<particle1->rot_mat2local*particle1->rot_mat2global<<endl;
	    //get the coordinates of the current point in the global system
	    //p = getSurface( phi );  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	    //std::cout<<"phi1= "<<phi(0)<<"\t phi2= "<<phi(1)<<std::endl;
	    Mat32r d_p2phi = A->rot_mat2global*A->Derivate_P2Phi(phi);//to global system
        //test P2Phi
        Vector3r p10 = A->getSurface(Vector2r(phi(0),phi(1)));
        Vector3r p11 = A->getSurface(Vector2r(phi(0)+1e-5,phi(1)));
        Vector3r p12 = A->getSurface(Vector2r(phi(0),phi(1)+1e-5));
        Mat32r P12Phi_d;
        for(int i=0;i<3;i++){
            P12Phi_d(i,0) = (p11(i)-p10(i))/1e-5;
        }
        for(int i=0;i<3;i++){
            P12Phi_d(i,1) = (p12(i)-p10(i))/1e-5;
        }
        cout<<"Approx. d_p12phi="<<A->rot_mat2global*P12Phi_d<<endl;
        cout<<"d_p12phi="<<d_p2phi<<endl;

	    Mat23r d_phi2n = A->Derivate_Phi2N(phi, contact1);
        //cout<<"d_phi2n="<<d_phi2n<<endl;
        Mat32r c2alpha;
	    c2alpha<<
			    -sin(para(0))*sin(para(1)), cos(para(0))*cos(para(1)),
			    cos(para(0))*sin(para(1)), sin(para(0))*cos(para(1)),
			    0., -sin(para(1));
	    Mat32r c2alpha1 = A->rot_mat2local*c2alpha;
        //cout<<"c2alpha ="<<c2alpha<<endl;
	    JacMat1 = d_p2phi*d_phi2n*c2alpha1;
        //test JacMat1

        Vector3r pa0 = SurfaceP(Vector2r(para(0),para(1)),A);
        Vector3r pa1 = SurfaceP(Vector2r(para(0)+1e-7,para(1)),A);
        Vector3r pa2 = SurfaceP(Vector2r(para(0),para(1)+1e-7),A);
        Mat32r JacMat1_d;
        for(int i=0;i<3;i++){
            JacMat1_d(i,0) = (pa1(i)-pa0(i))/1e-7;
        }
        for(int i=0;i<3;i++){
            JacMat1_d(i,1) = (pa2(i)-pa0(i))/1e-7;
        }
        cout<<"Jacobian mat1="<<JacMat1<<endl;
        cout<<"Approx. Jacobian mat1="<<JacMat1_d<<endl;

*/

        //test jacobian matrix ends
        //std::cerr << "watch 6"<< "\n";
	    //contactNormal = p1 - p2;
	    if(isnan(contactNormal[0])){
            //FIXME:For very flat shapes, the compatation is etremely expensive.
            //Not fatal errorr for only one iteration, so just pass it but with error log outputing.
            //FIXME: some issues may prompt out for some extreme shapes.
            //std::cerr << "id1="<<interaction->getId1()<<"id2="<<interaction->getId2()<<"\n";
            //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
            //scene->stopAtIter = scene->iter + 1;
            time_t tmNow = time(NULL);
            char tmp[64];
            strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&tmNow) );
            const char* filename = "SuperquadricsLaw_err.log";
            //sprintf(filename,"%d",gettid4());
            FILE * fin = fopen(filename,"a");
            fprintf(fin,"\n*******%s****NAN CONTACT NORMAL***\n",tmp);
            fprintf(fin,"LM %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
            fprintf(fin,"id1= %d id2= \t%d\n",interaction->getId1(), interaction->getId2());

            fprintf(fin,"point1\t%e\t%e\t%e point2\t%e\t%e\t%e\n",p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);

            fprintf(fin,"contact normal \t%e\t%e\t%e\n",contactNormal[0],contactNormal[1],contactNormal[2]);
            fprintf(fin,"para\t%e\t%e\n",para[0],para[1]);
            //fprintf(fin,"depth\t%e\n",depth);
            fclose(fin);
        }

        #ifdef _USE_LM_OPTIMIZATION
	        bang->contactPoint=center;
	        bang->PenetrationDepth=pDepth; //contactNormal.norm();


            //testing
            //SHOWN NORMAL BEGINS
            #define GEOM_GL1
            #ifdef GEOM_GL
            Vector3r contact(0.,0.,0.);
            Vector2r phi;
	        contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	        contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	        contact(2) = sin(para(1));

	        phi= A->Normal2Phi(A->rot_mat2local*contact);
            Vector3r n1,n2;
            n1 = A->getNormal(phi);
            n1.normalize();
            bang->point11 = A->rot_mat2global*n1*A->getr_max()+p1;

            phi= B->Normal2Phi(B->rot_mat2local*contact*(-1.0));
            n2 = B->getNormal(phi);
            n2.normalize();
            //cout<<"r_max="<<B->getr_max()<<endl;
            //cout<<"n2="<<n2.norm()<<endl;
            //Vector3r tmp1;
            //tmp1 = B->rot_mat2global*n2;
            //cout<<"tmp1="<<tmp1<<" magnitude: "<<tmp1.norm()<<endl;
	        bang->point22 = B->rot_mat2global*n2*B->getr_max()+p2;
            #endif
            //cout<<"length of the normal stick: "<<(bang->point22 - p2).norm()<<endl;
            //SHOW NORMAL ENDS
            #ifdef _DEBUG_LM_NM
            //output
            if(0){//if(isnan(p1(0)) || isnan(p2(0))){
                //cout<<"NAN ,p1"<<p1<<"p2="<<p2<<endl;//scene->stopAtIter = scene->iter + 1;
                FILE * fin = fopen("dis_func.dat","w");
                double dist;bool stop_flag;
                for(int i =0;i<40;i++){
                    for(int j =0;j<40;j++){
                        para(0) = Mathr::PI*0.05*0.5*(i-20)+1e-5;
                        para(1) = Mathr::PI*0.025*0.5*(j-20)+1e-5;
                        //Disfunction(para, center, A, B,position1, position2,p1,p2, stop_flag);
                        dist = DisContact(para, center, A, B,position1, position2,p1,p2, stop_flag);
                        //dist = (p2 - p1).norm();
                        fprintf(fin,"\t%e\t%e\t%e\n",para(0),para(1),dist);
                    }
                }
                para(0)=-0.0886504;para(1)= 0.497741;
                dist = DisContact(para, center, A, B,position1, position2,p1,p2, stop_flag);
                //dist = (p2 - p1).norm();
                fprintf(fin,"\t%e\t%e\t%e\n",para(0),para(1),dist);
                //Disfunction(para, center, A, B,position1, position2,p1,p2, stop_flag);
                //dist = (p2 - p1).norm();
                //fprintf(fin,"\t%e\t%e\t%e\n",para[0],para[1],dist);
		        fclose(fin);
            }
            #endif
        #else
	        bang->contactPoint=0.5*(p1+p2);
	        bang->PenetrationDepth=nm_info[1];//info[4];
            #ifdef _DEBUG_LM_NM
            Vector3r contact(0.,0.,0.);
            Vector2r phi;
	        contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
	        contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	        contact(2) = sin(para(1));
	        contact = globalContact*contact;//to the global sys.
	        phi= A->Normal2Phi(A->rot_mat2local*contact);
            Vector3r n1,n2;
            n1 = A->getNormal(phi);
            n1.normalize();
            bang->point11 = A->rot_mat2global*n1*A->getr_max()+p1;

            phi= B->Normal2Phi(B->rot_mat2local*contact*(-1.0));
            n2 = B->getNormal(phi);
            n2.normalize();
            //cout<<"r_max="<<B->getr_max()<<endl;
            //cout<<"n2="<<n2.norm()<<endl;
            //Vector3r tmp1;
            //tmp1 = B->rot_mat2global*n2;
            //cout<<"tmp1="<<tmp1<<" magnitude: "<<tmp1.norm()<<endl;
	        bang->point22 = B->rot_mat2global*n2*B->getr_max()+p2;
            //output
            if(1){//if(isnan(p1(0)) || isnan(p2(0))){
                //cout<<"NAN ,p1"<<p1<<"p2="<<p2<<endl;//scene->stopAtIter = scene->iter + 1;
                FILE * fin = fopen("dis_func_NM.dat","w");
                double dist;bool stop_flag;
                for(int i =0;i<40;i++){
                    for(int j =0;j<40;j++){
                        para[0] = Mathr::PI*0.05*(i-20)+1e-5;
                        para[1] = Mathr::PI*0.025*(j-20)+1e-5;
                        dist = Disfun(para, A, B, position1, position2,p1,p2,globalContact, touching);
                        //dist = (p2 - p1).norm();
                        fprintf(fin,"\t%e\t%e\t%e\n",para[0],para[1],dist);
                    }
                }
		        fclose(fin);
            }
            #endif
        #endif
        contactNormal.normalize();
	    /*if (!touching){
	            bang->contactAngle = para;
	            bang->PenetrationDepth=0;
	            bang->isShearNew = true;
	    }*/
	    //std::cerr << "watchpoint 3: " << "\n";
	}//non-spherical end
	bang->precompute(state1,state2,scene,interaction,contactNormal,bang->isShearNew,shift2);
	bang->normal= contactNormal;
	bang->contactAngle = para;


	/*
	FILE * fin = fopen("Interactions.dat","a");
	fprintf(fin,"************** IDS %d %d **************\n",interaction->id1, interaction->id2);
	fprintf(fin,"volume\t%e\n",volume);
	fprintf(fin,"centroid\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2]);
	PrintPolyhedron2File(PA,fin);
	PrintPolyhedron2File(PB,fin);
	PrintPolyhedron2File(Int,fin);
	fclose(fin);
	*/
    //std::cerr << "watch 6"<< "\n";
	return true;
}
//**********************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Superquadricss. */
//This is for Hertz-Mindlin model

bool Ig2_Superquadrics_Superquadrics_SuperquadricsGeom2::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){
    //std::cerr << "watch 1"<< "\n";
	//get polyhedras
	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
	Vector3r position1 = se31.position;
	Vector3r position2 = se32.position;
    //std::cerr << "watch 2"<< "\n";
	Superquadrics* A = static_cast<Superquadrics*>(cm1.get());
	Superquadrics* B = static_cast<Superquadrics*>(cm2.get());
        //cout<<"watch point"<<endl;
	bool isNew = !interaction->geom;
	//std::cerr << "SuperquadricsGeom:" <<se31.position[0]<<"\t"<<se31.position[1]<<"\t"<<se31.position[2] << "\n";
	//std::cerr << "SuperquadricsGeom:" <<se32.position[0]<<"\t"<<se32.position[1]<<"\t"<<se32.position[2] << "\n";

	//Vector3r p = A->getPosition();
	//std::cerr << "Gl1_Superquadrics:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";
        //A->rot_mat2local = (se31.orientation).conjugate().toRotationMatrix();//to particle's system
	//A->rot_mat2global = (se31.orientation).toRotationMatrix(); //to global system
	//B->rot_mat2local = (se32.orientation).conjugate().toRotationMatrix();//to particle's system
	//B->rot_mat2global = (se32.orientation).toRotationMatrix(); //to global system
	shared_ptr<SuperquadricsGeom2> bang;
	Vector2r dels(0.01,0.01);
	Vector2r para(0.0,0.0);
    //std::cerr << "watch 3"<< "\n";
	if (isNew) {
		// new interaction
		//cout<<"new intersection"<<endl;
		bang=shared_ptr<SuperquadricsGeom2>(new SuperquadricsGeom2());
		//bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
		dels = Vector2r(0.1,0.1);
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperquadricsGeom2>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
		//Vector3r v1 = state1.vel;
		//Vector3r v2 = state2.vel;//maybe a bug
		//cout<<"vvv1"<<v1(0)<<" "<<v1(1)<<" "<<v1(2)<<"  "<<v1.isZero(0)<<endl;
		//cout<<"vvv2"<<v2(0)<<" "<<v2(1)<<" "<<v2(2)<<"  "<<v2.isZero(0)<<endl;

		//if (v1.isZero(0)&&v2.isZero(0)){
		//        bang->PenetrationDepth=0;

               //         return true;
		//}
		dels = Vector2r(0.001,0.001);
	}
	Vector3r contactNormal(0,0,0), center(0,0,0);
	Vector3r p1,p2;
	//if the two particles are spherical
	if (A->isSphere && B->isSphere){
	        //particles are spherical
			//std::cerr << "spherical particles--"<< "\n";
	        Vector3r dist = position2 -position1;
	        double r_sum = A->getr_max()+B->getr_max();
	        double depth = r_sum - dist.norm();
	        if (depth < 1E-18){//no contact
                        bang->PenetrationDepth=0;
                        bang->isShearNew = true;
                        return true;
                }
                //compute contact geometric quantities
                dist.normalize();
                p1 = position1 + A->getr_max()*dist;
                p2 = position2 - B->getr_max()*dist;
                contactNormal = p1 - p2;
	        contactNormal.normalize();
	        bang->point1 = p1;
                bang->point2 = p2;

	        bang->contactPoint=0.5*(p1+p2);
	        bang->PenetrationDepth = depth;
			bang->curvatures_sum = r_sum/(A->getr_max() * B->getr_max());
			bang->isSphere = true;
	}else{
	//contact detection using bounding spheres
        double r_sum = A->getr_max()+B->getr_max();

        if ((position1 - position2).norm() > r_sum){//no contact
                bang->PenetrationDepth=0;
                bang->isShearNew = true;
                return true;
        }
        Matrix3r globalContact = Matrix3r::Identity();
        //cout<<"watch point2"<<endl;
	globalContact = CRotationMat(position1, position2);
	//cout<<"globalContact="<<globalContact<<endl;
	//std::cerr << "watchpoint 1: " << "\n";
	//para should be noted
	para = bang->contactAngle;

	//nelder-mead simplex method
	//Vector2r pmin;
	bool touching(true);
	double nm_info[2];

	//std::cerr << "watch 4"<< "\n";

	//CheckDC(para, A, B, globalContact);
	//para = neldermead(Disfunction, para, dels,nm_info, A,B,rot_mat1,rot_mat2,position1,position2,p1,p2,globalContact,center,touching);//not used
	//para =neldermead(Disfunction, para, dels,nm_info, A,B,position1,position2,p1,p2,globalContact,center,touching);//
	para =neldermead_simplex(para, dels,nm_info, A,B,position1,position2,p1,p2,globalContact,touching);//
        //CheckDC(para, A, B, globalContact);
    //para = pmin;

        //std::cerr << "touching: " <<touching<< "\n";
    //std::cerr << "watch 4.1"<< "\n";
	if (!touching){
	        bang->contactAngle = para;
	        bang->PenetrationDepth=0;
	        bang->isShearNew = true;
	        return true;}
	//for testing the two closest points
	//std::cerr << "watch 5"<< "\n";
	contactNormal = p1 - p2;
	contactNormal.normalize();
	//cout<<"contact normal"<<contactNormal(0)<<" "<<contactNormal(1)<<" "<<contactNormal(2)<<endl;
	//Real norm=contactNormal.norm();  // normal is unit vector now
	//Vector3r contactNormal1 = contactNormal/norm;
	//cout<<"contact normal111"<<contactNormal1(0)<<" "<<contactNormal1(1)<<" "<<contactNormal1(2)<<endl;
	//std::cerr << "watchpoint 2: " << "\n";

	//if((se32.position-centroid).dot(normal)<0) normal*=-1;
        bang->point1 = p1;
        bang->point2 = p2;
        //testing
        //bang->point11 = A->N_alpha12(para, globalContact,1)*0.5+p1;
	//bang->point22 = B->N_alpha12(para, globalContact,-1)*0.5+p2;
	bang->contactPoint=0.5*(p1+p2);
	bang->PenetrationDepth=nm_info[1];//info[4];

	//std::cerr << "watchpoint 3: " << "\n";

    //compute contact geometry for Hertz-Mindlin model
    //begin
    	//contact
	    Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);

	    contact(0) = cos(para(0))*sin(para(1));  //the base vectors are e_i, i=1,2,3
	    contact(1) = sin(para(0))*sin(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
	    contact(2) = cos(para(1));			   //e_i is pointing in the direction from the first particle center to the second particle center.

	    //std::cout<<"contact normal in its local sys: "<<contact<<std::endl;
	    contact = globalContact*contact;	   //to global

	    //Particle 1
	    contact1 = A->rot_mat2local*contact; //to particle's local system
	    //Particle 2
	    contact2 = B->rot_mat2local*contact; //to particle's local system


        Vector3r CurvatureDir_max1, CurvatureDir_max2,CurvatureDir_min1, CurvatureDir_min2,CurvatureDir1,CurvatureDir2;
        bool curvatureMaxflag1,curvatureMaxflag2;
        Vector2r curvatures1 = A->Curvatures(A->Normal2Phi( contact1),CurvatureDir_max1,CurvatureDir_min1,curvatureMaxflag1);
        Vector2r curvatures2 = B->Curvatures(B->Normal2Phi( (-1)*contact2),CurvatureDir_max2,CurvatureDir_min2,curvatureMaxflag2);//caution: particle 1 (default) or particle 2? if sign = -1, particle 2
        //using which curvature?
        if(curvatureMaxflag1 && curvatureMaxflag2){//using the maximum curvature
            CurvatureDir1 = CurvatureDir_max1;
            CurvatureDir2 = CurvatureDir_max2;
        }else{//using the other one
            CurvatureDir1 = CurvatureDir_min1;
            CurvatureDir2 = CurvatureDir_min2;
        }
        //calculate the angle between the two maximum principal curvatures from both particles.Note that, the sign has been considered
        Vector3r NormDir = contactNormal;

        Vector3r dir1 =Vector3r(CurvatureDir1(0)*NormDir(0),CurvatureDir1(1)*NormDir(1),CurvatureDir1(2)*NormDir(2));
        Vector3r dir2 =Vector3r(CurvatureDir2(0)*NormDir(0),CurvatureDir2(1)*NormDir(1),CurvatureDir2(2)*NormDir(2));
        //cout<<"dir1:"<<dir1<<"\t dir2"<<dir2<<"normDir"<<NormDir<<endl;
        dir1 -= CurvatureDir1;
        dir2 -= CurvatureDir2;
        //cout<<"direction:"<<dir1.norm()<<'\t'<<dir2.norm()<<endl;
        double angle = acos(dir1.dot(dir2)/(dir1.norm()*dir2.norm()));
        if(isnan(angle)){//this case is very rare. We force angle equal to zero for robust running.
            angle = 0.0;
        }
	    //cout<<"angle between the two maximum principal curvatures=========="<<angle*180./M_PIl<<endl;
        //test curvatures
        bang->curvatures1 = curvatures1;
        bang->curvatures2 = curvatures2;
        bang->curvatureAngle = angle;
	    //SuperqHerz(fabs(curvatues1(0)),fabs(curvatues1(1)),fabs(curvatues2(0)),fabs(curvatues2(1)),angle);
        bang->HertzMindlin(curvatures1, curvatures2, angle);

    //end
	}//non-spherical end
	bang->precompute(state1,state2,scene,interaction,contactNormal,bang->isShearNew,shift2);
	//bang->normal= contactNormal;
	bang->contactAngle = para;

    //std::cerr << "watch 6"<< "\n";
	return true;
}

//**********************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Polyhedron and Wall. */

bool Ig2_Wall_Superquadrics_SuperquadricsGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){

	const int& PA(cm1->cast<Wall>().axis);
	const int& sense(cm1->cast<Wall>().sense);
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	Superquadrics* B = static_cast<Superquadrics*>(cm2.get());

	bool isNew = !interaction->geom;
	shared_ptr<SuperquadricsGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<SuperquadricsGeom>(new SuperquadricsGeom());
		bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperquadricsGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}

	Vector3r normal = Vector3r(0,0,0);
	normal[PA] = sense;
	Vector3r wall_pos = se31.position;
	Vector3r B_pos = se32.position;
    // is the particle spherical?
    if (B->isSphere){
        //the particle is spherical
        Vector3r p1, point2;
        p1 = B_pos;
        p1[PA] = wall_pos[PA];//projection of the particle center on the wall
        point2 = B_pos - normal*B->getr_max();

        double depth = B->getr_max() -(B_pos - p1).norm();
        if (depth < 0.0){//no touching between the wall and the particle
            //cout<<"watch dog2"<<endl;
            bang->PenetrationDepth= 0;
            bang->isShearNew = true;
            return true;
        }

        bang->point1 = p1;
        bang->point2 = point2;
        bang->contactPoint= 0.5*(point2 + p1);
        bang->PenetrationDepth= depth;
    }else{

	    Matrix3r rot_mat1 = B->rot_mat2local;//(se32.orientation).conjugate().toRotationMatrix();//to particle's system
	    Matrix3r rot_mat2 = B->rot_mat2global;//(se32.orientation).toRotationMatrix();//to global system


	    //cout<<"wall_pos"<<wall_pos<<endl;
	    //cout<<"B_pos"<<B_pos<<endl;

	    Vector2r phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
	    Vector3r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
	    //check whether point2 is at the negative side of the wall.
	    Vector3r p1;
        p1 = point2;
        p1[PA] = wall_pos[PA];
	    Vector3r d = point2 - p1;	//the distance vector
	    //d[PA] -= wall_pos[PA];
	    //cout<<"distance v="<<d<<endl;
	    if (normal[PA]*d[PA]>=0.0){//(normal.dot(d) < 0.) {//no touching between the wall and the particle
		    //cout<<"watch dog2"<<endl;
		    bang->PenetrationDepth= 0;
		    bang->isShearNew = true;
		    return true;
	    }
	    //cout<<"watch dog1"<<endl;


        //cout<<"distance v="<<d.norm()<<endl;
        //cout<<"d--"<<d(0)<<" "<<d(1)<<" "<<d(2)<<endl;

        bang->point1 = p1;
        bang->point2 = point2;
	    bang->contactPoint= 0.5*(point2 + p1);
	    bang->PenetrationDepth= d.norm();
        bang->point11 = p1;
        bang->point22 = point2;
	}//non-spherical end
	//std::cerr << "watchpoint 3: " << "\n";
	bang->precompute(state1,state2,scene,interaction,normal,bang->isShearNew,shift2);
	bang->normal= normal;
	//bang->contactAngle = para;

	return true;
}
//**********************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Polyhedron and Wall. */

bool Ig2_Wall_Superquadrics_SuperquadricsGeom2::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){

	const int& PA(cm1->cast<Wall>().axis);
	const int& sense(cm1->cast<Wall>().sense);
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	Superquadrics* B = static_cast<Superquadrics*>(cm2.get());

	bool isNew = !interaction->geom;
	shared_ptr<SuperquadricsGeom2> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<SuperquadricsGeom2>(new SuperquadricsGeom2());
		bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperquadricsGeom2>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}

	Vector3r normal = Vector3r(0,0,0);
	normal[PA] = sense;
	Vector3r wall_pos = se31.position;
	Vector3r B_pos = se32.position;
        // is the particle spherical?
        if (B->isSphere){
                //the particle is spherical
                Vector3r p1, point2;
                p1 = B_pos;
                p1[PA] = wall_pos[PA];//projection of the particle center on the wall
                point2 = B_pos - normal*B->getr_max();

                double depth = B->getr_max() -(B_pos - p1).norm();
                if (depth < 0.0){//no touching between the wall and the particle
		        //cout<<"watch dog2"<<endl;
		        bang->PenetrationDepth= 0;
		        bang->isShearNew = true;
		        return true;
	        }


                bang->point1 = p1;
                bang->point2 = point2;
	        bang->contactPoint= 0.5*(point2 + p1);
	        bang->PenetrationDepth= depth;
        }else{

	Matrix3r rot_mat1 = B->rot_mat2local;//(se32.orientation).conjugate().toRotationMatrix();//to particle's system
	Matrix3r rot_mat2 = B->rot_mat2global;//(se32.orientation).toRotationMatrix();//to global system


	//cout<<"wall_pos"<<wall_pos<<endl;
	//cout<<"B_pos"<<B_pos<<endl;

	Vector2r phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
	Vector3r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
	//check whether point2 is at the negative side of the wall.
	Vector3r p1;
        p1 = point2;
        p1[PA] = wall_pos[PA];
	Vector3r d = point2 - p1;	//the distance vector
	//d[PA] -= wall_pos[PA];
	//cout<<"distance v="<<d<<endl;
	if (normal[PA]*d[PA]>=0.0){//(normal.dot(d) < 0.) {//no touching between the wall and the particle
		//cout<<"watch dog2"<<endl;
		bang->PenetrationDepth= 0;
		bang->isShearNew = true;
		return true;
	}
	//cout<<"watch dog1"<<endl;


        //cout<<"distance v="<<d.norm()<<endl;
        //cout<<"d--"<<d(0)<<" "<<d(1)<<" "<<d(2)<<endl;

        bang->point1 = p1;
        bang->point2 = point2;
	bang->contactPoint= 0.5*(point2 + p1);
	bang->PenetrationDepth= d.norm();
	}//non-spherical end
	//std::cerr << "watchpoint 3: " << "\n";
	bang->precompute(state1,state2,scene,interaction,normal,bang->isShearNew,shift2);
	bang->normal= normal;
	//bang->contactAngle = para;

	return true;
}
//**********************************************************************************
/*! Create Superquadrics (collision geometry) from colliding Polyhedron and Facet. */

bool Ig2_Facet_Superquadrics_SuperquadricsGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){


	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
	Facet*   A = static_cast<Facet*>(cm1.get());
	Superquadrics* B = static_cast<Superquadrics*>(cm2.get());

	bool isNew = !interaction->geom;

	shared_ptr<SuperquadricsGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<SuperquadricsGeom>(new SuperquadricsGeom());
		//bang->sep_plane.assign(3,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperquadricsGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}
/*
	//find intersection Superquadrics
	Polyhedron Int;
	Int = Polyhedron_Polyhedron_intersection(PA,PB,ToCGALPoint(bang->contactPoint),ToCGALPoint(se31.position),ToCGALPoint(se32.position), bang->sep_plane);

	//volume and centroid of intersection
	double volume;
	Vector3r centroid;
	P_volume_centroid(Int, &volume, &centroid);
 	if(isnan(volume) || volume<=1E-25 || volume > B->GetVolume()) {bang->equivalentPenetrationDepth=0; return true;}
	if (!Is_inside_Polyhedron(PB, ToCGALPoint(centroid)))  {bang->equivalentPenetrationDepth=0; return true;}

	//find normal direction
	Vector3r normal = FindNormal(Int, PA, PB);
	if((se32.position-centroid).dot(normal)<0) normal*=-1;

	//calculate area of projection of Intersection into the normal plane
        double area = volume/1E-8;
	//double area = CalculateProjectionArea(Int, ToCGALVector(normal));
	//if(isnan(area) || area<=1E-20) {bang->equivalentPenetrationDepth=0; return true;}

	// store calculated stuff in bang; some is redundant
	bang->equivalentCrossSection=area;
	bang->contactPoint=centroid;
	bang->penetrationVolume=volume;
	bang->equivalentPenetrationDepth=0;
	//bang->precompute(state1,state2,scene,interaction,normal,bang->isShearNew,shift2);
	bang->normal=normal;
*/
	return true;
}
