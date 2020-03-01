#include "PolySuperellipsoid.hpp"
#include "PolySuperellipsoid_Ig2.hpp"
#include "GJK.hpp"
//#include "sudodem/lib/base/Math.hpp"
#include <cmath>
//#include <sys/syscall.h>
#define _USE_MATH_DEFINES
#define USE_GJK_no
#define STATISTICS_no
#define NUM_ITERATION_no //test efficiency used
#define _USE_LM_OPTIMIZATION
#define _DEBUG_LM_NM1//debug flag
#define GEOM_GL_no//show contact points,no
//LM: Levenberg-Marquardt method
//NM: Nerld-Mead simplex method

#include<time.h>

SUDODEM_PLUGIN(/* self-contained in hpp: */ (Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom)
		(Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom)
		(Ig2_Facet_PolySuperellipsoid_PolySuperellipsoidGeom)
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
/*****************************************************/
template <typename T>
void Disfunction(Vector2r para, T *particle1, T *particle2,
		     Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag)
{
	Vector3r d;
	//stop_flag = false;
	//contact
	Vector3r contact(0.,0.,0.),contact1(0.,0.,0.),contact2(0.,0.,0.);
  //Vector2r phi(0.,0.);
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
	//phi= particle1->Normal2Phi(contact1);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
	//get the coordinates of the current point in the global system
	//p1 = particle1->getSurfaceMC( phi );
	p1 = particle1->Normal2SurfaceMC(contact1);
	p1 = particle1->rot_mat2global*p1 + position1;
	//p1 = particle1->Normal2SurfaceMCgl(contact) + position1;
	//point 2 on Particle 2
	contact2 = particle2->rot_mat2local*contact; //to particle's local system
	//phi = particle2->Normal2Phi( (-1)*contact2);  //caution: particle 1 (default) or particle 2? if sign = -1, particle 2
    //cout<<"phi of particle 2 = "<<phi<<endl;
        //get the coordinates of the current point in the global system
	//p2 = particle2->getSurfaceMC( phi );
	p2 = particle2->Normal2SurfaceMC((-1)*contact2);
    //cout<<"point of par 2= "<<p2<<endl;
	p2 = particle2->rot_mat2global*p2 + position2;
	//p2 = particle2->Normal2SurfaceMCgl((-1)*contact) + position2;
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
	//center = 0.5*(p1 + p2);
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

//secent Jacobian matrix with approximation
template <typename T>
Mat32r JacfunctionApp(Vector2r para, Vector2r paraDelta, T *particle1, T *particle2,Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag)
{
    Vector3r p_f0,p_f1,p_f2;
    Vector2r para0;
    Disfunction(para,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f0 = p2 - p1;//f(p) = p2 - p1, and
    para0 = Vector2r(para(0)+paraDelta(0),para(1));
    Disfunction(para0,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f1 = p2 - p1;
    para0 = Vector2r(para(0),para(1)+paraDelta(1));
    Disfunction(para0,particle1, particle2,position1, position2, p1, p2, stop_flag);
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

template <typename T>
Mat32r JacfunctionApp2(Vector2r para, Vector2r paraDelta, T *particle1, T *particle2,Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag)
{
    Vector3r p_f0,p_f1,p_f2;
    Vector2r para0;
    //Disfunction(para,center,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f0 = p2 - p1;//f(p) = p2 - p1, and
    para0 = Vector2r(para(0)+paraDelta(0),para(1));
    Disfunction(para0,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f1 = p2 - p1;
    para0 = Vector2r(para(0),para(1)+paraDelta(1));
    Disfunction(para0,particle1, particle2,position1, position2, p1, p2, stop_flag);
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
  void (*func)(Vector2r p, T *particle1, T *particle2,
		   Vector3r position1, Vector3r position2,Vector3r &p1,Vector3r &p2, bool &stop_flag), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  //Mat32r (*jacf)(Vector2r p, T *particle1, T *particle2),  /* function to evaluate the Jacobian \part x / \part p */
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
Vector3r &p1,Vector3r &p2, bool &touching, Vector3r &contact
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

	//cout<<"Dp="<<Dp<<endl;
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
  	(*func)(p, particle1, particle2,position1, position2,p1,p2, stop_flag); nfev=1;
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
        jac0 = JacfunctionApp(p, Dp, particle1, particle2,position1, position2,p1,p2, stop_flag);nfev+=2;
				njev++;
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
			//cout<<"Dp="<<Dp<<"Dp_L2"<<Dp_L2<<endl;
        	//cout<<"pDp="<<pDp(0)<<"\t"<<pDp(1)<<endl;

			//std::cout<<"solved --- pDp-"<<pDp[0]<<"\t"<<pDp[1]<<std::endl;
        	//if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        	if(Dp_L2 < 1e-3){//if(Dp_L2<=eps2*p_L2){ //relative change in p is small, stop
					//if(Dp_L2 < 1e-5){//the relative change in p is small, and switch to use gjk
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

        	 (*func)(pDp, particle1, particle2,position1, position2, p1, p2, stop_flag); ++nfev; /* evaluate function at p + Dp */
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

//***********************************************************************************/
//GJK

template <typename T>
bool getPenetration(T *particle1, T *particle2, double a_margin, double b_margin,
										Vector3r position1, Vector3r position2,
										Vector3r& v, Vector3r& pa, Vector3r& pb,int &num_iter)
{
	SD_Scalar margin = min(a_margin,b_margin);
	int num_iterations;
	//cout<<"hybrid_depth:into "<<endl;
	if (margin > SD_Scalar(0.0))
	{
		SD_Scalar margin2 = margin * margin;

		DT_GJK gjk;

#ifdef STATISTICS
		num_iterations = 0;
#endif
		SD_Scalar dist2 = Mathr::MAX_REAL;
		//cout<<"using margin"<<endl;
		//cout<<"v="<<v.normalized()<<endl;
		do
		{
			//
			Vector3r erosion = v.normalized()*margin;
			#ifdef STATISTICS
			cout<<"inside angle0="<<v.dot(pb-pa)/(v.norm()*(pb-pa).norm())<<"erosion="<<erosion.norm()<<endl;
			#endif
			Vector3r  p = particle1->support(-v) + position1;//support point of object A
			Vector3r  q = particle2->support(v) + position2;//support point of object B
			#ifdef STATISTICS
			double angle1 = v.dot(p-q)/(v.norm()*(p-q).norm());//phi in [0, pi]
			cout<<"angle1 inside "<<angle1<<endl;
			#endif
			//cout<<"p="<<p<<"pa="<<pa<<endl;
			//cout<<"v0="<<v.normalized()<<endl;
			//cout<<"p="<<p<<"q="<<q<<endl;
			#ifdef STATISTICS
			++num_iterations;
			#endif
			#ifdef NUM_ITERATION
			++num_iter;
			#endif
			Vector3r w = p - q;//support point of Minkowski difference A-B
			// SD_Scalar vn = v.norm();
			//
			// if (vn > 0.0)
			// {
			//  vn = 1/ vn;
			//  w += v*vn*margin;
			// }
			Vector3r v1 = v;
			SD_Scalar delta = v.dot(w);
			if(delta > 0.0){//not contact
				#ifdef STATISTICS
				cout<<"nfuc1="<<num_iterations<<endl;
				#endif
				return false;
			}

			double angle = delta/(v.norm()*w.norm());//phi in [0, pi]
			//cout<<"inside angle"<<angle<<"wm"<<w.norm()/margin<<endl;
			if(fabs(angle) > 0.95 && w.norm()<2.0*margin)
			{
				pa = p;
				pb = q;
				#ifdef STATISTICS
				cout<<"nfuc22="<<num_iterations<<endl;
				#endif
				return true;

			}

			p += erosion;
			q -= erosion;
			w = p - q;
			/*if (delta > SD_Scalar(0.0) && delta * delta > dist2 * margin2)//the enlarged objects do not interact.
			{
				return false;
			}*/

			//cout<<"gjkw)="<<gjk.inSimplex(w)<<endl;
			//if (gjk.inSimplex(w) || dist2 - delta <= dist2 * 1e-3)
			if (gjk.inSimplex(w) || dist2 - v.dot(w) <= dist2 * 1e-3)
			//if (gjk.inSimplex(w) || fabs(v.normalized().dot(w)) >= 0.95*w.norm())
			{ //the objects only interact in the margin

				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= erosion;//v * (margin / s);
				pb += erosion;//v * (margin / s);
				#ifdef STATISTICS
				//w = pa - pb;
				//v = pb - pa;
				//cout<<"pa="<<pa<<"pb="<<pb<<endl;
				cout<<"gjk.inSimplex(w)="<<gjk.inSimplex(w)<<endl;

				cout<<"p="<<p<<"pa="<<pa+erosion<<"pb"<<pb-erosion<<endl;
				cout<<"nfuc2="<<num_iterations<<"p/m="<<(pa-pb).norm()/margin<<"angle"<<v.dot(w)/(v.norm()*w.norm())<<"v="<<v.normalized()<<endl;
				#endif
				return true;
			}

			gjk.addVertex(w, p, q);//add a vertex to the convex hull

      if (gjk.isAffinelyDependent())
      {
        gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= erosion;//v * (margin / s);
				pb += erosion;//v * (margin / s);
				#ifdef STATISTICS
				cout<<"nfuc3="<<num_iterations<<"p/m="<<(pa-pb).norm()/margin<<endl;
				#endif
				return true;
      }


			//cout<<"v1="<<v<<endl;
			if (!gjk.closest(v))
			{
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				erosion = v.normalized()*margin;
				pa -= erosion;//v * (margin / s);
				pb += erosion;//v * (margin / s);
				#ifdef STATISTICS
				cout<<"nfuc4="<<num_iterations<<"p/m="<<(pa-pb).norm()/margin<<endl;
				#endif
				return true;
			}
			//cout<<"v2="<<v<<endl;
			#ifdef SAFE_EXIT
			SD_Scalar prev_dist2 = dist2;
			#endif

			dist2 = v.squaredNorm();

			#ifdef SAFE_EXIT
			//if (prev_dist2 - dist2 <= (Mathr::EPSILON * prev_dist2))
			if (prev_dist2 - dist2 <= 1E-3 * prev_dist2)
			{
  			gjk.backup_closest(v);
				dist2 = v.squaredNorm();
				gjk.compute_points(pa, pb);
				//SD_Scalar s = sqrt(dist2);
				//assert(s > SD_Scalar(0.0));
				erosion = v.normalized()*margin;
				pa -= erosion;//v * (margin / s);
				pb += erosion;//v * (margin / s);
				#ifdef STATISTICS
				//cout<<"gjk.inSimplex(w)="<<gjk.inSimplex(w)<<endl;
				p = particle1->support(-v) + position1;
				cout<<"p="<<p<<"pa="<<pa<<endl;
				cout<<"angle"<<v.dot(w)/(v.norm()*w.norm())<<"v="<<v.normalized()<<endl;
				cout<<"angle1="<<v.dot(pb-pa)/(v.norm()*(pb-pa).norm())<<endl;
				cout<<"nfuc5="<<num_iterations<<"p/m="<<(pa-pb).norm()/margin<<endl;
				#endif
				return true;
			}
			#endif
		}
		while (!gjk.fullSimplex() && dist2 > 1e-16 * gjk.maxVertex());
	}
	//cout<<"second GJK phase"<<endl;
	// Second GJK phase. compute points on the boundary of the offset object
//	return penetration_depth(a, b, v, pa, pb);
	#ifdef STATISTICS
	cout<<"nfuc6="<<num_iterations<<"p/m="<<(pa-pb).norm()/margin<<endl;
	#endif
	return false;

}

//**********************************************************************************
/*! Create PolySuperellipsoid (collision geometry) from colliding PolySuperellipsoids. */

bool Ig2_PolySuperellipsoid_PolySuperellipsoid_PolySuperellipsoidGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){
    //std::cerr << "watch 1"<< "\n";
	//get polyhedras
	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
    //std::cerr << "watch 2"<< "\n";
	PolySuperellipsoid* A = static_cast<PolySuperellipsoid*>(cm1.get());
	PolySuperellipsoid* B = static_cast<PolySuperellipsoid*>(cm2.get());
  Vector3r position1 = se31.position;
  Vector3r position2 = se32.position + shift2;
        //cout<<"watch point"<<endl;
	bool isNew = !interaction->geom;
	//std::cerr << "pos1:" <<position1 << "\n";
	//std::cerr << "pos2:" <<position2 << "\n";

	//Vector3r p = A->getPosition();
	//std::cerr << "Gl1_PolySuperellipsoid:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";
        //A->rot_mat2local = (se31.orientation).conjugate().toRotationMatrix();//to particle's system
	//A->rot_mat2global = (se31.orientation).toRotationMatrix(); //to global system
	//B->rot_mat2local = (se32.orientation).conjugate().toRotationMatrix();//to particle's system
	//B->rot_mat2global = (se32.orientation).toRotationMatrix(); //to global system
	shared_ptr<PolySuperellipsoidGeom> bang;
	Vector2r dels(0.01,0.01);
	Vector2r para(1e-5,1e-5);
    //std::cerr << "watch 3"<< "\n";

	if (isNew) {
		// new interaction
		//cout<<"new intersection"<<endl;
		bang=shared_ptr<PolySuperellipsoidGeom>(new PolySuperellipsoidGeom());
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
						if(para(0)<0) para(0) += Mathr::TWO_PI;
            //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
            para(1) = atan2(d0(2),sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
            para(0) += 1e-5;para(1) += 1e-5;
            //cout<<"new para="<<para<<endl;
            //to do ...
        #endif
	}else{
		// use data from old interaction
    bang=SUDODEM_PTR_CAST<PolySuperellipsoidGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
    para = bang->contactAngle;
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
	    bool touching(true);
			#ifdef NUM_ITERATION
			int jac_num=0; int func = 0;
			#endif
	    double nm_info[2];
	    double opts[5], info[10];
        //cout<<"para="<<para<<endl;
	    /* optimization control parameters */
	    //opts[0]=0.01; opts[1]=1E-15; opts[2]=1E-5; opts[3]=1E-20;opts[4]=1e-6;
        opts[0]=1.0; opts[1]=1E-15; opts[2]=1E-3; opts[3]=1E-20;opts[4]=1e-6;


	    /* invoke the optimization function */
	    //lmder_z(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2, globalContact,center,touching);
			//test GJK
			#ifdef USE_GJK
			double a_margin(0.05e-3),b_margin(0.1e-3);
			Vector3r v(0.0,0.0,0.0);//separating axis vector
			if (isNew) {
				v = position1 - position2;
			}else{
				//v = position2 - position1;//cache v later v = p1 - p2;
				v(0) = -cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
				v(1) = -sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
				v(2) = -sin(para(1));//para(0):[0,2pi];para(1):[-pi/2,pi/2]
				p1 = bang->point1;
				p2 =	bang->point2;
				v = p2 - p1;
				cout<<"outside0 v="<<v.normalized()<<"v="<<v<<endl;
			}
			int num_iter=0;
			touching = getPenetration(A,B,a_margin,b_margin, position1, position2, v, p1, p2, num_iter);
			para(0) = atan2(-v(1),-v(0));//[-pi,pi]
			para(0) = (para(0)<0)?para(0)+Mathr::TWO_PI:para(0);//[0,2pi]
			para(1) = atan2(-v(2),sqrt(v(0)*v(0)+v(1)*v(1)));//[-pi/2,pi/2]
			//debug

			cout<<"outside v="<<v.normalized()<<"p2-p1="<<(p2-p1).normalized()<<endl;
			Vector3r p = A->support(-(p1-p2)) + position1;
			Vector3r p3 = B->support(p1-p2) + position2;
			cout<<"d1="<<(p3-p).norm()<<"d2="<<(p2-p1).norm()<<endl;
			Vector3r p0 = A->rot_mat2global*A->Normal2SurfaceMC(A->rot_mat2local*(-v)) + position1;
			cout<<"outside p"<<p<<"p1"<<p1<<"p2"<<p2<<endl;
			cout<<"outside p0"<<p0<<"p3"<<p3<<endl;
			double angle1 = (p1-p2).dot(p3-p)/((p1-p2).norm()*(p3-p).norm());//phi in [0, pi]
			cout<<"angle1 outside "<<angle1<<endl;

			//test GJK ends
			#else
        lmder2(Disfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,touching,contactNormal);
				#ifdef NUM_ITERATION
				jac_num += info[8];
				func += info[7];
				#endif
				//if(info[6]==2){
				//	printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);

				//}
				//need to check if the norm of p1 - p2 is large enough. If it is, then we do GJK otherwise not
				if(info[7]>20 && info[6]==2){//if(false){//info[6]==2){//too small the step is
					//double a_margin(0.05e-3),b_margin(0.05e-3);
					//double angle = (p2-p1).dot(contactNormal)/(p2-p1).norm();
					//cout<<"angle="<<angle<<endl;
					//cout<<"flat face..."<<interaction->id1<<" "<<interaction->id2<<endl;
					double margin = 0.1*min(A->getr_max(),B->getr_max());
					Vector3r v;
					v(0) = -cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
					v(1) = -sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
					v(2) = -sin(para(1));
					int gjk_num_iter = 0;
					getPenetration(A,B,margin,margin, position1, position2, v, p1, p2, gjk_num_iter);
					#ifdef NUM_ITERATION
					//jac_num += info[7];
					func += gjk_num_iter;
					#endif
					//v = p2 - p1;
					contactNormal = -v.normalized();
					para(0) = atan2(-v(1),-v(0));//[-pi,pi]
					para(0) = (para(0)<0)?para(0)+Mathr::TWO_PI:para(0);//[0,2pi]
					para(1) = atan2(-v(2),sqrt(v(0)*v(0)+v(1)*v(1)));//[-pi/2,pi/2]
				}
				#ifdef STATISTICS
				printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
				#endif

        if(isnan(para[0]) || isnan(para[1]) || 7==info[6]){
            //FIXME: .
            if(touching){
                cout<<"Para or d is NAN. RESET para."<<endl;
                Vector3r d0 = position2 - position1;
                //atan2(y,x): arc tangent of y/x
                para(0) = atan2(d0(1),d0(0));
								if(para(0)<0) para(0) += Mathr::TWO_PI;
                //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
                para(1) = atan2(d0(2),sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[-pi/2,pi/2]
                para(0) += 1e-5;para(1) += 1e-5;
                //may be not safe
                lmder2(Disfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,touching,contactNormal);
                //std::cerr << "id1="<<interaction->getId1()<<"id2="<<interaction->getId2()<<"\n";
                //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
                //scene->stopAtIter = scene->iter + 1;
								//
            }
        }
				//cout<<"touching!"<<touching<<endl;
				pDepth = (p1-p2).dot(contactNormal);
        if(touching){
            //check if the contact depth is not too small.
						//cout<<"p2p1"<<10.0*(p2-p1).norm()<<"rmax"<< A->getr_max()<<endl;
            if(10.0*pDepth > A->getr_max()){//if(10.0*(p2-p1).norm() > A->getr_max()){//the contact depth is tremendously large which might be at local minimia. The common contact depth should be 3 order of magnitude less than particle size.FIXME:the coefficient of 10.0 may be not the best try.
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
										if(para(0)<0) para(0) += Mathr::TWO_PI;
                    //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
                    para(1) = atan2(d0(2),sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[-pi/2,pi/2]
                    para(0) += 1e-5;para(1) += 1e-5;
                    //may be not safe
                    lmder2(Disfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,touching,contactNormal);
										pDepth = (p1-p2).dot(contactNormal);
										if(10.0*pDepth > A->getr_max()){//give up
											bang->contactAngle = para;
											bang->PenetrationDepth=0;
											bang->isShearNew = true;
											return true;
										}
										#ifdef NUM_ITERATION
										jac_num += info[8];
										func += info[7];
										#endif
										cout<<"recalc depth "<<interaction->id1<<" "<<interaction->id2<<endl;
										/*cout<<A->isInside(p2-position1)<<" "<<B->isInside(p1-position2)<<endl;
										Vector3r pp1 = A->rot_mat2local*(p2-position1)+A->getMassCenter();
										Vector6r rxyz = A->getrxyz();
										Vector2r eps = A->geteps();
										double eps1 = eps[0];double eps2 = eps[1];
								    double rr= pow( pow(fabs(pp1(0)/rxyz[(pp1[0]>0?0:1)]),2.0/eps1) + pow(fabs(pp1(1)/rxyz[(pp1[1]>0?2:3)]),2.0/eps1), eps1/eps2) + pow(fabs(pp1(2)/rxyz[(pp1[2]>0?4:5)]),2.0/eps2);
										cout <<"rr="<<rr<<endl;
										*/
										/*
										Vector3r v(0.0,0.0,0.0);
										Vector3r v2(0.0,0.0,0.0);
										v2 = position1 - position2;
										double margin = 0.1*min(A->getr_max(),B->getr_max());
										//double a_margin(0.1e-3),b_margin(0.1e-3);
										getPenetration(A,B,margin,margin, position1, position2, v2, p1, p2, touching);
										v = p1 - p2;
										para(0) = atan2(v(1),v(0));
										if(para(0)<0) para(0) += Mathr::TWO_PI;
										para(1) = atan2(v(2),sqrt(v(0)*v(0)+v(1)*v(1)));//[0,pi]
										contactNormal = -v2.normalized();
										pDepth = v.dot(contactNormal);
										*/
                }
            }
        }
				#endif //USE_GJK ends
	    //contactNormal = -contactNormal;//make contact normal along the direction from particle 1 to particle 2.
        //#ifdef _DEBUG_LM_NM
	    //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
        //#endif
	    //printf("Best fit parameters: %.7g %.7g \n", para(0), para(1));


        //std::cerr << "touching: " <<touching<< "\n";
        //std::cerr << "watch 4.1"<< "\n";
				/*if(!touching){
					//sprintf(filename,"%d",gettid4());
					FILE * fin = fopen("ruleout.dat","a");
					fprintf(fin,"%g %g\n", info[7], info[8]);//nfuc and jacfunc
					//fprintf(fin,"depth\t%e\n",depth);
					fclose(fin);
				}*/
			//output iterations for test
			#ifdef NUM_ITERATION
				if(!isNew){//record the second iteration
				//sprintf(filename,"%d",gettid4());
				FILE * fin = fopen("exactcomp.dat","a");
				fprintf(fin,"%d %d\n", func, jac_num);//nfuc and jacfunc
				//fprintf(fin,"depth\t%e\n",depth);
				fclose(fin);
			}
			#endif
	    if (!touching){
	            bang->contactAngle = para;
	            bang->PenetrationDepth=0;
	            bang->isShearNew = true;
	            return true;}

            bang->point1 = p1;
            bang->point2 = p2;

        //test jacobian matrix ends
        //std::cerr << "watch 6"<< "\n";
	    //contactNormal = p1 - p2;//using the outward normal
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
            const char* filename = "PolySuperellipsoidLaw_err.log";
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
					bang->contactPoint=0.5*(p1+p2);
					//correct the penetration
					//pDepth < 0?
					if(pDepth<0){
						//cout<<"error:pDepth<0 "<<interaction->id1<<" "<<interaction->id2<<endl;
						bang->contactAngle = para;
						bang->PenetrationDepth=0;
						bang->isShearNew = true;
						return true;
					}
	        bang->PenetrationDepth=pDepth;//(p1-p2).dot(contactNormal);//contactNormal.norm();


            //testing
            //SHOWN NORMAL BEGINS
            #ifdef GEOM_GL
            Vector3r contact(0.,0.,0.);
            Vector2r phi;
		        contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
		        contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
		        contact(2) = sin(para(1));
						Vector3r v2 = p2 - p1;

		        phi= A->Normal2Phi(A->rot_mat2local*contact);
            Vector3r n1,n2;
            n1 = A->getNormal(phi);
            n1.normalize();
						double angle = v2.dot(A->rot_mat2global*n1)/(v2.norm());//phi in [0, pi]
						//cout<<"angle outside "<<angle<<endl;
						//Vector3r pp = A->support(-v) + position1;
            bang->point11 = A->rot_mat2global*n1*A->getr_max()+p1;

            phi= B->Normal2Phi(B->rot_mat2local*contact*(-1.0));
            n2 = B->getNormal(phi);
            n2.normalize();
            //cout<<"r_max="<<B->getr_max()<<endl;
            //cout<<"n2="<<n2.norm()<<endl;
            //Vector3r tmp1;
            //tmp1 = B->rot_mat2global*n2;
            //cout<<"tmp1="<<tmp1<<" magnitude: "<<tmp1.norm()<<endl;
						//Vector3r pp2 = B->support(v) + position2;
		        bang->point22 = B->rot_mat2global*n2*B->getr_max()+p2;
						#else
						bang->point11 = p1;
						bang->point22 = p2;
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

	return true;
}
//**********************************************************************************
/*! Create PolySuperellipsoid (collision geometry) from colliding Polyhedron and Wall. */

bool Ig2_Wall_PolySuperellipsoid_PolySuperellipsoidGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){

	const int& PA(cm1->cast<Wall>().axis);
	const int& sense(cm1->cast<Wall>().sense);
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	PolySuperellipsoid* B = static_cast<PolySuperellipsoid*>(cm2.get());

	bool isNew = !interaction->geom;
	shared_ptr<PolySuperellipsoidGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<PolySuperellipsoidGeom>(new PolySuperellipsoidGeom());
		bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<PolySuperellipsoidGeom>(interaction->geom);
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
	    Vector3r point2 = rot_mat2*( B->getSurfaceMC(phi)) + B_pos;	//the closest point on particle B to the wall
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
/*! Create PolySuperellipsoid (collision geometry) from colliding Polyhedron and Facet. */

bool Ig2_Facet_PolySuperellipsoid_PolySuperellipsoidGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){


	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
	Facet*   A = static_cast<Facet*>(cm1.get());
	PolySuperellipsoid* B = static_cast<PolySuperellipsoid*>(cm2.get());

	bool isNew = !interaction->geom;

	shared_ptr<PolySuperellipsoidGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<PolySuperellipsoidGeom>(new PolySuperellipsoidGeom());
		//bang->sep_plane.assign(3,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<PolySuperellipsoidGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}
/*
	//find intersection PolySuperellipsoid
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
