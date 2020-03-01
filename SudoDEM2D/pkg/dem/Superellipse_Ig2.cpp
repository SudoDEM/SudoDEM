#include "Superellipse.hpp"
#include "Superellipse_Ig2.hpp"
#include <cmath>
//#include <sys/syscall.h>
#define _USE_MATH_DEFINES

#define _USE_LM_OPTIMIZATION
#define _DEBUG_LM_NM1
//LM: Levenberg-Marquardt method
//NM: Nerld-Mead simplex method

#include<time.h>

SUDODEM_PLUGIN(/* self-contained in hpp: */ (Ig2_Superellipse_Superellipse_SuperellipseGeom)
		(Ig2_Wall_Superellipse_SuperellipseGeom)
		(Ig2_Fwall_Superellipse_SuperellipseGeom)
		//(Ig2_Facet_Superellipse_SuperellipseGeom)
	   );
CREATE_LOGGER(Ig2_Fwall_Superellipse_SuperellipseGeom);
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
/*
Matrix3r CRotationMat(Vector2r position1, Vector2r position2)//positions of particle 1 and particle 2.
{
	//Matrix3r globalContact;
	Vector2r e1,n1,x1(1.,0,0);
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
*/


/*****************************************************/
template <typename T>
void Disfunction(Real para,Vector2r &center, T *particle1, T *particle2,
		     Vector2r position1, Vector2r position2,Vector2r &p1,Vector2r &p2, bool &stop_flag)
{
	Vector2r d;
	//stop_flag = false;
	//contact
	Vector2r contact(cos(para),sin(para)),contact1(0.,0.),contact2(0.,0.);
  double phi(0.0);
	/////////////////continue
	//point 1 on Particle 1

	//std::cout<<"contact normal in global sys: "<<contact<<std::endl;
  //no need to use vector dot or product.
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

//secent Jacobian matrix with approximation
template <typename T>
Vector2r JacfunctionApp(Real para, Real paraDelta, Vector2r &center, T *particle1, T *particle2,Vector2r position1, Vector2r position2,Vector2r &p1,Vector2r &p2, bool &stop_flag)
{
    Vector2r p_f0,p_f1,p_f2;
    Real para0;
    Disfunction(para,center,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f0 = p2 - p1;//f(p) = p2 - p1, and
    para0 = para + paraDelta;
    Disfunction(para0,center,particle1, particle2,position1, position2, p1, p2, stop_flag);
    p_f1 = p2 - p1;

	return (p_f1 - p_f0)/paraDelta;
}

template <typename T>	//return distance
void lmder2(
  void (*func)(Real p, Vector2r &center, T *particle1, T *particle2,
		   Vector2r position1, Vector2r position2,Vector2r &p1,Vector2r &p2, bool &stop_flag), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  //Vector2r (*jacf)(Real p, T *particle1, T *particle2),  /* function to evaluate the Jacobian \part x / \part p */
  Real &p,         /* I/O: initial parameter estimates. On output has the estimated solution */
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
Vector2r position1, Vector2r position2,
Vector2r &p1,Vector2r &p2, Vector2r &center, bool &touching, Vector2r &contact
)
{
	register int k;
	touching = true;
	Vector2r e, p_f(0.,0.), pDp_f(0.,0.),deltaF,hx(0.,0.);
	Real jacTf, diag_jacTjac(0);
	Vector2r jac;
	Real jacTjac, Dp(1e-7), pDp;
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
  //  cout<<"para0="<<p<<endl;
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
        //Vector2r contact(cos(p),sin(p));
				contact(0) = cos(p);
				contact(1) = sin(p);
        //double angle = acos(p_f.dot(contact)/(p_f.norm()*contact.norm()));//phi in [0, pi]
        double angle = p_f.dot(contact)/(p_f.norm());//phi in [0, pi]
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
        Vector2r jac0,jac1;
        jac0 = JacfunctionApp(p, Dp,center, particle1, particle2,position1, position2,p1,p2, stop_flag);nfev+=2;
        //jac0 = JacfunctionApp(p, Vector2r(1e-7,1e-7),center, particle1, particle2,position1, position2,p1,p2, stop_flag);nfev+=3;

        //jac1 = JacfunctionSim(p, particle1, particle2);
		//cout<<"Jacobian mat="<<jac<<endl;
        //cout<<"Approx. Jacobian mat="<<jac0<<endl;
        jac = jac0;
        //cout<<"Jacobian mat simplized ="<<jac1<<endl;
    	/* J^T J, J^T f */
	 	jacTjac = jac.squaredNorm();//jac.transpose()*jac;
	 	jacTf = jac.dot(-1.0*p_f);//jac.transpose()*(-1.0*p_f);//jac.transpose()*p_f;
 		//cout<<"jacTf="<<jacTf<<endl;
		//cout<<"p2"<<p<<endl;
	  	/* Compute ||J^T f||_inf and ||p||^2 */
		jacTf_inf = jacTf;
		p_L2 = fabs(p);
		//cout<<"jacTf_inf"<<jacTf_inf<<endl;
    	/* check for convergence */
    	/*if((jacTf_inf <= eps1)){
      		Dp_L2=0.0; // no increment for p in this case
      		stop=1;
      		break;
    	}*/
		deltaF = pDp_f - p_f;
		pDp_feL2 = fabs(pDp_f.norm() - p_f.norm());
		/*if(pDp_feL2 <=eps2*p_fL2)
		{
			stop=8;
			break;
		}*/
   		/* compute initial damping factor */
    	if(k==0){
			tmp = jacTjac;//find max diagonal element
      		mu=tau*tmp;
    	}
		diag_jacTjac = jacTjac;
    	/* determine increment using adaptive damping */
   	 	while(1){
      		/* augment normal equations */
	  		//jacTjac += mu*Matrix2d::Identity();
			jacTjac +=mu;
      		/* solve augmented equations */
			//cout<<"jacTjac= "<<jacTjac<<endl;
	  		//Dp = jacTjac.inverse()*jacTf;//jacTjac.inverse()*(-jacTf);//caution:sigular
        //    double n_a = jacTjac(0,0)*jacTjac(1,1) - jacTjac(0,1)*jacTjac(0,1);
        //    Dp(0) = (jacTf(0)*jacTjac(1,1) - jacTf(1)*jacTjac(0,1))/n_a;
        //    Dp(1) = (jacTf(1)*jacTjac(0,0) - jacTf(0)*jacTjac(0,1))/n_a;
        Dp = jacTf/jacTjac;
 			//cout<<"Dp="<<Dp(0)<<"\t"<<Dp(1)<<endl;

			Dp_L2 = fabs(Dp);

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
            pDp_L2 = fabs(pDp);
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


			dL = Dp*(mu*Dp + jacTf);//0.5*Dp.transpose()*(mu*Dp - jacTf);
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
			jacTjac = diag_jacTjac;//restore diagonal J^T J entries

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

//**********************************************************************************
/*! Create Superellipse (collision geometry) from colliding Superellipses. */

bool Ig2_Superellipse_Superellipse_SuperellipseGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){
    //std::cerr << "watch 1"<< "\n";
	const Se2r& se21=state1.se2;
	const Se2r& se22=state2.se2;
	Vector2r position1 = se21.position;
	Vector2r position2 = se22.position + shift2;
	//cout<<"shift2="<<shift2<<endl;
	//cout<<"pos1="<<position1<<"pos2="<<position2<<endl;
    //std::cerr << "watch 2"<< "\n";
	Superellipse* A = static_cast<Superellipse*>(cm1.get());
	Superellipse* B = static_cast<Superellipse*>(cm2.get());
        //cout<<"watch point"<<endl;
	bool isNew = !interaction->geom;
	//std::cerr << "pos1:" <<position1 << "\n";
	//std::cerr << "pos2:" <<position2 << "\n";

	//Vector2r p = A->getPosition();
	//std::cerr << "Gl1_Superellipse:" <<p(0)<<"\t"<<p(1)<<"\t"<<p(2) << "\n";
        //A->rot_mat2local = (se21.orientation).conjugate().toRotationMatrix();//to particle's system
	//A->rot_mat2global = (se21.orientation).toRotationMatrix(); //to global system
	//B->rot_mat2local = (se22.orientation).conjugate().toRotationMatrix();//to particle's system
	//B->rot_mat2global = (se22.orientation).toRotationMatrix(); //to global system
	shared_ptr<SuperellipseGeom> bang;
	Real dels(0.01), para(1e-5);
    //std::cerr << "watch 3"<< "\n";

	if (isNew) {
		// new interaction
		//cout<<"new intersection"<<endl;
		bang=shared_ptr<SuperellipseGeom>(new SuperellipseGeom());
		//bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector2r(0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
		dels = 0.1;
        para = bang->contactAngle;
        //cout<<"para0="<<para<<endl;
        #ifdef _USE_LM_OPTIMIZATION
            //cout<<"new geom!..."<<endl;
            Vector2r d0 = position2 - position1;
            //atan2(y,x): arc tangent of y/x
            para = atan2(d0(1),d0(0));
            //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
            //para(1) = atan(d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
            para += 1e-5;
            //cout<<"new para="<<para<<endl;
            //to do ...
        #endif
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperellipseGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
        para = bang->contactAngle;
        //cout<<"para="<<para<<endl;
		//Vector2r v1 = state1.vel;
		//Vector2r v2 = state2.vel;//maybe a bug
		//cout<<"vvv1"<<v1(0)<<" "<<v1(1)<<" "<<v1(2)<<"  "<<v1.isZero(0)<<endl;
		//cout<<"vvv2"<<v2(0)<<" "<<v2(1)<<" "<<v2(2)<<"  "<<v2.isZero(0)<<endl;

		//if (v1.isZero(0)&&v2.isZero(0)){
		//        bang->PenetrationDepth=0;

               //         return true;
		//}
		dels = 0.001;
	}
	Vector2r contactNormal(0,0), center(0,0);
	double pDepth(0.0);
	Vector2r p1,p2;
	//if the two particles are spherical
	if (A->isSphere && B->isSphere){
	        //particles are spherical
	        Vector2r dist = position2 -position1;
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


	    //Levenberg-Marquadt method
        //para(0) = 0.001;para(1)=0.001;
	    double opts[5], info[10];
        //cout<<"para="<<para<<endl;
	    /* optimization control parameters */
	    //opts[0]=0.01; opts[1]=1E-15; opts[2]=1E-5; opts[3]=1E-20;opts[4]=1e-6;
        opts[0]=1.0; opts[1]=1E-15; opts[2]=1E-3; opts[3]=1E-20;opts[4]=1e-6;


	    /* invoke the optimization function */
	    //lmder_z(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2, globalContact,center,touching);

        //lmder2(Disfunction, Jacfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching);
        lmder2(Disfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching, contactNormal);
        if(isnan(para) || 7==info[6]){
            //FIXME: .
            if(touching){
                cout<<"Para or d is NAN. RESET para."<<endl;
                Vector2r d0 = position2 - position1;
                //atan2(y,x): arc tangent of y/x
                para = atan2(d0(1),d0(0));
                //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
                //para(1) = atan(d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
                para += 1e-5;
                //may be not safe
                lmder2(Disfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching, contactNormal);
                //std::cerr << "id1="<<interaction->getId1()<<"id2="<<interaction->getId2()<<"\n";
                //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
                //scene->stopAtIter = scene->iter + 1;
            }
        }
        //cout<<"para1="<<para<<endl;
				pDepth = (p1-p2).dot(contactNormal);
        if(touching){
            //check if the contact depth is not too small.
            if(10.0*(p2-p1).norm() > A->getr_max()){//the contact depth is tremendously large which might be at local minimia. The common contact depth should be 3 order of magnitude less than particle size.FIXME:the coefficient of 10.0 may be not the best try.
                if(!(A->isInside(p2-position1)) || !(B->isInside(p1-position2))){
                    //FIXME:the inside-outside function may not be sufficiently accurate for tiny contact depth during computation.
                    //cout<<"sure? not inside A or B?"<<endl;
                    //cout<<"local minimia!"<<endl;
                    //cout<<"contact depth="<<(p2-p1).norm()<<endl;
                    //std::cerr << "id1="<<interaction->getId1()<<"id2="<<interaction->getId2()<<"\n";
                    //scene->stopAtIter = scene->iter + 1;
                    Vector2r d0 = position2 - position1;
                    //atan2(y,x): arc tangent of y/x
                    para = atan2(d0(1),d0(0));
                    //cout<<"costheta="<<d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)+d0(2)*d0(2))<<endl;
                    //para(1) = atan(d0(2)/sqrt(d0(0)*d0(0)+d0(1)*d0(1)));//[0,pi]
                    para += 1e-5;
                    //may be not safe
                    lmder2(Disfunction, para, 50, opts, info,A,B, position1, position2, p1, p2,center,touching, contactNormal);
										pDepth = (p1-p2).dot(contactNormal);
                }
            }
        }
	    //contactNormal = -contactNormal;//make contact normal along the direction from particle 1 to particle 2.
        //#ifdef _DEBUG_LM_NM
	    //printf("Levenberg-Marquardt returned in %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
        //#endif
	    //printf("Best fit parameters: %.7g \n", para);


        //std::cerr << "touching: " <<touching<< "\n";
        //std::cerr << "watch 4.1"<< "\n";

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
	    //Vector2r contactNormal1 = contactNormal/norm;
	    //cout<<"contact normal111"<<contactNormal1(0)<<" "<<contactNormal1(1)<<" "<<contactNormal1(2)<<endl;
	    //std::cerr << "watchpoint 2: " << "\n";

	    //if((se22.position-centroid).dot(normal)<0) normal*=-1;
            bang->point1 = p1;
            bang->point2 = p2;
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
            const char* filename = "SuperellipseLaw_err.log";
            //sprintf(filename,"%d",gettid4());
            FILE * fin = fopen(filename,"a");
            fprintf(fin,"\n*******%s****NAN CONTACT NORMAL***\n",tmp);
            fprintf(fin,"LM %g iter, reason %g, nfuc %g n jacfunc %g\n", info[5], info[6], info[7], info[8]);
            fprintf(fin,"id1= %d id2= \t%d\n",interaction->getId1(), interaction->getId2());

            fprintf(fin,"point1\t%e\t%e\t%e point2\t%e\t%e\t%e\n",p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);

            fprintf(fin,"contact normal \t%e\t%e\t%e\n",contactNormal[0],contactNormal[1],contactNormal[2]);
            fprintf(fin,"para\t%e\n",para);
            //fprintf(fin,"depth\t%e\n",depth);
            fclose(fin);
        }


	        bang->contactPoint=center;
	        bang->PenetrationDepth= pDepth;//contactNormal.norm();


            //testing
            //SHOWN NORMAL BEGINS
            #define GEOM_GL
            #ifdef GEOM_GL
            //cout<<"para="<<para<<" pos="<<position1<<endl;
            Vector2r contact(cos(para),sin(para));
            Real phi;
						//cout<<"p1 and p2="<<p1<<" "<<p2;
	        phi= A->Normal2Phi(A->rot_mat2local*contact);
            Vector2r n1,n2;
            //cout<<"phi="<<phi<<endl;
            n1 = A->getNormal(phi);
            //p1 = A->getSurface(phi) + position1;
            //bang->point1 = p1;
            n1.normalize();
            bang->point11 = A->rot_mat2global*n1*A->getr_max()+p1;
            //cout<<"n1="<<n1<<"rot1"<<A->rot_mat2local<<endl;
            //cout<<"surf="<<A->getSurface(phi)<<endl;
            phi= B->Normal2Phi(B->rot_mat2local*contact*(-1.0));
            n2 = B->getNormal(phi);
            n2.normalize();
            //cout<<"r_max="<<B->getr_max()<<endl;
            //cout<<"n2="<<n2.norm()<<endl;
            //Vector2r tmp1;
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
/*! Create Superellipse (collision geometry) from colliding Polyhedron and Wall. */

bool Ig2_Wall_Superellipse_SuperellipseGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){

	const int& PA(cm1->cast<Wall>().axis);
	const int& sense(cm1->cast<Wall>().sense);
	const Se2r& se21=state1.se2; const Se2r& se22=state2.se2;
	Superellipse* B = static_cast<Superellipse*>(cm2.get());

	bool isNew = !interaction->geom;
	shared_ptr<SuperellipseGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<SuperellipseGeom>(new SuperellipseGeom());
		bang->contactAngle = 0.0;
		bang->contactPoint = Vector2r(0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperellipseGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}

	Vector2r normal = Vector2r(0,0);
	normal[PA] = sense;
	Vector2r wall_pos = se21.position;
	Vector2r B_pos = se22.position;
  //cout<<"is Sphere? "<<B->isSphere<<endl;
    // is the particle spherical?
    if (B->isSphere){
        //the particle is spherical
        Vector2r p1, point2;
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

	    Matrix2r rot_mat1 = B->rot_mat2local;//(se22.orientation).conjugate().toRotationMatrix();//to particle's system
	    Matrix2r rot_mat2 = B->rot_mat2global;//(se22.orientation).toRotationMatrix();//to global system

      //cout<<"rot1="<<rot_mat1<<endl;
      //cout<<"rot2="<<rot_mat2<<endl;
	    //cout<<"wall_pos"<<wall_pos<<endl;
	    //cout<<"B_pos"<<B_pos<<endl;

	    Real phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
	    Vector2r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
	    //check whether point2 is at the negative side of the wall.
	    Vector2r p1;
        p1 = point2;
        p1[PA] = wall_pos[PA];
	    Vector2r d = point2 - p1;	//the distance vector
	    //d[PA] -= wall_pos[PA];
	    //cout<<"distance v="<<d<<endl;
	    if (normal[PA]*d[PA]>=0.0){//(normal.dot(d) < 0.) {//no touching between the wall and the particle
		    //cout<<"watch dog2"<<endl;
		    bang->PenetrationDepth= 0;
		    bang->isShearNew = true;
		    return true;
	    }

        bang->point1 = p1;
        bang->point2 = point2;
	    bang->contactPoint= 0.5*(point2 + p1);
	    bang->PenetrationDepth= d.norm();
        bang->point11 = p1 + d;
        bang->point22 = point2;

	}//non-spherical end
	//std::cerr << "watchpoint 3: " << "\n";
	bang->precompute(state1,state2,scene,interaction,normal,bang->isShearNew,shift2);
	bang->normal= normal;
	//bang->contactAngle = para;

	return true;
}
//**********************************************************************************
bool Ig2_Fwall_Superellipse_SuperellipseGeom::go(const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const State& state1,
							const State& state2,
							const Vector2r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& interaction)
{

	const Se2r& se21=state1.se2; const Se2r& se22=state2.se2;
	Fwall*  fwall = static_cast<Fwall*>(cm1.get());
	Superellipse* B = static_cast<Superellipse*>(cm2.get());

	bool isNew = !interaction->geom;
	shared_ptr<SuperellipseGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<SuperellipseGeom>(new SuperellipseGeom());
		bang->contactAngle = 0.0;
		bang->contactPoint = Vector2r(0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  	bang=SUDODEM_PTR_CAST<SuperellipseGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}

	Vector2r normal = fwall->normal;

	Vector2r fwall_pos = se21.position;
	Vector2r B_pos = se22.position + shift2;



	//cout<<"is Sphere? "<<B->isSphere<<endl;
		// is the particle spherical?
		if (B->isSphere){
				//the particle is spherical
				//Matrix2r fw_rot = se21.rotation.toRotationMatrix();
				//Matrix2r fw = fw_rot.transpose();
				// local orientation
				Vector2r cl = (B_pos - fwall_pos);  // "contact line" in Fwall-local coords

				//Vector2r normal = fwall->normal;
				Real L = normal.dot(cl);
				if (L<0) {normal=-normal; L=-L; }

				Real diskRadius = B->getr_max();
				Real penetrationDepth = diskRadius - L;
				if (penetrationDepth < 1E-18 && !interaction->isReal() && !force) { // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false
					//	TIMING_DELTAS_CHECKPOINT("Ig2_Fwall_Disk_ScGeom");
					bang->PenetrationDepth= 0;
					bang->isShearNew = true;
					return true;
				}
			  //projection of cl along vu of Fwall
				Vector2r cp = cl - L*normal;
				Real cp_l = cp.norm();
				if (cp_l > fwall->vl*0.5){
					penetrationDepth *= 0.5;
				}
				//normal = fw_rot*normal; // in global orientation
				bang->contactPoint = B_pos - (diskRadius-0.5*penetrationDepth)*normal;
				//bang->point1 = p1;
				//bang->point2 = point2;
				bang->PenetrationDepth= penetrationDepth;
		}else{

			Matrix2r rot_mat1 = B->rot_mat2local;//(se22.orientation).conjugate().toRotationMatrix();//to particle's system
			Matrix2r rot_mat2 = B->rot_mat2global;//(se22.orientation).toRotationMatrix();//to global system
			//to particle B's local coords
			//Matrix2r rot_mat1 = (se22.rotation).inverse().toRotationMatrix();//we do not rotate fwall
			Real phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
			Vector2r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
			//check whether point2 is at the negative side of the wall.
			Real depth = (point2 -fwall_pos).dot(normal);	//the distance vector

			if (depth > 0.0){//no touching between the wall and the particle
				//cout<<"watch dog2"<<endl;
				bang->PenetrationDepth= 0;
				bang->isShearNew = true;
				return true;
			}
			//point at the fwall
			Vector2r p1 = point2 - depth*normal;

			bang->point1 = p1;
			bang->point2 = point2;
			bang->contactPoint= 0.5*(point2 + p1);
			bang->PenetrationDepth= -depth;
				//bang->point11 = p1 + d;
				//bang->point22 = point2;

	}//non-spherical end
	bang->precompute(state1,state2,scene,interaction,normal,bang->isShearNew,shift2);
	bang->normal= normal;
	return true;
}

bool Ig2_Fwall_Superellipse_SuperellipseGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector2r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	//LOG_WARN("Swapped interaction order for "<<c->getId2()<<"&"<<c->getId1());
	return go(cm2,cm1,state2,state1,-shift2,force,c);
}
//**********************************************************************************
/*! Create Superellipse (collision geometry) from colliding Polyhedron and Facet. */
/*
bool Ig2_Facet_Superellipse_SuperellipseGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector2r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){


	const Se2r& se21=state1.se2;
	const Se2r& se22=state2.se2;
	Facet*   A = static_cast<Facet*>(cm1.get());
	Superellipse* B = static_cast<Superellipse*>(cm2.get());

	bool isNew = !interaction->geom;

	shared_ptr<SuperellipseGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<SuperellipseGeom>(new SuperellipseGeom());
		//bang->sep_plane.assign(3,0);
		bang->contactPoint = Vector2r(0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<SuperellipseGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}
	return true;
}
*/
