#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/common/Aabb.hpp>
#include"Superellipse.hpp"

#include<time.h>

#define _USE_MATH_DEFINES


SUDODEM_PLUGIN(/* self-contained in hpp: */ (Superellipse) (SuperellipseGeom) (Bo1_Superellipse_Aabb) (SuperellipsePhys)  (SuperellipseMat)
 									(Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys) (SuperellipseLaw)
	/* some code in cpp (this file): */
	#ifdef SUDODEM_OPENGL
		(Gl1_Superellipse) (Gl1_SuperellipseGeom) /*(Gl1_SuperellipsePhys)*/
	#endif
	);

//move to Math.hpp in the future
////////////////////////////////some auxiliary functions/////////////////////////////////////////////////////////////////
double cot(double x){return tan(M_PI_2l-x);}
int Sign(double f){ if(f<0) return -1; if(f>0) return 1; return 0; }
/////////////////////////gamma func----referring to NR3.0
double gammln(const double xx) {
	int j;
	double x,tmp,y,ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}
///////////////////////beta func
double beta(const double z, const double w) {
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}
//***************************************************************************************
/*Constractor*/

Superellipse::Superellipse(double x, double y, double ep1)
{       createIndex();
	rx=x;
	ry=y;
	eps=ep1;
  rxy = Vector2r(x,y);
  ref_rxy = Vector2r(x,y);
  init = false;//this is very important
	//initialize to calculate some basic geometric parameters.
	Initial();
	//get surface
	//Surface = getSurface();
}

/*sssssssssssssss*/
void Superellipse::Initial()
{       //cout<<"watch init"<<init<<endl;
        if (init) return;
        //
        rx = rxy(0);
        ry = rxy(1);

        //
        //cout<<"rxyz"<<rxyz(0)<<" "<<rxyz(1)<<" "<<rxyz(2)<<endl;
        //cout<<"eps"<<eps(0)<<" "<<eps(1)<<endl;
        r_max = (rx > ry) ? rx : ry;
        r_max *= 1.05;//for virtual contact
        if (eps < 1.0){r_max *= pow(2.0,0.5);}
	//The area and moments of inetia can be expressed in terms of Beta fuction.
	//The detailed derivations are shown in the reference, i.e., A. Jaklič, A. Leonardis,
	//F. Solina, Segmentation and Recovery of Superellipse, Springer Netherlands, Dordrecht, 2000.
	double a,beta1,beta2;//some coefficients
	a = rx*ry*eps;
	beta1 = beta(0.5*eps, 0.5*eps + 1.0);
	beta2 = beta(1.5*eps, 0.5*eps);

	Area = 2.*a*beta1;     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
	Inertia = 2*a*(pow(rx, 2.)+pow(ry, 2.))*beta2;	//I_zz


	Orientation = Rotationr::Identity();
	//rotation matrices
  rot_mat2local = (Orientation).inverse().toRotationMatrix();//to particle's system
	rot_mat2global = (Orientation).toRotationMatrix();//to global system

	//
	Position = Vector2r::Zero();
	//test
	//Quaternionr Rot(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	//Rot.normalize();
	//Orientation =Rot;
	//test_end
	ContactNorm = Vector2r(0.,0.);
	CurrentP = Vector2r(0.,0.);

        //initialization done
	init = true;
}


//****************************************************************************************
/* Destructor */

Superellipse::~Superellipse(){}
SuperellipseGeom::~SuperellipseGeom(){}

//****************************************************************************************
bool Superellipse::isInside(Vector2r p)//check point p is inside the surface using the inside-outside function
{
    //double alpha1,alpha2;
    p = rot_mat2local*p;
		//cout<<"!(1>a) ="<<!(1>a)<<endl;
		//cout<<"~(1>a) ="<<~(1>a)<<endl;
    return 1 > pow(fabs(p(0)/rx),2.0/eps) + pow(fabs(p(1)/ry),2.0/eps);
}
//****************************************************************************************
Vector2r Superellipse::getSurface(Real phi) const//in local coordinates system
{
	Vector2r Surf;
	Surf(0) = Sign(cos(phi))*rx*pow(fabs(cos(phi)), eps);
	Surf(1) = Sign(sin(phi))*ry*pow(fabs(sin(phi)), eps);

	return Surf;
}

Vector2r Superellipse::getNormal(Real phi)
{
	Vector2r n;
	n(0) = Sign(cos(phi)) /rx *pow(fabs(cos(phi)), 2.-eps);
	n(1) = Sign(sin(phi)) /ry *pow(fabs(sin(phi)), 2.-eps);
	return n;
}


Real Superellipse::Normal2Phi(Vector2r n)
{       //n.normalize();
    //cout<<"contact normal: "<<n<<endl;
	//Vector2r phi;
	/*
	phi(0) = atan(Sign(n(1))/Sign(n(0))*pow((fabs(ry*n(1)/rx/n(0))),1/(2.-eps1)));
	if (fabs(rx*n(0))>fabs(ry*n(1)))
	{
		phi(1) = atan(Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(cos(phi(0))),2.-eps1) / fabs(rx*n(0))), 1/(2.-eps2) ) );
	}
	else{
		phi(1) = atan(Sign(n(2)) * pow(fabs(rz*n(2)* pow(fabs(sin(phi(0))),2.-eps1) / fabs(ry*n(1))), 1/(2.-eps2) ) );
	}
	*/
	//use atan2(y,x), which returns arctan(y/x).
	return atan2(Sign(n(1))*pow((fabs(ry*n(1))),1/(2.-eps)), Sign(n(0))*pow((fabs(rx*n(0))),1/(2.-eps))); //atan2 returns in [-pi,pi]//atan(x) returns in [-pi/2, pi/2]
}
/////////////////////////////////////////////////////////////////////////////////////////////////
////the following fuctions aim to calculating the derivates of Points to Phi.
////////////////////////////////////////
//
/*Mat32r Superellipse::Derivate_P2Phi(Vector2r phi)  //3*2 matrix
{	Vector2r P;
	P(0) = Sign(cos(phi(0)))*rx*pow(fabs(cos(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	P(1) = Sign(sin(phi(0)))*ry*pow(fabs(sin(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
	P(2) = Sign(sin(phi(1)))*rz*pow(fabs(sin(phi(1))), eps2);
	Mat32r p2phi;
	p2phi <<
			P(0)*eps1*fabs(tan(phi(0))), P(0)*eps2*fabs(tan(phi(1))),
			P(1)*eps1*fabs(cot(phi(0))), P(1)*eps2*fabs(tan(phi(1))),
			P(2)*0.					  , P(2)*eps2*fabs(cot(phi(1)));

	return p2phi;
}*/


//*************************************************************************************
//****************************************************************************************
/* AaBb overlap checker  */

void Bo1_Superellipse_Aabb::go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se2r& se2, const Body*){

	Superellipse* t=static_cast<Superellipse*>(ig.get());
	if (!t->IsInitialized()){
	        t->Initial();
					//FIXME:Directly changing rotation in python (i.e., something like O.bodies[i].state.ori = xxx) is not allowed temporarily, because this action would not update rotation matrices:rot_mat2local and rot_mat2global.
	        t->rot_mat2local = (se2.rotation).inverse().toRotationMatrix();//to particle's system
	        t->rot_mat2global = (se2.rotation).toRotationMatrix(); //to global system
	}
	double r_max = t->getr_max();
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());

  /*
  //use the tightly closed Box
  Matrix2r rot_mat1 = t->rot_mat2local;//(se22.orientation).conjugate().toRotationMatrix();//to particle's system
  Matrix2r rot_mat2 = t->rot_mat2global;//(se22.orientation).toRotationMatrix();//to global system

  Real phi = t->Normal2Phi(rot_mat1*Vector2r(1,0));//phi in particle's local system
  Vector2r point1 = rot_mat2*( t->getSurface(phi));
  phi = t->Normal2Phi(rot_mat1*Vector2r(0,1));//phi in particle's local system
  Vector2r point2 = rot_mat2*( t->getSurface(phi));
  aabb->min=se2.position - Vector2r(point1[0],point2[1]);
  aabb->max=se2.position + Vector2r(point1[0],point2[1]);
  */
  //use a fixed aabb
	Vector2r halfsize(r_max,r_max);
	//the aabb should be optimized in future.
	if(!scene->isPeriodic){
		aabb->min=se2.position-halfsize; aabb->max=se2.position+halfsize;
		return;
	}
	if(scene->cell->hasShear()) {//debug for 2d
	Vector2r refHalfSize(halfsize);
	const Real& _cos=scene->cell->getCos();
		//cerr<<"cos["<<i<<"]"<<cos[i]<<" ";
		//halfsize+=.5*refHalfSize*(1/_cos-1);
		halfsize = refHalfSize/_cos;
	}
	//cerr<<" || "<<halfSize<<endl;
	aabb->min = scene->cell->unshearPt(se2.position)-halfsize;
	aabb->max = scene->cell->unshearPt(se2.position)+halfsize;


}

//**********************************************************************************
/* Plotting */

#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	bool Gl1_Superellipse::wire;

	int  Gl1_Superellipse::Slices;
	int  Gl1_Superellipse::preSlices=5;

	//int Gl1_Superellipse::glStripedDiskList=-1;
	//int Gl1_Superellipse::glGlutDiskList=-1;

	void Gl1_Superellipse::initSolidGlList(Superellipse* t) {

		//Generate the list. Only once for each qtView, or more if quality is modified.
		glDeleteLists(t->m_glSolidListNum,1);
	  t->m_glSolidListNum = glGenLists(1);
	  glNewList(t->m_glSolidListNum,GL_COMPILE);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		// render the disk now
		glBegin(GL_TRIANGLE_FAN);
		glVertex2f(t->getrxy()[0],0);
	  double del_ang = Mathr::TWO_PI/Slices;
		float angle=0.0;
		for (int i=1 ;i<Slices;i++)
		{		angle = del_ang*i;
		    glVertex2v(t->getSurface(angle));
		}
		glVertex2f(t->getrxy()[0],0);

		glEnd();
		glEndList();

	}

	void Gl1_Superellipse::initWireGlList(Superellipse* t){
		//Generate the "no-stripes" display list, each time quality is modified
		  glDeleteLists(t->m_glWireListNum,1);
		  t->m_glWireListNum = glGenLists(1);
			glNewList(t->m_glWireListNum,GL_COMPILE);
			glEnable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);
			glBegin(GL_LINE_LOOP);
			glVertex2f(t->getrxy()[0],0);
		  double del_ang = Mathr::TWO_PI/Slices;
			float angle=0.0;
			for (int i=1 ;i<Slices;i++)
			{		angle = del_ang*i;
			    glVertex2v(t->getSurface(angle));
			}
			//glVertex2f(t->getrxy()[0],0);

			glEnd();
		glEndList();
	}
	void Gl1_Superellipse::go(const shared_ptr<Shape>& cm, const shared_ptr<State>&state,bool wire2,const GLViewInfo&)
	{
		glMaterialv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,Vector3r(cm->color[0],cm->color[1],cm->color[2]));
		glColor3v(cm->color);
		Superellipse* t=static_cast<Superellipse*>(cm.get());
		if (!t->IsInitialized()){
	                t->Initial();
	                t->rot_mat2local = state->ori.inverse().toRotationMatrix();//to particle's system
	                t->rot_mat2global = state->ori.toRotationMatrix(); //to global system
	        }

        bool somethingChanged = (t->m_GL_slices!=Slices|| glIsList(t->m_glSolidListNum)!=GL_TRUE);
		if (somethingChanged) {
            //std::cout<<"haha preslices"<<t->m_GL_slices<<"slices"<<slices<<std::endl;
            initSolidGlList(t); initWireGlList(t);t->m_GL_slices=Slices;}
        if(wire||wire2)
        {
            glCallList(t->m_glWireListNum);//glCallList(t->m_glListNum);
        }else{
            glCallList(t->m_glSolidListNum);
        }
	}
	//
	void Gl1_SuperellipseGeom::go(const shared_ptr<IGeom>& ig, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool){
	        SuperellipseGeom* pg=static_cast<SuperellipseGeom*>(ig.get());
	        Vector2r point1,point2;
          //cout<<"p1,p11,p2,p11"<<pg->point1<<pg->point2<<pg->point11<<pg->point22<<endl;
	        point1 = (!scene->isPeriodic ? pg->point1 : scene->cell->wrapShearedPt(pg->point1));
	        point2 = (!scene->isPeriodic ? pg->point2 : scene->cell->wrapShearedPt(pg->point2));
	        //testing

	        Vector2r point11,point22;
	        point11 = (!scene->isPeriodic ? pg->point11 : scene->cell->wrapShearedPt(pg->point11));
	        point22 = (!scene->isPeriodic ? pg->point22 : scene->cell->wrapShearedPt(pg->point22));

          /*
          Vector2r p1=Body::byId(i->getId1(),scene)->state->pos;
          const Vector2r& size=scene->cell->getSize();
          Vector2r shift2(i->cellDist[0]*size[0],i->cellDist[1]*size[1]);
          // in sheared cell, apply shear on the mutual position as well
          shift2=scene->cell->shearPt(shift2);
          Vector2r rel=Body::byId(i->getId2(),scene)->state->pos+shift2-p1;
          if(scene->isPeriodic) p1=scene->cell->wrapShearedPt(p1);
          glBegin(GL_LINES); glVertex2v(p1);glVertex2v(Vector2r(p1+rel));glEnd();
          */


	        //draw
	        glColor3f(0.5, 1., 1.0);
          glPointSize(10.0f);
          //the two closest points
          glBegin(GL_POINTS);
          glVertex2v(point1);
          glVertex2v(point2);
          glEnd();

          //Line connecting the two closest points
          glLineStipple (1, 0x0F0F);
          glBegin(GL_LINES);glLineWidth (3.0);

          glVertex2v(point1); glVertex2v(point2);
          //
          glVertex2v(point1); glVertex2v(point11);
          glVertex2v(point2); glVertex2v(point22);

          glEnd();
	}

#endif
//**********************************************************************************
//!Precompute data needed for rotating tangent vectors attached to the interaction

void SuperellipseGeom::precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector2r&
currentNormal, bool isNew, const Vector2r& shift){
	if(!isNew) {
		//rotAngle = scene->dt*0.5*(rbp1.angVel + rbp2.angVel);
    preNormal = normal;
	}
	else preNormal = currentNormal;
	//Update contact normal
	normal=currentNormal;

	//Precompute shear increment
	Vector2r c1x = (contactPoint - rbp1.pos);
	Vector2r c2x = (contactPoint - rbp2.pos - shift2);
	//std::cout<<"precompute---"<<c1x<<" "<<c2x<<std::endl;
	Vector2r del = Vector2r(-normal[1],normal[0]);
	Vector2r relativeVelocity = rbp2.vel - rbp1.vel - rbp2.angVel*c2x.norm()*del  - rbp1.angVel*c1x.norm()*del;
	//Vector2r relativeVelocity = (rbp2.vel+rbp2.angVel.cross(c2x)) - (rbp1.vel+rbp1.angVel.cross(c1x));
  //cout<<"relv1="<<relativeVelocity<<" del="<<del<<endl;
	Vector2r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(c->cellDist) : Vector2r::Zero();
	relativeVelocity += shiftVel;
	//keep the normal and shear parts
	relativeVn = normal.dot(relativeVelocity)*normal;
	relativeVs = relativeVelocity-relativeVn;
  //cout<<"shiftVel="<<shiftVel<<" relv2="<<relativeVelocity<<" revn="<<relativeVn<<" revs="<<relativeVs<<endl;
	shearInc = relativeVs*scene->dt;
  shift2 = shift;
}

//**********************************************************************************
Vector2r& SuperellipseGeom::rotate(Vector2r& shearForce) const {
	// approximated rotations
	/*double sin_theta = sin(rotAngle);
	double cos_theta = cos(rotAngle);
	Matrix2r rot;
	rot << cos_theta, -sin_theta,
				 sin_theta, cos_theta;
	shearForce = rot*shearForce;
 */
 Vector2r tangent = Vector2r(-normal[1],normal[0]);
 Vector2r preTangent = Vector2r(-preNormal[1],preNormal[0]);//tangent is left-hand side of normal
 shearForce = shearForce.dot(preTangent)*tangent;
 return shearForce;
}


//**********************************************************************************
/* Material law, physics */

void Ip2_SuperellipseMat_SuperellipseMat_SuperellipsePhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<SuperellipseMat>& mat1 = SUDODEM_PTR_CAST<SuperellipseMat>(b1);
	const shared_ptr<SuperellipseMat>& mat2 = SUDODEM_PTR_CAST<SuperellipseMat>(b2);
	interaction->phys = shared_ptr<SuperellipsePhys>(new SuperellipsePhys());
	const shared_ptr<SuperellipsePhys>& contactPhysics = SUDODEM_PTR_CAST<SuperellipsePhys>(interaction->phys);
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
Real SuperellipseLaw::getPlasticDissipation() {return (Real) plasticDissipation;}
void SuperellipseLaw::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}
Real SuperellipseLaw::elasticEnergy()
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
// Apply forces on superellipse in collision based on geometric configuration
bool SuperellipseLaw::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){

		if (!I->geom) {return true;}
		const shared_ptr<SuperellipseGeom>& contactGeom(SUDODEM_PTR_DYN_CAST<SuperellipseGeom>(I->geom));
		if(!contactGeom) {return true;}
		const Body::id_t idA=I->getId1(), idB=I->getId2();
		const shared_ptr<Body>& A=Body::byId(idA), B=Body::byId(idB);

                State* de1 = Body::byId(idA,scene)->state.get();
	        State* de2 = Body::byId(idB,scene)->state.get();

		SuperellipsePhys* phys = dynamic_cast<SuperellipsePhys*>(I->phys.get());

    Matrix2r cellHsize; Vector2r ab_min,ab_max;
    ab_min = B->bound->min;ab_max = B->bound->max;
    if(scene->isPeriodic){
      cellHsize=scene->cell->hSize;
      Vector2r shift2=cellHsize*I->cellDist.cast<Real>();
      ab_min += shift2;
      ab_max += shift2;
    }
		//erase the interaction when aAbB shows separation, otherwise keep it to be able to store previous separating plane for fast detection of separation
		if (A->bound->min[0] >= ab_max[0] || ab_min[0] >= A->bound->max[0] || A->bound->min[1] >= ab_max[1] || ab_min[1] >= A->bound->max[1])  {
		        phys->normalForce = Vector2r(0.,0.); phys->shearForce = Vector2r(0.,0.);
			scene->interactions->requestErase(I);
			return false;
		}

		//zero penetration depth means no interaction force
		if(!(contactGeom->PenetrationDepth > 1E-18) ) {
            phys->normalForce = Vector2r(0.,0.); phys->shearForce = Vector2r(0.,0.);
            return true;
        }
		Vector2r normalForce=contactGeom->normal*contactGeom->PenetrationDepth*phys->kn;
		Vector2r shearForce1=Vector2r::Zero();//zhswee deprecated SuperellipseLaw.shearForce
		//shear force: in case the polyhdras are separated and come to contact again, one should not use the previous shear force
		//std::cerr << "p1:"<<gettid4()<<shearForce[0]<< "\n";
		///////////////////////
		bool useDamping=(phys->betan > 1e-5 || phys->betas > 1e-5);//using viscous damping?
		// tangential and normal stiffness coefficients, recomputed from betan,betas at every step
        Real cn=0, cs=0;
        Vector2r normalViscous = Vector2r::Zero();
        Vector2r shearViscous = Vector2r::Zero();
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
	        Vector2r normTemp = normalForce - normalViscous; // temporary normal force
	        // viscous force should not exceed the value of current normal force, i.e. no attraction force should be permitted if particles are non-adhesive
	        //std::cerr << "nF1:"<<normalForce.norm()<< "\n";
	        if (normTemp.dot(contactGeom->normal)<0.0){

				normalViscous = normalForce; // normal viscous force is such that the total applied force is null - it is necessary to compute energy correctly!
				normalForce = Vector2r::Zero();
			}
			else{normalForce -= normalViscous;}
	                //std::cerr << "nF2:"<<normalForce.norm()<< "\n";
        }
        //if(isnan(phys->shearForce[0])){std::cerr << "NAN, phys->shearForce:"<<shearForce1<< "\n";}
		if (contactGeom->isShearNew){
			//shearForce = Vector2r::Zero();
			shearForce1 = Vector2r(0.,0.);
		}else{
			//shearForce = contactGeom->rotate(shearForce);

			shearForce1 = contactGeom->rotate(phys->shearForce);
			//std::cerr << "p2:"<<gettid4() <<shearForce[0]<< "\n";
        }
		const Vector2r& shearDisp = contactGeom->shearInc;
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
			//almost the same with additional Vector2r instatinated for energy tracing,
			//duplicated block to make sure there is no cost for the instanciation of the vector when traceEnergy==false
			if(shearForce.squaredNorm() > maxFs){
				Real ratio = sqrt(maxFs) / shearForce.norm();
				Vector2r trialForce=shearForce;//store prev force for definition of plastic slip
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
		Vector2r T = (B->state->pos-contactGeom->contactPoint).cross(F);
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
          const char* filename = "SuperellipseLaw_err.log";
		  //sprintf(filename,"%d",gettid4());
		  FILE * fin = fopen(filename,"a");
          fprintf(fin,"\n*******%s*******\n",tmp);
          Vector2r normal = contactGeom->normal;
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
          normalForce = Vector2r::Zero();
          shearForce1 = Vector2r::Zero();
		}
		Vector2r F = -normalForce-shearForce1+shearViscous;	//the normal force acting on particle A (or 1) is equal to the normal force.
		if (contactGeom->PenetrationDepth != contactGeom->PenetrationDepth) exit(1);
		scene->forces.addForce (idA,F);
		scene->forces.addForce (idB, -F);
		//
		Vector2r l1 = contactGeom->contactPoint - A->state->pos;
		double ang_F = atan2(F[1],F[0]);
		double theta1 = ang_F - atan2(l1[1],l1[0]);
		Vector2r l2 = contactGeom->contactPoint - B->state->pos - contactGeom->shift2;
		double theta2 = ang_F + Mathr::PI - atan2(l2[1],l2[0]);

		scene->forces.addTorque(idA, l1.norm()*F.norm()*sin(theta1));
		scene->forces.addTorque(idB, l2.norm()*F.norm()*sin(theta2));
		//scene->forces.addTorque(idA, -(A->state->pos-contactGeom->contactPoint).cross(F));
		//scene->forces.addTorque(idB, (B->state->pos-contactGeom->contactPoint).cross(F));
		//needed to be able to acces interaction forces in other parts of sudodem
		shearForce=shearForce1;
		phys->normalForce = normalForce;
		phys->shearForce = shearForce1;
		//may need to add phys->shearViscous and phys->normalViscous
		return true;
}
