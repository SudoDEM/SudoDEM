//@2016 Sway Zhao,  zhswee@gmail.com
//


#include<sudodem/pkg/dem/Superquadrics.hpp>
#include<sudodem/pkg/dem/PolySuperellipsoid.hpp>

#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/lib/pyutil/doc_opts.hpp>
#include<cmath>

#include<numpy/ndarrayobject.h>
#include <boost/concept_check.hpp>
#include <chrono> //test function running time

//quasi-random numbers for Monte Carlo //boost 1.70, random
#ifdef BOOST170
#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#endif

using namespace std::chrono;
namespace py = boost::python;
//***********************************************************************************

int Sign(double f){ if(f<0) return -1; if(f>0) return 1; return 0; }

//***********************************************************************************
double triangleArea(Vector3r p1,Vector3r p2,Vector3r p3){
    Vector3r v12 = p2 -p1;
    Vector3r v13 = p3 - p1;
    Vector3r vcross = v12.cross(v13);
    return 0.5*vcross.norm();

}
double polygonArea(Vector3r p1,Vector3r p2,Vector3r p3,Vector3r p4){
   //calculate surface area
    double area = 0.0;
    area = triangleArea(p1,p2,p3) + triangleArea(p1,p3,p4);

    return area;

}
double getSurfaceArea(double rx,double ry,double rz, double eps1, double eps2,int w, int h){
////
    std::vector<Vector3r> vertices;//discretized vertices of a particle
    double area = 0.0;
    int i,j;
    double a=0.0,b=0.0,phi0,phi1;
    double hStep=M_PI/(h-1);
    double wStep=2*M_PI/w;

    Vector3r Surf;
    //surface discretization
    for(a=0.0,i=0;i<h;i++,a+=hStep){
      for(b=0.0,j=0;j<w;j++,b+=wStep)
      {
        phi0 = b;
        phi1 = a-M_PI_2l;

        //get surface point
        Surf(0) = Sign(cos(phi0))*rx*pow(fabs(cos(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
        Surf(1) = Sign(sin(phi0))*ry*pow(fabs(sin(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
        Surf(2) = Sign(sin(phi1))*rz*pow(fabs(sin(phi1)), eps2);
	    vertices.push_back(Surf);
      }
    }

    //polygons of surface slices

	for(i=0;i<h-1;i++)
	{
	    for(j=0;j<w-1;j++)
	    {

            area += polygonArea(vertices[i*w+j],vertices[i*w+j+1],vertices[(i+1)*w+j+1],vertices[(i+1)*w+j]);
	    }

        area += polygonArea(vertices[i*w+j],vertices[i*w],vertices[(i+1)*w],vertices[(i+1)*w+j]);

	}
    return area;
}

//get the surface area of super-ellipsoid by ID
double getSurfArea(int id,int w, int h){
    const shared_ptr<Body>& b=Body::byId(id);
    Vector3r rxyz;
    Vector2r eps;
    if (w<10){w=10;}
    if(h<10){h=10;}
    double area(0.0);
	if (b->shape->getClassName()=="Superquadrics"){
		const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);

    	rxyz = A->getrxyz();
        eps = A->geteps();
		area = getSurfaceArea(rxyz(0),rxyz(1),rxyz(2),eps(0),eps(1),w, h);

    }
    return area;
}

//particle pairs with given penetration depth and contact normal
std::vector<Vector3r> PositionPB(Vector6r rxyzA, Vector2r epsA,Vector6r rxyzB, Vector2r epsB,Vector2r phi, double pd,double q0, double q1, double q2, double q3){
  //Particle A:without rotation locating at the origin
  PolySuperellipsoid A = PolySuperellipsoid(epsA,rxyzA);
  Vector3r normal = A.getNormal(phi);//also contact normal
  normal.normalize();
  Vector3r PA = A.getSurfaceMC(phi);//point A
  Vector3r PB = PA + (-normal*pd);
  //Particle B: with random rotation
  PolySuperellipsoid B = PolySuperellipsoid(epsB,rxyzB);
  Quaternionr Ori;
  Ori.w() =q0;
  Ori.x() =q1;
  Ori.y() =q2;
  Ori.z() =q3;
  B.setOrientation(Ori);
  Vector3r PB0 = B.Normal2SurfaceMC(-B.rot_mat2local*normal);//point B at B's local coordinate
  Vector3r posB = PB - B.rot_mat2global*PB0;//position of particle B (located by mass center)
  std::vector<Vector3r> data;
  data.push_back(0.5*(PA+PB));//contact point
  data.push_back(posB);
  return data;
}
//time consumptopm of the support function
double EST_Support(Vector6r rxyzA, Vector2r epsA,Vector2r para, int num_exe,double q0, double q1, double q2, double q3){
  PolySuperellipsoid A = PolySuperellipsoid(epsA,rxyzA);
  Vector3r contact;
  contact(0) = cos(para(0))*cos(para(1));  //the base vectors are e_i, i=1,2,3
  contact(1) = sin(para(0))*cos(para(1));  //e_1 = (s^2 - s^1)/||s^2 - s^1||
  contact(2) = sin(para(1));
  Vector3r normal = A.getNormal(para);//also contact normal
  Vector3r PA = A.getSurfaceMC(para);//point A
  Quaternionr Ori;
  Ori.w() =q0;
  Ori.x() =q1;
  Ori.y() =q2;
  Ori.z() =q3;
  A.setOrientation(Ori);
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  for(int i=0;i<num_exe;i++){
    //A.support(contact);
    (normal - PA).norm();
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration2 = duration_cast<nanoseconds>( t2 - t1 ).count();
  return duration2;
}
//new Superquadrics

void Superquadrics_test(Vector6r rxyz, Vector2r eps, Vector2r phi){
  PolySuperellipsoid A = PolySuperellipsoid(eps,rxyz);
  Vector3r normal = A.getNormal(phi);


  double x,y,z,a,alpha;
  normal.normalize();
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  Vector2r phi1 = A.Normal2Phi(normal);
  Vector3r pos = A.getSurface(phi1);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration1 = duration_cast<nanoseconds>( t2 - t1 ).count();
  t1 = high_resolution_clock::now();
  Vector3r pos3 = A.Normal2SurfaceMC(normal);
  t2 = high_resolution_clock::now();
  auto duration2 = duration_cast<nanoseconds>( t2 - t1 ).count();
  cout <<"duration1="<< duration1<<", "<<duration2<<endl;
  cout<<"phi="<<phi<<"normal="<<normal<<endl;
  cout<<"pos3="<<pos3<<endl;
}
#ifdef BOOST170
std::vector<double> MonteCarlo_test(Vector6r rxyz, Vector2r eps,long int num_max = 1000000){
  PolySuperellipsoid A = PolySuperellipsoid(eps,rxyz);
  double v1 = A.getVolume();
  //cout<<"computing volume="<<v1<<endl;
  double x_l, y_l,z_l;
  x_l = rxyz(0)+rxyz(1);
  y_l = rxyz(2)+rxyz(3);
  z_l = rxyz(4)+rxyz(5);
  long int num = 0,count = 0;//1 million points,num_max=1000000
  double x,y,z;
  Vector3r center = Vector3r::Zero();
  Vector3r inertia = Vector3r::Zero();
  //pseudo-random numbers
  /*
  for(;num<=num_max;num++){
    x = x_l*double(rand())/RAND_MAX-rxyz(1);
    y = y_l*double(rand())/RAND_MAX-rxyz(3);
    z = z_l*double(rand())/RAND_MAX-rxyz(5);
    if(A.isInside(Vector3r(x,y,z)-A.getMassCenter())){
      count++;
      center += Vector3r(x,y,z);
      inertia += Vector3r(y*y+z*z,x*x+z*z,x*x+y*y);
    }
  }*/
  //quasi-random numbers
  static const std::size_t dimension = 3;
  // Create a generator
  typedef boost::variate_generator<boost::random::sobol&, boost::uniform_01<double> > quasi_random_gen_t;
  // Initialize the engine to draw randomness out of thin air
  boost::random::sobol engine(dimension);
  // Glue the engine and the distribution together
  quasi_random_gen_t gen(engine, boost::uniform_01<double>());

  std::vector<double> sample(dimension);

  // At this point you can use std::generate, generate member f-n, etc.
  //std::generate(sample.begin(), sample.end(), gen);
  //engine.generate(sample.begin(), sample.end());

  for(;num<=num_max;num++){
    std::generate(sample.begin(), sample.end(), gen);
    //cout<<"quasi-random number="<<sample[0]<<" "<<sample[1]<<" "<<sample[2]<<endl;
    x = x_l*sample[0]-rxyz(1);
    y = y_l*sample[1]-rxyz(3);
    z = z_l*sample[2]-rxyz(5);
    if(A.isInside(Vector3r(x,y,z)-A.getMassCenter())){
      count++;
      center += Vector3r(x,y,z);
      inertia += Vector3r(y*y+z*z,x*x+z*z,x*x+y*y);
    }
  }


  double mv = double(count)/double(num_max)*x_l*y_l*z_l;
  double alpha = x_l*y_l*z_l/double(num_max);
  center /= double(count);
  inertia *= alpha;
  inertia -= Vector3r(center(1)*center(1)+center(2)*center(2),center(0)*center(0)+center(2)*center(2),center(1)*center(1)+center(0)*center(0))*mv;
  //cout<<"Monte Carlo Volume="<<mv<<endl;
  //cout<<"center1="<<A.getMassCenter()<<" center2="<<center<<endl;
  //cout<<"inertia1="<<A.getInertia()<<" inertia2="<<inertia<<endl;
  Vector3r inertia1 = A.getInertia();
  //return
  std::vector<double> data;
  data.push_back(v1);//analytical
  data.push_back(mv);//MC
  for(int i=0;i<3;i++){
    data.push_back(inertia1[i]);//analytical
    data.push_back(inertia[i]);//MC
  }
  return data;
}
#endif
shared_ptr<Body> NewSuperquadrics(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,bool rotate){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superquadrics>(new Superquadrics(x,y,z,ep1,ep2));
	Superquadrics* A = static_cast<Superquadrics*>(body->shape.get());
    A->isSphere = false;//non-spherical by default.
	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getVolume();
	body->state->inertia =A->getInertia()*body->material->density;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//
	if(rotate){
		double heading = double(rand())/RAND_MAX*2.0*Mathr::PI;//rotate around z
		double attitude = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double bank = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double c1 = cos(heading/2);
    double s1 = sin(heading/2);
    double c2 = cos(attitude/2);
    double s2 = sin(attitude/2);
    double c3 = cos(bank/2);
    double s3 = sin(bank/2);
    double c1c2 = c1*c2;
    double s1s2 = s1*s2;
    Quaternionr Ori;


    Ori.w() =c1c2*c3 - s1s2*s3;
    Ori.x() =c1c2*s3 + s1s2*c3;
    Ori.y() =s1*c2*c3 + c1*s2*s3;
    Ori.z() =c1*s2*c3 - s1*c2*s3;
    A->setOrientation(A->getOrientation()+Ori);
	 }
	body->state->ori = A->getOrientation();
	return body;
}

shared_ptr<Body> NewSuperquadrics2(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,bool rotate, bool isSphere){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superquadrics>(new Superquadrics(x,y,z,ep1,ep2));
	Superquadrics* A = static_cast<Superquadrics*>(body->shape.get());

  A->isSphere = isSphere;
	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getVolume();
	body->state->inertia =A->getInertia()*body->material->density;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//
	if(rotate){
		double heading = double(rand())/RAND_MAX*2.0*Mathr::PI;//rotate around z
		double attitude = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double bank = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double c1 = cos(heading/2);
    double s1 = sin(heading/2);
    double c2 = cos(attitude/2);
    double s2 = sin(attitude/2);
    double c3 = cos(bank/2);
    double s3 = sin(bank/2);
    double c1c2 = c1*c2;
    double s1s2 = s1*s2;
    Quaternionr Ori;


    Ori.w() =c1c2*c3 - s1s2*s3;
    Ori.x() =c1c2*s3 + s1s2*c3;
    Ori.y() =s1*c2*c3 + c1*s2*s3;
    Ori.z() =c1*s2*c3 - s1*c2*s3;
    A->setOrientation(A->getOrientation()+Ori);
	 }
	body->state->ori = A->getOrientation();
	return body;
}
shared_ptr<Body> NewSuperquadrics_rot(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,double q0, double q1, double q2, double q3){//with a specified quatanion
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superquadrics>(new Superquadrics(x,y,z,ep1,ep2));
	Superquadrics* A = static_cast<Superquadrics*>(body->shape.get());

	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getVolume();
	body->state->inertia =A->getInertia()*body->material->density;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//

  Quaternionr Ori;


  Ori.w() = q0;
  Ori.x() = q1;
  Ori.y() = q2;
  Ori.z() = q3;
  A->setOrientation(Ori);

	body->state->ori = A->getOrientation();
	return body;
}

shared_ptr<Body> NewSuperquadrics_rot2(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,double q0, double q1, double q2, double q3, bool isSphere){//with a specified quatanion
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superquadrics>(new Superquadrics(x,y,z,ep1,ep2));
	Superquadrics* A = static_cast<Superquadrics*>(body->shape.get());
    A->isSphere = isSphere;
	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getVolume();
	body->state->inertia =A->getInertia()*body->material->density;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//

  Quaternionr Ori;


  Ori.w() = q0;
  Ori.x() = q1;
  Ori.y() = q2;
  Ori.z() = q3;
  A->setOrientation(Ori);

	body->state->ori = A->getOrientation();
	return body;
}
/*
shared_ptr<Body> NewSuperquadrics(Vector3r rxyz, Vector2r eps, shared_ptr<Material> mat,bool rotate){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superquadrics>(new Superquadrics(rxyz,eps));
	Superquadrics* A = static_cast<Superquadrics*>(body->shape.get());
        body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getVolume();
	body->state->inertia =A->getInertia()*body->material->density;

	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//
	if(rotate){
		double heading = double(rand())/RAND_MAX*2.0*Mathr::PI;//rotate around z
		double attitude = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double bank = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double c1 = cos(heading/2);
                double s1 = sin(heading/2);
                double c2 = cos(attitude/2);
                double s2 = sin(attitude/2);
                double c3 = cos(bank/2);
                double s3 = sin(bank/2);
                double c1c2 = c1*c2;
                double s1s2 = s1*s2;
                Quaternionr Ori;


                Ori.w() =c1c2*c3 - s1s2*s3;
                Ori.x() =c1c2*s3 + s1s2*c3;
	        Ori.y() =s1*c2*c3 + c1*s2*s3;
	        Ori.z() =c1*s2*c3 - s1*c2*s3;
	        A->setOrientation(A->getOrientation()+Ori);
	 }
	body->state->ori = A->getOrientation();
	return body;
}*/
shared_ptr<Body> NewPolySuperellipsoid(Vector2r eps,Vector6r rxyz, shared_ptr<Material> mat,bool rotate, bool isSphere, double inertiaScale = 1.0){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<PolySuperellipsoid>(new PolySuperellipsoid(eps,rxyz));
	PolySuperellipsoid* A = static_cast<PolySuperellipsoid*>(body->shape.get());

  A->isSphere = isSphere;
	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getVolume();
	body->state->inertia =A->getInertia()*body->material->density*inertiaScale;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//
	if(rotate){
		double heading = double(rand())/RAND_MAX*2.0*Mathr::PI;//rotate around z
		double attitude = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double bank = double(rand())/RAND_MAX*2.0*Mathr::PI; //y
		double c1 = cos(heading/2);
    double s1 = sin(heading/2);
    double c2 = cos(attitude/2);
    double s2 = sin(attitude/2);
    double c3 = cos(bank/2);
    double s3 = sin(bank/2);
    double c1c2 = c1*c2;
    double s1s2 = s1*s2;
    Quaternionr Ori;


    Ori.w() =c1c2*c3 - s1s2*s3;
    Ori.x() =c1c2*s3 + s1s2*c3;
    Ori.y() =s1*c2*c3 + c1*s2*s3;
    Ori.z() =c1*s2*c3 - s1*c2*s3;
    A->setOrientation(A->getOrientation()+Ori);
	 }
	body->state->ori = A->getOrientation();
	return body;
}
shared_ptr<Body> NewPolySuperellipsoid_rot(Vector2r eps, Vector6r rxyz, shared_ptr<Material> mat,double q0, double q1, double q2, double q3, bool isSphere, double inertiaScale = 1.0){//with a specified quatanion
	shared_ptr<Body> body = NewPolySuperellipsoid(eps,rxyz, mat,false, isSphere, inertiaScale);
  Quaternionr Ori;
  Ori.w() = q0;
  Ori.x() = q1;
  Ori.y() = q2;
  Ori.z() = q3;
  PolySuperellipsoid* A = static_cast<PolySuperellipsoid*>(body->shape.get());
  A->setOrientation(Ori);

	body->state->ori = A->getOrientation();
	return body;
}
//**********************************************************************************
//output some basic data for post-processing, e.g., POV-ray
void outputPOV(const char* filename){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write
	//writing file head of a POV file
	fprintf(fin,"#version 3.6\n");
	fprintf(fin,"//using the command: povray **.pov +A0.01 +W1600 +H1200\n");
        fprintf(fin,"// Right-handed coordinate system where the z-axis points upwards\n");
        //writing camera property
        fprintf(fin,"//camera {\n");
        fprintf(fin,"//	location <4,-2,1>//<8,-10,7>\n");
        fprintf(fin,"//	sky z\n");
        fprintf(fin,"//	right 0.15*x*image_width/image_height\n");
        fprintf(fin,"//	up 0.15*z\n");
        fprintf(fin,"//	look_at <0.2,0.04,0.04>\n");
        fprintf(fin,"//}\n");

        fprintf(fin,"camera{ //orthographic angle 45\n");
        fprintf(fin,"        location <0.6,0.6,0.5>*50\n");
        fprintf(fin,"        sky z\n");
        fprintf(fin,"        look_at  <0,0,0.05>*10\n");
        fprintf(fin,"        right x*image_width/image_height\n");
        fprintf(fin,"        translate <-0.05,0.01,-0.2>*50\n");
        fprintf(fin,"      }\n");
        //writing background
        fprintf(fin,"// White background\n");
        fprintf(fin,"background{rgb 1}\n");
        //writing lights
        fprintf(fin,"// Two lights with slightly different colors\n");
        fprintf(fin,"light_source{<4,8,5>*10 color rgb <1,1,1>}\n");
        fprintf(fin,"light_source{<12,-6>*10 color rgb <1,1,1>}\n");

	//fprintf(fin,"volume\t%e\n",volume);
	//fprintf(fin,"centroid\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2]);

        //output superquadrics
        Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	//particlesVolume = 0;
	Quaternionr quaternion;
	double q0,q1,q2,q3,phi,theta,psi;
	Vector3r pos, color,rxyz;
	Vector2r eps;
	for ( ; bi!=biEnd; ++bi )
	{
		const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		if (b->shape->getClassName()=="Superquadrics"){
			const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());
			State* state=b->state.get();
		        //Euler Angles from Quaternion
		        quaternion = state->ori;
		        q0 = quaternion.w();
		        q1 = quaternion.x();
		        q2 = quaternion.y();
		        q3 = quaternion.z();
		        phi = atan2(2.0*(q0*q1 + q2*q3), 1.0 - 2.0*(q1*q1 + q2*q2))/Mathr::PI*180.0;
		        theta = asin(2.0*(q0*q2 - q3*q1))/Mathr::PI*180.0;
		        psi = atan2(2.0*(q0*q3 + q1*q2), 1.0 - 2.0*(q2*q2 + q3*q3))/Mathr::PI*180.0;
		        pos = state->pos;
		        color = A->color;
		        rxyz = A->getrxyz();
		        eps = A->geteps();
		        //output
		        fprintf(fin,"superellipsoid{ <%f,%f>\n",eps(0),eps(1));
		        fprintf(fin, "  texture{ pigment{ color rgb<%f,%f,%f>}\n",color(0),color(1),color(2));
		        fprintf(fin, "           finish { phong 1}\n");
		        fprintf(fin, "         } // end of texture\n");
		        fprintf(fin, "  scale <%f,%f,%f>\n",rxyz(0),rxyz(1),rxyz(2));
		        fprintf(fin, "  rotate<%f,%f,%f>\n",phi,theta,psi);
		        fprintf(fin, "  translate<%f,%f,%f>\n",pos(0),pos(1),pos(2));
		        fprintf(fin, "}\n");
		}
		if (b->shape->getClassName()=="Sphere"){
			const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b->shape);
			State* state=b->state.get();
		        pos = state->pos;
		        color = A->color;
		        //output
		        fprintf(fin,"sphere{ <%f,%f,%f>,%f\n",pos(0),pos(1),pos(2),A->radius);
		        fprintf(fin, "  texture{ pigment{ color rgb<%f,%f,%f>}\n",color(0),color(1),color(2));
		        fprintf(fin, "           finish { phong 1}}\n");
		        fprintf(fin, "}\n");
		}

	}


	fclose(fin);


}
//**********************************************************************************
//output some basic data for post-processing, e.g., POV-ray
void PolySuperellipsoidPOV(const char* filename,std::vector<int> ids, double scale){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write
	//writing file head of a POV file
	//fprintf(fin,"#version 3.6\n");
  //fprintf(fin, "#include \"colors.inc\"");
  //fprintf(fin, "#declare singleColor = %s;\n", singleColor?"true":"false");
	//fprintf(fin,"//using the command: povray **.pov +A0.01 +W1600 +H1200\n");
  //fprintf(fin,"// Right-handed coordinate system where the z-axis points upwards\n");
  //writing camera property
  /*
  fprintf(fin,"//camera {\n");
  fprintf(fin,"//	location <4,-2,1>//<8,-10,7>\n");
  fprintf(fin,"//	sky z\n");
  fprintf(fin,"//	right 0.15*x*image_width/image_height\n");
  fprintf(fin,"//	up 0.15*z\n");
  fprintf(fin,"//	look_at <0.2,0.04,0.04>\n");
  fprintf(fin,"//}\n");
  */
  //writing background
  // White background
  //fprintf(fin,"background{rgb 1}\n");
  //writing lights
  //fprintf(fin,"// Two lights with slightly different colors\n");
  //fprintf(fin,"light_source{<4,8,5>*10 color rgb <1,1,1>}\n");
  //fprintf(fin,"light_source{<12,-6>*10 color rgb <1,1,1>}\n");
  //fprintf(fin, "#declare transparency=%.1f\n", 0.3);
  //define myremain to clip an object for an octant
  fprintf(fin,"#macro myremain(sign1,sign2,sign3)\n");//remain the part identified by sign1,sign2 and sign3
  fprintf(fin,"	clipped_by{plane{-sign1*x,0}}\n");
  fprintf(fin,"	clipped_by{plane{-sign2*y,0}}\n");
  fprintf(fin,"	clipped_by{plane{-sign3*z,0}}\n");
  fprintf(fin,"#end\n");
  //fprintf(fin, "#declare Random_3 = seed (%d);\n#declare Random_4 = seed (%d);\n#declare Random_5 = seed (%d);\n", 1432,7242,9912);
  //random color
  /*
  fprintf(fin, "#macro qtom(w,x,y,z)\n");
  fprintf(fin, "sqw = w*q\n");
  fprintf(fin, "sqx = x*x\n");
  fprintf(fin, "sqy = y*y\n");
  fprintf(fin, "sqz = z*z\n");

  fprintf(fin, "m00 = ( sqx - sqy - sqz + sqw)\n");
  fprintf(fin, "m11 = (-sqx + sqy - sqz + sqw)\n");
  fprintf(fin, "m22 = (-sqx - sqy + sqz + sqw)\n");

  fprintf(fin, "tmp1 = x*y\n");
  fprintf(fin, "tmp2 = z*w\n");
  fprintf(fin, "m10 = 2.0 * (tmp1 + tmp2)\n");
  fprintf(fin, "m01 = 2.0 * (tmp1 - tmp2)\n");

  fprintf(fin, "tmp1 = x*z\n");
  fprintf(fin, "tmp2 = y*w\n");
  fprintf(fin, "m20 = 2.0 * (tmp1 - tmp2)\n");
  fprintf(fin, "m02 = 2.0 * (tmp1 + tmp2)\n");
  fprintf(fin, "tmp1 = y*z\n");
  fprintf(fin, "tmp2 = x*w\n");
  fprintf(fin, "m21 = 2.0 * (tmp1 + tmp2)\n");
  fprintf(fin, "m12 = 2.0 * (tmp1 - tmp2)\n");
  fprintf(fin, "matrix <m00,m01,m02,m10,m11,m12,m20,m21,m22>\n");
  fprintf(fin, "#end\n");
  */
  fprintf(fin,"#macro randomColor(r,g,b)\n");
  fprintf(fin, "texture{ pigment{ color rgb<r,g,b> transmit para_trans} finish { phong para_phong }}\n");//transmit for transparency, higher for more transparent.
  fprintf(fin,"#end\n");
  fprintf(fin,"#macro octantSuper(ep1,ep2,rx,ry,rz,sign1,sign2,sign3)\n");
  fprintf(fin,"    superellipsoid{ <ep1,ep2> \n#if (singleColor = false)\n randomColor(rand(Random_r),rand(Random_g),rand(Random_b)) \n#end\nmyremain(sign1,sign2,sign3) scale <rx,ry,rz>}\n#end\n");


  fprintf(fin,"#macro polySuperEllipsoid(ep1,ep2,a1,a2,b1,b2,c1,c2,t1,t2,t3,r1,r2,r3,r,g,b)\nunion{\n");

  for(int k=0;k<2;k++){
    for(int j=0;j<2;j++){
      for(int i=0;i<2;i++)
        //print each octant
        fprintf(fin,"octantSuper(ep1,ep2,%s,%s,%s,%d,%d,%d)\n",i>0?"a1":"a2",j>0?"b1":"b2",k>0?"c1":"c2",(2*i-1),(2*j-1),(2*k-1));
    }
  }
  fprintf(fin,"rotate<r1,r2,r3>\ntranslate<t1,t2,t3>\n#if (singleColor = true)\nrandomColor(r,g,b)\n#end\n}\n#end\n");//translate should be following rotate.


	//fprintf(fin,"volume\t%e\n",volume);
	//fprintf(fin,"centroid\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2]);

  //output superquadrics
  Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	//particlesVolume = 0;
	Quaternionr quaternion;
	double q0,q1,q2,q3,phi,theta,psi;
	Vector3r pos, color,translateV; Vector6r rxyz;
	Vector2r eps;
  std::vector<int> ids_out = ids;
  if (!ids_out.size()){//output all particles by default
    for ( ; bi!=biEnd; ++bi ) ids_out.push_back((*bi)->id);
  }
	for (std::vector<int>::iterator id = ids_out.begin(); id!=ids_out.end(); ++id)
	{
		const shared_ptr<Body>& b=Body::byId(*id);//const shared_ptr<Body>& b = *bi;
		if (b->shape->getClassName()=="PolySuperellipsoid"){
			const shared_ptr<PolySuperellipsoid>& A = SUDODEM_PTR_CAST<PolySuperellipsoid> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());
			State* state=b->state.get();
		        //Euler Angles from Quaternion
		        quaternion = state->ori;
            quaternion.normalize();
		        q0 = quaternion.w();
		        q1 = quaternion.x();
		        q2 = quaternion.y();
		        q3 = quaternion.z();
		        phi = atan2(2.0*(q0*q1 + q2*q3), 1.0 - 2.0*(q1*q1 + q2*q2))/Mathr::PI*180.0;
		        theta = asin(2.0*(q0*q2 - q3*q1))/Mathr::PI*180.0;
		        psi = atan2(2.0*(q0*q3 + q1*q2), 1.0 - 2.0*(q2*q2 + q3*q3))/Mathr::PI*180.0;
            //testing
            //Roll pitch and yaw in Radians
          /*float roll = phi/180.0*Mathr::PI, pitch = theta/180.0*Mathr::PI, yaw = psi/180.0*Mathr::PI;
            Quaternionr q;
            q = AngleAxisr(yaw, Vector3r::UnitZ())
                * AngleAxisr(pitch, Vector3r::UnitY())
                * AngleAxisr(roll, Vector3r::UnitX());
            q.normalize();
            std::cout << "Quaternion1" << std::endl << quaternion.coeffs() << std::endl;
            std::cout << "Quaternion2" << std::endl << q.coeffs() << std::endl;
            Vector3r euler1 = q.toRotationMatrix().eulerAngles(0,1,2);
            std::cout << "roll, pitch, yaw"<< std::endl << euler1 << std::endl;
            Vector3r euler = quaternion.toRotationMatrix().eulerAngles(0,1,2);
            euler = euler/Mathr::PI*180.0;
            std::cout << "Euler from quaternion in roll, pitch, yaw"<< std::endl << euler << std::endl;
            */
		        pos = state->pos;
		        color = A->color;
		        rxyz = A->getrxyz();
		        eps = A->geteps();
		        //output
            translateV = pos - quaternion.toRotationMatrix()*A->getMassCenter();
            rxyz *= scale;
            translateV *= scale;
            fprintf(fin,"polySuperEllipsoid(%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%f,%f,%f,%f,%f,%f)\n",
            eps(0),eps(1),rxyz(0),rxyz(1),rxyz(2),rxyz(3),rxyz(4),rxyz(5),translateV(0),translateV(1),translateV(2),phi,theta,psi,color(0),color(1),color(2));
		}
		if (b->shape->getClassName()=="Sphere"){
			const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b->shape);
			State* state=b->state.get();
		        pos = state->pos;
		        color = A->color;
		        //output
		        fprintf(fin,"sphere{ <%f,%f,%f>,%f\n",pos(0),pos(1),pos(2),A->radius);
		        fprintf(fin, "  texture{ pigment{ color rgb<%f,%f,%f>}\n",color(0),color(1),color(2));
		        fprintf(fin, "           finish { phong 1}}\n");
		        fprintf(fin, "}\n");
		}

	}


	fclose(fin);


}
void outputVTK(const char* filename,int slices){
    std::cout << "writing Polydata into a VTK file" << std::endl;
    std::ofstream f;
    f.open(filename);
    //writing the file head
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Sway data processing" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;

    //output superquadrics
    Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	//particlesVolume = 0;
	Quaternionr quaternion;
	Vector3r pos, color,rxyz;
	Vector2r eps;

	int w=2*slices,h=slices;
    //point counter
    unsigned int point_counter = 0;
    vector<Vector3r> vertices;//discretized vertices of a particle

	for ( ; bi!=biEnd; ++bi )
	{
		const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		if (b->shape->getClassName()=="Superquadrics"){
            point_counter ++;
			const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());
			State* state=b->state.get();
            ////
		    int i,j;
		    double a=0.0,b=0.0;
		    double hStep=M_PI/(h-1);
		    double wStep=2*M_PI/w;

            Matrix3r Rot = (state->ori).toRotationMatrix();
            pos = state->pos;
            rxyz = A->getrxyz();
		    eps = A->geteps();
	        //surface discretization
		    for(a=0.0,i=0;i<h;i++,a+=hStep)
		      for(b=0.0,j=0;j<w;j++,b+=wStep)
		      {
                Vector3r Surf = A->getSurface(Vector2r(b,a-M_PI_2l));

                /*
                std::cout<<rxyz(0)<<" "<<rxyz(1)<<" "<<rxyz(2)<<" "<<eps(0)<<" "<<eps(1)<<" "<<phi(0)<<" "<<phi(1)<<std::endl;
                Surf(0) = Sign(cos(phi(0)))*rxyz(0)*pow(fabs(cos(phi(0))), eps(0))*pow(fabs(cos(phi(1))), eps(1));
	            Surf(1) = Sign(sin(phi(0)))*rxyz(1)*pow(fabs(sin(phi(0))), eps(0))*pow(fabs(cos(phi(1))), eps(1));
	            Surf(2) = Sign(sin(phi(1)))*rxyz(2)*pow(fabs(sin(phi(1))), eps(1));
                std::cout<<a<<" "<<b<<" "<<Surf(0)<<" "<<Surf(1)<<" "<<Surf(2)<<std::endl;
                */
                Surf = Rot*Surf + pos;
                //std::cout<<a<<" "<<b<<" "<<Surf[0]<<" "<<Surf[1]<<" "<<Surf[2]<<std::endl;
			    vertices.push_back(Surf);
		      }
		    }


	}
    //output data to the vtk file
    f << "POINTS   "<<vertices.size() <<" float" << std::endl;//index is from 0 in VTK
    //f << "0   0   0" << std::endl;//putting a point to take the place.
    for ( vector<Vector3r>::iterator it = vertices.begin(); it != vertices.end(); ++ it)
    {

        f <<  std::setprecision(12) << (*it)[0] << " " << std::setprecision(12) << (*it)[1] << " " << std::setprecision(12) << (*it)[2]<< std::endl;

    }
    //polygons
    f << "POLYGONS " << point_counter*(h-1)*w <<" "<<point_counter*(h-1)*w*5<< std::endl;
    unsigned long start_index = 0;
    for(unsigned long kk = 0; kk < point_counter; kk++ ){
        start_index = kk*w*h;
        //polygons of surface slices
        int i,j;
		for(i=0;i<h-1;i++)
		{
		    for(j=0;j<w-1;j++)
		    {
                f << "4 " << start_index+i*w+j << " "<< start_index+i*w+j+1 << " "<< start_index+(i+1)*w+j+1 << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		    }
		    f << "4 " << start_index+i*w+j << " "<< start_index+i*w << " "<< start_index+(i+1)*w << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		}
    }
	f.close();


}
//**********************************************************************************
Vector3r getSurface_PolySuperEllipsoid2(Vector2r phi,double R,Vector6r eps, Vector6r pShape){
  //vector6r eps:eps1_1,eps1_2,eps2_1,eps2_2,eps2_3,eps2_4

  assert(phi[0]>=0);
  double x,y,z,eps1,eps2;
  x = R*Sign(cos(phi(0)));
  y = R*Sign(sin(phi(0)));
  z = R*Sign(sin(phi(1)));
  eps1 = (z>0?eps[0]:eps[1]);
  eps2 = (x>0?(y>0?eps[2]:eps[5]):(y>0?eps[3]:eps[4]));
  x *= pShape[(x>0?0:1)]*pow(fabs(cos(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
  y *= pShape[(y>0?2:3)]*pow(fabs(sin(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
  z *= pShape[(z>0?4:5)]*pow(fabs(sin(phi(1))), eps2);
  return Vector3r(x,y,z);
}
Vector3r getSurface_PolySuperEllipsoid(Vector2r phi,double R,Vector2r eps, Vector6r pShape){
  assert(phi[0]>=0);
  double x,y,z;
  x = R*Sign(cos(phi(0)))*pow(fabs(cos(phi(0))), eps[0])*pow(fabs(cos(phi(1))), eps[1]);
  y = R*Sign(sin(phi(0)))*pow(fabs(sin(phi(0))), eps[0])*pow(fabs(cos(phi(1))), eps[1]);
  z = R*Sign(sin(phi(1)))*pow(fabs(sin(phi(1))), eps[1]);
  x *= pShape[(x>0?0:1)];
  y *= pShape[(y>0?2:3)];
  z *= pShape[(z>0?4:5)];
  return Vector3r(x,y,z);
}
void outputVTK_PolySuperEllipsoid2(const char* filename,double R,Vector6r eps,Vector6r pShape,double slices){
    //a poly-ellipsoid is controled by one size parameter and six shape parameters
    //Vector6r pshape:alpha1, alpha2, beta1, beta2,gamma1,gamma2
    std::cout << "writing Polydata into a VTK file" << std::endl;
    std::ofstream f;
    f.open(filename);
    //writing the file head
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Sway data processing" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;

//	Quaternionr quaternion;
//	Vector3r pos, color,rxyz;

	int w=2*slices,h=slices;
    //point counter
    unsigned int point_counter = 0;
    vector<Vector3r> vertices;//discretized vertices of a particle

//	for ( ; bi!=biEnd; ++bi )
//	{
		//const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		//if (b->shape->getClassName()=="Superquadrics"){
      point_counter ++;
		//	const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
		//	State* state=b->state.get();
            ////
		    int i,j;
		    double a=0.0,b=0.0;
		    double hStep=M_PI/(h-1);
		    double wStep=2*M_PI/w;

            //Matrix3r Rot = (state->ori).toRotationMatrix();
            //pos = state->pos;
	        //surface discretization
		    for(a=0.0,i=0;i<h;i++,a+=hStep)
		      for(b=0.0,j=0;j<w;j++,b+=wStep)
		      {
                Vector3r Surf = getSurface_PolySuperEllipsoid2(Vector2r(b, a-M_PI_2l), R,eps,pShape);

                //Surf = Rot*Surf + pos;
                //std::cout<<a<<" "<<b<<" "<<Surf[0]<<" "<<Surf[1]<<" "<<Surf[2]<<std::endl;
			    vertices.push_back(Surf);
		      }
		    //}


//	}
    //output data to the vtk file
    f << "POINTS   "<<vertices.size() <<" float" << std::endl;//index is from 0 in VTK
    //f << "0   0   0" << std::endl;//putting a point to take the place.
    for ( vector<Vector3r>::iterator it = vertices.begin(); it != vertices.end(); ++ it)
    {

        f <<  std::setprecision(12) << (*it)[0] << " " << std::setprecision(12) << (*it)[1] << " " << std::setprecision(12) << (*it)[2]<< std::endl;

    }
    //polygons
    f << "POLYGONS " << point_counter*(h-1)*w <<" "<<point_counter*(h-1)*w*5<< std::endl;
    unsigned long start_index = 0;
    for(unsigned long kk = 0; kk < point_counter; kk++ ){
        start_index = kk*w*h;
        //polygons of surface slices
        int i,j;
		for(i=0;i<h-1;i++)
		{
		    for(j=0;j<w-1;j++)
		    {
                f << "4 " << start_index+i*w+j << " "<< start_index+i*w+j+1 << " "<< start_index+(i+1)*w+j+1 << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		    }
		    f << "4 " << start_index+i*w+j << " "<< start_index+i*w << " "<< start_index+(i+1)*w << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		}
    }
	f.close();


}
void outputVTK_PolySuperEllipsoid(const char* filename,double R,Vector2r eps,Vector6r pShape,double slices){
    //a poly-ellipsoid is controled by one size parameter and six shape parameters
    //Vector6r pshape:alpha1, alpha2, beta1, beta2,gamma1,gamma2
    std::cout << "writing Polydata into a VTK file" << std::endl;
    std::ofstream f;
    f.open(filename);
    //writing the file head
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Sway data processing" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;

//	Quaternionr quaternion;
//	Vector3r pos, color,rxyz;

	int w=2*slices,h=slices;
    //point counter
    unsigned int point_counter = 0;
    vector<Vector3r> vertices;//discretized vertices of a particle

//	for ( ; bi!=biEnd; ++bi )
//	{
		//const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		//if (b->shape->getClassName()=="Superquadrics"){
      point_counter ++;
		//	const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
		//	State* state=b->state.get();
            ////
		    int i,j;
		    double a=0.0,b=0.0;
		    double hStep=M_PI/(h-1);
		    double wStep=2*M_PI/w;

            //Matrix3r Rot = (state->ori).toRotationMatrix();
            //pos = state->pos;
	        //surface discretization
		    for(a=0.0,i=0;i<h;i++,a+=hStep)
		      for(b=0.0,j=0;j<w;j++,b+=wStep)
		      {
                Vector3r Surf = getSurface_PolySuperEllipsoid(Vector2r(b, a-M_PI_2l), R,eps,pShape);

                //Surf = Rot*Surf + pos;
                //std::cout<<a<<" "<<b<<" "<<Surf[0]<<" "<<Surf[1]<<" "<<Surf[2]<<std::endl;
			    vertices.push_back(Surf);
		      }
		    //}


//	}
    //output data to the vtk file
    f << "POINTS   "<<vertices.size() <<" float" << std::endl;//index is from 0 in VTK
    //f << "0   0   0" << std::endl;//putting a point to take the place.
    for ( vector<Vector3r>::iterator it = vertices.begin(); it != vertices.end(); ++ it)
    {

        f <<  std::setprecision(12) << (*it)[0] << " " << std::setprecision(12) << (*it)[1] << " " << std::setprecision(12) << (*it)[2]<< std::endl;

    }
    //polygons
    f << "POLYGONS " << point_counter*(h-1)*w <<" "<<point_counter*(h-1)*w*5<< std::endl;
    unsigned long start_index = 0;
    for(unsigned long kk = 0; kk < point_counter; kk++ ){
        start_index = kk*w*h;
        //polygons of surface slices
        int i,j;
		for(i=0;i<h-1;i++)
		{
		    for(j=0;j<w-1;j++)
		    {
                f << "4 " << start_index+i*w+j << " "<< start_index+i*w+j+1 << " "<< start_index+(i+1)*w+j+1 << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		    }
		    f << "4 " << start_index+i*w+j << " "<< start_index+i*w << " "<< start_index+(i+1)*w << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		}
    }
	f.close();


}
//**********************************************************************************
Vector3r getSurface_PolyEllipsoid(double theta, double phi,double R,Vector6r pShape){
  //theta:[0,2pi)
  //phi:[--.5pi,0.5pi]
  assert(theta>=0);
  double x,y,z;
  x = R*cos(theta)*cos(phi);
  y = R*sin(theta)*cos(phi);
  z = R*sin(phi);
  x *= pShape[(x>0?0:1)];
  y *= pShape[(y>0?2:3)];
  z *= pShape[(z>0?4:5)];
  return Vector3r(x,y,z);
}
void outputVTK_polyEllipsoid(const char* filename,double R,Vector6r pShape,double slices){
    //a poly-ellipsoid is controled by one size parameter and six shape parameters
    //Vector6r pshape:alpha1, alpha2, beta1, beta2,gamma1,gamma2
    std::cout << "writing Polydata into a VTK file" << std::endl;
    std::ofstream f;
    f.open(filename);
    //writing the file head
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Sway data processing" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;

//	Quaternionr quaternion;
//	Vector3r pos, color,rxyz;

	int w=2*slices,h=slices;
    //point counter
    unsigned int point_counter = 0;
    vector<Vector3r> vertices;//discretized vertices of a particle

//	for ( ; bi!=biEnd; ++bi )
//	{
		//const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		//if (b->shape->getClassName()=="Superquadrics"){
      point_counter ++;
		//	const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
		//	State* state=b->state.get();
            ////
		    int i,j;
		    double a=0.0,b=0.0;
		    double hStep=M_PI/(h-1);
		    double wStep=2*M_PI/w;

            //Matrix3r Rot = (state->ori).toRotationMatrix();
            //pos = state->pos;
	        //surface discretization
		    for(a=0.0,i=0;i<h;i++,a+=hStep)
		      for(b=0.0,j=0;j<w;j++,b+=wStep)
		      {
                Vector3r Surf = getSurface_PolyEllipsoid(b, a-M_PI_2l, R,pShape);

                //Surf = Rot*Surf + pos;
                //std::cout<<a<<" "<<b<<" "<<Surf[0]<<" "<<Surf[1]<<" "<<Surf[2]<<std::endl;
			    vertices.push_back(Surf);
		      }
		    //}


//	}
    //output data to the vtk file
    f << "POINTS   "<<vertices.size() <<" float" << std::endl;//index is from 0 in VTK
    //f << "0   0   0" << std::endl;//putting a point to take the place.
    for ( vector<Vector3r>::iterator it = vertices.begin(); it != vertices.end(); ++ it)
    {

        f <<  std::setprecision(12) << (*it)[0] << " " << std::setprecision(12) << (*it)[1] << " " << std::setprecision(12) << (*it)[2]<< std::endl;

    }
    //polygons
    f << "POLYGONS " << point_counter*(h-1)*w <<" "<<point_counter*(h-1)*w*5<< std::endl;
    unsigned long start_index = 0;
    for(unsigned long kk = 0; kk < point_counter; kk++ ){
        start_index = kk*w*h;
        //polygons of surface slices
        int i,j;
		for(i=0;i<h-1;i++)
		{
		    for(j=0;j<w-1;j++)
		    {
                f << "4 " << start_index+i*w+j << " "<< start_index+i*w+j+1 << " "<< start_index+(i+1)*w+j+1 << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		    }
		    f << "4 " << start_index+i*w+j << " "<< start_index+i*w << " "<< start_index+(i+1)*w << " "<< start_index+(i+1)*w+j <<" "<<std::endl;
		}
    }
	f.close();


}
//**********************************************************************************
//output geometrical info of an assembly
void outputParticles(const char* filename){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write


	//fprintf(fin,"volume\t%e\n",volume);
	fprintf(fin,"#For poly-superellipsoids(15 columns:rx\try\trz\teps1\teps2\tx\ty\tz\tq0\tq1\tq2\tq3); for superellipsoids (12 columns: rx1,rx2\try1,ry2\trz1,rz2\teps1\teps2\tx\ty\tz\tq0\tq1\tq2\tq3); for spheres (4 columns: r,x,y,z)\n");

        //output superquadrics
        Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	//particlesVolume = 0;
	Quaternionr quaternion;
	double q0,q1,q2,q3,phi,theta,psi;
	Vector3r pos, color;
	Vector2r eps;
	for ( ; bi!=biEnd; ++bi )
	{
		const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		if (b->shape->getClassName()=="Superquadrics"){
			const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());
			State* state=b->state.get();
      //Euler Angles from Quaternion
      quaternion = state->ori;
      q0 = quaternion.w();
      q1 = quaternion.x();
      q2 = quaternion.y();
      q3 = quaternion.z();

      pos = state->pos;
  	  Vector3r rxyz = A->getrxyz();
      eps = A->geteps();
      //output
      fprintf(fin,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",rxyz(0),rxyz(1),rxyz(2),eps(0),eps(1),pos(0),pos(1),pos(2),q0,q1,q2,q3);
		}else if(b->shape->getClassName()=="PolySuperellipsoid"){
      const shared_ptr<PolySuperellipsoid>& A = SUDODEM_PTR_CAST<PolySuperellipsoid> (b->shape);
      State* state=b->state.get();
      quaternion = state->ori;
      q0 = quaternion.w();
      q1 = quaternion.x();
      q2 = quaternion.y();
      q3 = quaternion.z();

      pos = state->pos;
      Vector6r rxyz = A->getrxyz();
      eps = A->geteps();
      //output
      fprintf(fin,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",rxyz(0),rxyz(1),rxyz(2),rxyz(3),rxyz(4),rxyz(5),eps(0),eps(1),pos(0),pos(1),pos(2),q0,q1,q2,q3);
    }else if(b->shape->getClassName()=="Sphere"){
      const shared_ptr<Sphere>& A = SUDODEM_PTR_CAST<Sphere> (b->shape);
      State* state=b->state.get();
      pos = state->pos;
      fprintf(fin, "%e\t%e\t%e\t%e\n", A->radius,pos(0),pos(1),pos(2));
    }
	}
	fclose(fin);
}
//output particles by a list of ids
void outputParticlesbyIds(const char* filename, std::vector<int> ids){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write


	//fprintf(fin,"volume\t%e\n",volume);
	fprintf(fin,"#rx\try\trz\teps1\teps2\tx\ty\tz\tq0\tq1\tq2\tq3\n");

    //output superquadrics
    Scene* scene=Omega::instance().getScene().get();

	//particlesVolume = 0;
	Quaternionr quaternion;
	double q0,q1,q2,q3,phi,theta,psi;
	Vector3r pos, color,rxyz;
	Vector2r eps;
	for (std::vector<int>::iterator id = ids.begin(); id!=ids.end(); ++id)
	{
        const shared_ptr<Body>& b=Body::byId(*id);
		if (b->shape->getClassName()=="Superquadrics"){
			const shared_ptr<Superquadrics>& A = SUDODEM_PTR_CAST<Superquadrics> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());
			State* state=b->state.get();
		        //Euler Angles from Quaternion
		        quaternion = state->ori;
		        q0 = quaternion.w();
		        q1 = quaternion.x();
		        q2 = quaternion.y();
		        q3 = quaternion.z();

		        pos = state->pos;
		    	rxyz = A->getrxyz();
		        eps = A->geteps();
		        //output
                fprintf(fin,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",rxyz(0),rxyz(1),rxyz(2),eps(0),eps(1),pos(0),pos(1),pos(2),q0,q1,q2,q3);
		}

	}


	fclose(fin);


}

//output positions of walls
void outputWalls(const char* filename){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write

    Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	Vector3r pos;
	for ( ; bi!=biEnd; ++bi )
	{
		const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		if (b->shape->getClassName()=="Wall"){

			State* state=b->state.get();
	        pos = state->pos;
	        //output
            fprintf(fin,"%e\t%e\t%e\n",pos(0),pos(1),pos(2));
		}

	}


	fclose(fin);


}

//**********************************************************************************//


BOOST_PYTHON_MODULE(_superquadrics_utils){
	// http://numpy.scipy.org/numpydoc/numpy-13.html mentions this must be done in module init, otherwise we will crash
	import_array();

	SUDODEM_SET_DOCSTRING_OPTS;
    py::def("Superquadrics_test",Superquadrics_test,"test the computation of points on a superquadric given a normal vector.");
    py::def("PositionPB",PositionPB,"generate a pair of particles for test.");
    py::def("EST_Support",EST_Support,"Estimate the time consumption of the support function for a particle. For test only.");
    //py::def("MonteCarlo_test",MonteCarlo_test,(py::args("num_max")=1000000),"test the computation of particle geometric quantities (volume, mass center, and moment of inertia) by MonteCarlo simulation.");
    py::def("outputWalls", outputWalls, "output positions of walls.");
    py::def("getSurfArea", getSurfArea, "getSurfArea(id,w,h): get the surface area of a super-ellipsoid by a given particle id with resolution w and h (e.g., w=10,h=10). Larger w and h will yield more accurate results.");
    py::def("outputParticles", outputParticles, "output particles.");
    py::def("outputParticlesbyIds", outputParticlesbyIds, "output particles by a list of ids.");
	  py::def("NewSuperquadrics", NewSuperquadrics, "Generate a Superquadric.");
    py::def("NewPolySuperellipsoid", NewPolySuperellipsoid,(py::args("inertiaScale")=1.0), "Generate a PolySuperellipsoid.");
    py::def("NewPolySuperellipsoid_rot", NewPolySuperellipsoid_rot, (py::args("inertiaScale")=1.0),"Generate a PolySuperellipsoid with specified quaternion components: q0(w),q1(x),q2(y),q3(z)..");
	  py::def("NewSuperquadrics2", NewSuperquadrics2, "Generate a Superquadric with isSphere option.");
	  py::def("NewSuperquadrics_rot", NewSuperquadrics_rot, "Generate a Super-ellipsoid with specified quaternion components: q0(w),q1(x),q2(y),q3(z).");
    py::def("NewSuperquadrics_rot2", NewSuperquadrics_rot2, "Generate a Super-ellipsoid with specified quaternion components: q0(w),q1(x),q2(y),q3(z) with isSphere option.");
	  py::def("outputPOV",outputPOV,"output Superquadrics into a POV file for post-processing using POV-ray.");
    py::def("PolySuperellipsoidPOV",PolySuperellipsoidPOV,(py::args("ids")=NULL,py::args("scale")=1.0),"output PolySuperellipsoids into a POV file for post-processing using POV-ray.");
    py::def("outputVTK",outputVTK,"output Superquadrics into a VTK file for post-processing using Paraview.");
    py::def("outputVTK_polyEllipsoid",outputVTK_polyEllipsoid,"output Poly-ellipsoid into a VTK file for post-processing using Paraview.[Only for test.]");
    py::def("outputVTK_polySuperEllipsoid",outputVTK_PolySuperEllipsoid,"output Poly-Superellipsoid into a VTK file for post-processing using Paraview.[Only for test.]");
    py::def("outputVTK_polySuperEllipsoid2",outputVTK_PolySuperEllipsoid2,"output Poly-Superellipsoid into a VTK file for post-processing using Paraview.[Only for test.]");
}
