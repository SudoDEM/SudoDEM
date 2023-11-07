//@2016 Sway Zhao,  zhswee@gmail.com
//
#include"sudodem/pkg/dem/GJKParticle.hpp"

#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/lib/pyutil/doc_opts.hpp>
#include<sudodem/lib/voro++/voro++.hh>

#include<cmath>

#include<numpy/ndarrayobject.h>
#include <boost/concept_check.hpp>

namespace py = boost::python;
using namespace voro;
//***********************************************************************************

int Sign(double f){ if(f<0) return -1; if(f>0) return 1; return 0; }
//***********************************************************************************
void ConvexHull(std::vector<Vector3r> verts,std::vector<Vector3r>& CH_vertices,T_MultiIndexBuf& CH_facetIndexBuf )
{   //std::cout<<"construct the convex hull..."<<std::endl;
	//std::cout<<"Input points:"<<std::endl;
	//for(int i=0;i<verts.size();i++){
	//	std::cout<<verts[i]<<std::endl;
	//}
	int curlong, totlong, exitcode;
    unsigned int count = verts.size();
    assert(count);
    facetT *facet;
    vertexT *vertex;
    vertexT **vertexp;

    std::vector<coordT> array;
	T_IndexBuf index;
    unsigned int i;
	//std::cout<<"size of m_vertices2="<<m_vertices2.size()<<std::endl;
	CH_vertices.clear();
    for (i = 0; i < count; i++)
	{
	array.push_back(verts[i][0]);
	array.push_back(verts[i][1]);
	array.push_back(verts[i][2]);
	//std::cout<<"x=="<<verts[i][0]<<std::endl;
    CH_vertices.push_back(verts[i]);
	index.push_back(i);
    }
	//std::cout<<"Vertices from the convex hull:"<<std::endl;
	//for(int i=0;i<m_vertices2.size();i++){
	//	std::cout<<m_vertices2[i]<<std::endl;
	//}
    double* a = &array[0];
    qh_init_A(stdin, stdout, stderr, 0, NULL);
    if ((exitcode = setjmp(qh errexit)))
	{
		exit(exitcode);
	}
    qh_initflags(options);
    qh_init_B(a, array.size()/3, 3, False);
    qh_qhull();
    qh_check_output();

    FORALLfacets
	{
		setT *vertices = qh_facet3vertex(facet);//get vertices of the current facet

		T_IndexBuf  facetIndices;

		FOREACHvertex_(vertices)//taverse the vertices
		{
			facetIndices.push_back(index[qh_pointid(vertex->point)]);
		}
		CH_facetIndexBuf.push_back(facetIndices);
    }

    qh NOerrexit = True;
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
}
//
void randomGeometry(int num_vertices,double m_radius,std::vector<Vector3r>&m_vertices){//generate vertices of a polyhedral particle based on a unit sphere
	SD_Scalar theta,phi;
	double x,y,z;
	for(int i=0;i<num_vertices;i++){
		theta = Mathr::TWO_PI*Mathr::UnitRandom();
		phi = Mathr::TWO_PI*0.5*Mathr::UnitRandom();
		x = m_radius*sin(phi)*cos(theta);
		y = m_radius*sin(phi)*sin(theta);
		z = m_radius*cos(phi);
		m_vertices.push_back(Vector3r(x,y,z));
	}
}

//**************************************************************************
/* Generator of randomly shaped polyhedron based on Voronoi tessellation*/

std::vector<Vector3r> GenerateRandomGeometry(){
	int seed = time(NULL);//Seed for random generator
	srand(seed);

	vector<Vector3r> nuclei;
	nuclei.push_back(Vector3r(5.,5.,5.));
	Vector3r trial;
	int iter = 0;
	bool isOK;
	//fill box 5x5x5 with randomly located nuclei with restricted minimal mutual distance 0.75 which gives approximate mean mutual distance 1;
	double dist_min2 = 0.75*0.75;
	while(iter<500){
		isOK = true;
		iter++;
		trial = Vector3r(double(rand())/RAND_MAX*5.+2.5,double(rand())/RAND_MAX*5.+2.5,double(rand())/RAND_MAX*5.+2.5);
		for(int i=0;i< (int) nuclei.size();i++) {
			isOK =  (nuclei[i] - trial).squaredNorm() > dist_min2;
			if (!isOK) break;
		}

		if(isOK){
			iter = 0;
			nuclei.push_back(trial);
		}
	}
	/*cout<<"random points"<<endl;
	for(std::vector<Vector3r>::iterator v=nuclei.begin();v!=nuclei.end();v++){
		cout<<*v<<endl;
	}*/
	std::vector<Vector3r> vert;
	//perform Voronoi tessellation
  //nuclei.erase(nuclei.begin());
	const double xmin = 2.5;
	const double xmax = 7.5;
	const double ymin = 2.5;
	const double ymax = 7.5;
	const double zmin = 2.5;
	const double zmax = 7.5;
	bool xpbc = false;
	bool ypbc = false;
	bool zpbc = false;
	int nx(1), ny(1), nz(1);
	container con(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, xpbc, ypbc, zpbc, 8);
	int id=0;
	for(std::vector<Vector3r>::iterator it = nuclei.begin();it!=nuclei.end();it++){
		con.put(id++, (*it)(0),(*it)(1),(*it)(2));
	}
	// loop over all voronoi cells
	 c_loop_all cla(con);
	 if(cla.start())
	 {
			 //std::cout << "started\n" << std::flush;
			 do
			 {
					 voronoicell_neighbor c;

					 if(con.compute_cell(c,cla))
					 {

							 //std::cout << "computed"  << std::endl;
							 double xc = 0;
							 double yc = 0;
							 double zc = 0;
							 // Get the position of the current particle under consideration
							 cla.pos(xc,yc,zc);

							 unsigned int id = cla.pid();
							 if(id == 0){//select the first vell enclosing the origin
								 //cell volume
								 //c.volume();
								// std::vector<int> f; // list of face vertices (bracketed, as ID)
								//c.face_vertices(f);

								std::vector<double> vertices;   // all vertices for this cell
								c.vertices(xc,yc,zc, vertices);

								//std::vector<int> w; // neighbors of faces
								//c.neighbors(w);

								/*for(std::vector<double>::iterator v=vertices.begin();v!=vertices.end();v++){
									vert.push_back(Vector3r(*v,*(++v),*(++v)));
								}*/
								//cout<<"test v="<<endl;
								for(int i=0;i<vertices.size()/3;i+=3){
									vert.push_back(Vector3r(vertices[i],vertices[i+1],vertices[i+2]));
									//cout<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[i+2]<<endl;
								}
								break;
							 }
						}
				} while (cla.inc());
		}


	//resize and rotate the voronoi cell
	/*
	Quaternionr Rot(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	Rot.normalize();
	for(int i=0; i< (int) v.size();i++) {
		v[i] = Rot*(Vector3r(v[i][0]*size[0],v[i][1]*size[1],v[i][2]*size[2]));
	}

	//to avoid patological cases (that should not be present, but CGAL works somehow unpredicable)
	if (v.size() < 8) {
		cout << "wrong " << v.size() << endl;
		v.clear();
		seed = rand();
		GenerateRandomGeometry();
	}
	*/
	/*cout<<"vertices= "<<endl;
	for(std::vector<Vector3r>::iterator v=vert.begin();v!=vert.end();v++){
		cout<<*v<<endl;
	}*/
	return vert;
}


void gen_randomPolyhedron(std::vector<Vector3r>&m_vertices, Vector3r extent){
	/*m_vertices:
	  extent:
		*/
	if(m_vertices.size()==0){
		//int num_vertices = 8;
		/*#define NON_RANDOM_VER
		#ifdef NON_RANDOM_VER
		for (int i = 0; i != num_vertices; ++i)
		{
			m_vertices.push_back(Vector3r(0.,0.,0.));
		}
		m_vertices[0] = Vector3r(0.,0.,0.);
		m_vertices[1] = Vector3r(m_radius,0.,0.);
		m_vertices[2] = Vector3r(m_radius,m_radius,0.);
		m_vertices[3] = Vector3r(0.,m_radius,0.);
		m_vertices[4] = Vector3r(0.,0.,m_radius);
		m_vertices[5] = Vector3r(m_radius,0.,m_radius);
		m_vertices[6] = Vector3r(m_radius,m_radius,m_radius);
		m_vertices[7] = Vector3r(0.,m_radius,m_radius);
		//m_vertices[8] = Vector3r(0.5,0.5,1.5);
		#else*/
		//randomGeometry(num_vertices,m_radius,m_vertices);
		//#endif
		m_vertices = GenerateRandomGeometry();
		for(std::vector<Vector3r>::iterator v=m_vertices.begin(); v!=m_vertices.end();v++) {
			*v = Vector3r((*v)(0)*extent(0),(*v)(1)*extent(1),(*v)(2)*extent(2));
		}
	}

}
//new GJKParticle

shared_ptr<Body> NewGJKParticle(vector<Vector3r> V,int shapeType, SD_Scalar radius, SD_Scalar height,shared_ptr<Material> mat,bool rotate){
	shared_ptr<Body> body(new Body);
	body->material=mat;
    if(1==shapeType){//polyhedron
        gen_randomPolyhedron(V,Vector3r(radius,radius,radius));
    }
	body->shape=shared_ptr<GJKParticle>(new GJKParticle(V,shapeType,radius,height));
	GJKParticle* A = static_cast<GJKParticle*>(body->shape.get());
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
shared_ptr<Body> NewGJKParticle2(vector<Vector3r> V,int shapeType, SD_Scalar radius, SD_Scalar margin, SD_Scalar height,shared_ptr<Material> mat,bool rotate){
	shared_ptr<Body> body(new Body);
	body->material=mat;
    /*if(1==shapeType){//polyhedron
        gen_randomPolyhedron(V,Vector3r(1.,1.,1.)*radius);
    }*/
	body->shape=shared_ptr<GJKParticle>(new GJKParticle(V,shapeType,radius,height,margin));
	GJKParticle* A = static_cast<GJKParticle*>(body->shape.get());
	A->isSphere = false;//non-spherical by default.
  //cout<<"shape pointer="<<A->m_shape<<"count="<<A->m_shape.use_count()<<endl;
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
//Spheres
shared_ptr<Body> GJKSphere(SD_Scalar radius, SD_Scalar margin, shared_ptr<Material> mat){
  vector<Vector3r> V;
  return NewGJKParticle2(V,0,radius, margin,0,mat,false);
}
//Polyhedron
shared_ptr<Body> GJKPolyhedron(vector<Vector3r> V,Vector3r extent, SD_Scalar margin,shared_ptr<Material> mat,bool rotate){
  gen_randomPolyhedron(V,extent);
	/*cout<<"vertices= "<<endl;
	for(std::vector<Vector3r>::iterator v=V.begin();v!=V.end();v++){
		cout<<*v<<endl;
	}*/

  return NewGJKParticle2(V,1,0,margin,0,mat,rotate);
}
//Cones
shared_ptr<Body> GJKCone(SD_Scalar radius, SD_Scalar height, SD_Scalar margin, shared_ptr<Material> mat,bool rotate){
	vector<Vector3r> V;
  return NewGJKParticle2(V,2,radius, margin,height,mat,rotate);
}
//Cylinders
shared_ptr<Body> GJKCylinder(SD_Scalar radius, SD_Scalar height, SD_Scalar margin,shared_ptr<Material> mat,bool rotate){
	vector<Vector3r> V;
  return NewGJKParticle2(V,3,radius, margin,height,mat,rotate);
}
//cuboids
shared_ptr<Body> GJKCuboid(Vector3r extent, SD_Scalar margin,shared_ptr<Material> mat,bool rotate){
	vector<Vector3r> V;
	V.push_back(extent);
  return NewGJKParticle2(V, 4, 0, margin,0,mat,rotate);
}
/*
shared_ptr<Body> NewGJKParticle2(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,bool rotate, bool isSphere){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<GJKParticle>(new GJKParticle(x,y,z,ep1,ep2));
	GJKParticle* A = static_cast<GJKParticle*>(body->shape.get());

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
shared_ptr<Body> NewGJKParticle_rot(double x, double y, double z, double ep1, double ep2, shared_ptr<Material> mat,double q0, double q1, double q2, double q3){//with a specified quatanion
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<GJKParticle>(new GJKParticle(x,y,z,ep1,ep2));
	GJKParticle* A = static_cast<GJKParticle*>(body->shape.get());

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
        fprintf(fin,"light_source{<4,8,5>*10 color rgb <0.77,0.75,0.75>}\n");
        fprintf(fin,"light_source{<12,-6>*10 color rgb <0.43,0.45,0.45>}\n");

	//fprintf(fin,"volume\t%e\n",volume);
	//fprintf(fin,"centroid\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2]);

        //output GJKParticle
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

		if (b->shape->getClassName()=="GJKParticle"){
			const shared_ptr<GJKParticle>& A = SUDODEM_PTR_CAST<GJKParticle> (b->shape);
			//particlesVolume += A->getVolume();
			//GJKParticle* A = static_cast<GJKParticle*>(b->shape.get());
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

    //output GJKParticle
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

		if (b->shape->getClassName()=="GJKParticle"){
            point_counter ++;
			const shared_ptr<GJKParticle>& A = SUDODEM_PTR_CAST<GJKParticle> (b->shape);
			//particlesVolume += A->getVolume();
			//GJKParticle* A = static_cast<GJKParticle*>(b->shape.get());
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
//output geometrical info of an assembly
void outputParticles(const char* filename){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write


	//fprintf(fin,"volume\t%e\n",volume);
	fprintf(fin,"#rx\try\trz\teps1\teps2\tx\ty\tz\tq0\tq1\tq2\tq3\n");

        //output GJKParticle
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

		if (b->shape->getClassName()=="GJKParticle"){
			const shared_ptr<GJKParticle>& A = SUDODEM_PTR_CAST<GJKParticle> (b->shape);
			//particlesVolume += A->getVolume();
			//GJKParticle* A = static_cast<GJKParticle*>(b->shape.get());
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
*/
//**********************************************************************************//


BOOST_PYTHON_MODULE(_gjkparticle_utils){
	// http://numpy.scipy.org/numpydoc/numpy-13.html mentions this must be done in module init, otherwise we will crash
	//import_array();

	SUDODEM_SET_DOCSTRING_OPTS;
	py::def("NewGJKParticle", NewGJKParticle, "Generate a GJKParticle.");
	py::def("NewGJKParticle2", NewGJKParticle2, "Generate a GJKParticle with margin.");
	py::def("GJKSphere", GJKSphere, "Generate a GJK Sphere with margin.");
	py::def("GJKPolyhedron", GJKPolyhedron, "Generate a GJK Polyhedron with margin.");
	py::def("GJKCone", GJKCone, "Generate a GJK Cone with margin.");
	py::def("GJKCylinder", GJKCylinder, "Generate a GJK Cylinder with margin.");
	py::def("GJKCuboid", GJKCuboid, "Generate a GJK Cuboid with margin.");

}
    //py::def("outputWalls", outputWalls, "output positions of walls.");
    //py::def("outputParticles", outputParticles, "output particles.");
//py::def("NewGJKParticle2", NewGJKParticle2, "Generate a Superquadric with isSphere option.");
	//py::def("NewGJKParticle_rot", NewGJKParticle_rot, "Generate a Superquadric with specified quaternion components: q0(w),q1(x),q2(y),q3(z).");
	//py::def("outputPOV",outputPOV,"output GJKParticle into a POV file for post-processing using POV-ray.");
    //py::def("outputVTK",outputVTK,"output GJKParticle into a VTK file for post-processing using Paraview.");
