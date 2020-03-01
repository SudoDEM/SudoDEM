#include<sudodem/core/Interaction.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/State.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/pkg/common/Aabb.hpp>

#include"GJKParticle.hpp"

#include<time.h>

#define _USE_MATH_DEFINES
//#define USE_MARGIN

SUDODEM_PLUGIN(/* self-contained in hpp: */ (GJKParticle) (GJKParticleGeom) (Bo1_GJKParticle_Aabb) (GJKParticlePhys) (GJKParticleMat) (Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys) (GJKParticleLaw)
	/* some code in cpp (this file): */
	#ifdef SUDODEM_OPENGL
		(Gl1_GJKParticle) (Gl1_GJKParticleGeom) /*(Gl1_GJKParticlePhys)*/
	#endif
	);


// Runtime

void DT_SetAccuracy(DT_Scalar max_error)
{
	if (max_error > SD_Scalar(0.0))
	{
		DT_Accuracy::setAccuracy(SD_Scalar(max_error));
	}
}

void DT_SetTolerance(DT_Scalar tol_error)
{
	if (tol_error > SD_Scalar(0.0))
	{
		DT_Accuracy::setTolerance(SD_Scalar(tol_error));
	}
}

SD_Scalar GJKParticle::supportH(const Vector3r& v) const {
	//testflag++;
	//if(!m_shape) cout<<"m_shape2"<<m_shape2<<" m_init="<<m_init<<"m_shapeTpye="<<m_shapeType<<"testflag="<<testflag<<endl;
	return m_shape->supportH(v);
}

Vector3r GJKParticle::support(const Vector3r &gn) const{
	Vector3r n = rot_mat2local*gn;
	//cout<<"gn"<<gn<<endl;
	//cout<<"mat"<<rot_mat2local<<endl;
	Vector3r sp;// = SUDODEM_PTR_CAST<DT_Convex>(m_shape)->support(n);//for polyhedron, we need rotate all vertices to the global
	switch(m_shapeType){
		case 1:
			sp = SUDODEM_PTR_CAST<DT_Polyhedron>(m_shape)->support(n);
			break;
		default:
			sp = SUDODEM_PTR_CAST<DT_Convex>(m_shape)->support(n);
	}
	//cout<<"sp"<<sp<<endl;
	return rot_mat2global*sp;
}

//**************************************************************
/*Constractor*/
//constructor
/*sssssssssssssss*/
///////////////////////////
void GJKParticle::ConvexHull(std::vector<Vector3r> verts )
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
	m_vertices2.clear();
    for (i = 0; i < count; i++)
	{
	array.push_back(verts[i][0]);
	array.push_back(verts[i][1]);
	array.push_back(verts[i][2]);
	//std::cout<<"x=="<<verts[i][0]<<std::endl;
    m_vertices2.push_back(verts[i]);
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
		m_facetIndexBuf.push_back(facetIndices);
    }

    qh NOerrexit = True;
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
}
void GJKParticle::Initial()
{
	if (m_init) return;

	switch(m_shapeType){
	case 0://sphere
		{

		#if defined(USE_MARGIN)
      	m_shape = shared_ptr<DT_Point>(new DT_Point(Vector3r(0.0f, 0.0f, 0.0f)));
		#else
      	m_shape = shared_ptr<DT_Sphere>(new DT_Sphere(m_radius));
		#endif
		//m_radius(radius)
		m_volume = 2./3.*Mathr::TWO_PI*pow(m_radius,3);
		m_centroid = Vector3r(0.,0.,0.);
		double inert = 0.4*m_volume*m_radius*m_radius;
		m_inertia = Vector3r(inert,inert,inert);
		m_orientation = Quaternionr::Identity();
		//Sphere_facets(10);
		break;
		}
	case 1://polyhedron
		{

		//create the convex hull
		ConvexHull(m_vertices);
		/*for(int i=0;i<m_vertices2.size();i++){
			std::cout<<m_vertices2[i]<<std::endl;
		}*/
		Polyhedron_volume_centroid();
		Polyhedron_inertia();
		//cout<<"m_shape"<<m_shape<<"use count="<<m_shape.use_count()<<"testflag"<<testflag<<endl;
		m_shape = shared_ptr<DT_Polyhedron>(new DT_Polyhedron(m_vertices2));

		//shared_ptr<DT_Polyhedron> shape(new DT_Polyhedron(m_vertices2));
		//cout<<"shape"<<shape<<"use count="<<shape.use_count()<<endl;
		//m_shape = SUDODEM_PTR_CAST<DT_Shape>(shape);
		//testflag = 20;
		//cout<<"m_shape"<<m_shape<<"use count="<<m_shape.use_count()<<endl;
		//if(!m_shape){cout<<"m_shape is a null pointer!"<<endl;}
		//m_facetIndexBuf = reinterpret_cast<DT_Polyhedron *>(m_shape)->m_facetIndexBuf;
		//m_vertices2 = reinterpret_cast<DT_Polyhedron *>(m_shape)->m_vertices2;
		//Polyhedron_facets(10);
        //data for opengl plotting
        m_facetVertices = m_vertices2;
        //std::cout<<"m_vertices2="<<m_vertices2<<std::endl;
        for(unsigned int i=0;i<m_facetIndexBuf.size();i++){
            T_IndexBuf Vindices = m_facetIndexBuf[i];
            //std::cerr<<"facet vertices= "<<Vindices[0]<<" "<<Vindices[1]<<" "<<Vindices[2]<<"\n";
            for(int tri=0; tri < (int) Vindices.size()-2; tri++){
                //triangular facet
                facetTri.push_back(Vindices[0]);
                //
                facetTri.push_back(Vindices[tri+2]);
                facetTri.push_back(Vindices[tri+1]);
            }

        }
		break;
		}
	case 2://cone
		{
		m_shape = shared_ptr<DT_Cone>(new DT_Cone(m_radius, m_height));
		m_displayList = 0;
		m_volume = Mathr::TWO_PI*m_radius*m_radius*m_height/6.0;
		m_centroid = Vector3r::Zero();//FIXME:not the real mass center, but the half height. The real mass center is below the top of 0.75*h.
		double inert_z = 0.3*m_volume*m_radius*m_radius;
		//double inert_x = 0.6*m_volume*m_height*m_height+0.5*inert_z;//about the apex
		//double inert_x = 3.0/80.0*m_volume*m_height*m_height+0.5*inert_z;//about the mass centroid
		double inert_x = 0.1*m_volume*m_height*m_height + 0.5*inert_z;//about the half height
		m_inertia = Vector3r(inert_x,inert_x,inert_z);
		m_orientation = Quaternionr::Identity();
		//
		//Cone_facets(10);
		break;
		}
	case 3://cylinder
		{
      	m_shape = shared_ptr<DT_Cylinder>(new DT_Cylinder(m_radius, m_height));
		m_displayList = 0;
		m_volume = Mathr::TWO_PI*m_radius*m_radius*m_height*0.5;
		m_centroid =  Vector3r::Zero();
		double inert_z = 0.5*m_volume*m_radius*m_radius;
		double inert_x = 0.5*inert_z+m_volume*m_height*m_height/12.0;
		m_inertia = Vector3r(inert_x,inert_x,inert_z);
		m_orientation = Quaternionr::Identity();
		//Cylinder_facets(10);
		break;
		}
		case 4://cube
			{
			Vector3r extend = m_vertices[0];//the first item is used for the extend of a box
			double x(extend[0]),y(extend[1]),z(extend[2]);
			m_shape = shared_ptr<DT_Box>(new DT_Box(extend));
			m_displayList = 0;
			m_volume = 8.0*x*y*z;
			m_centroid =  Vector3r::Zero();
			double inert_x = m_volume*(y*y + z*z)/3.0;
			double inert_y = m_volume*(x*x + z*z)/3.0;
			double inert_z = m_volume*(y*y + x*x)/3.0;
			m_inertia = Vector3r(inert_x,inert_y,inert_z);
			m_orientation = Quaternionr::Identity();
			//
			m_facetVertices.clear();
			//bottom
			m_facetVertices.push_back(Vector3r(-x,-y,-z));
			m_facetVertices.push_back(Vector3r( x,-y,-z));
			m_facetVertices.push_back(Vector3r( x, y,-z));
			m_facetVertices.push_back(Vector3r(-x, y,-z));
			//top
			m_facetVertices.push_back(Vector3r(-x,-y, z));
			m_facetVertices.push_back(Vector3r( x,-y, z));
			m_facetVertices.push_back(Vector3r( x, y, z));
			m_facetVertices.push_back(Vector3r(-x, y, z));
			break;
			}
	}
    //GetSurfaceTriangulation(10);
	//m_margin = 0.0f;
	m_xform.setIdentity();
	//setBBox();
	//std::cout<<"creating an object!"<<std::endl;
	//
	//calc_volume_centroid();
	//calc_inertia();
    //initialization done
		rot_mat2local = (m_orientation).conjugate().toRotationMatrix();//to particle's system
		rot_mat2global = (m_orientation).toRotationMatrix();//to global system
	m_shape2 = m_shape;
	m_init = true;
}


/*
void GJKParticle::randomGeometry(int num_vertices){//generate vertices of a polyhedral particle based on a unit sphere
	SD_Scalar theta,phi;
	double x,y,z;
	for(int i=0;i<num_vertices;i++){
		theta = Mathr::TWO_PI*MT_random();
		phi = Mathr::TWO_PI*0.5*MT_random();
		x = m_radius*MT_sin(phi)*MT_cos(theta);
		y = m_radius*MT_sin(phi)*MT_sin(theta);
		z = m_radius*MT_cos(phi);
		m_vertices.push_back(Vector3r(x,y,z));
	}


}
*/


void GJKParticle::Polyhedron_volume_centroid(){
	m_centroid = Vector3r(0.,0.,0.);
	m_volume = 0.0;

	Vector3r basepoint = m_vertices2[0];//pick one vertex of the convex hull as a base point
	Vector3r A,B,C;

	double vtetra;

	//compute centroid and volume
	//std::cout<<"facet num="<<m_facetIndexBuf.size()<<std::endl;
	for(unsigned int i=0;i<m_facetIndexBuf.size();i++){
		T_IndexBuf Vindices = m_facetIndexBuf[i];
		A = m_vertices2[Vindices[0]];
		B = m_vertices2[Vindices[1]];
		//std::cout<<"vertex num="<<Vindices.size()<<std::endl;
		for(unsigned int j = 2; j<Vindices.size();j++){
			C = m_vertices2[Vindices[j]];
			//std::cout<<"ABC"<<A-C<<B-C<<basepoint-C<<std::endl;
			vtetra = std::fabs((basepoint-C).dot((A-C).cross(B-C)))/6.;
			m_volume += vtetra;
			//std::cout<<"vtetra"<<vtetra<<std::endl;
			m_centroid += (basepoint+A+B+C) / 4. * vtetra;
			B = C;
		}
	}
	m_centroid /= m_volume;
	//convert vertices of particle into its local coordinate system.
	for(int i=0; i< (int) m_vertices2.size();i++) {
		m_vertices2[i] -=  m_centroid;
		m_vertices[i] -=  m_centroid;//FIXME:here no vertices are swept out when constructing the convex hull, i.e., m_vertices = m_vertices2.
	}
	m_position = m_centroid;
	//std::cout<<"centroid "<<m_centroid<<" volume "<<m_volume<<std::endl;
}
/*! Calculates tetrahedron inertia relative to the origin (0,0,0), with unit density (scales linearly).

See article F. Tonon, "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates", http://www.scipub.org/fulltext/jms2/jms2118-11.pdf
*/

//Centroid MUST be [0,0,0]
Matrix3r GJKParticle::TetraInertiaTensor(Vector3r av,Vector3r bv,Vector3r cv,Vector3r dv){
	#define x1 av[0]
	#define y1 av[1]
	#define z1 av[2]
	#define x2 bv[0]
	#define y2 bv[1]
	#define z2 bv[2]
	#define x3 cv[0]
	#define y3 cv[1]
	#define z3 cv[2]
	#define x4 dv[0]
	#define y4 dv[1]
	#define z4 dv[2]

	// Jacobian of transformation to the reference 4hedron
	double detJ=(x2-x1)*(y3-y1)*(z4-z1)+(x3-x1)*(y4-y1)*(z2-z1)+(x4-x1)*(y2-y1)*(z3-z1)
		-(x2-x1)*(y4-y1)*(z3-z1)-(x3-x1)*(y2-y1)*(z4-z1)-(x4-x1)*(y3-y1)*(z2-z1);
	detJ=fabs(detJ);
	double a=detJ*(y1*y1+y1*y2+y2*y2+y1*y3+y2*y3+
		y3*y3+y1*y4+y2*y4+y3*y4+y4*y4+z1*z1+z1*z2+
		z2*z2+z1*z3+z2*z3+z3*z3+z1*z4+z2*z4+z3*z4+z4*z4)/60.;
	double b=detJ*(x1*x1+x1*x2+x2*x2+x1*x3+x2*x3+x3*x3+
		x1*x4+x2*x4+x3*x4+x4*x4+z1*z1+z1*z2+z2*z2+z1*z3+
		z2*z3+z3*z3+z1*z4+z2*z4+z3*z4+z4*z4)/60.;
	double c=detJ*(x1*x1+x1*x2+x2*x2+x1*x3+x2*x3+x3*x3+x1*x4+
		x2*x4+x3*x4+x4*x4+y1*y1+y1*y2+y2*y2+y1*y3+
		y2*y3+y3*y3+y1*y4+y2*y4+y3*y4+y4*y4)/60.;
	// a' in the article etc.
	double a__=detJ*(2*y1*z1+y2*z1+y3*z1+y4*z1+y1*z2+
		2*y2*z2+y3*z2+y4*z2+y1*z3+y2*z3+2*y3*z3+
		y4*z3+y1*z4+y2*z4+y3*z4+2*y4*z4)/120.;
	double b__=detJ*(2*x1*z1+x2*z1+x3*z1+x4*z1+x1*z2+
		2*x2*z2+x3*z2+x4*z2+x1*z3+x2*z3+2*x3*z3+
		x4*z3+x1*z4+x2*z4+x3*z4+2*x4*z4)/120.;
	double c__=detJ*(2*x1*y1+x2*y1+x3*y1+x4*y1+x1*y2+
		2*x2*y2+x3*y2+x4*y2+x1*y3+x2*y3+2*x3*y3+
		x4*y3+x1*y4+x2*y4+x3*y4+2*x4*y4)/120.;

	Matrix3r ret; ret<<
		a   , -c__, -b__,
		-c__, b   , -a__,
		-b__, -a__, c   ;
	return ret;

	#undef x1
	#undef y1
	#undef z1
	#undef x2
	#undef y2
	#undef z2
	#undef x3
	#undef y3
	#undef z3
	#undef x4
	#undef y4
	#undef z4
}
void GJKParticle::Polyhedron_inertia(){
	//compute inertia
	Vector3r inertia;
	Quaternionr orientation;
	orientation = Quaternionr::Identity();

	Vector3r origin = Vector3r::Zero();//m_position;
	double vtetra;//volume of tetrahedron
	Vector3r ctetra;//centroid of tetrahedron
	Matrix3r Itet1, Itet2;
	Matrix3r inertia_tensor(Matrix3r::Zero());
	Vector3r A,B,C;
	for(unsigned int i=0;i<m_facetIndexBuf.size();i++){//for(int i=0; i<(int) faceTri.size(); i+=3){
		T_IndexBuf Vindices = m_facetIndexBuf[i];
		A = m_vertices2[Vindices[0]];
		B = m_vertices2[Vindices[1]];
		//std::cout<<"vertex num="<<Vindices.size()<<std::endl;
		for(unsigned int j = 2; j<Vindices.size();j++){
			C = m_vertices2[Vindices[j]];
			vtetra = std::fabs((origin-C).dot((A-C).cross(B-C)))/6.;

			ctetra = (origin+A+B+C) / 4.;
			Itet1 = TetraInertiaTensor(origin-ctetra, A-ctetra, B-ctetra, C-ctetra);
			ctetra = ctetra-origin;
			Itet2<<
				ctetra[1]*ctetra[1]+ctetra[2]*ctetra[2], -ctetra[0]*ctetra[1], -ctetra[0]*ctetra[2],
				-ctetra[0]*ctetra[1], ctetra[0]*ctetra[0]+ctetra[2]*ctetra[2], -ctetra[2]*ctetra[1],
				-ctetra[0]*ctetra[2], -ctetra[2]*ctetra[1], ctetra[1]*ctetra[1]+ctetra[0]*ctetra[0];
			inertia_tensor = inertia_tensor + Itet1 + Itet2*vtetra;
			B = C;
		}
	}
	//
	//std::cout<<"inertia tensor="<<inertia_tensor<<std::endl;
	//std::cout<<"iner="<<fabs(inertia_tensor(0,1))<<" "<<fabs(inertia_tensor(0,2))<<" "<<fabs(inertia_tensor(1,2))<<std::endl;
	double inertia_tensor_trace = inertia_tensor(0,0) + inertia_tensor(1,1) + inertia_tensor(2,2);
	if(fabs(inertia_tensor(0,1))+fabs(inertia_tensor(0,2))+fabs(inertia_tensor(1,2)) < 1E-5*inertia_tensor_trace){//FIXME: put a reasonable value of the relative tolerance if possible
		// no need to rotate, inertia already diagonal
		orientation = Quaternionr::Identity();
		inertia = Vector3r(inertia_tensor(0,0),inertia_tensor(1,1),inertia_tensor(2,2));
		//std::cout<<"inertia1="<<inertia.sum()<<std::endl;
	}else{
		// calculate eigenvectors of I
		Vector3r rot;
		Matrix3r I_rot(Matrix3r::Zero()), I_new(Matrix3r::Zero());
		matrixEigenDecomposition(inertia_tensor,I_rot,I_new);
		//std::cout<<"I_rot1>>"<<I_rot<<std::endl;
		//std::cout<<"I_new>>"<<I_new<<std::endl;
		// I_rot = eigenvectors of inertia_tensors in columns
		// I_new = eigenvalues on diagonal
		// set positove direction of vectors - otherwise reloading does not work
		Matrix3r sign(Matrix3r::Zero());
		double max_v_signed = I_rot(0,0);
		double max_v = fabs(I_rot(0,0));
		if (max_v < fabs(I_rot(1,0))) {max_v_signed = I_rot(1,0); max_v = fabs(I_rot(1,0));}
		if (max_v < fabs(I_rot(2,0))) {max_v_signed = I_rot(2,0); max_v = fabs(I_rot(2,0));}
		sign(0,0) = max_v_signed/max_v;
		//std::cout<<"max_v>>"<<max_v<<std::endl;
		max_v_signed = I_rot(0,1);
		max_v = fabs(I_rot(0,1));
		if (max_v < fabs(I_rot(1,1))) {max_v_signed = I_rot(1,1); max_v = fabs(I_rot(1,1));}
		if (max_v < fabs(I_rot(2,1))) {max_v_signed = I_rot(2,1); max_v = fabs(I_rot(2,1));}
		sign(1,1) = max_v_signed/max_v;
		sign(2,2) = 1.;
		//std::cout<<"max_v>>"<<max_v<<std::endl;
		I_rot = I_rot*sign;
		//std::cout<<"sign>>"<<sign<<std::endl;
		// force the eigenvectors to be right-hand oriented
		Vector3r third = (I_rot.col(0)).cross(I_rot.col(1));
		I_rot(0,2) = third[0];
		I_rot(1,2) = third[1];
		I_rot(2,2) = third[2];


		inertia = Vector3r(I_new(0,0),I_new(1,1),I_new(2,2));
		//std::cout<<"I_rot>>"<<I_rot<<std::endl;

		orientation = Quaternionr(I_rot);
		//std::cout<<"inertia2="<<inertia.sum()<<std::endl;
		//std::cout<<"orientation>>"<<orientation<<std::endl;
		//rotate the voronoi cell so that x - is maximal inertia axis and z - is minimal inertia axis
		orientation.normalize();  //not needed
		//std::cout<<"orientation>>"<<orientation<<std::endl;
		Matrix3r rot_mat = (orientation.conjugate()).toRotationMatrix();
		//std::cout<<"rot_mat>>"<<rot_mat<<std::endl;
		for(int i=0; i< (int) m_vertices2.size();i++) {
			m_vertices2[i] =  rot_mat*m_vertices2[i];
			//m_vertices2[i] =  orientation.conjugate()*m_vertices2[i];
			//m_vertices[i] =  orientation.conjugate()*m_vertices[i];//FIXME:it is assumed that m_vertices is equal to m_vertices2.
		}
		//std::cout<<"after rotated to its local coordinate systems"<<std::endl;
		//for(int i=0;i<m_vertices2.size();i++){
		//	std::cout<<m_vertices2[i]<<std::endl;
		//}
		//rotate also the CGAL structure Polyhedron
		//Matrix3d rot_mat = (orientation.conjugate()).toRotationMatrix();
		//Transformation t_rot(rot_mat(0,0),rot_mat(0,1),rot_mat(0,2),rot_mat(1,0),rot_mat(1,1),rot_mat(1,2),rot_mat(2,0),rot_mat(2,1),rot_mat(2,2),1.);
		//std::transform( P.points_begin(), P.points_end(), P.points_begin(), t_rot);

	}

	m_inertia = inertia;
	//m_orientation = orientation;
	setOrientation(orientation);
}


//****************************************************************************************
/* Destructor */
GJKParticle::~GJKParticle()
{
//	delete m_shape;
}
GJKParticleGeom::~GJKParticleGeom(){}


//*************************************************************************************
//****************************************************************************************
/* AaBb overlap checker  */

void Bo1_GJKParticle_Aabb::go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*){
  //if(!ig){cout<<"the pointer of ig is null at aabb_go."<<endl;}
	GJKParticle* t=static_cast<GJKParticle*>(ig.get());
	//if(!t){cout<<"the pointer of GJKParticle is null at aabb_go."<<endl;}
	if (!t->IsInitialized()){ std::cout<<"not initialized yet!"<<std::endl; t->Initial();
		t->rot_mat2local = (se3.orientation).conjugate().toRotationMatrix();//to particle's system
		t->rot_mat2global = (se3.orientation).toRotationMatrix(); //to global system
	}
	//cout<<"position aabb1--"<<se3.position<<endl;
	//t->setSe3(se3);
	//cout<<"orientation aabb2--"<<se3.orientation<<endl;
	//t->setOrientation(se3.orientation);


	/*if (!t->IsInitialized()){
	        t->Initial();
	        t->rot_mat2local = (se3.orientation).conjugate().toRotationMatrix();//to particle's system
	        t->rot_mat2global = (se3.orientation).toRotationMatrix(); //to global system
	}*/
	//double r_max = t->getr_max();
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());

        //Vector3r mincoords(-p_x(0),-p_y(1),-p_z(2)), maxcoords(p_x(0),p_y(1),p_z(2));
	//Vector3r mincoords(, maxcoords(r_max,r_max,r_max);
	//the aabb should be optimized in future.
	//return SUDODEM_PTR_CAST<DT_Convex>(m_shape)->supportH(v);
	Matrix3r rot = se3.orientation.toRotationMatrix();
	//cout<<"rot="<<rot<<endl;
	Vector3r min(- t->supportH(-rot.row(0)) ,
               - t->supportH(-rot.row(1)),
		       - t->supportH(-rot.row(2)) );
	Vector3r max( t->supportH( rot.row(0)) ,
               + t->supportH( rot.row(1)),
               + t->supportH( rot.row(2)) );

	//cout<<"min max1 "<<min<<max<<endl;
	min += se3.position;
	max += se3.position;
	//min += se3.position-t->getCentroid();
	//max += se3.position-t->getCentroid();
	//cout<<"min max"<<min<<max<<endl;
	//Vector3r max = t->getBBox().getMax();
	double margin = t->getMargin();
	aabb->min=Vector3r(min[0]-margin,min[1]-margin,min[2]-margin);//se3.position+mincoords;
	aabb->max=Vector3r(max[0]+margin,max[1]+margin,max[2]+margin);//se3.position+maxcoords;
}

//**********************************************************************************
/* Plotting */

#ifdef SUDODEM_OPENGL
	#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
	bool Gl1_GJKParticle::wire;
	int Gl1_GJKParticle::slices;

    void Gl1_GJKParticle::initSolidGlList(GJKParticle* t) {

	    //Generate the list. Only once for each qtView, or more if quality is modified.
	    glDeleteLists(t->m_glSolidListNum,1);
	    t->m_glSolidListNum = glGenLists(1);
	    glNewList(t->m_glSolidListNum,GL_COMPILE);
	    glEnable(GL_LIGHTING);
	    glShadeModel(GL_SMOOTH);//glShadeModel(GL_FLAT);
        std::vector<Vector3r> facetVertices;
	    // render the sphere now
        switch(t->getShapeType()){
	    case 0://sphere
		    {


	            //to do...
	            int w=2*slices,h=slices;

                ////
                int i,j;
                double a=0.0,b=0.0,c=0.0;
                double hStep=M_PI/(h-1);
                double wStep=2*M_PI/w;
                double x,y,z;

                //surface discretization
                for(a=0.0,i=0;i<h;i++,a+=hStep)
	            {
                  for(b=0.0,j=0;j<w;j++,b+=wStep)
                  {	c = a-M_PI_2l;

		            x = t->m_radius*cos(c)*cos(b);
		            y = t->m_radius*cos(c)*sin(b);
                    z = t->m_radius*sin(c);

	                facetVertices.push_back(Vector3r(x,y,z));
                  }
                }
	            //FIXME:the two polar points should be just single points.
                //polygons of surface slices
                Vector3r vert;
                glBegin(GL_TRIANGLE_STRIP);
	            for(i=0;i<h-1;i++)
	            {   //Use TRIANGLE_STRIP for faster display of adjacent facets
		            //glBegin(GL_TRIANGLE_STRIP);
	                for(j=0;j<w-1;j++)
	                {
                        vert = facetVertices[(i+1)*w+j];vert.normalize();
                        glNormal3v(vert);glVertex3v(facetVertices[(i+1)*w+j]);
                        vert = facetVertices[i*w+j];vert.normalize();
                        glNormal3v(vert);glVertex3v(facetVertices[i*w+j]);
		            }
                        vert = facetVertices[(i+1)*w];vert.normalize();
                        glNormal3v(vert);glVertex3v(facetVertices[(i+1)*w]);
                        vert = facetVertices[i*w];vert.normalize();
                        glNormal3v(vert);glVertex3v(facetVertices[i*w]);

                    //glEnd();
	            }
                glEnd();
                //std::cout<<"call the makelist function!"<<std::endl;


		    break;
		    }
	    case 1://polyhedron
		    {	//Polyhedron_facets(slices);
	            //facetVertices = t->m_vertices2;
                vector<int> facetTri = t->facetTri;

                facetVertices = t->m_facetVertices;
								//
								double margin = t->getMargin();
								//cout<<"margin: "<<margin<<endl;
								/*if(margin > 0.0){//add the margin to the outside of the particle
									for(std::vector<Vector3r>::iterator it = facetVertices.begin();it!=facetVertices.end();it++){
										Vector3r tmp = *it;
										*it = tmp.normalized()*margin+tmp;
										//cout<<"origin vertice: "<<tmp<<endl;
									}
									//for(std::vector<Vector3r>::iterator it = facetVertices.begin();it!=facetVertices.end();it++){
									//	cout<<"updated vertice: "<<*it<<endl;
									//}
								}*/
                //glDisable(GL_CULL_FACE);
                Vector3r faceCenter, n;
                glBegin(GL_TRIANGLES);
                #define __ONEFACE(a,b,c) n=(facetVertices[b]-facetVertices[a]).cross(facetVertices[c]-facetVertices[a]); n.normalize(); faceCenter=(facetVertices[a]+facetVertices[b]+facetVertices[c])/3.; if((faceCenter).dot(n)<0)n=-n; glNormal3v(n); glVertex3v(facetVertices[a]); glVertex3v(facetVertices[b]); glVertex3v(facetVertices[c]);
				    for(int tri=0; tri < (int) facetTri.size(); tri+=3) {__ONEFACE(facetTri[tri],facetTri[tri+1],facetTri[tri+2]);}
			    #undef __ONEFACE

                glEnd();

		    break;
		    }
	    case 2://cone
		    {

	            Vector3r top_vert = Vector3r(0.0,0.0,0.5*t->m_height);
                Vector3r bottom_center = Vector3r(0.0,0.0,-0.5*t->m_height);
                Vector3r vert;
	            double x,y,z;
	            z = -0.5*t->m_height;
	            float theta = Mathr::TWO_PI/float(slices);
	            //save vertices on the base
	            for(int ii = 0; ii<slices; ii++)
	            {	x = t->m_radius*cos(ii*theta);
		            y = t->m_radius*sin(ii*theta);
		            facetVertices.push_back(Vector3r(x,y,z));
	            }
	            //save the top vertex
	            //facetVertices.push_back(top_vert);
	            //adjancy of facets
	            //side facets
                glBegin(GL_TRIANGLE_FAN);
                glVertex3v(top_vert);   //the top vertice
	            for(int i=0;i<slices;i++){//number of side facets = num_segments
                    vert = facetVertices[i];vert.normalize();
		            glNormal3v(vert);glVertex3v(facetVertices[i]);
	            }
                vert = facetVertices[0];vert.normalize();
                glNormal3v(vert);glVertex3v(facetVertices[0]);
                glEnd();
                //the bottom face
                glBegin(GL_POLYGON);
                //glVertex3v(bottom_center);   //the position of bottom center
	            for(int i=slices-1;i>=0;i--){//number of side facets = num_segments
                    vert = facetVertices[i];vert.normalize();
		            glNormal3v(vert);glVertex3v(facetVertices[i]);
	            }
                //vert = facetVertices[slices-1];vert.normalize();
                //glNormal3v(vert);glVertex3v(facetVertices[slices-1]);
                glEnd();

			    //Cone_facets(slices);
		    break;
		    }
	    case 3://cylinder
		    {
	            double x,y,z1,z2;
	            z1 = 0.5*t->m_height;
                z2 = -0.5*t->m_height;
	            float theta = Mathr::TWO_PI/float(slices);
	            //save vertices on the top and bottom facets
	            for(int ii = 0; ii<slices; ii++)
	            {	x = t->m_radius*cos(ii*theta);
		            y = t->m_radius*sin(ii*theta);
		            facetVertices.push_back(Vector3r(x,y,z1));//top
                    facetVertices.push_back(Vector3r(x,y,z2));//bottom
	            }
	            Vector3r vert;
                glBegin(GL_TRIANGLE_STRIP);
	            //adjancy of facets
	            //side facets
	            for(int i=0;i<facetVertices.size();i++){//
                    vert = facetVertices[i];vert.normalize();
		            glNormal3v(vert);glVertex3v(facetVertices[i]);
	            }
                vert = facetVertices[0];vert.normalize();
                glNormal3v(vert);glVertex3v(facetVertices[0]);
                vert = facetVertices[1];vert.normalize();
                glNormal3v(vert);glVertex3v(facetVertices[1]);
                glEnd();
                //top and bottom faces
                glBegin(GL_POLYGON);
	            for(int i=0;i<slices;i++){//
                    //vert = facetVertices[i];vert.normalize();
		            /*glNormal3v(vert);*/glVertex3v(facetVertices[i*2]);
	            }
                //vert = facetVertices[slices-1];vert.normalize();
                //glNormal3v(vert);glVertex3v(facetVertices[slices-1]);
                glEnd();
                glBegin(GL_POLYGON);

	            for(int i=slices-1;i>=0;i--){//
                    //vert = facetVertices[i];vert.normalize();
		            /*glNormal3v(vert);*/glVertex3v(facetVertices[i*2+1]);
	            }
                glEnd();
		    break;
		    }
				case 4://cube
				 {
					facetVertices = t->m_facetVertices;
					glBegin(GL_TRIANGLE_STRIP);
				//adjancy of facets
				//side facets
				Vector3r vert1,vert2;
				for(int i=0;i<4;i++){//
					vert1 = facetVertices[i+4];vert1.normalize();
					glNormal3v(vert1);glVertex3v(facetVertices[i+4]);
					vert2 = facetVertices[i];vert2.normalize();
					glNormal3v(vert2);glVertex3v(facetVertices[i]);
				}
				vert1 = facetVertices[4];vert1.normalize();
				glNormal3v(vert1);glVertex3v(facetVertices[4]);
				vert1 = facetVertices[0];vert1.normalize();
				glNormal3v(vert1);glVertex3v(facetVertices[0]);
				glEnd();
				glBegin(GL_POLYGON);
				for(int i=0;i<4;i++){//
							//vert = facetVertices[i];vert.normalize();
					/*glNormal3v(vert);*/glVertex3v(facetVertices[i+4]);
				}
					//vert = facetVertices[slices-1];vert.normalize();
					//glNormal3v(vert);glVertex3v(facetVertices[slices-1]);
				glEnd();
				glBegin(GL_POLYGON);
				for(int i=3;i>=0;i--){//
							//vert = facetVertices[i];vert.normalize();
					/*glNormal3v(vert);*/glVertex3v(facetVertices[i]);
				}
					//vert = facetVertices[slices-1];vert.normalize();
					//glNormal3v(vert);glVertex3v(facetVertices[slices-1]);
				glEnd();
				 break;
				 }
	    }
	    glEndList();

    }

    void Gl1_GJKParticle::initWireGlList(GJKParticle* t){
	    //Generate the "no-stripes" display list, each time quality is modified
	    glDeleteLists(t->m_glWireListNum,1);
	    t->m_glWireListNum = glGenLists(1);
	    glNewList(t->m_glWireListNum,GL_COMPILE);
	    glDisable(GL_LIGHTING);
	    //	glShadeModel(GL_SMOOTH);

        std::vector<Vector3r> facetVertices;
        switch(t->getShapeType()){
	    case 0://sphere
		    {

            	        //to do...
	            int w=2*slices,h=slices;

                ////
                int i,j;
                double a=0.0,b=0.0,c=0.0;
                double hStep=M_PI/(h-1);
                double wStep=2*M_PI/w;
                double x,y,z;

                //surface discretization
                for(a=0.0,i=0;i<h;i++,a+=hStep)
	            {
                  for(b=0.0,j=0;j<w;j++,b+=wStep)
                  {	c = a-M_PI_2l;

		            x = t->m_radius*cos(c)*cos(b);
		            y = t->m_radius*cos(c)*sin(b);
                    z = t->m_radius*sin(c);

	                facetVertices.push_back(Vector3r(x,y,z));
                  }
                }
	            //FIXME:the two polar points should be just single points.
                //polygons of surface slices
                glBegin(GL_LINE_STRIP);
                //draw the latitudes first
	            for(i=1;i<h-1;i++)
	            {
		            //glBegin(GL_LINE_LOOP);
	                for(j=0;j<w;j++)
	                {	glVertex3v(facetVertices[i*w+j]);

		            }
                    glVertex3v(facetVertices[i*w]);
                    //glEnd();
	            }
                glVertex3v(facetVertices[(h-1)*w]);//the top pole

                //downwards the bottom pole at the opposite side
                for(i=h-2;i>=0;i--)
                {
                    glVertex3v(facetVertices[i*w+slices]);////FIXME: it is assumed that w = 2*slices; the opposite column is the "slices".
                }

                //draw the rest longitudes
                for(j=1;j<slices;j++)
                {
                    for(i=1;i<h;i++)
                    {
                        glVertex3v(facetVertices[i*w+j]);
                    }
                    for(i=h-1;i>=0;i--)
                    {
                        glVertex3v(facetVertices[i*w+j+slices]);
                    }

                }
                glEnd();
                /*for(j=0;i<w-1;j++)
	            {
		            glBegin(GL_LINE_STRIP);
	                for(i=0;i<h;i++)
	                {	glVertex3v(facetVertices[i*w+j]);

		            }
                    glEnd();
	            }*/
                //std::cout<<"call the makelist function!"<<std::endl;
                //glutSolidSphere(t->m_radius,slices,slices*2);
                //glutWireSphere(t->m_radius,slices,slices*2);
                //std::cout<<"call the makelist function!"<<std::endl;

		    break;
		    }
	    case 1://polyhedron
		    {
                vector<int> facetTri = t->facetTri;
                facetVertices = t->m_facetVertices;
                glBegin(GL_LINES);//FIXME:duplicate edges are redrawn.
                for(int tri=0; tri < (int) facetTri.size(); tri+=3){
                        //triangular facet
                        glVertex3v(facetVertices[ facetTri[tri]]);
                        glVertex3v(facetVertices[ facetTri[tri+1]]);
                        glVertex3v(facetVertices[ facetTri[tri+2]]);
                }
                glEnd();
		    break;
		    }
	    case 2://cone
		    {
                //glutWireCone(t->m_radius,t->m_height,slices, slices);//slices, stacks
			    //Cone_facets(slices);
	            Vector3r top_vert = Vector3r(0.0,0.0,0.5*t->m_height);
                Vector3r bottom_center = Vector3r(0.0,0.0,-0.5*t->m_height);
                Vector3r vert;
	            double x,y,z;
	            z = -0.5*t->m_height;
	            float theta = Mathr::TWO_PI/float(slices);
	            //save vertices on the base
	            for(int ii = 0; ii<slices; ii++)
	            {	x = t->m_radius*cos(ii*theta);
		            y = t->m_radius*sin(ii*theta);
		            facetVertices.push_back(Vector3r(x,y,z));
	            }
	            //save the top vertex
	            //facetVertices.push_back(top_vert);
	            //adjancy of facets
	            //side facets
                glBegin(GL_LINES);

	            for(int i=0;i<slices;i++){//number of side facets = num_segments
                    glVertex3v(top_vert);   //the top vertice
                    glVertex3v(facetVertices[i]);
	            }
                glEnd();
                //the bottom face
                glBegin(GL_LINE_LOOP);
                //glVertex3v(bottom_center);   //the position of bottom center
	            for(int i=slices-1;i>=0;i--){//number of side facets = num_segments
                    glVertex3v(facetVertices[i]);
	            }
                //vert = facetVertices[slices-1];vert.normalize();
                //glNormal3v(vert);glVertex3v(facetVertices[slices-1]);
                glEnd();
		    break;
		    }
	    case 3://cylinder
		    {
	            double x,y,z1,z2;
	            z1 = 0.5*t->m_height;
                z2 = -0.5*t->m_height;
	            float theta = Mathr::TWO_PI/float(slices);
	            //save vertices on the top and bottom facets
	            for(int ii = 0; ii<slices; ii++)
	            {	x = t->m_radius*cos(ii*theta);
		            y = t->m_radius*sin(ii*theta);
		            facetVertices.push_back(Vector3r(x,y,z1));//top
                    facetVertices.push_back(Vector3r(x,y,z2));//bottom
	            }
	            Vector3r vert;
                glBegin(GL_LINES);
	            //adjancy of facets
	            //side facets
	            for(int i=0;i<facetVertices.size();i++){//
		            glVertex3v(facetVertices[i]);
	            }
                glVertex3v(facetVertices[0]);
                glVertex3v(facetVertices[1]);
                glEnd();
                //top and bottom faces
                glBegin(GL_LINE_LOOP);
	            for(int i=0;i<slices;i++){//
                    //vert = facetVertices[i];vert.normalize();
		            /*glNormal3v(vert);*/glVertex3v(facetVertices[i*2]);
	            }
                //vert = facetVertices[slices-1];vert.normalize();
                //glNormal3v(vert);glVertex3v(facetVertices[slices-1]);
                glEnd();
                glBegin(GL_LINE_LOOP);

	            for(int i=slices-1;i>=0;i--){//
                    //vert = facetVertices[i];vert.normalize();
		            /*glNormal3v(vert);*/glVertex3v(facetVertices[i*2+1]);
	            }
                glEnd();
		    break;
		    }
				case 4://cube
				 {
					facetVertices = t->m_facetVertices;
					glBegin(GL_LINE_LOOP);
					for(int i=0;i<4;i++){
						glVertex3v(facetVertices[i]);
					}
					glVertex3v(facetVertices[0]);
					for(int i=4;i<8;i++){
						glVertex3v(facetVertices[i]);
					}
					glVertex3v(facetVertices[4]);
					glEnd();
					glBegin(GL_LINES);
					for(int i=1;i<4;i++){glVertex3v(facetVertices[i]);glVertex3v(facetVertices[i+4]);}
					glEnd();
				 break;
				 }
	    }
	    glEndList();
    }






	void Gl1_GJKParticle::go(const shared_ptr<Shape>& cm, const shared_ptr<State>&state,bool wire2,const GLViewInfo&)
	{
		glMaterialv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,Vector3r(cm->color[0],cm->color[1],cm->color[2]));
		glColor3v(cm->color);
		GJKParticle* t=static_cast<GJKParticle*>(cm.get());
		//
            if (!t->IsInitialized()){
	            t->Initial();
	        }
		//if (t->getGLslices()!=slices){// || !t->IsGLInit()) {
        // vector<int> facetTri;
       /* if(pre_slices!=slices)
        {   std::cout<<"haha preslices"<<pre_slices<<"slices"<<slices<<std::endl;
            pre_slices = slices;
            t->GetSurfaceTriangulation(slices);
        }*/
        bool somethingChanged = (t->m_GL_slices!=slices|| glIsList(t->m_glSolidListNum)!=GL_TRUE);
		if (somethingChanged) {
            //std::cout<<"haha preslices"<<t->m_GL_slices<<"slices"<<slices<<std::endl;
            initSolidGlList(t); initWireGlList(t);t->m_GL_slices=slices;}
        if(!wire)
        {
            glCallList(t->m_glSolidListNum);
        }else{
            //std::cout<<"glSolidListNum="<<t->m_glSolidListNum<<"   glWireListNum="<<t->m_glWireListNum<<std::endl;
            glCallList(t->m_glWireListNum);//glCallList(t->m_glListNum);
        }
	}

	void Gl1_GJKParticleGeom::go(const shared_ptr<IGeom>& ig, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool){
	        GJKParticleGeom* pg=static_cast<GJKParticleGeom*>(ig.get());
	        Vector3r point1,point2;
	        point1 = pg->point1;
	        point2 = pg->point2;
	        //testing

	        Vector3r point11,point22;
	        point11 = pg->point11;
	        point22 = pg->point22;

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

void GJKParticleGeom::precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r&
currentNormal, bool isNew, const Vector3r& shift2){

	if(!isNew) {
		orthonormal_axis = normal.cross(currentNormal);
		Real angle = scene->dt*0.5*currentNormal.dot(rbp1.angVel + rbp2.angVel);
		twist_axis = angle*currentNormal;}
	else twist_axis=orthonormal_axis=Vector3r::Zero();
	//Update contact normal
	normal=currentNormal;
	//Precompute shear increment
	Vector3r c1x = (contactPoint - rbp1.pos);
	Vector3r c2x = (contactPoint - rbp2.pos + shift2);
	Vector3r relativeVelocity = (rbp2.vel+rbp2.angVel.cross(c2x)) - (rbp1.vel+rbp1.angVel.cross(c1x));
	//keep the shear part only
	relativeVn = normal.dot(relativeVelocity)*normal;
	relativeVs = relativeVelocity-relativeVn;
	shearInc = relativeVs*scene->dt;
}

//**********************************************************************************
Vector3r& GJKParticleGeom::rotate(Vector3r& shearForce) const {
	// approximated rotations

	//std::cerr << "rto:"<<gettid4()<<orthonormal_axis[0]<<orthonormal_axis[1]<<orthonormal_axis[2]<< "\n";
	//std::cerr << "rtt:"<<gettid4()<<twist_axis[0]<<twist_axis[1]<<twist_axis[2]<< "\n";
	//std::cerr << "rtn:"<<gettid4()<<normal[0]<<normal[1]<<normal[2]<< "\n";
	//char filename[20];
	//sprintf(filename,"%d",gettid4());
	//FILE * fin = fopen(filename,"a");
	if(shearForce.squaredNorm()<1E-18){return shearForce;}
	//fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
	//fprintf(fin,"shear_force1\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	//fprintf(fin,"orthonormal_axis\t%e\t%e\t%e\n",orthonormal_axis[0],orthonormal_axis[1],orthonormal_axis[2]);
	//fprintf(fin,"twist_axis\t%e\t%e\t%e\n",twist_axis[0],twist_axis[1],twist_axis[2]);
	//fprintf(fin,"normal\t%e\t%e\t%e\n",normal[0],normal[1],normal[2]);
	shearForce -= shearForce.cross(orthonormal_axis);
	//fprintf(fin,"shear_force2\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	shearForce -= shearForce.cross(twist_axis);
	//fprintf(fin,"shear_force3\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	//NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	shearForce -= normal.dot(shearForce)*normal;
	//fprintf(fin,"shear_force4\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
	//fclose(fin);
	return shearForce;
}


//**********************************************************************************
/* Material law, physics */

void Ip2_GJKParticleMat_GJKParticleMat_GJKParticlePhys::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return;
	const shared_ptr<GJKParticleMat>& mat1 = SUDODEM_PTR_CAST<GJKParticleMat>(b1);
	const shared_ptr<GJKParticleMat>& mat2 = SUDODEM_PTR_CAST<GJKParticleMat>(b2);
	interaction->phys = shared_ptr<GJKParticlePhys>(new GJKParticlePhys());
	const shared_ptr<GJKParticlePhys>& contactPhysics = SUDODEM_PTR_CAST<GJKParticlePhys>(interaction->phys);
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
Real GJKParticleLaw::getPlasticDissipation() {return (Real) plasticDissipation;}
void GJKParticleLaw::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}
Real GJKParticleLaw::elasticEnergy()
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
// Apply forces on GJKParticle in collision based on geometric configuration
bool GJKParticleLaw::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I){

		if (!I->geom) {return true;}
		const shared_ptr<GJKParticleGeom>& contactGeom(SUDODEM_PTR_DYN_CAST<GJKParticleGeom>(I->geom));
		if(!contactGeom) {return true;}
		const Body::id_t idA=I->getId1(), idB=I->getId2();
		const shared_ptr<Body>& A=Body::byId(idA), B=Body::byId(idB);

        State* de1 = Body::byId(idA,scene)->state.get();
	    State* de2 = Body::byId(idB,scene)->state.get();

		GJKParticlePhys* phys = dynamic_cast<GJKParticlePhys*>(I->phys.get());

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
		Vector3r normalForce=contactGeom->normal*contactGeom->PenetrationDepth*phys->kn;
		Vector3r shearForce1=Vector3r::Zero();//zhswee deprecated GJKParticleLaw.shearForce
		//shear force: in case the polyhdras are separated and come to contact again, one should not use the previous shear force
		//std::cerr << "p1:"<<gettid4()<<shearForce[0]<< "\n";
		///////////////////////
		bool useDamping=(phys->betan > 1e-5 || phys->betas > 1e-5);//using viscous damping?
		//FIXME:BUGS WHEN USING VISCOUS DAMPING. Particles are likely to bound and fly away suddently.
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
		if (contactGeom->isShearNew)
		{
			//shearForce = Vector3r::Zero();
			shearForce1 = Vector3r(0.,0.,0.);
		}else{
			//shearForce = contactGeom->rotate(shearForce);
			//shearForce1 = contactGeom->rotate(phys->shearForce);
			shearForce1 = phys->shearForce;
			shearForce1 = contactGeom->rotate(shearForce1);
		}
		if(isnan(shearForce1[0])){
			std::cerr << "shearF:"<<phys->shearForce(0)<<" "<<phys->shearForce(1)<<" "<<phys->shearForce(2)<< "\n";
            std::cerr<<"isShearNeew: "<<contactGeom->isShearNew<<"\n";
            std::cerr<<"shearForce1: "<<shearForce1(0)<<" "<<shearForce1(1)<<" "<<shearForce1(2)<<"\n";
			//std::cerr << "relativeVs:"<<contactGeom->relativeVs(0)<<" "<<contactGeom->relativeVs(1)<<" "<<contactGeom->relativeVs(2)<< "\n";
		}
		const Vector3r& shearDisp = contactGeom->shearInc;
		shearForce1 -= phys->ks*shearDisp;

        //std::cerr << "p3:"<<gettid4() <<shearDisp[0]<< "\n";
		Real maxFs = normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);

		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if( shearForce1.squaredNorm() > maxFs ){
			Real ratio = sqrt(maxFs) / shearForce1.norm();
			shearForce1 *= ratio;

		}
		else if (useDamping){ // add current contact damping if we do not slide and if damping is requested
		        shearViscous = cs*contactGeom->relativeVs; // get shear viscous component
		        //std::cerr << "sF1:"<<shearForce1(0)<<" "<<shearForce1(1)<<" "<<shearForce1(2)<< "\n";
		        //std::cerr << "sF1:"<<shearForce1.norm()<< "\n";
		        //shearForce1 -= shearViscous;
		        //std::cerr << "sF2:"<<shearForce1.norm()<< "\n";
		        if(isnan(shearForce1[0])){
					std::cerr << "shearViscous:"<<shearViscous(0)<<" "<<shearViscous(1)<<" "<<shearViscous(2)<< "\n";
					std::cerr << "relativeVs:"<<contactGeom->relativeVs(0)<<" "<<contactGeom->relativeVs(1)<<" "<<contactGeom->relativeVs(2)<< "\n";
				}
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

		if(isnan(shearForce1[0])||isnan(normalForce[0])){
		  //char filename[20];
          time_t tmNow = time(NULL);
          char tmp[64];
          strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&tmNow) );
		  const char* filename = "JGKParticleLaw_err.log";
		  //sprintf(filename,"%d",gettid4());
		  FILE * fin = fopen(filename,"a");
          fprintf(fin,"\n*******%s*******\n",tmp);
          fprintf(fin,"isShearNew\t%d\n",contactGeom->isShearNew);
          Vector3r normal = contactGeom->normal;
          Real depth = contactGeom->PenetrationDepth;
		  //fprintf(fin,"normal_force\t%e\t%e\t%e\n",normalForce[0],normalForce[1],normalForce[2]);
		  //fprintf(fin,"shear_force\t%e\t%e\t%e\n",shearForce[0],shearForce[1],shearForce[2]);
		  fprintf(fin,"shear_force1\t%e\t%e\t%e\n",shearForce1[0],shearForce1[1],shearForce1[2]);
          fprintf(fin,"shearViscous\t%e\t%e\t%e\n",shearViscous[0],shearViscous[1],shearViscous[2]);
          fprintf(fin,"relativeVs\t%e\t%e\t%e\n",contactGeom->relativeVs[0],contactGeom->relativeVs[1],contactGeom->relativeVs[2]);
		  //fprintf(fin,"F\t%e\t%e\t%e\n",F[0],F[1],F[2]);
		  fprintf(fin,"ks\t%e maxFs\t%e sfn\t%e\n",phys->ks,maxFs,phys->shearForce.squaredNorm());
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
          contactGeom->isShearNew = true;
		}
		Vector3r F = -normalForce-shearForce1 + shearViscous;	//the normal force acting on particle A (or 1) is equal to the normal force.
		if (contactGeom->PenetrationDepth != contactGeom->PenetrationDepth) exit(1);
		scene->forces.addForce (idA,F);
		scene->forces.addForce (idB, -F);
		scene->forces.addTorque(idA, -(A->state->pos-contactGeom->contactPoint).cross(F));
		scene->forces.addTorque(idB, (B->state->pos-contactGeom->contactPoint).cross(F));
		//needed to be able to acces interaction forces in other parts of sudodem
		//shearForce=shearForce1;
		phys->normalForce = normalForce;
		phys->shearForce = shearForce1;
		//may need to add phys->shearViscous and phys->normalViscous
		return true;
}
