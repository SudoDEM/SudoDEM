#include "GJKParticle_Ig2.hpp"
#include <cmath>
//#include <sys/syscall.h>
#define _USE_MATH_DEFINES


SUDODEM_PLUGIN(/* self-contained in hpp: */ (Ig2_GJKParticle_GJKParticle_GJKParticleGeom)
		(Ig2_Wall_GJKParticle_GJKParticleGeom)
		(Ig2_Facet_GJKParticle_GJKParticleGeom)
	   );

//**********************************************************************************
/*some auxiliary functions*/

//////////////////////////////////////////////
//#define LM_INFO_SZ    	 10
// #define LM_ERROR         -1
// #define LM_INIT_MU    	 1E-03
// #define LM_STOP_THRESH	 1E-17
//
// #define EPSILON       1E-12
// #define ONE_THIRD     0.3333333333333334 /* 1.0/3.0 */
// #define LM_REAL_MIN -DBL_MAX
// #define LM_CNST(x) (x)
// #define LM_REAL_MAX DBL_MAX
// #define LM_FINITE finite
// #define LM_ISINF(x) isinf(x)

//using namespace std;
//////////////////////////////////////////////
/*
bool penetration_depth(const GJKParticle *objA, const GJKParticle* objB, Vector3r& v, Vector3r& pa, Vector3r& pb)
{

return hybrid_penetration_depth(DT_Transform(objA->m_xform, static_cast<const DT_Convex&>(*(objA->m_shape))), objA->m_margin,
									DT_Transform(objB->m_xform, static_cast<const DT_Convex&>(*(objB->m_shape))), objB->m_margin, v, pa, pb);

}
bool penetration_depth(const GJKParticle *objA,const Matrix3r& rot_matA,  Vector3r& positionA,const GJKParticle* objB,const Matrix3r& rot_matB, Vector3r& positionB, Vector3r& v, Vector3r& pa, Vector3r& pb)
{
return hybrid_penetration_depth(DT_TransM2(rot_matA,positionA, (objA->m_shape)), objA->m_margin,
									DT_TransM2(rot_matB, positionB,(objB->m_shape)), objB->m_margin, v, pa, pb);

}
bool intersect(const GJKParticle* a, const GJKParticle* b, Vector3r& v)
{
	DT_Transform ta(a->m_xform, (const DT_Convex&)*(a->m_shape));
	DT_Transform tb(b->m_xform, (const DT_Convex&)*(b->m_shape));
    return intersect((a->m_margin > SD_Scalar(0.0) ? static_cast<const DT_Convex&>(DT_Minkowski(ta, DT_Sphere(a->m_margin))) : static_cast<const DT_Convex&>(ta)),
			         (b->m_margin > SD_Scalar(0.0) ? static_cast<const DT_Convex&>(DT_Minkowski(tb, DT_Sphere(b->m_margin))) : static_cast<const DT_Convex&>(tb)), v);

}


bool common_point(const GJKParticle* a, const GJKParticle* b, Vector3r& v, Vector3r& pa, Vector3r& pb)
{

	DT_Transform ta(a->m_xform, (const DT_Convex&)*(a->m_shape));
	DT_Transform tb(b->m_xform, (const DT_Convex&)*(b->m_shape));
    return common_point((a->m_margin > SD_Scalar(0.0) ? static_cast<const DT_Convex&>(DT_Minkowski(ta, DT_Sphere(a->m_margin))) : static_cast<const DT_Convex&>(ta)),
			         (b->m_margin > SD_Scalar(0.0) ? static_cast<const DT_Convex&>(DT_Minkowski(tb, DT_Sphere(b->m_margin))) : static_cast<const DT_Convex&>(tb)), v, pa, pb);

}


SD_Scalar closest_points(const GJKParticle* a, const GJKParticle* b,
						 Vector3r& pa, Vector3r& pb)
{
	DT_Transform ta(a->m_xform, (const DT_Convex&)*(a->m_shape));
	DT_Transform tb(b->m_xform, (const DT_Convex&)*(b->m_shape));
    return closest_points((a->m_margin > SD_Scalar(0.0) ? static_cast<const DT_Convex&>(DT_Minkowski(ta, DT_Sphere(a->m_margin))) : static_cast<const DT_Convex&>(ta)),
						  (b->m_margin > SD_Scalar(0.0) ? static_cast<const DT_Convex&>(DT_Minkowski(tb, DT_Sphere(b->m_margin))) : static_cast<const DT_Convex&>(tb)), Mathr::MAX_REAL, pa, pb);

}
*/
//
bool penetration_depth(const DT_Convex& a, const DT_Convex& b,
                       Vector3r& v, Vector3r& pa, Vector3r& pb)
{
	DT_GJK gjk;

	#ifdef STATISTICS
    num_iterations = 0;
	#endif

	SD_Scalar dist2 = Mathr::MAX_REAL;

  do
	{
		Vector3r  p = a.support(-v);
		Vector3r  q = b.support(v);
		Vector3r w = p - q;

    if (v.dot(w) > SD_Scalar(0.0))//the origin is outside the Minkowski difference A-B, thus no contacting
		{
			return false;
		}
		//add the vertex to the convex hull
		gjk.addVertex(w, p, q);

    if (gjk.isAffinelyDependent())
    {
			#ifdef STATISTICS
            ++num_irregularities;
			#endif

      return false;
    }


		#ifdef STATISTICS
		        ++num_iterations;
		#endif
    if (!gjk.closest(v))
		{
			#ifdef STATISTICS
			            ++num_irregularities;
			#endif
      return false;
    }

		#ifdef SAFE_EXIT
				SD_Scalar prev_dist2 = dist2;
		#endif

		dist2 = v.squaredNorm();

		#ifdef SAFE_EXIT
			if (prev_dist2 - dist2 <= Mathr::EPSILON * prev_dist2)
			{
	 			return false;
			}
		#endif
    }
    while (!gjk.fullSimplex() && dist2 > DT_Accuracy::tol_error * gjk.maxVertex());

    SolvePenDepth SPD;
	return SPD.penDepth(gjk, a, b, v, pa, pb);

}
template <typename T>
bool getPenetration(T *particle1, T *particle2,	Vector3r position1, Vector3r position2,
										Vector3r& v, Vector3r& pa, Vector3r& pb)
{
	SD_Scalar margin = min(particle1->getMargin(),particle1->getMargin());
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
		//v.normalize();
		do
		{
			//
			Vector3r erosion = v.normalized()*margin;
			Vector3r  p = particle1->support(-v) + position1;//support point of object A
			Vector3r  q = particle2->support(v) + position2;//support point of object B
			//cout<<"v0="<<v.normalized()<<endl;
			//cout<<"p="<<p<<"q="<<q<<endl;
			#ifdef STATISTICS
			++num_iterations;
			#endif
			Vector3r w = p - q;//support point of Minkowski difference A-B
			// SD_Scalar vn = v.norm();
			//
			// if (vn > 0.0)
			// {
			//  vn = 1/ vn;
			//  w += v*vn*margin;
			// }

			SD_Scalar delta = v.dot(w);

			if(delta > 0.0){//not contact//&& delta * delta > dist2 * margin2
				#ifdef STATISTICS
				cout<<"nfuc1="<<num_iterations<<endl;
				#endif
				return false;
			}

			double angle = delta/(v.norm()*w.norm());//phi in [0, pi]
			//cout<<"angle"<<angle<<"wm"<<w.norm()/margin<<endl;
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



			if (!gjk.closest(v))
			{
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= erosion;//v * (margin / s);
				pb += erosion;//v * (margin / s);
				#ifdef STATISTICS
				cout<<"nfuc4="<<num_iterations<<"p/m="<<(pa-pb).norm()/margin<<endl;
				#endif
				return true;
			}

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
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= erosion;//v * (margin / s);
				pb += erosion;//v * (margin / s);
				#ifdef STATISTICS
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

template <typename T>
bool getPenetration2(T *particle1, T *particle2,	Vector3r position1, Vector3r position2,
										Vector3r& v, Vector3r& pa, Vector3r& pb)
{
	SD_Scalar margin = min(particle1->getMargin(),particle1->getMargin());
	//cout<<"margin="<<margin<<endl;
	//cout<<"hybrid_depth:into "<<endl;
	if (margin > SD_Scalar(0.0))
	{
		SD_Scalar margin2 = 4.0*margin * margin;

		DT_GJK gjk;

#ifdef STATISTICS
		num_iterations = 0;
#endif
		SD_Scalar dist2 = Mathr::MAX_REAL;
		//cout<<"using margin"<<endl;
		do
		{
			//
			//cout<<"v"<<v<<"p1"<<particle1->support(-v)<<"p2"<<particle2->support(v)<<endl;
			Vector3r  p = particle1->support(-v) + position1;//support point of object A
			Vector3r  q = particle2->support(v) + position2;//support point of object B


			Vector3r w = p - q;//support point of Minkowski difference A-B
			//cout<<"p"<<p<<"q"<<q<<"w"<<w<<endl;
			// SD_Scalar vn = v.norm();
			//
			// if (vn > 0.0)
			// {
			//  vn = 1/ vn;
			//  w += v*vn*margin;
			// }

			SD_Scalar delta = v.dot(w);

			if (delta > SD_Scalar(0.0) && delta * delta > dist2 * margin2)//the enlarged objects do not interact.
			{
				return false;
			}

			if (gjk.inSimplex(w) || dist2 - delta <= dist2 * 1e-6)
			{ //the objects only interact in the margin
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= v * (margin / s);
				pb += v * (margin / s);
				return true;
			}

			gjk.addVertex(w, p, q);//add a vertex to the convex hull

      if (gjk.isAffinelyDependent())
      {
				#ifdef STATISTICS
                ++num_irregularities;
				#endif
        gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= v * (margin / s);
				pb += v * (margin / s);
				return true;
      }


			#ifdef STATISTICS
			++num_iterations;
			#endif
			if (!gjk.closest(v))
			{
				#ifdef STATISTICS
				++num_irregularities;
				#endif
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= v * (margin / s);
				pb += v * (margin / s);
				return true;
			}

			#ifdef SAFE_EXIT
			SD_Scalar prev_dist2 = dist2;
			#endif

			dist2 = v.squaredNorm();

			#ifdef SAFE_EXIT
			if (prev_dist2 - dist2 <= 1e-10 * prev_dist2)
			{
  			gjk.backup_closest(v);
				dist2 = v.squaredNorm();
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= v * (margin / s);
				pb += v * (margin / s);
				return true;
			}
			#endif
		}
		while (!gjk.fullSimplex() && dist2 > 1e-16 * gjk.maxVertex());

	}
	//cout<<"second GJK phase"<<endl;
	// Second GJK phase. compute points on the boundary of the offset object
//	return penetration_depth(a, b, v, pa, pb);
 cout<<"second phase is not allowed!"<<endl;
 return false;

}
//**********************************************************************************
/*! Create GJKParticle (collision geometry) from colliding GJKParticles. */

bool Ig2_GJKParticle_GJKParticle_GJKParticleGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){
    //std::cerr << "watch 1"<< "\n";
	//get polyhedras
	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
	Vector3r position1 = se31.position;
	Vector3r position2 = se32.position;
    //std::cerr << "watch 2"<< "\n";
	GJKParticle* A = static_cast<GJKParticle*>(cm1.get());
	GJKParticle* B = static_cast<GJKParticle*>(cm2.get());
        //cout<<"watch point"<<endl;
	bool isNew = !interaction->geom;

	shared_ptr<GJKParticleGeom> bang;
	//Vector2r dels(0.01,0.01);
	Vector2r para(0.0,0.0);
	Vector3r v(0.0,0.0,0.0);//separating axis vector
    v = position1 -position2;
    //std::cerr << "watch 3"<< "\n";
	if (isNew) {
		// new interaction
		//cout<<"new intersection"<<endl;
		bang=shared_ptr<GJKParticleGeom>(new GJKParticleGeom());
		//bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
		//dels = Vector2r(0.1,0.1);
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<GJKParticleGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
        v = bang->contactSA;
		//Vector3r v1 = state1.vel;
		//Vector3r v2 = state2.vel;//maybe a bug
		//cout<<"vvv1"<<v1(0)<<" "<<v1(1)<<" "<<v1(2)<<"  "<<v1.isZero(0)<<endl;
		//cout<<"vvv2"<<v2(0)<<" "<<v2(1)<<" "<<v2(2)<<"  "<<v2.isZero(0)<<endl;

		//if (v1.isZero(0)&&v2.isZero(0)){
		//        bang->PenetrationDepth=0;

               //         return true;
		//}
		//dels = Vector2r(0.001,0.001);
	}
	Vector3r contactNormal(0,0,0), center(0,0,0);
	Vector3r p1(0.0,0.0,0.0),p2(0.0,0.0,0.0);

	//if the two particles are spherical
	if (!A->getShapeType() && !B->getShapeType()){//shapeType 0 is for sphere
	        //particles are spherical
	        Vector3r dist = position2 -position1;
	        double r_sum = A->getRadius()+B->getRadius();
	        double depth = r_sum - dist.norm();
	        if (depth < 1E-18){//no contact
                        bang->PenetrationDepth=0;
                        bang->isShearNew = true;
                        return true;
                }
                //compute contact geometric quantities
                dist.normalize();
                p1 = position1 + A->getRadius()*dist;
                p2 = position2 - B->getRadius()*dist;
                contactNormal = p1 - p2;
	        contactNormal.normalize();
	        bang->point1 = p1;
                bang->point2 = p2;

	        bang->contactPoint=0.5*(p1+p2);
	        bang->PenetrationDepth = depth;
	}else{
	//contact detection using bounding spheres
     /*   double r_sum = A->getRadius()+B->getRadius();

        if ((position1 - position2).norm() > r_sum){//no contact
                bang->PenetrationDepth=0;
                bang->isShearNew = true;
                return true;
        }
	*/
        Matrix3r globalContact = Matrix3r::Identity();

	//v = bang->contactSA;
	bool touching(true);
	//A->setPosition(position1);
	//B->setPosition(position2);
	//A->setOrientation(se31.orientation);
	//B->setOrientation(se32.orientation);
	//touching = penetration_depth(A, B, v, p1, p2);
	//Matrix3r rot_matA = (se31.orientation).toRotationMatrix();//to global system
	//Matrix3r rot_matB = (se32.orientation).toRotationMatrix();//to global system
	//Vector3r v1 = v;
	//touching = penetration_depth(A,rot_matA,position1,B,rot_matB,position2,v, p1,p2);
	//cout<<"v1"<<v<<endl;
	touching = getPenetration2(A,B, position1, position2,v, p1, p2);
	//cout<<"v2="<<v<<endl;
	if (!touching){
	        //bang->contactAngle = para;
			bang->contactSA = v;
	        bang->PenetrationDepth=0;
	        bang->isShearNew = true;
	        return true;}
	//for testing the two closest points
	//std::cerr << "watch 5"<< "\n";
	contactNormal = p1 - p2;
	//debugging
	/*bug fixed
	if(isinf(contactNormal.norm())){//zero contact normal:the bug arises when reloading data.
		//char filename[20];
		time_t tmNow = time(NULL);
		char tmp[64];
		strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&tmNow) );
		const char* filename = "JGKParticleLaw_err.log";
		//sprintf(filename,"%d",gettid4());
		FILE * fin = fopen(filename,"a");
		fprintf(fin,"\n**Error in Comp. of Contact Normal*****%s*******\n",tmp);

		fprintf(fin,"point1\t%e\t%e\t%e\n",p1[0],p1[1],p1[2]);
		fprintf(fin,"point2\t%e\t%e\t%e\n",p2[0],p2[1],p2[2]);
		fprintf(fin,"v\t%e\t%e\t%e\n",v[0],v[1],v[2]);
		fprintf(fin,"v1\t%e\t%e\t%e\n",v1[0],v1[1],v1[2]);

		fclose(fin);
		//throw runtime_error("P4 #");
				LOG_WARN("Error in computation of contact normal");
				scene->stopAtIter = scene->iter + 1;
	}
	*/
	//
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
	bang->PenetrationDepth=contactNormal.norm();
	/*if (!touching){
			bang->contactSA = v;
	        bang->PenetrationDepth=0;
	        bang->isShearNew = true;
	}*/
	//std::cerr << "watchpoint 3: " << "\n";
	}//non-spherical end
    contactNormal.normalize();
		bang->point11 = p1;//contactNormal*0.5+p1;
		bang->point22 = p2;//-contactNormal*0.5+p2;
	bang->precompute(state1,state2,scene,interaction,contactNormal,bang->isShearNew,shift2);
	//contactNormal.normalize();
	//bang->normal= contactNormal;
	bang->contactSA = v;
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
/*! Create GJKParticle (collision geometry) from colliding Polyhedron and Wall. */

bool Ig2_Wall_GJKParticle_GJKParticleGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){

	const int& PA(cm1->cast<Wall>().axis);
	const int& sense(cm1->cast<Wall>().sense);
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	GJKParticle* B = static_cast<GJKParticle*>(cm2.get());

	bool isNew = !interaction->geom;
	shared_ptr<GJKParticleGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<GJKParticleGeom>(new GJKParticleGeom());
		bang->contactAngle = Vector2r(0,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<GJKParticleGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}

	Vector3r normal = Vector3r(0,0,0);
	normal[PA] = sense;
	Vector3r wall_pos = se31.position;
	Vector3r B_pos = se32.position;
        // is the particle spherical?
        if (!B->getShapeType()){
                //the particle is spherical
                Vector3r p1, point2;
                p1 = B_pos;
                p1[PA] = wall_pos[PA];//projection of the particle center on the wall
                point2 = B_pos - normal*B->getRadius();

                double depth = B->getRadius() -(B_pos - p1).norm();
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

	//Matrix3r rot_mat1 = B->rot_mat2local;//(se32.orientation).conjugate().toRotationMatrix();//to particle's system
	//Matrix3r rot_mat2 = B->rot_mat2global;//(se32.orientation).toRotationMatrix();//to global system
	//Matrix3r rot1 = (se32.orientation).conjugate().toRotationMatrix();//to particle's system
	//Matrix3r rot2 = se32.orientation.toRotationMatrix();

	Vector3r point2 =  B->support(-normal);
  point2 += se32.position;

	//cout<<"wall_pos"<<wall_pos<<endl;
	//cout<<"B_pos"<<B_pos<<endl;

	//Vector2r phi = B->Normal2Phi(-rot_mat1*normal);//phi in particle's local system
	//Vector3r point2 = rot_mat2*( B->getSurface(phi)) + B_pos;	//the closest point on particle B to the wall
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
		bang->point11 = p1;//normal*0.5+p1;
		bang->point22 = point2;//-normal*0.5+point2;
	}//non-spherical end
	//std::cerr << "watchpoint 3: " << "\n";
	bang->precompute(state1,state2,scene,interaction,normal,bang->isShearNew,shift2);
	//bang->normal= normal;
	//bang->contactAngle = para;

	return true;
}

//**********************************************************************************
/*! Create GJKParticle (collision geometry) from colliding Polyhedron and Facet. */

bool Ig2_Facet_GJKParticle_GJKParticleGeom::go(const shared_ptr<Shape>& cm1,const shared_ptr<Shape>& cm2,const State& state1,const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& interaction){


	const Se3r& se31=state1.se3;
	const Se3r& se32=state2.se3;
	Facet*   A = static_cast<Facet*>(cm1.get());
	GJKParticle* B = static_cast<GJKParticle*>(cm2.get());

	bool isNew = !interaction->geom;

	shared_ptr<GJKParticleGeom> bang;
	if (isNew) {
		// new interaction
		bang=shared_ptr<GJKParticleGeom>(new GJKParticleGeom());
		//bang->sep_plane.assign(3,0);
		bang->contactPoint = Vector3r(0,0,0);
		bang->isShearNew = true;
		interaction->geom = bang;
	}else{
		// use data from old interaction
  		bang=SUDODEM_PTR_CAST<GJKParticleGeom>(interaction->geom);
		bang->isShearNew = bang->PenetrationDepth<=0;
	}
/*
	//find intersection GJKParticle
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
