#ifndef GJKPARTICLE_IG2_H
#define GJKPARTICLE_IG2_H
#pragma once
#include"GJKParticle.hpp"
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/Box.hpp>
#include"GJK.hpp"
//---------------------------------------------------------------------------
static const SD_Scalar rel_error = SD_Scalar(1.0e-3);

SD_Scalar DT_Accuracy::rel_error2 = rel_error * rel_error;
SD_Scalar DT_Accuracy::depth_tolerance = SD_Scalar(1.0) + SD_Scalar(2.0) * rel_error;
SD_Scalar DT_Accuracy::tol_error = Mathr::EPSILON;
class TriangleComp
{
public:

    bool operator()(const Triangle *face1, const Triangle *face2)
    {
        return face1->getDist2() > face2->getDist2();
    }
} triangleComp;

class SolvePenDepth
{
private:
    int       MaxSupportPoints = 100;
    int       MaxFacets         = 200;

    Vector3r  pBuf[100];
    Vector3r  qBuf[100];
    Vector3r  yBuf[100];


    Triangle *triangleHeap[200];
    int  num_triangles;
public:
    SolvePenDepth(){
    //int       MaxSupportPoints = 100;
    //int       MaxFacets         = 200;
    }
    inline void addCandidate(Triangle *triangle, SD_Scalar upper2)
//inline void SolvePenDepth::addCandidate(Triangle *triangle, SD_Scalar upper2)
    {
        if (triangle->isClosestInternal() && triangle->getDist2() <= upper2)
        {
            triangleHeap[num_triangles++] = triangle;
            std::push_heap(&triangleHeap[0], &triangleHeap[num_triangles], triangleComp);//the triangle with the shortest distance is at the top of the heap
    #ifdef DEBUG
            std::cout << " accepted" << std::endl;
    #endif
        }
        else
        {
    #ifdef DEBUG
            std::cout << " rejected, ";
            if (!triangle->isClosestInternal())
                {
                    std::cout << "closest point not internal";
                }
            else
            {
                std::cout << "triangle is further than upper bound";
            }
            std::cout << std::endl;
    #endif
        }
    }

    inline int originInTetrahedron(const Vector3r& p1, const Vector3r& p2,
                                   const Vector3r& p3, const Vector3r& p4)
    {
        Vector3r normal1 = (p2 - p1).cross(p3 - p1);
	    if ((normal1.dot(p1) > SD_Scalar(0.0)) == (normal1.dot(p4) > SD_Scalar(0.0)))
        {
            return 4;
        }

        Vector3r normal2 = (p4 - p2).cross(p3 - p2);
        if ((normal2.dot(p2) > SD_Scalar(0.0)) == (normal2.dot(p1) > SD_Scalar(0.0)))
        {
            return 1;
        }

        Vector3r normal3 = (p4 - p3).cross(p1 - p3);
        if ((normal3.dot(p3) > SD_Scalar(0.0)) == (normal3.dot(p2) > SD_Scalar(0.0)))
        {
            return 2;
        }

        Vector3r normal4 = (p2 - p4).cross(p1 - p4);
        if ((normal4.dot(p4) > SD_Scalar(0.0)) == (normal4.dot(p3) > SD_Scalar(0.0)))
        {
            return 3;
        }

        return 0;
    }

    bool penDepth(const DT_GJK& gjk, const DT_Convex& a, const DT_Convex& b,
                  Vector3r& v, Vector3r& pa, Vector3r& pb)
    {
        //get some space
        //const int       MaxSupportPoints = 100;
        //Vector3r  pBuf[MaxSupportPoints];
        //Vector3r  qBuf[MaxSupportPoints];
        //Vector3r yBuf[MaxSupportPoints];

	    int num_verts = gjk.getSimplex(pBuf, qBuf, yBuf);
        SD_Scalar tolerance = Mathr::EPSILON * gjk.maxVertex();

        num_triangles = 0;
        TriangleStore g_triangleStore;
        g_triangleStore.clear();

        switch (num_verts)
        {
        case 1:
            // Touching contact. Yes, we have a collision,
            // but no penetration.
            return false;
        case 2:
        {
            // We have a line segment inside the Minkowski sum containing the
            // origin. Blow it up by adding three additional support points.

            Vector3r dir  = yBuf[1] - yBuf[0];

            dir /= dir.norm();
		        Vector3r d_abs = dir.cwiseAbs();
		        int axis = d_abs[0] < d_abs[1] ? (d_abs[0] < d_abs[2] ? 0 : 2) : (d_abs[1] < d_abs[2] ? 1 : 2);

            //int        axis = dir.furthestAxis();

            static const SD_Scalar sin_60 = sqrt(SD_Scalar(3.0)) * SD_Scalar(0.5);

            Quaternionr rot(dir[0] * sin_60, dir[1] * sin_60, dir[2] * sin_60, SD_Scalar(0.5));
            Matrix3r rot_mat(rot);

            Vector3r aux1 = dir.cross(Vector3r(axis == 0, axis == 1, axis == 2));
            Vector3r aux2 = rot_mat * aux1;
            Vector3r aux3 = rot_mat * aux2;

            pBuf[2] = a.support(aux1);
            qBuf[2] = b.support(-aux1);
            yBuf[2] = pBuf[2] - qBuf[2];

            pBuf[3] = a.support(aux2);
            qBuf[3] = b.support(-aux2);
            yBuf[3] = pBuf[3] - qBuf[3];

            pBuf[4] = a.support(aux3);
            qBuf[4] = b.support(-aux3);
            yBuf[4] = pBuf[4] - qBuf[4];

            if (originInTetrahedron(yBuf[0], yBuf[2], yBuf[3], yBuf[4]) == 0)
            {
                pBuf[1] = pBuf[4];
                qBuf[1] = qBuf[4];
                yBuf[1] = yBuf[4];
            }
            else if (originInTetrahedron(yBuf[1], yBuf[2], yBuf[3], yBuf[4]) == 0)
            {
                pBuf[0] = pBuf[4];
                qBuf[0] = qBuf[4];
                yBuf[0] = yBuf[4];
            }
            else
            {
                // Origin not in initial polytope
                return false;
            }

            num_verts = 4;
        }
        // Fall through allowed!!
        case 4:
        {
            int bad_vertex = originInTetrahedron(yBuf[0], yBuf[1], yBuf[2], yBuf[3]);

            if (bad_vertex == 0)
            {
                Triangle *f0 = g_triangleStore.newTriangle(yBuf, 0, 1, 2);
                Triangle *f1 = g_triangleStore.newTriangle(yBuf, 0, 3, 1);
                Triangle *f2 = g_triangleStore.newTriangle(yBuf, 0, 2, 3);
                Triangle *f3 = g_triangleStore.newTriangle(yBuf, 1, 3, 2);

                if (!(f0 && f0->getDist2() > SD_Scalar(0.0) &&
                      f1 && f1->getDist2() > SD_Scalar(0.0) &&
                      f2 && f2->getDist2() > SD_Scalar(0.0) &&
                      f3 && f3->getDist2() > SD_Scalar(0.0)))
			    {
				    return false;
			    }

                link(Edge(f0, 0), Edge(f1, 2));
                link(Edge(f0, 1), Edge(f3, 2));
                link(Edge(f0, 2), Edge(f2, 0));
                link(Edge(f1, 0), Edge(f2, 2));
                link(Edge(f1, 1), Edge(f3, 0));
                link(Edge(f2, 1), Edge(f3, 1));

                addCandidate(f0, Mathr::MAX_REAL);
                addCandidate(f1, Mathr::MAX_REAL);
                addCandidate(f2, Mathr::MAX_REAL);
                addCandidate(f3, Mathr::MAX_REAL);
                break;
            }

            if (bad_vertex < 4)
            {
                pBuf[bad_vertex - 1] = pBuf[4];
                qBuf[bad_vertex - 1] = qBuf[4];
                yBuf[bad_vertex - 1] = yBuf[4];

            }

            num_verts = 3;

        }
        // Fall through allowed!!
        case 3:
        {
            // We have a triangle inside the Minkowski sum containing
            // the origin. First blow it up.

            Vector3r v1     = yBuf[1] - yBuf[0];
            Vector3r v2     = yBuf[2] - yBuf[0];
            Vector3r vv     = v1.cross(v2);

            pBuf[3] = a.support(vv);
            qBuf[3] = b.support(-vv);
            yBuf[3] = pBuf[3] - qBuf[3];
            pBuf[4] = a.support(-vv);
            qBuf[4] = b.support(vv);
            yBuf[4] = pBuf[4] - qBuf[4];

            Triangle* f0 = g_triangleStore.newTriangle(yBuf, 0, 1, 3);
            Triangle* f1 = g_triangleStore.newTriangle(yBuf, 1, 2, 3);
            Triangle* f2 = g_triangleStore.newTriangle(yBuf, 2, 0, 3);
            Triangle* f3 = g_triangleStore.newTriangle(yBuf, 0, 2, 4);
            Triangle* f4 = g_triangleStore.newTriangle(yBuf, 2, 1, 4);
            Triangle* f5 = g_triangleStore.newTriangle(yBuf, 1, 0, 4);

            if (!(f0 && f0->getDist2() > SD_Scalar(0.0) &&
                  f1 && f1->getDist2() > SD_Scalar(0.0) &&
                  f2 && f2->getDist2() > SD_Scalar(0.0) &&
                  f3 && f3->getDist2() > SD_Scalar(0.0) &&
                  f4 && f4->getDist2() > SD_Scalar(0.0) &&
                  f5 && f5->getDist2() > SD_Scalar(0.0)))
            {
                return false;
            }

            link(Edge(f0, 1), Edge(f1, 2));
            link(Edge(f1, 1), Edge(f2, 2));
            link(Edge(f2, 1), Edge(f0, 2));

            link(Edge(f0, 0), Edge(f5, 0));
            link(Edge(f1, 0), Edge(f4, 0));
            link(Edge(f2, 0), Edge(f3, 0));

            link(Edge(f3, 1), Edge(f4, 2));
            link(Edge(f4, 1), Edge(f5, 2));
            link(Edge(f5, 1), Edge(f3, 2));

            addCandidate(f0, Mathr::MAX_REAL);
            addCandidate(f1, Mathr::MAX_REAL);
            addCandidate(f2, Mathr::MAX_REAL);
            addCandidate(f3, Mathr::MAX_REAL);
            addCandidate(f4, Mathr::MAX_REAL);
            addCandidate(f5, Mathr::MAX_REAL);

            num_verts = 5;
        }
        break;
        }

        // We have a polytope inside the Minkowski sum containing
        // the origin.

        if (num_triangles == 0)
        {
            return false;
        }

        // at least one triangle on the heap.

        Triangle *triangle = 0;

        SD_Scalar upper_bound2 = Mathr::MAX_REAL;

        do
        {
            triangle = triangleHeap[0];
            std::pop_heap(&triangleHeap[0], &triangleHeap[num_triangles], triangleComp);
            --num_triangles;

            if (!triangle->isObsolete())
            {
                if (num_verts == MaxSupportPoints)
                {
    #ifdef DEBUG
                    std::cout << "Ouch, no convergence!!!" << std::endl;
    #endif
                    assert(false);
                    break;
                }

                pBuf[num_verts] = a.support( triangle->getClosest());
                qBuf[num_verts] = b.support(-triangle->getClosest());
                yBuf[num_verts] = pBuf[num_verts] - qBuf[num_verts];

                int index = num_verts++;
                SD_Scalar far_dist = yBuf[index].dot(triangle->getClosest());

                // Make sure the support mapping is OK.
                assert(far_dist > SD_Scalar(0.0));
                SD_Scalar far_dist2 = far_dist * far_dist / triangle->getDist2();
                Math_set_min(upper_bound2, far_dist2);

                SD_Scalar error = far_dist - triangle->getDist2();
                if (error <= Math_max(1e-6 * far_dist, tolerance)
    #if 1
                    || yBuf[index] == yBuf[(*triangle)[0]]
                    || yBuf[index] == yBuf[(*triangle)[1]]
                    || yBuf[index] == yBuf[(*triangle)[2]]
    #endif
                    )
                {
                    break;
                }


                // Compute the silhouette cast by the new vertex
                // Note that the new vertex is on the positive side
                // of the current triangle, so the current triangle will
                // not be in the convex hull. Start local search
                // from this triangle.

                int i = g_triangleStore.getFree();

                if (!triangle->silhouette(yBuf, index, g_triangleStore))
                {
                    break;
                }

                while (i != g_triangleStore.getFree())
                {
                    Triangle *newTriangle = &g_triangleStore[i];
                    //assert(triangle->getDist2() <= newTriangle->getDist2());

                    addCandidate(newTriangle, upper_bound2);

                    ++i;
                }
            }
        }
        while (num_triangles > 0 && triangleHeap[0]->getDist2() <= upper_bound2);

    #ifdef DEBUG
        std::cout << "#triangles left = " << num_triangles << std::endl;
    #endif

        v = triangle->getClosest();
        pa = triangle->getClosestPoint(pBuf);
        pb = triangle->getClosestPoint(qBuf);
        return true;
    }
};
/*
///////////////////////////////////
bool intersect(const DT_Convex& a, const DT_Convex& b, Vector3r& v)
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

        if (v.dot(w) > SD_Scalar(0.0))
		{
			return false;
		}

		gjk.addVertex(w);
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

    v = Vector3r::Zero();//.setValue(SD_Scalar(0.0), SD_Scalar(0.0), SD_Scalar(0.0));

    return true;
}




bool common_point(const DT_Convex& a, const DT_Convex& b,
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

        if (v.dot(w) > SD_Scalar(0.0))
		{
			return false;
		}

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

	gjk.compute_points(pa, pb);

    v = Vector3r::Zero();//.setValue(SD_Scalar(0.0), SD_Scalar(0.0), SD_Scalar(0.0));

    return true;
}

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
*/
/*
bool hybrid_penetration_depth(const DT_Convex& a, SD_Scalar a_margin,
							  const DT_Convex& b, SD_Scalar b_margin,
                              Vector3r& v, Vector3r& pa, Vector3r& pb)
{
	SD_Scalar margin = a_margin + b_margin;
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
		do
		{
			//

			Vector3r  p = a.support(-v);//support point of object A
			Vector3r  q = b.support(v);//support point of object B


			Vector3r w = p - q;//support point of Minkowski difference A-B
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

			if (gjk.inSimplex(w) || dist2 - delta <= dist2 * DT_Accuracy::rel_error2)
			{ //the objects only interact in the margin
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= v * (a_margin / s);
				pb += v * (b_margin / s);
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
				pa -= v * (a_margin / s);
				pb += v * (b_margin / s);
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
				pa -= v * (a_margin / s);
				pb += v * (b_margin / s);
				return true;
			}

			#ifdef SAFE_EXIT
			SD_Scalar prev_dist2 = dist2;
			#endif

			dist2 = v.squaredNorm();

			#ifdef SAFE_EXIT
			if (prev_dist2 - dist2 <= Mathr::EPSILON * prev_dist2)
			{
  			gjk.backup_closest(v);
				dist2 = v.squaredNorm();
				gjk.compute_points(pa, pb);
				SD_Scalar s = sqrt(dist2);
				assert(s > SD_Scalar(0.0));
				pa -= v * (a_margin / s);
				pb += v * (b_margin / s);
				return true;
			}
			#endif
		}
		while (!gjk.fullSimplex() && dist2 > DT_Accuracy::tol_error * gjk.maxVertex());

	}
	//cout<<"second GJK phase"<<endl;
	// Second GJK phase. compute points on the boundary of the offset object
//	return penetration_depth(a, b, v, pa, pb);

	return penetration_depth((a_margin > SD_Scalar(0.0) ?
							  static_cast<const DT_Convex&>(DT_Minkowski(a, DT_Sphere(a_margin))) :
							  static_cast<const DT_Convex&>(a)),
							 (b_margin > SD_Scalar(0.0) ?
							  static_cast<const DT_Convex&>(DT_Minkowski(b, DT_Sphere(b_margin))) :
							  static_cast<const DT_Convex&>(b)), v, pa, pb);

}

*/


/*
SD_Scalar closest_points(const DT_Convex& a, const DT_Convex& b, SD_Scalar max_dist2,
                         Vector3r& pa, Vector3r& pb)
{
	Vector3r v(SD_Scalar(0.0), SD_Scalar(0.0), SD_Scalar(0.0));

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

		SD_Scalar delta = v.dot(w);
		if (delta > SD_Scalar(0.0) && delta * delta > dist2 * max_dist2)
		{
			return Mathr::MAX_REAL;
		}

		if (gjk.inSimplex(w) || dist2 - delta <= dist2 * DT_Accuracy::rel_error2)
		{
            break;
		}

		gjk.addVertex(w, p, q);
        if (gjk.isAffinelyDependent())
        {
#ifdef STATISTICS
            ++num_irregularities;
#endif
            break;
        }

#ifdef STATISTICS
        ++num_iterations;
        if (num_iterations > 1000)
		{
			std::cout << "v: " << v << " w: " << w << std::endl;
		}
#endif
        if (!gjk.closest(v))
		{
#ifdef STATISTICS
            ++num_irregularities;
#endif
            break;
        }

#ifdef SAFE_EXIT
		SD_Scalar prev_dist2 = dist2;
#endif

		dist2 = v.squaredNorm();

#ifdef SAFE_EXIT
		if (prev_dist2 - dist2 <= Mathr::EPSILON * prev_dist2)
		{
            gjk.backup_closest(v);
            dist2 = v.squaredNorm();
			break;
		}
#endif
    }
    while (!gjk.fullSimplex() && dist2 > DT_Accuracy::tol_error * gjk.maxVertex());

	assert(!gjk.emptySimplex());

	if (dist2 <= max_dist2)
	{
		gjk.compute_points(pa, pb);
	}

	return dist2;
}
*/
////////////////////////////////////

//***************************************************************************
/*! Create GJKParticle (collision geometry) from colliding GJKParticles. */
class Ig2_GJKParticle_GJKParticle_GJKParticleGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_GJKParticle_GJKParticle_GJKParticleGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(GJKParticle,GJKParticle);
		DEFINE_FUNCTOR_ORDER_2D(GJKParticle,GJKParticle);
		SUDODEM_CLASS_BASE_DOC(Ig2_GJKParticle_GJKParticle_GJKParticleGeom,IGeomFunctor,"Create/update geometry of collision between 2 GJKParticles");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_GJKParticle_GJKParticle_GJKParticleGeom);

//***************************************************************************
/*! Create GJKParticle (collision geometry) from colliding Wall & GJKParticle. */
class Ig2_Wall_GJKParticle_GJKParticleGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Wall_GJKParticle_GJKParticleGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Wall,GJKParticle);
		DEFINE_FUNCTOR_ORDER_2D(Wall,GJKParticle);
		SUDODEM_CLASS_BASE_DOC(Ig2_Wall_GJKParticle_GJKParticleGeom,IGeomFunctor,"Create/update geometry of collision between Wall and GJKParticle");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Wall_GJKParticle_GJKParticleGeom);
//***************************************************************************
/*! Create GJKParticle (collision geometry) from colliding Facet & GJKParticle. */
class Ig2_Facet_GJKParticle_GJKParticleGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Facet_GJKParticle_GJKParticleGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Facet,GJKParticle);
		DEFINE_FUNCTOR_ORDER_2D(Facet,GJKParticle);
		SUDODEM_CLASS_BASE_DOC(Ig2_Facet_GJKParticle_GJKParticleGeom,IGeomFunctor,"Create/update geometry of collision between Facet and GJKParticle");
		DECLARE_LOGGER;
	private:
};
REGISTER_SERIALIZABLE(Ig2_Facet_GJKParticle_GJKParticleGeom);
#endif
