#ifndef GJKPARTICLE_SHAPES_H
#define GJKPARTICLE_SHAPES_H

#include <vector>
#include <algorithm>
#include <new>

#include<sudodem/lib/base/Math.hpp>

typedef MT::Transform<SD_Scalar> MT_Transform;
typedef DT_Array<unsigned int, unsigned int> DT_IndexArray;

//#define DEBUG
//#define SAFE_EXIT

#ifdef STATISTICS
int num_iterations = 0;
int num_irregularities = 0;
#endif

enum DT_ShapeType {
    COMPLEX,
    CONVEX
};

class DT_Convex{
public:
  virtual ~DT_Convex() {}
	virtual DT_ShapeType getType() const { return CONVEX; }

	virtual SD_Scalar supportH(const Vector3r& v) const {	return v.dot(support(v)); }
  virtual Vector3r support(const Vector3r& v) const = 0;
protected:
	DT_Convex() {}
};

//Point
class DT_Point : public DT_Convex {
public:
    DT_Point(const Vector3r& point) : m_point(point) {}
	SD_Scalar supportH(const Vector3r& v) const
	{
		return v.dot(m_point);
	}

	Vector3r support(const Vector3r&) const
	{
		return m_point;
	}
private:
	Vector3r m_point;
};

//LineSegment
class DT_LineSegment : public DT_Convex {
public:
    DT_LineSegment(const Vector3r& end1, const Vector3r& end2) :
	   m_end1(end1),
	   m_end2(end2) {}

	SD_Scalar supportH(const Vector3r& v) const
	{
		return Math_max(v.dot(m_end1),v.dot(m_end2));
	}

	Vector3r support(const Vector3r& v) const
	{
		return v.dot(m_end1) > v.dot(m_end2) ?	m_end1 : m_end2;
	}
	const Vector3r& getEnd1() const { return m_end1; }
	const Vector3r& getEnd2() const { return m_end2; }

private:
	Vector3r m_end1;//point at end1
	Vector3r m_end2;//point at end2
};

//TriEdge
class Edge;
class Triangle;
class TriangleStore;

typedef unsigned short Index_t;

bool link(const Edge& edge0, const Edge& edge1);
void half_link(const Edge& edge0, const Edge& edge1);


class Edge {
private:
    Triangle *m_triangle;
    int       m_index;

public:
    Edge() {}
    Edge(Triangle *triangle, int index)
      : m_triangle(triangle),
        m_index(index)
	{}

    bool silhouette(const Vector3r *verts, Index_t index, TriangleStore& triangleStore) const;

    Triangle *triangle() const { return m_triangle; }
    int       index()    const { return m_index; }

	Index_t  getEnd1() const;//end point 1
	Index_t  getEnd2() const;//end point 2

};




class Triangle {
private:
    Index_t     m_indices[3];//indices of three vertices
    Edge        m_adjEdges[3];//three edges

    bool        m_obsolete;

    SD_Scalar   m_det;
    SD_Scalar   m_lambda1;
    SD_Scalar   m_lambda2;
    Vector3r  m_closest;
    SD_Scalar   m_dist2;

public:
    Triangle() {}
    Triangle(Index_t i0, Index_t i1, Index_t i2)
	  :	m_obsolete(false)//construct a triangle with the indices of vertices
    {
		m_indices[0] = i0;
		m_indices[1] = i1;
		m_indices[2] = i2;
    }

    Index_t operator[](int i) const { return m_indices[i]; }

	const Edge& getAdjEdge(int i) const { return m_adjEdges[i]; }

    void setObsolete(bool obsolete)
	{
#ifdef DEBUG
		std::cout << "Triangle " <<  m_indices[0] << ' ' << m_indices[1] << ' ' << m_indices[2] << " obsolete" << std::endl;
#endif
		m_obsolete = obsolete;
	}

	bool isObsolete() const { return m_obsolete; }


    bool computeClosest(const Vector3r *verts);

    const Vector3r& getClosest() const { return m_closest; }

    bool isClosestInternal() const
	{
		return m_lambda1 >= SD_Scalar(0.0) &&
			   m_lambda2 >= SD_Scalar(0.0) &&
			   m_lambda1 + m_lambda2 <= m_det;
    }

	bool isVisibleFrom(const Vector3r *verts, Index_t index) const
	{
		Vector3r lever = verts[index] - m_closest;
		return m_closest.dot(lever) > SD_Scalar(0.0);
	}

    SD_Scalar getDist2() const { return m_dist2; }

    Vector3r getClosestPoint(const Vector3r *points) const
	{
		const Vector3r& p0 = points[m_indices[0]];

		return p0 + (m_lambda1 * (points[m_indices[1]] - p0) +
					 m_lambda2 * (points[m_indices[2]] - p0)) / m_det;
    }

    bool silhouette(const Vector3r *verts, Index_t index, TriangleStore& triangleStore);

	friend bool link(const Edge& edge0, const Edge& edge1);
	friend void half_link(const Edge& edge0, const Edge& edge1);
};


static const int MaxTriangles = 200;

class TriangleStore
{
private:
	Triangle m_triangles[MaxTriangles];//alocate mem
	int      m_free;
public:
	TriangleStore()
	  : m_free(0)
	{}

	void clear() { m_free = 0; }

	int getFree() const { return m_free; }

	Triangle& operator[](int i) { return m_triangles[i]; }
	Triangle& last() { return m_triangles[m_free - 1]; }

	void setFree(int backup) { m_free = backup; }


	Triangle *newTriangle(const Vector3r *verts, Index_t i0, Index_t i1, Index_t i2)
	{
		Triangle *newTriangle = 0;
		if (m_free != MaxTriangles)
		{
			newTriangle = &m_triangles[m_free++];
			new (newTriangle) Triangle(i0, i1, i2);
			if (!newTriangle->computeClosest(verts))
			{
				--m_free;
				newTriangle = 0;
			}
		}

		return newTriangle;
	}
};

//extern TriangleStore g_triangleStore;
inline int circ_next(int i) { return (i + 1) % 3; }
inline int circ_prev(int i) { return (i + 2) % 3; }

inline Index_t Edge::getEnd1() const
{
	return (*m_triangle)[m_index];
}

inline Index_t Edge::getEnd2() const
{
	return (*m_triangle)[circ_next(m_index)];
}

///////////////////////////
///////////////////Triangle
class DT_Triangle : public DT_Convex {
public:
	DT_Triangle(Vector3r* verts)
	{
		m_vertices[0] = verts[0];
		m_vertices[1] = verts[1];
		m_vertices[2] = verts[2];
	}
	DT_Triangle(Vector3r a, Vector3r b,Vector3r c)
	{
		m_vertices[0] = a;
		m_vertices[1] = b;
		m_vertices[2] = c;
	}
	SD_Scalar supportH(const Vector3r& v) const
	{
		return Math_max(Math_max(v.dot((*this)[0]), v.dot((*this)[1])), v.dot((*this)[2]));
	}
	int maxAxis(const Vector3r& v) const
	{
		return v[0] < v[1] ? (v[1] < v[2] ? 2 : 1) : (v[0] < v[2] ? 2 : 0);
	}

	Vector3r support(const Vector3r& v) const
	{
		Vector3r dots(v.dot((*this)[0]), v.dot((*this)[1]), v.dot((*this)[2]));

		return (*this)[maxAxis(dots)];
	}
    Vector3r operator[](int i) const { return m_vertices[i]; }

private:
	Vector3r	m_vertices[3];
};
////////////////////////////////
/////////Sphere
class DT_Sphere : public DT_Convex {
public:
   DT_Sphere(SD_Scalar radius) : m_radius(radius) {}

	SD_Scalar supportH(const Vector3r& v) const
	{
		return m_radius * v.norm();
	}

	Vector3r support(const Vector3r& v) const
	{
	   SD_Scalar s = v.norm();

		if (s > SD_Scalar(0.0))
		{
			s = m_radius / s;
			return Vector3r(v[0] * s, v[1] * s, v[2] * s);
		}
		else
		{
			return Vector3r(m_radius, SD_Scalar(0.0), SD_Scalar(0.0));
		}
	}
protected:
    SD_Scalar m_radius;
};



//////////////////box
class DT_Box : public DT_Convex {
public:
    DT_Box(SD_Scalar x, SD_Scalar y, SD_Scalar z) :
        m_extent(x, y, z)
	{}

    DT_Box(const Vector3r& e) :
		m_extent(e)
	{}

	SD_Scalar supportH(const Vector3r& v) const
	{
		return v.cwiseAbs().dot(m_extent);
	}

	Vector3r support(const Vector3r& v) const
	{
    return Vector3r(v[0] < 0 ? -m_extent[0] : (v[0] > 0 ? m_extent[0]: 0),
                    v[1] < 0 ? -m_extent[1] : (v[1] > 0 ? m_extent[1]: 0),
                    v[2] < 0 ? -m_extent[2] : (v[2] > 0 ? m_extent[2]: 0));
	}
    const Vector3r& getExtent() const { return m_extent; }

protected:
	typedef unsigned int T_Outcode;

	T_Outcode outcode(const Vector3r& p) const
	{
		return (p[0] < -m_extent[0] ? 0x01 : 0x0) |
			   (p[0] >  m_extent[0] ? 0x02 : 0x0) |
			   (p[1] < -m_extent[1] ? 0x04 : 0x0) |
			   (p[1] >  m_extent[1] ? 0x08 : 0x0) |
			   (p[2] < -m_extent[2] ? 0x10 : 0x0) |
			   (p[2] >  m_extent[2] ? 0x20 : 0x0);
	}

    Vector3r m_extent;
};


//////Hull
class DT_Hull : public DT_Convex {
public:
	DT_Hull(const DT_Convex& lchild, const DT_Convex& rchild) :
		m_lchild(lchild),
		m_rchild(rchild)
	{}

	virtual SD_Scalar supportH(const Vector3r& v) const
	{
		return Math_max(m_lchild.supportH(v), m_rchild.supportH(v));
	}

	virtual Vector3r support(const Vector3r& v) const
	{
		Vector3r lpnt = m_lchild.support(v);
		Vector3r rpnt = m_rchild.support(v);
		return v.dot(lpnt) > v.dot(rpnt) ? lpnt : rpnt;
	}

private:
	const DT_Convex& m_lchild;
	const DT_Convex& m_rchild;
};


////////cone
class DT_Cone : public DT_Convex {
public:
	DT_Cone(SD_Scalar r, SD_Scalar h) :
        bottomRadius(r),
        halfHeight(h * SD_Scalar(0.5)),
        sinAngle(r / sqrt(r * r + h * h))
    {}

	Vector3r support(const Vector3r& v) const
	{
		SD_Scalar v_len = v.norm();//it might be not necessary to get the norm due to a unit vector applied.

		if (v[2] > v_len * sinAngle)//the apex is at the z-axia and the base is at the xoy plane.
		{
			return Vector3r(0.0,0.0, halfHeight);
		}
		else
		{
		    SD_Scalar s = sqrt(v[0] * v[0] + v[1] * v[1]);
		    if (s != SD_Scalar(0.0))
			{
		        SD_Scalar d = bottomRadius / s;
		        return Vector3r(v[0] * d, v[1] * d, -halfHeight);
		    }
		    else
			{
				return Vector3r(0.0, 0.0, -halfHeight);
			}
		}
	}
protected:
    SD_Scalar bottomRadius;
    SD_Scalar halfHeight;
    SD_Scalar sinAngle;
};

////////////////cylinder
class DT_Cylinder : public DT_Convex {
public:
    DT_Cylinder(SD_Scalar r, SD_Scalar h) :
        radius(r),
        halfHeight(h * SD_Scalar(0.5)) {}

	Vector3r support(const Vector3r& v) const
	{
		SD_Scalar s = sqrt(v[0] * v[0] + v[1] * v[1]);
		if (s != SD_Scalar(0.0))
		{
		    SD_Scalar d = radius / s;
		    return Vector3r(v[0] * d, v[1] * d, v[2] < 0.0 ? -halfHeight : halfHeight);
		}
		else
		{
		    return Vector3r(0.0, 0.0, v[2] < 0.0 ? -halfHeight : halfHeight);
		}
	}
protected:
    SD_Scalar radius;
    SD_Scalar halfHeight;
};
//////////////////Minkowski
class DT_Minkowski : public DT_Convex {
public:
	DT_Minkowski(const DT_Convex& lchild, const DT_Convex& rchild)
     : m_lchild(lchild),
       m_rchild(rchild)
   {}

	virtual SD_Scalar supportH(const Vector3r& v) const
	{
		return m_lchild.supportH(v) + m_rchild.supportH(v);
	}

	virtual Vector3r support(const Vector3r& v) const
	{
		return m_lchild.support(v) + (Vector3r)m_rchild.support(v);
	}

private:
	const DT_Convex& m_lchild;
	const DT_Convex& m_rchild;
};
///////////Transform

class DT_Transform : public DT_Convex {
public:
	DT_Transform(const MT_Transform& xform, const DT_Convex& child) :
		m_xform(xform),
		m_child(child)
	{}

	virtual SD_Scalar supportH(const Vector3r& v) const
	{
		return m_child.supportH(m_xform.getBasis().transpose()*v) +
			   v.dot(m_xform.getOrigin());
	}

	virtual Vector3r support(const Vector3r& v) const
	{
		return m_xform(m_child.support(m_xform.getBasis().transpose()*v));
	}

private:
	const MT_Transform& m_xform;
	const DT_Convex&    m_child;
};
///////////////
class DT_TransM : public DT_Convex {
public:
	DT_TransM(const Matrix3r& TM, const Vector3r& origin, const DT_Convex& child) :
		m_TM(TM),
		m_origin(origin),
		m_child(child)
	{}

	virtual SD_Scalar supportH(const Vector3r& v) const
	{
		return m_child.supportH(m_TM.transpose()*v) +
			   v.dot(m_origin);
	}

	virtual Vector3r support(const Vector3r& v) const
	{
		Vector3r x = m_child.support(m_TM.transpose()*v);
		return Vector3r(m_TM.row(0).dot(x) + m_origin[0],
			   m_TM.row(1).dot(x) + m_origin[1],
			   m_TM.row(2).dot(x) + m_origin[2]);
	}

private:
	//const MT_Transform& m_xform;
	const DT_Convex&    m_child;
	const Matrix3r&	m_TM;
	const Vector3r&		m_origin;
};
class DT_TransM2 : public DT_Convex {
public:
	DT_TransM2(const Matrix3r& TM, const Vector3r& origin, const shared_ptr<DT_Convex>& child) :
		m_TM(TM),
		m_origin(origin),
		m_child(child)
	{}

	virtual SD_Scalar supportH(const Vector3r& v) const
	{
		return m_child->supportH(m_TM.transpose()*v) +
			   v.dot(m_origin);
	}

	virtual Vector3r support(const Vector3r& v) const
	{
		Vector3r x = m_child->support(m_TM.transpose()*v);
		return Vector3r(m_TM.row(0).dot(x) + m_origin[0],
			   m_TM.row(1).dot(x) + m_origin[1],
			   m_TM.row(2).dot(x) + m_origin[2]);
	}

private:
	//const MT_Transform& m_xform;
	const shared_ptr<DT_Convex>    m_child;
	const Matrix3r&	m_TM;
	const Vector3r&		m_origin;
};
/////////////////Polyhedron
#ifdef HAVE_CONFIG_H
# include "config.h"
# if HAVE_QHULL_QHULL_A_H
#  define QHULL
# endif
#endif

#define QHULL
typedef std::vector<unsigned int> T_IndexBuf;
typedef std::vector<T_IndexBuf> T_MultiIndexBuf;

static char options[] = "qhull Qts i Tv";

#define DK_HIERARCHY
#ifdef QHULL
extern "C"{
#include <sudodem/lib/qhull/qhull_a.h>
}

class DT_Polyhedron : public DT_Convex {
	typedef DT_Array<DT_IndexArray> T_MultiIndexArray;
public:
	/*DT_Polyhedron()
		: m_verts(0),
		  m_cobound(0)
	{}*/

	DT_Polyhedron(std::vector<Vector3r> vertices);
	virtual ~DT_Polyhedron();

    virtual SD_Scalar supportH(const Vector3r& v) const;
    virtual Vector3r support(const Vector3r& v) const;

	const Vector3r& operator[](int i) const { return m_verts[i]; }
    unsigned int numVerts() const { return m_count; }
    void adjacency(unsigned int, const std::vector<Vector3r>& , const char *);
	T_MultiIndexBuf adjacency_graph(unsigned int, const std::vector<Vector3r >&, const char *);
	T_MultiIndexBuf simplex_adjacency_graph(unsigned int, const char *);

    T_MultiIndexBuf       m_facetIndexBuf;//vertex ids of each facet
    std::vector<Vector3r > m_vertices2;
private:
	unsigned int              m_count;
	//Vector3r			 *m_verts;
  std::vector<Vector3r > m_verts;
	T_MultiIndexArray    *m_cobound;

  unsigned int            m_start_vertex;
	mutable unsigned int      m_curr_vertex;
};



#else //define Polytope


	class DT_Polytope : public DT_Convex {
	public:
		DT_Polytope() {}
		/*DT_Polytope(const DT_VertexBase *base, unsigned int count, const unsigned int *indices)
		  : m_base(base),
			m_index(count, indices)
		{}*/
		DT_Polytope(std::vector<Vector3r> vertices): m_vertices(vertices)
		{}

		virtual SD_Scalar supportH(const Vector3r& v) const;
		virtual Vector3r support(const Vector3r& v) const;

		Vector3r operator[](int i) const { assert(0 <= i && i < m_vertices.size()); return m_vertices[i]; }
		unsigned int numVerts() const { return m_index.size(); }

	protected:
		//const DT_VertexBase *m_base;
		//DT_IndexArray        m_index;
		std::vector<Vector3r> m_vertices;
	};

	SD_Scalar DT_Polytope::supportH(const Vector3r& v) const
	{
		SD_Scalar h = (*this)[0].dot(v), d;
		unsigned int i;
		for (i = 1; i < numVerts(); ++i)
		{
		    if ((d = (*this)[i].dot(v)) > h)
			{
				h = d;
			}
		}
		return h;
	}

	Vector3r DT_Polytope::support(const Vector3r& v) const
	{
		int c = 0;
		SD_Scalar h = (*this)[0].dot(v), d;
		unsigned int i;
		for (i = 1; i < numVerts(); ++i)
		{
		    if ((d = (*this)[i].dot(v)) > h)
			{
				c = i;
				h = d;
			}
		}
		return (*this)[c];
	}

typedef DT_Polytope DT_Polyhedron;

#endif
////////////////////
#endif
