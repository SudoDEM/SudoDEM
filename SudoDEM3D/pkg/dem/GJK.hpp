#ifndef GJK_H
#define GJK_H
#pragma once
#include "GJKParticle_shapes.h"

//---------------------------------------------------------------------------
//-----------------GJK algorithm-----------------------------
//#define USE_BACKUP_PROCEDURE
//#define JOHNSON_ROBUST
#define FAST_CLOSEST
#define UNROLL_LOOPS
#define SAFE_EXIT

//static const SD_Scalar rel_error = SD_Scalar(1.0e-3);

//SD_Scalar DT_Accuracy::rel_error2 = rel_error * rel_error;
//SD_Scalar DT_Accuracy::depth_tolerance = SD_Scalar(1.0) + SD_Scalar(2.0) * rel_error;
//SD_Scalar DT_Accuracy::tol_error = Mathr::EPSILON;

////////////////////

class DT_GJK {
private:
	inline static bool subseteq(unsigned int a, unsigned int b) { return (a & b) == a; }//a is a subset of b
	inline static bool contains(unsigned int a, unsigned int b) { return (a & b) != 0x0; }//a contains b (b has only one non-empty bit)

public:
	DT_GJK() :
		m_bits(0x0),//m_bits identifies the simplex
		m_all_bits(0x0)
	{}

	bool emptySimplex() const { return m_bits == 0x0; }
	bool fullSimplex() const { return m_bits == 0xf; }

	void reset()
	{
		m_bits = 0x0;
		m_all_bits = 0x0;
	}

	bool inSimplex(const Vector3r& w) const//w is in the simplex
	{ int i;
		unsigned int bit;
		for (i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
            if (contains(m_all_bits, bit) && w == m_y[i])//m_all_bits has a bit equal to bit.
			{
				return true;
			}
		}
		return false;
	}

    bool isAffinelyDependent() const
	{
        SD_Scalar sum(0.0);
		//loop the vertices of the simplex
		int i;
		unsigned int bit;
		for ( i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
      if (contains(m_all_bits, bit))//the current point is a vertex of the simplex
      {
          sum += m_det[m_all_bits][i];
			}
		}
		return sum <= SD_Scalar(0.0);
	}

	//add a vertex to the simplex (convex hull)
	void addVertex(const Vector3r& w)
	{
		assert(!fullSimplex());
		m_last = 0;
    m_last_bit = 0x1;
    while (contains(m_bits, m_last_bit))//place w at the new vertex of the simplex
		{
			++m_last;
			m_last_bit <<= 1;
		}
		m_y[m_last] = w;
		m_ylen2[m_last] = w.squaredNorm();
    m_all_bits = m_bits | m_last_bit;//m_all_bits contains the whole hull

		update_cache();
		compute_det();
	}

	void addVertex(const Vector3r& w, const Vector3r& p, const Vector3r& q)
	{
		addVertex(w);
		m_p[m_last] = p;//supporting points p and q
		m_q[m_last] = q;
	}
	//get the vertices of the simplex and the corresponding supporing points p and q.
	int getSimplex(Vector3r *pBuf, Vector3r *qBuf, Vector3r *yBuf) const
	{
		int num_verts = 0;
		//loop all vertices of the simplex
		int i;
		unsigned int bit;
		for (i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
			if (contains(m_bits, bit))//the point is a vertex of the simplex
			{
				pBuf[num_verts] = m_p[i];
				qBuf[num_verts] = m_q[i];
				yBuf[num_verts] = m_y[i];

				#ifdef DEBUG
								std::cout << "Point " << i << " = " << m_y[i] << std::endl;
				#endif

				++num_verts;
			}
		}
		return num_verts;//return simplex dimension
    }

	void compute_points(Vector3r& p1, Vector3r& p2)
	{
		SD_Scalar sum = SD_Scalar(0.0);
		p1 = Vector3r::Zero();
		p2 = Vector3r::Zero();
		//loop the simplex
		int i;
		unsigned int bit;
		for ( i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
			if (contains(m_bits, bit))
			{
				sum += m_det[m_bits][i];
				p1 += m_p[i] * m_det[m_bits][i];
				p2 += m_q[i] * m_det[m_bits][i];
			}
		}

		assert(sum > SD_Scalar(0.0));//here may be failed , Sway when using 1000 ellipsoids
		SD_Scalar s = SD_Scalar(1.0) / sum;
		p1 *= s;
		p2 *= s;
	}

	bool closest(Vector3r& v)
	{
#ifdef FAST_CLOSEST
		unsigned int s;
		for (s = m_bits; s != 0x0; --s)
		{
			if (subseteq(s, m_bits) && valid(s | m_last_bit))
			{
				m_bits = s | m_last_bit;
				v = compute_vector(m_bits);
				return true;
			}
		}
		if (valid(m_last_bit))
		{
			m_bits = m_last_bit;
			m_maxlen2 = m_ylen2[m_last];
			v = m_y[m_last];
			return true;
		}
#else
		unsigned int s;
		for (s = m_all_bits; s != 0x0; --s)
		{
			if (subseteq(s, m_all_bits) && valid(s))
			{
				m_bits = s;
				v = compute_vector(m_bits);
				return true;
			}
		}
#endif

		// Original GJK calls the backup procedure at this point.
#ifdef USE_BACKUP_PROCEDURE
        backup_closest(v);
#endif
		return false;
	}

	void backup_closest(Vector3r& v)
	{
		SD_Scalar min_dist2 = Mathr::MAX_REAL;

      unsigned int s;
		for (s = m_all_bits; s != 0x0; --s)
		{
			if (subseteq(s, m_all_bits) && proper(s))
			{
				Vector3r u = compute_vector(s);
				SD_Scalar dist2 = u.squaredNorm();
				if (dist2 < min_dist2)
				{
					min_dist2 = dist2;
					m_bits = s;
					v = u;
				}
			}
		}
	}

	SD_Scalar maxVertex() const { return m_maxlen2; }


private:
	//void update_cache();
	void compute_det();
	inline void update_cache()
	{
		int i;
		unsigned int bit;
	  for (i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
	    if (contains(m_bits, bit))
			{
				m_edge[i][m_last] = m_y[i] - m_y[m_last];
				m_edge[m_last][i] = -m_edge[i][m_last];

	#ifdef JOHNSON_ROBUST
				m_norm[i][m_last] = m_norm[m_last][i] = m_edge[i][m_last].squaredNorm();
	#endif

			}
		}
	}

	bool valid(unsigned int s)
	{
		int i;
		unsigned int bit;
		for (i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
			if (contains(m_all_bits, bit))
			{
				if (contains(s, bit))
				{
					if (m_det[s][i] <= SD_Scalar(0.0))
					{
						return false;
					}
				}
				else if (m_det[s | bit][i] > SD_Scalar(0.0))
				{
					return false;
				}
			}
		}
		return true;
	}

	bool proper(unsigned int s)
	{
		int i;
		unsigned int bit;
		for (i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
			if (contains(s, bit) && m_det[s][i] <= SD_Scalar(0.0))
			{
				return false;
			}
		}
		return true;
	}

	Vector3r compute_vector(unsigned int s)
	{
    Vector3r v(0.0,0.0,0.0);
		m_maxlen2 = SD_Scalar(0.0);

    SD_Scalar sum = SD_Scalar(0.0);
	  //loop the simplex
		int i;
		unsigned int bit;
		for ( i = 0, bit = 0x1; i < 4; ++i, bit <<= 1)
		{
			if (contains(s, bit))
			{
				sum += m_det[s][i];
				//Math_set_max(m_maxlen2, m_ylen2[i]);
				if(m_maxlen2 < m_ylen2[i])
				{
					m_maxlen2 = m_ylen2[i];
				}
				v += m_y[i] * m_det[s][i];
			}
		}

		assert(sum > SD_Scalar(0.0));

		return v / sum;
	}

private:
	SD_Scalar	m_det[16][4]; // cached sub-determinants
  Vector3r	m_edge[4][4];

  #ifdef JOHNSON_ROBUST
    SD_Scalar	m_norm[4][4];
  #endif

	Vector3r	m_p[4];    // support points of object A in local coordinates
	Vector3r	m_q[4];    // support points of object B in local coordinates
	Vector3r	m_y[4];   // support points of A - B in world coordinates
	SD_Scalar	m_ylen2[4];   // Squared lengths support points y

	SD_Scalar	m_maxlen2; // Maximum squared length to a vertex of the current
	                      // simplex

	unsigned int		m_bits;      // identifies current simplex
	unsigned int		m_last;      // identifies last found support point
	unsigned int		m_last_bit;  // m_last_bit == 0x1 << last
	unsigned int		m_all_bits;  // m_all_bits == m_bits  | m_last_bit
};


#ifdef JOHNSON_ROBUST

inline void DT_GJK::compute_det()
{
    m_det[m_last_bit][m_last] = SD_Scalar(1.0);

  if (m_bits != 0x0)//not null
  {
	int i;
	unsigned int si;
  for (i = 0, si = 0x1; i < 4; ++i, si <<= 1)//traverse
	{
        if (contains(m_bits, si))
		{
            unsigned int s2 = si | m_last_bit;
            m_det[s2][i] = m_edge[m_last][i].dot(m_y[m_last]);
            m_det[s2][m_last] = m_edge[i][m_last].dot(m_y[i]);

			int j;
			unsigned int sj;
            for (j = 0, sj = 0x1; j < i; ++j, sj <<= 1)
			{
                if (contains(m_bits, sj))
				{
					int k;
                    unsigned int s3 = sj | s2;

					k = m_norm[i][j] < m_norm[m_last][j] ? i : m_last;
                    m_det[s3][j] = m_det[s2][i] * m_edge[k][j].dot(m_y[i]) +
                                   m_det[s2][m_last] * m_edge[k][j].dot(m_y[m_last]);
					k = m_norm[j][i] < m_norm[m_last][i] ? j : m_last;
                    m_det[s3][i] = m_det[sj|m_last_bit][j] * m_edge[k][i].dot(m_y[j]) +
                                   m_det[sj|m_last_bit][m_last] * m_edge[k][i].dot(m_y[m_last]);
					k = m_norm[i][m_last] < m_norm[j][m_last] ? i : j;
                    m_det[s3][m_last] = m_det[sj|si][j] * m_edge[k][m_last].dot(m_y[j]) +
                                        m_det[sj|si][i] * m_edge[k][m_last].dot(m_y[i]);
                }
            }
        }
    }

    if (m_all_bits == 0xf)
	{
		int k;

		k = m_norm[1][0] < m_norm[2][0] ? (m_norm[1][0] < m_norm[3][0] ? 1 : 3) : (m_norm[2][0] < m_norm[3][0] ? 2 : 3);

        m_det[0xf][0] = m_det[0xe][1] * m_edge[k][0].dot(m_y[1]) +
                        m_det[0xe][2] * m_edge[k][0].dot(m_y[2]) +
                        m_det[0xe][3] * m_edge[k][0].dot(m_y[3]);

		k = m_norm[0][1] < m_norm[2][1] ? (m_norm[0][1] < m_norm[3][1] ? 0 : 3) : (m_norm[2][1] < m_norm[3][1] ? 2 : 3);

        m_det[0xf][1] = m_det[0xd][0] * m_edge[k][1].dot(m_y[0]) +
                        m_det[0xd][2] * m_edge[k][1].dot(m_y[2]) +
                        m_det[0xd][3] * m_edge[k][1].dot(m_y[3]);

		k = m_norm[0][2] < m_norm[1][2] ? (m_norm[0][2] < m_norm[3][2] ? 0 : 3) : (m_norm[1][2] < m_norm[3][2] ? 1 : 3);

        m_det[0xf][2] = m_det[0xb][0] * m_edge[k][2].dot(m_y[0]) +
                        m_det[0xb][1] * m_edge[k][2].dot(m_y[1]) +
                        m_det[0xb][3] * m_edge[k][2].dot(m_y[3]);

		k = m_norm[0][3] < m_norm[1][3] ? (m_norm[0][3] < m_norm[2][3] ? 0 : 2) : (m_norm[1][3] < m_norm[2][3] ? 1 : 2);

        m_det[0xf][3] = m_det[0x7][0] * m_edge[k][3].dot(m_y[0]) +
                        m_det[0x7][1] * m_edge[k][3].dot(m_y[1]) +
                        m_det[0x7][2] * m_edge[k][3].dot(m_y[2]);
    }
  }
}

#else

inline void DT_GJK::compute_det()
{
    m_det[m_last_bit][m_last] = SD_Scalar(1.0);

    if (m_bits != 0x0)
    {

#ifdef UNROLL_LOOPS

    if (contains(m_bits, 0x1))
    {
        unsigned int s2 = 0x1 | m_last_bit;
        m_det[s2][0] = m_edge[m_last][0].dot(m_y[m_last]);
        m_det[s2][m_last] = m_edge[0][m_last].dot(m_y[0]);
    }

    if (contains(m_bits, 0x2))
    {
        unsigned int s2 = 0x2 | m_last_bit;
        m_det[s2][1] = m_edge[m_last][1].dot(m_y[m_last]);
        m_det[s2][m_last] = m_edge[1][m_last].dot(m_y[1]);

        if (contains(m_bits, 0x1))
        {
            unsigned int s3 = 0x1 | s2;
            unsigned int t2 = 0x1 | m_last_bit;
            m_det[s3][0] = m_det[s2][1] * m_edge[1][0].dot(m_y[1]) + m_det[s2][m_last] * m_edge[1][0].dot(m_y[m_last]);
            m_det[s3][1] = m_det[t2][0] * m_edge[0][1].dot(m_y[0]) + m_det[t2][m_last] * m_edge[0][1].dot(m_y[m_last]);
            m_det[s3][m_last] = m_det[0x3][0] * m_edge[0][m_last].dot(m_y[0]) + m_det[0x3][1] * m_edge[0][m_last].dot(m_y[1]);
        }
    }

    if (contains(m_bits, 0x4))
    {
        unsigned int s2 = 0x4 | m_last_bit;
        m_det[s2][2] = m_edge[m_last][2].dot(m_y[m_last]);
        m_det[s2][m_last] = m_edge[2][m_last].dot(m_y[2]);

        if (contains(m_bits, 0x1))
        {
            unsigned int s3 = 0x1 | s2;
            unsigned int t2 = 0x1 | m_last_bit;
            m_det[s3][0] = m_det[s2][2] * m_edge[2][0].dot(m_y[2]) + m_det[s2][m_last] * m_edge[2][0].dot(m_y[m_last]);
            m_det[s3][2] = m_det[t2][0] * m_edge[0][2].dot(m_y[0]) + m_det[t2][m_last] * m_edge[0][2].dot(m_y[m_last]);
            m_det[s3][m_last] = m_det[0x5][0] * m_edge[0][m_last].dot(m_y[0]) + m_det[0x5][2] * m_edge[0][m_last].dot(m_y[2]);
        }

        if (contains(m_bits, 0x2))
        {
            unsigned int s3 = 0x2 | s2;
            unsigned int t2 = 0x2 | m_last_bit;
            m_det[s3][1] = m_det[s2][2] * m_edge[2][1].dot(m_y[2]) + m_det[s2][m_last] * m_edge[2][1].dot(m_y[m_last]);
            m_det[s3][2] = m_det[t2][1] * m_edge[1][2].dot(m_y[1]) + m_det[t2][m_last] * m_edge[1][2].dot(m_y[m_last]);
            m_det[s3][m_last] = m_det[0x6][1] * m_edge[1][m_last].dot(m_y[1]) + m_det[0x6][2] * m_edge[1][m_last].dot(m_y[2]);
        }
    }

    if (contains(m_bits, 0x8))
    {
        unsigned int s2 = 0x8 | m_last_bit;
        m_det[s2][3] = m_edge[m_last][3].dot(m_y[m_last]);
        m_det[s2][m_last] = m_edge[3][m_last].dot(m_y[3]);

        if (contains(m_bits, 0x1))
        {
            unsigned int s3 = 0x1 | s2;
            unsigned int t2 = 0x1 | m_last_bit;
            m_det[s3][0] = m_det[s2][3] * m_edge[3][0].dot(m_y[3]) + m_det[s2][m_last] * m_edge[3][0].dot(m_y[m_last]);
            m_det[s3][3] = m_det[t2][0] * m_edge[0][3].dot(m_y[0]) + m_det[t2][m_last] * m_edge[0][3].dot(m_y[m_last]);
            m_det[s3][m_last] = m_det[0x9][0] * m_edge[0][m_last].dot(m_y[0]) + m_det[0x9][3] * m_edge[0][m_last].dot(m_y[3]);
        }

        if (contains(m_bits, 0x2))
        {
            unsigned int s3 = 0x2 | s2;
            unsigned int t2 = 0x2 | m_last_bit;
            m_det[s3][1] = m_det[s2][3] * m_edge[3][1].dot(m_y[3]) + m_det[s2][m_last] * m_edge[3][1].dot(m_y[m_last]);
            m_det[s3][3] = m_det[t2][1] * m_edge[1][3].dot(m_y[1]) + m_det[t2][m_last] * m_edge[1][3].dot(m_y[m_last]);
            m_det[s3][m_last] = m_det[0xa][1] * m_edge[1][m_last].dot(m_y[1]) + m_det[0xa][3] * m_edge[1][m_last].dot(m_y[3]);
        }

        if (contains(m_bits, 0x4))
        {
            unsigned int s3 = 0x4 | s2;
            unsigned int t2 = 0x4 | m_last_bit;
            m_det[s3][2] = m_det[s2][3] * m_edge[3][2].dot(m_y[3]) + m_det[s2][m_last] * m_edge[3][2].dot(m_y[m_last]);
            m_det[s3][3] = m_det[t2][2] * m_edge[2][3].dot(m_y[2]) + m_det[t2][m_last] * m_edge[2][3].dot(m_y[m_last]);
            m_det[s3][m_last] = m_det[0xc][2] * m_edge[2][m_last].dot(m_y[2]) + m_det[0xc][3] * m_edge[2][m_last].dot(m_y[3]);
        }

    }


#else

	int i;
	unsigned int si;
    for (i = 0, si = 0x1; i < 4; ++i, si <<= 1)
	{
        if (contains(m_bits, si))
		{
            unsigned int s2 = si | m_last_bit;
            m_det[s2][i] = m_edge[m_last][i].dot(m_y[m_last]);
            m_det[s2][m_last] = m_edge[i][m_last].dot(m_y[i]);

			int j;
			unsigned int sj;
            for (j = 0, sj = 0x1; j < i; ++j, sj <<= 1)
			{
                if (contains(m_bits, sj))
				{
                    unsigned int s3 = sj | s2;
                    m_det[s3][j] = m_det[s2][i] * m_edge[i][j].dot(m_y[i]) +
                                   m_det[s2][m_last] * m_edge[i][j].dot(m_y[m_last]);
                    m_det[s3][i] = m_det[sj|m_last_bit][j] * m_edge[j][i].dot(m_y[j]) +
                                   m_det[sj|m_last_bit][m_last] * m_edge[j][i].dot(m_y[m_last]);
                    m_det[s3][m_last] = m_det[sj|si][j] * m_edge[j][m_last].dot(m_y[j]) +
                                        m_det[sj|si][i] * m_edge[j][m_last].dot(m_y[i]);
                }
            }
        }
    }

#endif

    if (m_all_bits == 0xf)
	{
        m_det[0xf][0] = m_det[0xe][1] * m_edge[1][0].dot(m_y[1]) +
                        m_det[0xe][2] * m_edge[1][0].dot(m_y[2]) +
                        m_det[0xe][3] * m_edge[1][0].dot(m_y[3]);
        m_det[0xf][1] = m_det[0xd][0] * m_edge[0][1].dot(m_y[0]) +
                        m_det[0xd][2] * m_edge[0][1].dot(m_y[2]) +
                        m_det[0xd][3] * m_edge[0][1].dot(m_y[3]);
        m_det[0xf][2] = m_det[0xb][0] * m_edge[0][2].dot(m_y[0]) +
                        m_det[0xb][1] * m_edge[0][2].dot(m_y[1]) +
                        m_det[0xb][3] * m_edge[0][2].dot(m_y[3]);
        m_det[0xf][3] = m_det[0x7][0] * m_edge[0][3].dot(m_y[0]) +
                        m_det[0x7][1] * m_edge[0][3].dot(m_y[1]) +
                        m_det[0x7][2] * m_edge[0][3].dot(m_y[2]);
    }
}
}

#endif
//////////////////////
////////////////////////pendepth
//const int       MaxSupportPoints = 100;
//const int       MaxFacets         = 200;

//static Vector3r  pBuf[MaxSupportPoints];
//static Vector3r  qBuf[MaxSupportPoints];
//static Vector3r yBuf[MaxSupportPoints];


//static Triangle *triangleHeap[MaxFacets];
//static int  num_triangles;


#endif
