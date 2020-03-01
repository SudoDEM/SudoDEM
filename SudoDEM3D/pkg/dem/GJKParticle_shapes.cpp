//#pragma once

#include <vector>
#include <algorithm>
#include <new>

#include "GJKParticle_shapes.h"

///////////////TriEdge

//TriangleStore g_triangleStore;
bool link(const Edge& edge0, const Edge& edge1)
{
	bool ok = edge0.getEnd1() == edge1.getEnd2() && edge0.getEnd2() == edge1.getEnd1();

    //assert(ok);

    if (ok)
    {
	edge0.triangle()->m_adjEdges[edge0.index()] = edge1;
    edge1.triangle()->m_adjEdges[edge1.index()] = edge0;
	}

    return ok;
}

void half_link(const Edge& edge0, const Edge& edge1)
{
	assert(edge0.getEnd1() == edge1.getEnd2() && edge0.getEnd2() == edge1.getEnd1());

	edge0.triangle()->m_adjEdges[edge0.index()] = edge1;
}


bool Edge::silhouette(const Vector3r *verts, Index_t index, TriangleStore& triangleStore) const
{
    if (!m_triangle->isObsolete())
    {
		if (!m_triangle->isVisibleFrom(verts, index))
        {
			Triangle *triangle = triangleStore.newTriangle(verts, index, getEnd2(), getEnd1());

			if (triangle)
			{
				half_link(Edge(triangle, 1), *this);
				return true;
			}

			return false;
		}
        else
        {
            m_triangle->setObsolete(true); // Triangle is visible

			int backup = triangleStore.getFree();

            if (!m_triangle->getAdjEdge(circ_next(m_index)).silhouette(verts, index, triangleStore))
			{
				m_triangle->setObsolete(false);

				Triangle *triangle = triangleStore.newTriangle(verts, index, getEnd2(), getEnd1());

				if (triangle)
				{
					half_link(Edge(triangle, 1), *this);
					return true;
				}

				return false;
			}
			else if (!m_triangle->getAdjEdge(circ_prev(m_index)).silhouette(verts, index, triangleStore))
			{
				m_triangle->setObsolete(false);

				triangleStore.setFree(backup);

				Triangle *triangle = triangleStore.newTriangle(verts, index, getEnd2(), getEnd1());

				if (triangle)
				{
					half_link(Edge(triangle, 1), *this);
					return true;
				}

				return false;
			}
        }
    }

	return true;
}


bool Triangle::computeClosest(const Vector3r *verts)
{
    const Vector3r& p0 = verts[m_indices[0]];//vertice 1
		//vectors v1 and v2
    Vector3r v1 = verts[m_indices[1]] - p0;
    Vector3r v2 = verts[m_indices[2]] - p0;
    SD_Scalar v1dv1 = v1.squaredNorm();
    SD_Scalar v1dv2 = v1.dot(v2);
    SD_Scalar v2dv2 = v2.squaredNorm();
    SD_Scalar p0dv1 = p0.dot(v1);
    SD_Scalar p0dv2 = p0.dot(v2);

    m_det = v1dv1 * v2dv2 - v1dv2 * v1dv2; // non-negative
    m_lambda1 = p0dv2 * v1dv2 - p0dv1 * v2dv2;
    m_lambda2 = p0dv1 * v1dv2 - p0dv2 * v1dv1;

    if (m_det > SD_Scalar(0.0))
	{
		m_closest = p0 + (m_lambda1 * v1 + m_lambda2 * v2) / m_det;
		m_dist2 = m_closest.squaredNorm();

		return true;
    }

    return false;
}


bool Triangle::silhouette(const Vector3r *verts, Index_t index, TriangleStore& triangleStore)
{
	//assert(isVisibleFrom(verts, index));

	int first = triangleStore.getFree();

	setObsolete(true);

	bool result = m_adjEdges[0].silhouette(verts, index, triangleStore) &&
	              m_adjEdges[1].silhouette(verts, index, triangleStore) &&
	              m_adjEdges[2].silhouette(verts, index, triangleStore);

	if (result)
	{
		int i, j;
		for (i = first, j = triangleStore.getFree()-1; i != triangleStore.getFree(); j = i++)
		{
			Triangle *triangle = &triangleStore[i];
			half_link(triangle->getAdjEdge(1), Edge(triangle, 1));
            if (!link(Edge(triangle, 0), Edge(&triangleStore[j], 2)))
            {
                return false;
            }
		}
	}

	return result;
}
///////////////////////////
void DT_Polyhedron::adjacency(unsigned int count, const std::vector<Vector3r>& verts, const char *flags)
{
	int curlong, totlong, exitcode;

    facetT *facet;
    vertexT *vertex;
    vertexT **vertexp;

    std::vector<coordT> array;
	T_IndexBuf index;
    unsigned int i;
    Vector3r tmp;
	//coordT tmp1 [3] = {0,0,0};
    for (i = 0; i != count; ++i)
	{
		if (flags == 0 || flags[i])
		{   tmp = verts[i];
			//tmp1[0] = tmp[0];
			//tmp1[1] = tmp[1];
			//tmp1[2] = tmp[2];
            //array.push_back(&tmp1[0]);
			array.push_back(verts[i][0]);
			array.push_back(verts[i][1]);
			array.push_back(verts[i][2]);
            m_vertices2.push_back(verts[i]);
			index.push_back(i);
		}
    }
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

    //T_IndexBuf *indexBuf = new T_IndexBuf[count];//array of vectors
    //T_MultiIndexBuf facetIndexBuf;
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

T_MultiIndexBuf DT_Polyhedron::adjacency_graph(unsigned int count, const std::vector<Vector3r>& verts, const char *flags)
{
	int curlong, totlong, exitcode;

    facetT *facet;
    vertexT *vertex;
    vertexT **vertexp;

    std::vector<coordT> array;
	T_IndexBuf index;
    unsigned int i;
	//coordT tmp1[3];
    for (i = 0; i != count; ++i)
	{
		if (flags == 0 || flags[i])
		{

			/*tmp1[0] = verts[i][0];
			tmp1[1] = verts[i][1];
			tmp1[2] = verts[i][2];
            array.push_back(&tmp1[0]);*/
			array.push_back(verts[i][0]);
			array.push_back(verts[i][1]);
			array.push_back(verts[i][2]);
			index.push_back(i);//index for point id.
		}
    }
    //m_vertices2 = array;
	double* a = &array[0];
    qh_init_A(stdin, stdout, stderr, 0, NULL);
    if ((exitcode = setjmp(qh errexit)))
	{
		exit(exitcode);
	}
    qh_initflags(options);
    qh_init_B(a, array.size()/3, 3, False);
    qh_qhull();//quick hull calculation
    qh_check_output();

    //T_IndexBuf *indexBuf = new T_IndexBuf[count];//array of vectors
		T_MultiIndexBuf indexBuf(count);
    //T_MultiIndexBuf facetIndexBuf;
    FORALLfacets//loop all facets
	{
		setT *vertices = qh_facet3vertex(facet);//get vertices of the current facet

		T_IndexBuf  facetIndices;

		FOREACHvertex_(vertices)//taverse the vertices
		{
			facetIndices.push_back(index[qh_pointid(vertex->point)]);//push back point ids of the vertices of the current facet
            //std::cout<<"out point id"<<qh_pointid(vertex->point)<<"corrd="<<vertex->point[0]<<" "<<vertex->point[1]<<" "<<vertex->point[2]<<std::endl;
			//index[qh_pointid(vertex->point)] is equal to qh_pointid(vertex->point).
			//std::cout<<"debug:point_id"<<qh_pointid(vertex->point)<<" "<<index[qh_pointid(vertex->point)]<<std::endl;
		}
		//m_facetIndexBuf.push_back(facetIndices);
		int i, j;
		for (i = 0, j = facetIndices.size()-1; i < (int)facetIndices.size(); j = i++)
		{
			indexBuf[facetIndices[j]].push_back(facetIndices[i]);//adding the neighboring vertices of a given vertex to a group (vector). The group id in indexBuf is equal to the vetex id.
		}
    }

    //indexBuf consists of connection information of each vertex with others, e.g., indexBuf[i] is a vector of ids of vertices connected with the vertex of id i.

    qh NOerrexit = True;
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

	return indexBuf;
}



T_MultiIndexBuf DT_Polyhedron::simplex_adjacency_graph(unsigned int count, const char *flags)
{
	//T_IndexBuf *indexBuf = new T_IndexBuf[count];
  T_MultiIndexBuf indexBuf(count);
	unsigned int index[4];

	unsigned int k = 0;
	unsigned int i;
	for (i = 0; i != count; ++i)
	{
		if (flags == 0 || flags[i])
		{
			index[k++] = i;
		}
	}

	assert(k <= 4);

	for (i = 0; i != k; ++i)
	{
		unsigned int j;
		for (j = 0; j != k; ++j)
		{
			if (i != j)
			{
				indexBuf[index[i]].push_back(index[j]);
			}
		}
	}

	return indexBuf;
}

#ifdef DK_HIERARCHY
//prune the duplicate vertices on the lower (finer) layers, yielding a slimmed-down multi-layer vertex adjacent graph.
void prune(unsigned int count, T_MultiIndexBuf *cobound)
{
	for (unsigned int i = 0; i != count; ++i)//the first dimension: traverse all vertices
	{
		assert(cobound[i].size());
		//traverse each layer of the ith vertex
		for (unsigned int j = 0; j != cobound[i].size() - 1; ++j)
		{
			T_IndexBuf::iterator it = cobound[i][j].begin();
			while (it != cobound[i][j].end())//traverse an adjancy graph
			{
				T_IndexBuf::iterator jt =
					std::find(cobound[i][j+1].begin(), cobound[i][j+1].end(), *it);

				if (jt != cobound[i][j+1].end())//find a *it, i.e., a duplicate vertex exsits on the upper layer.
				{		//delete the duplicate vertex on the lower layer (current cobound)
            if (it != cobound[i][j].end() - 1)//'it' is not the last element
            {
					    std::swap(*it, cobound[i][j].back());
					    cobound[i][j].pop_back();
            }
            else
            {
                cobound[i][j].pop_back();
                break;
            }
				}
				else
				{
					++it;
				}
			}
		}
	}
}

#endif

DT_Polyhedron::DT_Polyhedron(std::vector<Vector3r> vertices)
{	unsigned int count = vertices.size();
	assert(count);

	std::vector<Vector3r> vertexBuf(vertices);
	//unsigned int i;
	/*
	for (unsigned int i = 0; i < count; i++)
	{	vertexBuf.push_back(vertices[i]);//copy points to vertexBuf
	}*/
    //get neighboring connection information
    //adjacency(count, &vertexBuf[0], 0);

	//T_IndexBuf *indexBuf = count > 4 ? adjacency_graph(count, &vertexBuf[0], 0) : simplex_adjacency_graph(count, 0);
	T_MultiIndexBuf indexBuf = count > 4 ? adjacency_graph(count, vertexBuf, 0) : simplex_adjacency_graph(count, 0);

	std::vector<Vector3r> pointBuf;

	for (unsigned int i = 0; i < count; i++)
	{
		if (!indexBuf[i].empty())//sweep the ineffective points those are enclosed in the HULL.
		{
			pointBuf.push_back(vertexBuf[i]);
		}
	}
  indexBuf.clear();
	//delete [] indexBuf;
    //pointBuf consists of all vertices of the hull
	m_count = pointBuf.size();
	m_verts = pointBuf;
	//m_verts = new Vector3r[m_count];
	//std::copy(pointBuf.begin(), pointBuf.end(), &m_verts[0]);

	T_MultiIndexBuf *cobound = new T_MultiIndexBuf[m_count];
    char *flags = new char[m_count];
		//flag is used to denote an adjacent vertex of P is whether or not an adjacent vertex of P on a coarse layer.
	std::fill(&flags[0], &flags[m_count], 1);//set each element of flags to 1

	unsigned int num_layers = 0;
	unsigned int layer_count = m_count;
	//construct the Dobkin-Kirkpatrick Hierarchical multi-layer
	while (layer_count > 4)//vertex number of the hull is greater than 4
	{
		//T_IndexBuf *indexBuf = adjacency_graph(m_count, m_verts, flags);//new adjacency with only hull points
		T_MultiIndexBuf indexBuf = adjacency_graph(m_count, m_verts, flags);//new adjacency with only hull points

		unsigned int i;
		for (i = 0; i != m_count; ++i)
		{
			if (flags[i])
			{
				assert(!indexBuf[i].empty());
				cobound[i].push_back(indexBuf[i]);//at the first iteration, each indexBuf is push into cobound respectively.
			}
		}

		++num_layers;
		indexBuf.clear();
		//delete [] indexBuf;

		std::fill(&flags[0], &flags[m_count], 0);//set each element of flag to 0

		for (i = 0; i != m_count; ++i)
		{
			if (cobound[i].size() == num_layers)
			{
				T_IndexBuf& curr_cobound = cobound[i].back();
				if (!flags[i] && curr_cobound.size() <= 8)
				{
					for (unsigned int j  = 0; j != curr_cobound.size(); ++j)
					{
						flags[curr_cobound[j]] = 1;//current cobound
					}
				}
			}
		}

		layer_count = 0;

		for (i = 0; i != m_count; ++i)
		{
			if (flags[i])//the vertex with flag = 1 is included in layer_count
			{
				++layer_count;
			}
		}
	}

	indexBuf = simplex_adjacency_graph(m_count, flags);

	for (unsigned int i = 0; i < m_count; i++)
	{
		if (flags[i])
		{
			assert(!indexBuf[i].empty());
			cobound[i].push_back(indexBuf[i]);
		}
	}

	++num_layers;
	indexBuf.clear();
	//delete [] indexBuf;
	delete [] flags;



#ifdef DK_HIERARCHY
	prune(m_count, cobound);
#endif

	m_cobound = new T_MultiIndexArray[m_count];

	for (unsigned int i = 0; i < m_count; i++)
	{
		new (&m_cobound[i]) T_MultiIndexArray(cobound[i].size());

		unsigned int j;
		for (j = 0; j != cobound[i].size(); ++j)
		{
			new (&m_cobound[i][j]) DT_IndexArray(cobound[i][j].size(), &cobound[i][j][0]);
		}
	}

	delete [] cobound;

	m_start_vertex = 0;
	while (m_cobound[m_start_vertex].size() != num_layers)
	{
		++m_start_vertex;
		assert(m_start_vertex < m_count);
	}

	m_curr_vertex = m_start_vertex;
}


DT_Polyhedron::~DT_Polyhedron()
{
	//delete [] m_verts;
    delete [] m_cobound;
}

#ifdef DK_HIERARCHY

SD_Scalar DT_Polyhedron::supportH(const Vector3r& v) const
{
    //m_curr_vertex = m_start_vertex;
		//use curr_vertex to replace m_curr_vertex
    unsigned int curr_vertex = m_start_vertex;
    //SD_Scalar d = (*this)[curr_vertex].dot(v);
		SD_Scalar d = m_verts[curr_vertex].dot(v);
    SD_Scalar h = d;

	for (int curr_layer = m_cobound[m_start_vertex].size(); curr_layer != 0; --curr_layer)
	{
		const DT_IndexArray& curr_cobound = m_cobound[curr_vertex][curr_layer-1];

		for (unsigned int i = 0; i != curr_cobound.size(); ++i)
		{
			//d = (*this)[curr_cobound[i]].dot(v);
			d = m_verts[curr_cobound[i]].dot(v);
			if (d > h)
			{
				curr_vertex = curr_cobound[i];
				h = d;
			}
		}
	}

    return h;
}

Vector3r DT_Polyhedron::support(const Vector3r& v) const
{
	//m_curr_vertex = m_start_vertex;
    //use curr_vertex to replace m_curr_vertex for OpenMP
    unsigned int curr_vertex = m_start_vertex;
    //SD_Scalar d = (*this)[curr_vertex].dot(v);
		SD_Scalar d = m_verts[curr_vertex].dot(v);
    SD_Scalar h = d;
	int curr_layer;
	for (curr_layer = m_cobound[m_start_vertex].size(); curr_layer != 0; --curr_layer)
	{
		const DT_IndexArray& curr_cobound = m_cobound[curr_vertex][curr_layer-1];
        unsigned int i;
		for (i = 0; i != curr_cobound.size(); ++i)
		{
			//d = (*this)[curr_cobound[i]].dot(v);
			d = m_verts[curr_cobound[i]].dot(v);
			if (d > h)
			{
				curr_vertex = curr_cobound[i];
				h = d;
			}
		}
	}

    //return (*this)[curr_vertex];
		return m_verts[curr_vertex];
}

#else

SD_Scalar DT_Polyhedron::supportH(const Vector3r& v) const
{
    int last_vertex = -1;
    SD_Scalar d = (*this)[m_curr_vertex].dot(v);
    SD_Scalar h = d;

	for (;;)
	{
        DT_IndexArray& curr_cobound = m_cobound[m_curr_vertex][0];
        int i = 0, n = curr_cobound.size();
        while (i != n &&
               (curr_cobound[i] == last_vertex ||
				(d = (*this)[curr_cobound[i]].dot(v)) - h <= fabs(h) * Mathr::EPSILON))
		{
            ++i;
		}

        if (i == n)
		{
			break;
		}

        last_vertex = m_curr_vertex;
        m_curr_vertex = curr_cobound[i];
        h = d;
    }
    return h;
}

Vector3r DT_Polyhedron::support(const Vector3r& v) const
{
	int last_vertex = -1;
    SD_Scalar d = (*this)[m_curr_vertex].dot(v);
    SD_Scalar h = d;

    for (;;)
	{
        DT_IndexArray& curr_cobound = m_cobound[m_curr_vertex][0];
        int i = 0, n = curr_cobound.size();
        while (i != n &&
               (curr_cobound[i] == last_vertex ||
				(d = (*this)[curr_cobound[i]].dot(v)) - h <= fabs(h) * Mathr::EPSILON))
		{
            ++i;
		}

        if (i == n)
		{
			break;
		}

		last_vertex = m_curr_vertex;
        m_curr_vertex = curr_cobound[i];
        h = d;
    }
    return (*this)[m_curr_vertex];
}

#endif
////////////////////
