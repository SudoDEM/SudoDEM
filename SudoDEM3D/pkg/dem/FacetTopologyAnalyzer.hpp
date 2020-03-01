// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/core/Interaction.hpp>
/*! Initializer for filling adjacency geometry data for facets.
 *
 * Common vertices and common edges are identified and mutual angle between facet faces
 * is written to Facet instances.
 * If facets don't move with respect to each other, this must be done only at the beginning.
 */
class FacetTopologyAnalyzer: public GlobalEngine{
	struct VertexData{
		VertexData(Body::id_t _id, int _vertexNo, Vector3r _pos, Real _coord): id(_id), vertexNo(_vertexNo), coord(_coord), pos(_pos){index=3*id+vertexNo; isLowestIndex=true; vertexId=-1;}
		//! Facet (body id) that we represent
		Body::id_t id;
		//! vertex number within this Facet
		int vertexNo;
		//! projected coordinate along projectionAxis
		Real coord;
		//! global coordinates
		Vector3r pos;
		//! index of this vertex (canonical ordering)
		size_t index;
		//! is this vertex the first one in the "identity" graph?
		bool isLowestIndex;
		//! vertices that are "identical" with this one, but have higer indices
		vector<shared_ptr<VertexData> > nextIdentical;
		//! id of vertex, once they have id's assigned
		long vertexId;
	};
	struct VertexComparator{
		bool operator()(const shared_ptr<VertexData>& v1, const shared_ptr<VertexData>& v2){return v1->coord<v2->coord;}
	};
	struct VertexIndexComparator{
		bool operator()(const shared_ptr<VertexData>& v1, const shared_ptr<VertexData>& v2){return v1->index<v2->index;}
	};
	struct FacetTopology{
		FacetTopology(Body::id_t _id): id(_id){vertices[0]=vertices[1]=vertices[2]=-1;}
		//! integrized vertices
		long vertices[3];
		//! facet id, for back reference
		Body::id_t id;
		long minVertex(){return min(vertices[0],min(vertices[1],vertices[2]));}
		long maxVertex(){return max(vertices[0],max(vertices[1],vertices[2]));}
		struct MinVertexComparator{
			bool operator()(const shared_ptr<FacetTopology>& t1, const shared_ptr<FacetTopology>& t2){ return t1->minVertex()<t2->minVertex();}
		};
	};
	public:
		void action();

	SUDODEM_CLASS_BASE_DOC_ATTRS(FacetTopologyAnalyzer,GlobalEngine,"Initializer for filling adjacency geometry data for facets.\n\nCommon vertices and common edges are identified and mutual angle between facet faces is written to Facet instances. If facets don't move with respect to each other, this must be done only at the beginng.",
		((Vector3r,projectionAxis,Vector3r::UnitX(),,"Axis along which to do the initial vertex sort"))
		((Real,relTolerance,1e-4,,"maximum distance of 'identical' vertices, relative to minimum facet size"))
		((long,commonEdgesFound,0,,"how many common edges were identified during last run. |yupdate|"))
		((long,commonVerticesFound,0,,"how many common vertices were identified during last run. |yupdate|"))
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(FacetTopologyAnalyzer);



