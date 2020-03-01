#include"FacetTopologyAnalyzer.hpp"
#include<sudodem/pkg/common/Facet.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Body.hpp>

CREATE_LOGGER(FacetTopologyAnalyzer);
SUDODEM_PLUGIN((FacetTopologyAnalyzer));
#ifndef FACET_TOPO
void FacetTopologyAnalyzer::action(){
	throw runtime_error("FACET_TOPO was not enabled in Facet.hpp at compile-time. Do not use FacetTopologyAnalyzer or recompile.");
}
#else
void FacetTopologyAnalyzer::action(){
	commonEdgesFound=0;
	LOG_DEBUG("Projection axis for analysis is "<<projectionAxis);
	vector<shared_ptr<VertexData> > vv;
	// minimum facet edge length (tolerance scale)
	Real minSqLen=numeric_limits<Real>::infinity();
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		shared_ptr<Facet> f=SUDODEM_PTR_DYN_CAST<Facet>(b->shape);
		if(!f) continue;
		const Vector3r& pos=b->state->pos;
		for(size_t i=0; i<3; i++){
			vv.push_back(shared_ptr<VertexData>(new VertexData(b->getId(),i,f->vertices[i]+pos,(f->vertices[i]+pos).dot(projectionAxis))));
			minSqLen=min(minSqLen,(f->vertices[i]-f->vertices[(i+1)%3]).squaredNorm());
		}
	}
	LOG_DEBUG("Added data for "<<vv.size()<<" vertices ("<<vv.size()/3.<<" facets).");
	Real tolerance=sqrt(minSqLen)*relTolerance;
	LOG_DEBUG("Absolute tolerance is "<<tolerance);
	sort(vv.begin(),vv.end(),VertexComparator());
	size_t nVertices=vv.size(), j;
	for(size_t i=0; i<nVertices; i++){
		j=i;
		while(++j<nVertices && (vv[j]->coord-vv[i]->coord)<=tolerance){
			shared_ptr<VertexData> &vi=vv[i], &vj=vv[j];
			if(std::abs(vi->pos[0]-vj->pos[0])<=tolerance &&
				std::abs(vi->pos[1]-vj->pos[1])<=tolerance &&
				std::abs(vi->pos[2]-vj->pos[2])<=tolerance &&
				(vi->pos-vj->pos).squaredNorm()<=tolerance*tolerance){
				// OK, these two vertices touch
				// LOG_TRACE("Vertices "<<vi->id<<"/"<<vi->vertexNo<<" and "<<vj->id<<"/"<<vj->vertexNo<<" close enough.");
				// add vertex to the nextIndetical of the one that has lower index; the one that is added will have isLowestIndex=false
				if(vi->index<vj->index){ vi->nextIdentical.push_back(vj); vj->isLowestIndex=false; }
				else{                    vj->nextIdentical.push_back(vi); vi->isLowestIndex=false; }
			}
		}
	}
	// identity chains start always at lower indices, this way we get all of them
	sort(vv.begin(),vv.end(),VertexIndexComparator());
	int maxVertexId=0;
	FOREACH(shared_ptr<VertexData>& v, vv){
		if(v->vertexId<0){
			assert(v->isLowestIndex);
			v->vertexId=maxVertexId++;
		}
		FOREACH(shared_ptr<VertexData>& vNext, v->nextIdentical){
			vNext->vertexId=v->vertexId;
		}
		if(v->vertexId>=0) continue; // already assigned
	}
	LOG_DEBUG("Found "<<maxVertexId<<" unique vertices.");
	commonVerticesFound=maxVertexId;
	// add FacetTopology for all facets; index within the topo array is the body id
	vector<shared_ptr<FacetTopology> > topo(scene->bodies->size()); // initialized with the default ctor
	FOREACH(shared_ptr<VertexData>& v, vv){
		if(!topo[v->id]) topo[v->id]=shared_ptr<FacetTopology>(new FacetTopology(v->id));
		topo[v->id]->vertices[v->vertexNo]=v->vertexId;
	}
	// make sure all facets have their vertex id's assigned
	// add non-empty ones to topo1 that will be used for adjacency search afterwards
	vector<shared_ptr<FacetTopology> > topo1;
	for(size_t i=0; i<topo.size(); i++){
		shared_ptr<FacetTopology> t=topo[i];
		if(!t) continue;
		if(t->vertices[0]<0 || t->vertices[1]<0 || t->vertices[2]<0){
			LOG_FATAL("Facet #"<<i<<": some vertex has no integrized vertexId assigned!!");
			LOG_FATAL("Vertices are: "<<t->vertices[0]<<","<<t->vertices[1]<<","<<t->vertices[2]);
			throw logic_error("Facet vertex has no integrized vertex assigned?!");
		}
		topo1.push_back(t);
	}
	std::sort(topo1.begin(),topo1.end(),FacetTopology::MinVertexComparator());
	size_t nTopo=topo1.size();
	for(size_t i=0; i<nTopo; i++){
		size_t j=i;
		const shared_ptr<FacetTopology>& ti(topo1[i]);
		long tiMaxVertex=ti->maxVertex();
		while(++j<nTopo){
			const shared_ptr<FacetTopology> &tj(topo1[j]);
			/* since facets are sorted by their min vertex id,
				we know that it is safe to skip all the rest
				as soon as max vertex of ti one is smaller than min vertex of tj, as i<=j */
			if(tj->minVertex()>tiMaxVertex) break;
			vector<size_t> vvv; // array of common vertices
			for(size_t k=0; k<3; k++){
				if     (ti->vertices[k]==tj->vertices[0]) vvv.push_back(ti->vertices[k]);
				else if(ti->vertices[k]==tj->vertices[1]) vvv.push_back(ti->vertices[k]);
				else if(ti->vertices[k]==tj->vertices[2]) vvv.push_back(ti->vertices[k]);
			}
			if(vvv.size()<2) continue;
			assert(vvv.size()!=3 /* same coords? nonsense*/ ); assert(vvv.size()==2);
			// from here, we know ti and tj are adjacent
			vector<int> edge(2,0); int &ei(edge[0]),&ej(edge[1]);
			// identify what edge are we at, for both facets
			for(int k=0; k<2; k++){
				for(edge[k]=0; edge[k]<3; edge[k]++){
					const shared_ptr<FacetTopology>& tt( k==0 ? ti : tj);
					size_t v1=tt->vertices[edge[k]],v2=tt->vertices[(edge[k]+1)%3];
					if((vvv[0]==v1 && vvv[1]==v2) || (vvv[0]==v2 && vvv[1]==v1)) break;
					if(edge[k]==2){
						LOG_FATAL("No edge identified for 2 vertices "<<vvv[0]<<","<<vvv[1]<<" (facet #"<<tt->id<<" is: "<<tt->vertices[0]<<","<<tt->vertices[1]<<","<<tt->vertices[2]<<")");
						throw logic_error("No edge identified for given vertices.");
					}
				}
			}
			// add adjacency information to the facet itself
			Facet *f1=SUDODEM_CAST<Facet*>((*scene->bodies)[ti->id]->shape.get()), *f2=SUDODEM_CAST<Facet*>((*scene->bodies)[tj->id]->shape.get());
			f1->edgeAdjIds[ei]=ti->id; f2->edgeAdjIds[ej]=tj->id;
			// normals are in the sense of vertex rotation (right-hand rule); therefore, if vertices of the adjacent edge are opposite on each facet, normals are in the same direction
			bool invNormals=(ti->vertices[ei]==tj->vertices[ej]);
			assert(
				( invNormals && (ti->vertices[ ei     ]==tj->vertices[ej]) && (ti->vertices[(ei+1)%3]==tj->vertices[(ej+1)%3]) ) ||
				(!invNormals && (ti->vertices[(ei+1)%3]==tj->vertices[ej]) && (ti->vertices[ ei     ]==tj->vertices[(ej+1)%3]) ));
			// angle between normals
			const shared_ptr<Body>& b1=Body::byId(ti->id,scene); const shared_ptr<Body>& b2=Body::byId(tj->id,scene);
			Vector3r n1g=b1->state->ori*f1->normal, n2g=b2->state->ori*f2->normal;
			//TRVAR2(n1g,n2g);
			Vector3r contEdge1g=b1->state->ori*(f1->vertices[(ei+1)%3]-f1->vertices[ei]); // vector of the edge of contact in global coords
			Quaternionr q12; q12.setFromTwoVectors(n1g,(invNormals?-1.:1.)*n2g); AngleAxisr aa12(q12); Real halfAngle=.5*aa12.angle();
			assert(halfAngle>=0 && halfAngle<=Mathr::HALF_PI);
			if(aa12.axis().dot(contEdge1g)<0 /* convex contact from the side of +n1 */ ) halfAngle*=-1.;
			f1->edgeAdjHalfAngle[ei]=halfAngle;
			f2->edgeAdjHalfAngle[ej]=(invNormals ? -halfAngle : halfAngle);
			commonEdgesFound++;
			LOG_TRACE("Found adjacent #"<<ti->id<<"+#"<<tj->id<<"; common edges "<<ei<<"+"<<ej<<", normals "<<(invNormals?"inversed":"congruent")<<", halfAngle "<<halfAngle<<")");
		}
	}
}
#endif


