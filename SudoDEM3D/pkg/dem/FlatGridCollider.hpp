// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#include<sudodem/pkg/common/Collider.hpp>
class NewtonIntegrator;
class FlatGridCollider: public Collider{
	struct Grid{
		typedef std::vector<Body::id_t> idVector;
		Vector3i size;
		Vector3r mn, mx;
		Real step;
		// convert point into its integral coordinates (can be outside grid, use fitGrid to coords inside)
		Vector3i pt2int(const Vector3r& pt){ Vector3i ret; for(int i=0;i<3;i++)ret[i]=floor((pt[i]-mn[1])/step); return ret; }
		std::vector<idVector> data;
		// force integral coordinate inside (0,sz-1)
		int fit(int i, int sz) const { return max(0,min(i,sz-1)); }
		Vector3i fitGrid(const Vector3i& v){ return Vector3i(fit(v[0],size[0]),fit(v[1],size[1]),fit(v[2],size[2])); }
		// linearize x,y,z → n in data vector
		size_t lin(int x, int y, int z) const { return fit(x,size[0])+size[0]*fit(y,size[1])+size[0]*size[1]*fit(z,size[2]); }
		// return vector of ids at (x,y,z)
		idVector& operator()(int x, int y, int z){ return data[lin(x,y,z)];}
		idVector& operator()(const Vector3i& v){ return data[lin(v[0],v[1],v[2])];}
		const idVector& operator()(int x, int y, int z) const { return data[lin(x,y,z)];}
	};
	Grid grid;
	int sphereIdx, facetIdx, wallIdx, boxIdx;
	// needed for maxVelSq
	shared_ptr<NewtonIntegrator> newton;
	// track maximum distance of fastest body, at every step in isActivated
	Real fastestBodyMaxDist;
	void initIndices();
	void updateGrid();
	void updateBodyCells(const shared_ptr<Body>& b);
	void updateCollisions();
	virtual void action();
	virtual bool isActivated();
	DECLARE_LOGGER;

	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(FlatGridCollider,Collider,"Non-optimized grid collider, storing grid as dense flat array. Each body is assigned to (possibly multiple) cells, which are arranged in regular grid between *aabbMin* and *aabbMax*, with cell size *step* (same in all directions). Bodies outsize (*aabbMin*, *aabbMax*) are handled gracefully, assigned to closest cells (this will create spurious potential interactions). *verletDist* determines how much is each body enlarged to avoid collision detection at every step.\n\n.. note::\n\tThis collider keeps all cells in linear memory array, therefore will be memory-inefficient for sparse simulations.\n\n.. warning::\n\tobjects :yref:`Body::bound` are not used, :yref:`BoundFunctors<BoundFunctor>` are not used either: assigning cells to bodies is hard-coded internally. Currently handles :yref:`Shapes<Shape>` are: :yref:`Sphere`.\n\n.. note::\n\tPeriodic boundary is not handled (yet).\n\n",
		((Real,verletDist,0,,"Length by which enlarge space occupied by each particle; avoids running collision detection at every step."))
		((Vector3r,aabbMin,Vector3r::Zero(),,"Lower corner of grid."))
		((Vector3r,aabbMax,Vector3r::Zero(),,"Upper corner of grid (approximate, might be rouded up to *minStep*."))
		((Real,step,0,,"Step in the grid (cell size)")),
		initIndices();
	);
};
REGISTER_SERIALIZABLE(FlatGridCollider);
