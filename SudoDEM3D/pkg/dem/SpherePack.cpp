// © 2009 Václav Šmilauer <eudoxos@arcig.cz>

#include<sudodem/pkg/dem/SpherePack.hpp>

#include<sudodem/core/Omega.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/dem/Shop.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include<sudodem/core/Timing.hpp>

// not a serializable in the sense of SUDODEM_PLUGIN

CREATE_LOGGER(SpherePack);

namespace py=boost::python;


void SpherePack::fromList(const py::list& l){
	pack.clear();
	size_t len=py::len(l);
	for(size_t i=0; i<len; i++){
		const py::tuple& t=py::extract<py::tuple>(l[i]);
		py::extract<Vector3r> vec(t[0]);
		if(vec.check()) { pack.push_back(Sph(vec(),py::extract<double>(t[1]),(py::len(t)>2?py::extract<int>(t[2]):-1))); continue; }
		PyErr_SetString(PyExc_TypeError, "List elements must be (Vector3, float) or (Vector3, float, int)!");
		py::throw_error_already_set();
	}
};

void SpherePack::fromLists(const vector<Vector3r>& centers, const vector<Real>& radii){
	pack.clear();
	if(centers.size()!=radii.size()) throw std::invalid_argument(("The same number of centers and radii must be given (is "+boost::lexical_cast<string>(centers.size())+", "+boost::lexical_cast<string>(radii.size())+")").c_str());
	size_t l=centers.size();
	for(size_t i=0; i<l; i++){
		add(centers[i],radii[i]);
	}
	cellSize=Vector3r::Zero();
}

py::list SpherePack::toList() const {
	py::list ret;
	FOREACH(const Sph& s, pack) ret.append(s.asTuple());
	return ret;
};

void SpherePack::fromFile(string file) {
	typedef boost::tuple<Vector3r,Real,int> tupleVector3rRealInt;
	vector<tupleVector3rRealInt> ss;
	Vector3r mn,mx;
	ss=Shop::loadSpheresFromFile(file,mn,mx,&cellSize);
	pack.clear();
	FOREACH(const tupleVector3rRealInt& s, ss) pack.push_back(Sph(boost::get<0>(s),boost::get<1>(s),boost::get<2>(s)));
}

void SpherePack::toFile(const string fname) const {
	ofstream f(fname.c_str());
	if(!f.good()) throw runtime_error("Unable to open file `"+fname+"'");
	if(cellSize!=Vector3r::Zero()){ f<<"##PERIODIC:: "<<cellSize[0]<<" "<<cellSize[1]<<" "<<cellSize[2]<<endl; }
	FOREACH(const Sph& s, pack){
		//if(s.clumpId>=0) throw std::invalid_argument("SpherePack with clumps cannot be (currently) exported to a text file.");
		f<<s.c[0]<<" "<<s.c[1]<<" "<<s.c[2]<<" "<<s.r<<" "<<s.clumpId<<endl;
	}
	f.close();
};

void SpherePack::fromSimulation() {
	pack.clear();
	Scene* scene=Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		shared_ptr<Sphere> intSph=SUDODEM_PTR_DYN_CAST<Sphere>(b->shape);
		if(!intSph) continue;
		pack.push_back(Sph(b->state->pos,intSph->radius,(b->isClumpMember()?b->clumpId:-1)));
	}
	if(scene->isPeriodic) {
		cellSize=scene->cell->getSize();
		isPeriodic = true;
	}
}

long SpherePack::makeCloud(Vector3r mn, Vector3r mx, Real rMean, Real rRelFuzz, int num, bool periodic, Real porosity, const vector<Real>& psdSizes, const vector<Real>& psdCumm, bool distributeMass, int seed, Matrix3r hSize){
	isPeriodic = periodic;
	static boost::minstd_rand randGen(seed!=0?seed:(int)TimingInfo::getNow(/* get the number even if timing is disabled globally */ true));
	static boost::variate_generator<boost::minstd_rand&, boost::uniform_real<Real> > rnd(randGen, boost::uniform_real<Real>(0,1));
	vector<Real> psdRadii; // holds plain radii (rather than diameters), scaled down in some situations to get the target number
	vector<Real> psdCumm2; // psdCumm but dimensionally transformed to match mass distribution
	Vector3r size;
	bool hSizeFound =(hSize!=Matrix3r::Zero());//is hSize passed to the function?
	size=mx-mn;
	if (!hSizeFound) {hSize=size.asDiagonal();}
	if (hSizeFound && !periodic) LOG_WARN("hSize can be defined only for periodic cells.");
	Real volume=hSize.determinant();
	Matrix3r invHsize =hSize.inverse();
	Real area=std::abs(size[0]*size[2]+size[0]*size[1]+size[1]*size[2]);//2 terms will be null if one coordinate is 0, the other is the area
	if (!volume) {
		if (hSizeFound) throw invalid_argument("The period defined by hSize has null length in at least one direction, this is not supported. Define flat boxes via min-max and keep hSize undefined if you want a 2D packing.");
		else LOG_WARN("The volume of the min-max box is null, we will assume that the packing is 2D. If it is not what you want then you defined wrong input values; check that min and max corners are defined correctly.");}
	int mode=-1; bool err=false;
	// determine the way we generate radii
	if(porosity<=0) {LOG_WARN("porosity must be >0, changing it for you. It will be ineffective if rMean>0."); porosity=0.5;}
	//If rMean is not defined, then in will be defined in RDIST_NUM
	if(rMean>0) mode=RDIST_RMEAN;
	else if(num>0 && psdSizes.size()==0) {
		mode=RDIST_NUM;
		// the term (1+rRelFuzz²) comes from the mean volume for uniform distribution : Vmean = 4/3*pi*Rmean*(1+rRelFuzz²)
		if (volume) rMean=pow(volume*(1-porosity)/(Mathr::PI*(4/3.)*(1+rRelFuzz*rRelFuzz)*num),1/3.);
		else {//The volume is null, we will generate a 2D packing with the following rMean
			if (!area) throw invalid_argument("The box defined has null volume AND null surface. Define at least maxCorner of the box, or hSize if periodic.");
			rMean=pow(area*(1-porosity)/(Mathr::PI*(1+rRelFuzz*rRelFuzz)*num),0.5);}
	}
	// transform sizes and cummulated fractions values in something convenient for the generation process
	if(psdSizes.size()>0){
		err=(mode>=0); mode=RDIST_PSD;
		if(psdSizes.size()!=psdCumm.size()) throw invalid_argument(("SpherePack.makeCloud: psdSizes and psdCumm must have same dimensions ("+boost::lexical_cast<string>(psdSizes.size())+"!="+boost::lexical_cast<string>(psdCumm.size())).c_str());
		if(psdSizes.size()<=1) throw invalid_argument("SpherePack.makeCloud: psdSizes must have at least 2 items");
		if((*psdCumm.begin())!=0. && (*psdCumm.rbegin())!=1.) throw invalid_argument("SpherePack.makeCloud: first and last items of psdCumm *must* be exactly 0 and 1.");
		psdRadii.reserve(psdSizes.size());
		for(size_t i=0; i<psdSizes.size(); i++) {
			psdRadii.push_back(/* radius, not diameter */ .5*psdSizes[i]);
			if(distributeMass) {
				//psdCumm2 is first obtained by integrating the number of particles over the volumic PSD (dN/dSize = totV*(dPassing/dSize)*1/(4/3πr³)). The total cumulated number will be the number of spheres in volume*(1-porosity), it is used to decide if the PSD will be scaled down. psdCumm2 is normalized below in order to fit in [0,1]. (Bruno C.)
				if (i==0) psdCumm2.push_back(0);
				else psdCumm2.push_back(psdCumm2[i-1] + 3.0* (volume?volume:(area*psdSizes[psdSizes.size()-1])) *(1-porosity)/Mathr::PI*(psdCumm[i]-psdCumm[i-1])/(psdSizes[i]-psdSizes[i-1])*(pow(psdSizes[i-1],-2)-pow(psdSizes[i],-2)));
			}
			LOG_DEBUG(i<<". "<<psdRadii[i]<<", cdf="<<psdCumm[i]<<", cdf2="<<(distributeMass?boost::lexical_cast<string>(psdCumm2[i]):string("--")));
			// check monotonicity
			if(i>0 && (psdSizes[i-1]>psdSizes[i] || psdCumm[i-1]>psdCumm[i])) throw invalid_argument("SpherePack:makeCloud: psdSizes and psdCumm must be both non-decreasing.");
		}
		// check the consistency between sizes, num, and poro if all three are imposed. If target number will not fit in (1-poro)*volume, scale down particles sizes
		if (num>1){
			appliedPsdScaling=1;
			if(distributeMass) {
				if (psdCumm2[psdSizes.size()-1]<num) appliedPsdScaling=pow(psdCumm2[psdSizes.size()-1]/num,1./3.);
			} else {
				double totVol=0;
				for(size_t i=1; i<psdSizes.size(); i++) totVol+= 4/3*Mathr::PI*(psdCumm[i]-psdCumm[i-1])*num*
					pow(0.5*(psdSizes[i]+psdSizes[i-1]),3)*(1+pow(0.5*(psdSizes[i]-psdSizes[i-1]),2));
				Real volumeRatio = totVol/((1-porosity)*(volume?volume:(area*psdSizes[psdSizes.size()-1])));
				if (volumeRatio>1) appliedPsdScaling=pow(volumeRatio,-1./3.);
			}
			if (appliedPsdScaling<1) for(size_t i=0; i<psdSizes.size(); i++) psdRadii[i]*=appliedPsdScaling;
		}
		//Normalize psdCumm2 so it's between 0 and 1
		if(distributeMass) for(size_t i=1; i<psdSizes.size(); i++) psdCumm2[i]/=psdCumm2[psdSizes.size()-1];
	}
	if(err || mode<0) throw invalid_argument("SpherePack.makeCloud: at least one of rMean, porosity, psdSizes & psdCumm arguments must be specified. rMean can't be combined with psdSizes.");
	// adjust uniform distribution parameters with distributeMass; rMean has the meaning (dimensionally) of _volume_
	const int maxTry=1000;
	if(periodic && volume && !hSizeFound)(cellSize=size);
	Real r=0;
	for(int i=0; (i<num) || (num<0); i++) {
		Real norm, rand;
		//Determine radius of the next sphere that will be placed in space. If (num>0), generate radii the deterministic way, in decreasing order, else radii are stochastic since we don't know what the final number will be
		if (num>0) rand = ((Real)num-(Real)i+0.5)/((Real)num+1.);
		else rand = rnd();
		int t;
		switch(mode){
			case RDIST_RMEAN:
			//FIXME : r is never defined, it will be zero at first iteration, but it will have values in the next ones.
			//I don't understand why it apparently works. Some magic?
			case RDIST_NUM:
				if(distributeMass) r=pow3Interp(rand,rMean*(1-rRelFuzz),rMean*(1+rRelFuzz));
				else r=rMean*(2*(rand-.5)*rRelFuzz+1); // uniform distribution in rMean*(1±rRelFuzz)
				break;
			case RDIST_PSD:
				if(distributeMass){
					int piece=psdGetPiece(rand,psdCumm2,norm);
					r=pow3Interp(norm,psdRadii[piece],psdRadii[piece+1]);
				} else {
					int piece=psdGetPiece(rand,psdCumm,norm);
					r=psdRadii[piece]+norm*(psdRadii[piece+1]-psdRadii[piece]);}
		}
		// try to put the sphere into a free spot
		for(t=0; t<maxTry; ++t){
			Vector3r c;
			if(!periodic) { for(int axis=0; axis<3; axis++) c[axis]=mn[axis]+(size[axis]?(size[axis]-2*r)*rnd()+r:0);}//we handle 2D with the special case size[axis]==0
			else { 	for(int axis=0; axis<3; axis++) c[axis]=rnd();//coordinates in [0,1]
				c=mn+hSize*c;}//coordinates in reference frame (inside the base cell)
			size_t packSize=pack.size(); bool overlap=false;
			if(!periodic) for(size_t j=0;j<packSize;j++) {if(pow(pack[j].r+r,2)>=(pack[j].c-c).squaredNorm()) {overlap=true; break;}}
			else {
				for(size_t j=0; j<packSize; j++){
					Vector3r dr=Vector3r::Zero();
					if (!hSizeFound) {//The box is axis-aligned, use the wrap methods
						for(int axis=0; axis<3; axis++) dr[axis]=size[axis]? min(cellWrapRel(c[axis],pack[j].c[axis],pack[j].c[axis]+size[axis]),cellWrapRel(pack[j].c[axis],c[axis],c[axis]+size[axis])) : 0;
					} else {//not aligned, find closest neighbor in a cube of size 1, then transform distance to cartesian coordinates
						Vector3r c1c2=invHsize*(pack[j].c-c);
						for(int axis=0; axis<3; axis++){
							if (std::abs(c1c2[axis])<std::abs(c1c2[axis] - Mathr::Sign(c1c2[axis]))) dr[axis]=c1c2[axis];
							else dr[axis] = c1c2[axis] - Mathr::Sign(c1c2[axis]);}
						dr=hSize*dr;//now in cartesian coordinates
					}
					if(pow(pack[j].r+r,2)>= dr.squaredNorm()){ overlap=true; break; }
				}
			}
			if(!overlap) { pack.push_back(Sph(c,r)); break; }
		}
		if (t==maxTry) {
			if(num>0) {
				if (mode!=RDIST_RMEAN) {
					//if rMean is not imposed, then we call makeCloud recursively, scaling the PSD down until the target num is obtained
					Real nextPoro = porosity+(1-porosity)/10.;
					LOG_WARN("Exceeded "<<maxTry<<" tries to insert non-overlapping sphere to packing. Only "<<i<<" spheres were added, although you requested "<<num<<". Trying again with porosity "<<nextPoro<<". The size distribution is being scaled down");
					pack.clear();
					return makeCloud(mn, mx, -1., rRelFuzz, num, periodic, nextPoro, psdSizes, psdCumm, distributeMass,seed,hSizeFound?hSize:Matrix3r::Zero());}
				else LOG_WARN("Exceeded "<<maxTry<<" tries to insert non-overlapping sphere to packing. Only "<<i<<" spheres were added, although you requested "<<num<<".");
			}
			return i;}
	}
	if (appliedPsdScaling<1) LOG_WARN("The size distribution has been scaled down by a factor pack.appliedPsdScaling="<<appliedPsdScaling);
	return pack.size();
}

void SpherePack::cellFill(Vector3r vol){
	Vector3i count;
	for(int i=0; i<3; i++) count[i]=(int)(ceil(vol[i]/cellSize[i]));
	LOG_DEBUG("Filling volume "<<vol<<" with cell "<<cellSize<<", repeat counts are "<<count);
	cellRepeat(count);
}

void SpherePack::cellRepeat(Vector3i count){
	if(cellSize==Vector3r::Zero()){ throw std::runtime_error("cellRepeat cannot be used on non-periodic packing."); }
	if(count[0]<=0 || count[1]<=0 || count[2]<=0){ throw std::invalid_argument("Repeat count components must be positive."); }
	size_t origSize=pack.size();
	pack.reserve(origSize*count[0]*count[1]*count[2]);
	for(int i=0; i<count[0]; i++){
		for(int j=0; j<count[1]; j++){
			for(int k=0; k<count[2]; k++){
				if((i==0) && (j==0) && (k==0)) continue; // original cell
				Vector3r off(cellSize[0]*i,cellSize[1]*j,cellSize[2]*k);
				for(size_t l=0; l<origSize; l++){
					const Sph& s=pack[l]; pack.push_back(Sph(s.c+off,s.r));
				}
			}
		}
	}
	cellSize=Vector3r(cellSize[0]*count[0],cellSize[1]*count[1],cellSize[2]*count[2]);
}

int SpherePack::psdGetPiece(Real x, const vector<Real>& cumm, Real& norm){
	int sz=cumm.size(); int i=0;
	while(i<sz && cumm[i]<=x) i++; // upper interval limit index
	if((i==sz-1) && cumm[i]<=x){ i=sz-2; norm=1.; return i;}
	i--; // lower interval limit intex
	norm=(x-cumm[i])/(cumm[i+1]-cumm[i]);
	//LOG_TRACE("count="<<sz<<", x="<<x<<", piece="<<i<<" in "<<cumm[i]<<"…"<<cumm[i+1]<<", norm="<<norm);
	return i;
}

py::tuple SpherePack::psd(int bins, bool mass) const {
	if(pack.size()==0) return py::make_tuple(py::list(),py::list()); // empty packing
	// find extrema
	Real minD=std::numeric_limits<Real>::infinity(); Real maxD=-minD;
	// volume, but divided by π*4/3
	Real vol=0; long N=pack.size();
	FOREACH(const Sph& s, pack){ maxD=max(2*s.r,maxD); minD=min(2*s.r,minD); vol+=pow(s.r,3); }
	if(minD==maxD){ minD-=.5; maxD+=.5; } // emulates what numpy.histogram does
	// create bins and bin edges
	vector<Real> hist(bins,0); vector<Real> cumm(bins+1,0); /* cummulative values compute from hist at the end */
	vector<Real> edges(bins+1); for(int i=0; i<=bins; i++){ edges[i]=minD+i*(maxD-minD)/bins; }
	// weight each grain by its "volume" relative to overall "volume"
	FOREACH(const Sph& s, pack){
		int bin=int(bins*(2*s.r-minD)/(maxD-minD)); bin=min(bin,bins-1); // to make sure
		if (mass) hist[bin]+=pow(s.r,3)/vol; else hist[bin]+=1./N;
	}
	for(int i=0; i<bins; i++) cumm[i+1]=min((Real)1.,cumm[i]+hist[i]); // cumm[i+1] is OK, cumm.size()==bins+1
	return py::make_tuple(edges,cumm);
}

/* possible enhacement: proportions parameter, so that the domain is not cube, but box with sides having given proportions */
long SpherePack::particleSD2(const vector<Real>& radii, const vector<Real>& passing, int numSph, bool periodic, Real cloudPorosity, int seed){
	typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
	typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorXr;

	int dim=radii.size()+1;
	if(passing.size()!=radii.size()) throw std::invalid_argument("SpherePack.particleSD2: radii and passing must have the same length.");
	MatrixXr M=MatrixXr::Zero(dim,dim);
	VectorXr rhs=VectorXr::Zero(dim);
	/*

	We know percentages for each fraction (Δpi) and their radii (ri), and want to find
	the number of sphere for each fraction Ni and total solid volume Vs. For each fraction,
	we know that the volume is equal to Ni*(4/3*πri³), which must be equal to Vs*Δpi (Δpi is
	relative solid volume of the i-th fraction).

	The last equation says that total number of particles (sum of fractions) is equal to N,
	which is the total number of particles requested by the user.

	   N1     N2     N3    Vs       rhs

	4/3πr₁³   0      0     -Δp₁   | 0
	  0     4/3πr₂³  0     -Δp₂   | 0
	  0       0    4/3πr₃³ -Δp₃   | 0
     1       1      1      0     | N

	*/
	for(int i=0; i<dim-1; i++){
		M(i,i)=(4/3.)*Mathr::PI*pow(radii[i],3);
		M(i,dim-1)=-(passing[i]-(i>0?passing[i-1]:0))/100.;
		M(dim-1,i)=1;
	}
	rhs[dim-1]=numSph;
	// NumsVs=M^-1*rhs: number of spheres and volume of solids
	VectorXr NumsVs(dim); NumsVs=M.inverse()*rhs;
	Real Vs=NumsVs[dim-1]; // total volume of solids
	Real Vtot=Vs/(1-cloudPorosity); // total volume of cell containing the packing
	Vector3r cellSize=pow(Vtot,1/3.)*Vector3r().Ones(); // for now, assume always cubic sample
	Real rMean=pow(Vs/(numSph*(4/3.)*Mathr::PI),1/3.); // make rMean such that particleSD will compute the right Vs (called a bit confusingly Vtot anyway) inversely
	// cerr<<"Vs="<<Vs<<", Vtot="<<Vtot<<", rMean="<<rMean<<endl;
	// cerr<<"cellSize="<<cellSize<<", rMean="<<rMean<<", numSph="<<numSph<<endl;
	return particleSD(Vector3r::Zero(),cellSize,rMean,periodic,"",numSph,radii,passing,false);
};

// TODO: header, python wrapper, default params

// Discrete particle size distribution
long SpherePack::particleSD(Vector3r mn, Vector3r mx, Real rMean, bool periodic, string name, int numSph, const vector<Real>& radii, const vector<Real>& passing, bool passingIsNotPercentageButCount, int seed){
	vector<Real> numbers;
	if(!passingIsNotPercentageButCount){
		Real Vtot=numSph*4./3.*Mathr::PI*pow(rMean,3.); // total volume of the packing (computed with rMean)

		// calculate number of spheres necessary per each radius to match the wanted psd
		// passing has to contain increasing values
		for (size_t i=0; i<radii.size(); i++){
			Real volS=4./3.*Mathr::PI*pow(radii[i],3.);
			if (i==0) {numbers.push_back(passing[i]/100.*Vtot/volS);}
			else {numbers.push_back((passing[i]-passing[i-1])/100.*Vtot/volS);} //
			cout<<"fraction #"<<i<<" ("<<passing[i]<<"%, r="<<radii[i]<<"): "<<numbers[i]<<" spheres, fraction/cloud volumes "<<volS<<"/"<<Vtot<<endl;
		}
	} else {
		FOREACH(Real p, passing) numbers.push_back(p);
	}

	static boost::minstd_rand randGen(seed!=0?seed:(int)TimingInfo::getNow(true));
	static boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > rnd(randGen, boost::uniform_real<>(0,1));

	const int maxTry=1000;
	Vector3r size=mx-mn;
	if(periodic)(cellSize=size);
	for (int ii=(int)radii.size()-1; ii>=0; ii--){
		Real r=radii[ii]; // select radius
		for(int i=0; i<numbers[ii]; i++) { // place as many spheres as required by the psd for the selected radius into the free spot
			int t;
			for(t=0; t<maxTry; ++t){
				Vector3r c;
				if(!periodic) { for(int axis=0; axis<3; axis++) c[axis]=mn[axis]+r+(size[axis]-2*r)*rnd(); }
				else { for(int axis=0; axis<3; axis++) c[axis]=mn[axis]+size[axis]*rnd(); }
				size_t packSize=pack.size(); bool overlap=false;
				if(!periodic){
					for(size_t j=0; j<packSize; j++){ if(pow(pack[j].r+r,2) >= (pack[j].c-c).squaredNorm()) { overlap=true; break; } }
				} else {
					for(size_t j=0; j<packSize; j++){
						Vector3r dr;
						for(int axis=0; axis<3; axis++) dr[axis]=min(cellWrapRel(c[axis],pack[j].c[axis],pack[j].c[axis]+size[axis]),cellWrapRel(pack[j].c[axis],c[axis],c[axis]+size[axis]));
						if(pow(pack[j].r+r,2)>= dr.squaredNorm()){ overlap=true; break; }
					}
				}
				if(!overlap) { pack.push_back(Sph(c,r)); break; }
			}
			if (t==maxTry) {
				if(numbers[ii]>0) LOG_WARN("Exceeded "<<maxTry<<" tries to insert non-overlapping sphere to packing. Only "<<i<<" spheres were added, although you requested "<<numbers[ii]<<" with radius "<<radii[ii]);
				return i;
			}
		}
	}
	return pack.size();
}

// 2d function
long SpherePack::particleSD_2d(Vector2r mn, Vector2r mx, Real rMean, bool periodic, string name, int numSph, const vector<Real>& radii, const vector<Real>& passing, bool passingIsNotPercentageButCount, int seed){
	vector<Real> numbers;
	if(!passingIsNotPercentageButCount){
		Real Vtot=numSph*4./3.*Mathr::PI*pow(rMean,3.); // total volume of the packing (computed with rMean)

		// calculate number of spheres necessary per each radius to match the wanted psd
		// passing has to contain increasing values
		for (size_t i=0; i<radii.size(); i++){
			Real volS=4./3.*Mathr::PI*pow(radii[i],3.);
			if (i==0) {numbers.push_back(passing[i]/100.*Vtot/volS);}
			else {numbers.push_back((passing[i]-passing[i-1])/100.*Vtot/volS);} //
			cout<<"fraction #"<<i<<" ("<<passing[i]<<"%, r="<<radii[i]<<"): "<<numbers[i]<<" spheres, fraction/cloud volumes "<<volS<<"/"<<Vtot<<endl;
		}
	} else {
		FOREACH(Real p, passing) numbers.push_back(p);
	}

	static boost::minstd_rand randGen(seed!=0?seed:(int)TimingInfo::getNow(true));
	static boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > rnd(randGen, boost::uniform_real<>(0,1));

	const int maxTry=1000;
	Vector2r size=mx-mn;
	//if(periodic)(cellSize=size); in this case, it must be defined in py script as cell needs the third dimension
	for (int ii=(int)radii.size()-1; ii>=0; ii--){
		Real r=radii[ii]; // select radius
		for(int i=0; i<numbers[ii]; i++) { // place as many spheres as required by the psd for the selected radius into the free spot
			int t;
			for(t=0; t<maxTry; ++t){
				Vector3r c;
				if(!periodic) { for(int axis=0; axis<2; axis++) c[axis]=mn[axis]+r+(size[axis]-2*r)*rnd(); }
				else { for(int axis=0; axis<2; axis++) c[axis]=mn[axis]+size[axis]*rnd(); }
				size_t packSize=pack.size(); bool overlap=false;
				if(!periodic){
					for(size_t j=0; j<packSize; j++){ if(pow(pack[j].r+r,2) >= (pack[j].c-c).squaredNorm()) { overlap=true; break; } }
				} else {
					for(size_t j=0; j<packSize; j++){
						Vector3r dr=Vector3r::Zero();
						for(int axis=0; axis<2; axis++) dr[axis]=min(cellWrapRel(c[axis],pack[j].c[axis],pack[j].c[axis]+size[axis]),cellWrapRel(pack[j].c[axis],c[axis],c[axis]+size[axis]));
						if(pow(pack[j].r+r,2)>= dr.squaredNorm()){ overlap=true; break; }
					}
				}
				if(!overlap) { pack.push_back(Sph(c,r)); break; }
			}
			if (t==maxTry) {
				if(numbers[ii]>0) LOG_WARN("Exceeded "<<maxTry<<" tries to insert non-overlapping sphere to packing. Only "<<i<<" spheres were added, although you requested "<<numbers[ii]<<" with radius "<<radii[ii]);
				return i;
			}
		}
	}
	return pack.size();
}

long SpherePack::makeClumpCloud(const Vector3r& mn, const Vector3r& mx, const vector<shared_ptr<SpherePack> >& _clumps, bool periodic, int num, int seed){
	// recenter given clumps and compute their margins
	vector<SpherePack> clumps; /* vector<Vector3r> margins; */ Vector3r boxMargins(Vector3r::Zero()); Real maxR=0;
	vector<Real> boundRad; // squared radii of bounding sphere for each clump
	FOREACH(const shared_ptr<SpherePack>& c, _clumps){
		SpherePack c2(*c);
		c2.translate(c2.midPt()); //recenter
		clumps.push_back(c2);
		Real r=0;
		FOREACH(const Sph& s, c2.pack) r=max(r,s.c.norm()+s.r);
		boundRad.push_back(r);
		Vector3r cMn,cMx; c2.aabb(cMn,cMx); // centered at zero now, this gives margin
		//margins.push_back(periodic?cMx:Vector3r::Zero());
		//boxMargins=boxMargins.cwise().max(cMx);
		FOREACH(const Sph& s, c2.pack) maxR=max(maxR,s.r); // keep track of maximum sphere radius
	}
	std::list<ClumpInfo> clumpInfos;
	Vector3r size=mx-mn;
	if(periodic)(cellSize=size);
	const int maxTry=200;
	int nGen=0; // number of clumps generated
	// random point coordinate generator, with non-zero margins if aperiodic
 	static boost::minstd_rand randGen(seed!=0?seed:(int)TimingInfo::getNow(/* get the number even if timing is disabled globally */ true));
 	static boost::variate_generator<boost::minstd_rand&, boost::uniform_real<Real> > rnd(randGen, boost::uniform_real<Real>(0,1));
	while(nGen<num || num<0){
		int clumpChoice=(int)(rnd()*(clumps.size()-1e-20));
		int tries=0;
		while(true){ // check for tries at the end
			Vector3r pos(0.,0.,0.);
			for(int i=0;i<3;i++){
				pos[i]=rnd()*(mx[i]-mn[i])+mn[i];
			}
			// TODO: check this random orientation is homogeneously distributed
			Quaternionr ori(rnd(),rnd(),rnd(),rnd()); ori.normalize();
			// copy the packing and rotate
			SpherePack C(clumps[clumpChoice]); C.rotateAroundOrigin(ori); C.translate(pos);
			const Real& rad(boundRad[clumpChoice]);
			ClumpInfo ci; // to be used later, but must be here because of goto's

			// find collisions
			// check against bounding spheres of other clumps, and only check individual spheres if there is overlap
			if(!periodic){
				// check overlap with box margins first
				if((pos+rad*Vector3r::Ones()).cwiseMax(mx)!=mx || (pos-rad*Vector3r::Ones()).cwiseMin(mn)!=mn){ FOREACH(const Sph& s, C.pack) if((s.c+s.r*Vector3r::Ones()).cwiseMax(mx)!=mx || (s.c-s.r*Vector3r::Ones()).cwiseMin(mn)!=mn) goto overlap; }
				// check overlaps with bounding spheres of other clumps
				FOREACH(const ClumpInfo& cInfo, clumpInfos){
					bool detailedCheck=false;
					// check overlaps between individual spheres and bounding sphere of the other clump
					if((pos-cInfo.center).squaredNorm()<pow(rad+cInfo.rad,2)){ FOREACH(const Sph& s, C.pack) if(pow(s.r+cInfo.rad,2)>(s.c-cInfo.center).squaredNorm()){ detailedCheck=true; break; }}
					// check sphere-by-sphere, since bounding spheres did overlap
					if(detailedCheck){ FOREACH(const Sph& s, C.pack) for(int id=cInfo.minId; id<=cInfo.maxId; id++) if((s.c-pack[id].c).squaredNorm()<pow(s.r+pack[id].r,2)) goto overlap;}
				}
			} else {
				FOREACH(const ClumpInfo& cInfo, clumpInfos){
					// bounding spheres overlap (in the periodic space)
					if(periPtDistSq(pos,cInfo.center)<pow(rad+cInfo.rad,2)){
						bool detailedCheck=false;
						// check spheres with bounding sphere of the other clump
						FOREACH(const Sph& s, C.pack) if(pow(s.r+cInfo.rad,2)>periPtDistSq(s.c,cInfo.center)){ detailedCheck=true; break; }
						// check sphere-by-sphere
						if(detailedCheck){ FOREACH(const Sph& s, C.pack) for(int id=cInfo.minId; id<=cInfo.maxId; id++) if(periPtDistSq(s.c,pack[id].c)<pow(s.r+pack[id].r,2)) goto overlap; }
					}
				}
			}

			#if 0
			// crude algorithm: check all spheres against all other spheres (slow!!)
			// use vtkPointLocator, add all existing points and check distance of r+maxRadius, then refine
			// for periodicity, duplicate points close than boxMargins to the respective boundary
			if(!periodic){
				for(size_t i=0; i<C.pack.size(); i++){
					for(size_t j=0; j<pack.size(); j++){
						const Vector3r& c(C.pack[i].c); const Real& r(C.pack[i].r);
						if(pow(r+pack[j].r,2)>=(c-pack[j].c).squaredNorm()) goto overlap;
						// check that we are not over the box boundary
						// this could be handled by adjusting the position random interval (by taking off the smallest radius in the clump)
						// but usually the margin band is relatively small and this does not make the code as hairy
						if((c+r*Vector3r::Ones()).cwise().max(mx)!=mx || (c-r*Vector3r::Ones()).cwise().min(mn)!=mn) goto overlap;
					}
				}
			}else{
				for(size_t i=0; i<C.pack.size(); i++){
					for(size_t j=0; j<pack.size(); j++){
						const Vector3r& c(C.pack[i].c); const Real& r(C.pack[i].r);
						Vector3r dr;
						for(int axis=0; axis<3; axis++) dr[axis]=min(cellWrapRel(c[axis],pack[j].c[axis],pack[j].c[axis]+size[axis]),cellWrapRel(pack[j].c[axis],c[axis],c[axis]+size[axis]));
						if(pow(pack[j].r+r,2)>= dr.squaredNorm()) goto overlap;
					}
				}
			}
			#endif

			// add the clump, if no collisions
			/*number clumps consecutively*/ ci.clumpId=nGen; ci.center=pos; ci.rad=rad; ci.minId=pack.size(); ci.maxId=pack.size()+C.pack.size();
			FOREACH(const Sph& s, C.pack){
				pack.push_back(Sph(s.c,s.r,ci.clumpId));
			}
			clumpInfos.push_back(ci);
			nGen++;
			//cerr<<"O";
			break; // break away from the try-loop

			overlap:
			//cerr<<".";
			if(tries++==maxTry){ // last loop
				if(num>0) LOG_WARN("Exceeded "<<maxTry<<" attempts to place non-overlapping clump. Only "<<nGen<<" clumps were added, although you requested "<<num);
				return nGen;
			}
		}
	}
	return nGen;
}

bool SpherePack::hasClumps() const { FOREACH(const Sph& s, pack){ if(s.clumpId>=0) return true; } return false; }
py::tuple SpherePack::getClumps() const{
	std::map<int,py::list> clumps;
	py::list standalone; size_t packSize=pack.size();
	for(size_t i=0; i<packSize; i++){
		const Sph& s(pack[i]);
		if(s.clumpId<0) { standalone.append(i); continue; }
		if(clumps.count(s.clumpId)==0) clumps[s.clumpId]=py::list();
		clumps[s.clumpId].append(i);
	}
	py::list clumpList;
	typedef std::pair<int,py::list> intListPair;
	FOREACH(const intListPair& c, clumps) clumpList.append(c.second);
	return py::make_tuple(standalone,clumpList);
}

