// Kovthanan …
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/core/GlobalEngine.hpp>
#include<sudodem/pkg/common/Collider.hpp>


class SpheresFactory: public GlobalEngine {
	shared_ptr<Collider> collider;
	protected:
		// Pick random position of a sphere. Should be override in derived engine.
		virtual void pickRandomPosition(Vector3r&/*picked position*/, Real/*sphere's radius*/);
		vector<Real> PSDCurMean;  //Current value of material in each bin
		vector<Real> PSDCurProc;  //Current value of material in each bin, in procents
		vector<Real> PSDNeedProc; //Need value of procent in each bin
		bool PSDuse;        //PSD or not
	public:
		virtual void action();
		struct SpherCoord{
			Vector3r c; Real r;
			SpherCoord(const Vector3r& _c, Real _r){ c=_c; r=_r;}
		};
	DECLARE_LOGGER;
	SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(SpheresFactory,GlobalEngine,"Engine for spitting spheres based on mass flow rate, particle size distribution etc. Initial velocity of particles is given by *vMin*, *vMax*, the *massFlowRate* determines how many particles to generate at each step. When *goalMass* is attained or positive *maxParticles* is reached, the engine does not produce particles anymore. Geometry of the region should be defined in a derived engine by overridden SpheresFactory::pickRandomPosition(). \n\nA sample script for this engine is in :ysrc:`scripts/spheresFactory.py`.",
		((Real,massFlowRate,NaN,,"Mass flow rate [kg/s]"))
		((Real,rMin,NaN,,"Minimum radius of generated spheres (uniform distribution)"))
		((Real,rMax,NaN,,"Maximum radius of generated spheres (uniform distribution)"))
		((Real,vMin,NaN,,"Minimum velocity norm of generated spheres (uniform distribution)"))
		((Real,vMax,NaN,,"Maximum velocity norm of generated spheres (uniform distribution)"))
		((Real,vAngle,NaN,,"Maximum angle by which the initial sphere velocity deviates from the normal."))
		((Vector3r,normal,Vector3r(NaN,NaN,NaN),,"Orientation of the region's geometry, direction of particle's velocites if normalVel is not set."))
		((Vector3r,normalVel,Vector3r(NaN,NaN,NaN),,"Direction of particle's velocites."))
		((int,materialId,-1,,"Shared material id to use for newly created spheres (can be negative to count from the end)"))
		((int,mask,-1,,"groupMask to apply for newly created spheres "))
		((Vector3r,color,Vector3r(-1,-1,-1),,"Use the color for newly created particles, if specified"))
		((vector<int>,ids,,,"ids of created bodies"))
		((Real,totalMass,0,,"Mass of spheres that was produced so far. |yupdate|"))
		((Real,totalVolume,0,,"Volume of spheres that was produced so far. |yupdate|"))
		((Real,goalMass,0,,"Total mass that should be attained at the end of the current step. |yupdate|"))
		((int,maxParticles,100,,"The number of particles at which to stop generating new ones regardless of massFlowRate. if maxParticles=-1 - this parameter is ignored ."))
		((Real,maxMass,-1,,"Maximal mass at which to stop generating new particles regardless of massFlowRate. if maxMass=-1 - this parameter is ignored."))
		((int,numParticles,0,,"Cummulative number of particles produces so far |yupdate|"))
		((int,maxAttempt,5000 ,,"Maximum number of attempts to position a new sphere randomly."))
		((bool,silent,false ,,"If true no complain about excessing maxAttempt but disable the factory (by set massFlowRate=0)."))
		((std::string,blockedDOFs,"" ,,"Blocked degress of freedom"))
		((vector<Real>,PSDsizes,,,"PSD-dispersion, sizes of cells, Diameter [m]"))
		((vector<Real>,PSDcum,,,"PSD-dispersion, cumulative procent meanings [-]"))
		((bool,PSDcalculateMass,true,,"PSD-Input is in mass (true), otherwise the number of particles will be considered."))
		((bool,stopIfFailed,true,,"If true, the SpheresFactory stops (sets massFlowRate=0), when maximal number of attempts to insert particle exceed."))
		((bool,exactDiam,true,,"If true, the particles only with the defined in PSDsizes diameters will be created. Otherwise the diameter will be randomly chosen in the range [PSDsizes[i-1]:PSDsizes[i]], in this case the length of PSDsizes should be  more on 1, than the length of PSDcum.")),
		PSDuse=false;
	);
};
REGISTER_SERIALIZABLE(SpheresFactory);

class CircularFactory: public SpheresFactory {
	protected:
		virtual void pickRandomPosition(Vector3r&, Real);
	public:
		virtual ~CircularFactory(){};
		DECLARE_LOGGER;
		SUDODEM_CLASS_BASE_DOC_ATTRS(CircularFactory,SpheresFactory,"Circular geometry of the SpheresFactory region. It can be disk (given by radius and center), or cylinder (given by radius, length and center).",
		((Real,radius,NaN,,"Radius of the region"))
		((Real,length,0,,"Length of the cylindrical region (0 by default)"))
		((Vector3r,center,Vector3r(NaN,NaN,NaN),,"Center of the region"))
	);
};
REGISTER_SERIALIZABLE(CircularFactory);

class BoxFactory: public SpheresFactory {
	protected:
		virtual void pickRandomPosition(Vector3r&, Real);
	public:
		virtual ~BoxFactory(){};
		DECLARE_LOGGER;
		SUDODEM_CLASS_BASE_DOC_ATTRS(BoxFactory,SpheresFactory,"Box geometry of the SpheresFactory region, given by extents and center",
		((Vector3r,extents,Vector3r(NaN,NaN,NaN),,"Extents of the region"))
		((Vector3r,center,Vector3r(NaN,NaN,NaN),,"Center of the region"))
	);
};
REGISTER_SERIALIZABLE(BoxFactory);
