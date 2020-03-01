// 2007 © Václav Šmilauer <eudoxos@arcig.cz>

#pragma once

#include<boost/lambda/lambda.hpp>

#include "sudodem/lib/base/Math.hpp"
#include "sudodem/lib/base/Logging.hpp"
#include "sudodem/core/Body.hpp"

#include<boost/function.hpp>


class Scene;
class Body;
class SimpleViscoelasticBodyParameters;
class ViscElMat;
class FrictMat;
class Interaction;

using boost::shared_ptr;
namespace py = boost::python;

/*! Miscillaneous utility functions which are believed to be generally useful.
 *
 * All data members are methods are static, no instance of Shop is created. It is not serializable either.
 */

class Shop{
	public:
		DECLARE_LOGGER;
/*
		//! create default disk, along with its bound etc.
		static shared_ptr<Body> disk(Vector2r center, Real radius, shared_ptr<Material> mat);
		//! create default box with everything needed

		//! return instance of default FrictMat
		static shared_ptr<FrictMat> defaultGranularMat();

		//! Return vector of pairs (center,radius) loaded from a file with numbers inside
		static vector<boost::tuple<Vector3r,Real,int> > loadDisksFromFile(const string& fname,Vector3r& minXYZ, Vector3r& maxXYZ, Vector3r* cellSize=NULL);

		//! Save disks in the current simulation into a text file
		static void saveDisksToFile(string fileName);

		//! Compute the total volume of disks
		static Real getDisksVolume(const shared_ptr<Scene>& rb=shared_ptr<Scene>(), int mask=-1);

		//! Compute the total mass of disks
		static Real getDisksMass(const shared_ptr<Scene>& rb=shared_ptr<Scene>(), int mask=-1);

		//! Compute porosity; volume must be given for aperiodic simulations
		static Real getPorosity(const shared_ptr<Scene>& rb=shared_ptr<Scene>(),Real volume=-1);

		//! Compute porosity by dividing given volume into a grid of voxels;
		static Real getVoxelPorosity(const shared_ptr<Scene>& rb=shared_ptr<Scene>(),int resolution=500,Vector3r start=Vector3r(0,0,0),Vector3r end=Vector3r(0,0,0));

		//! Estimate timestep based on P-wave propagation speed
		static Real PWaveTimeStep(const shared_ptr<Scene> rb=shared_ptr<Scene>());

		//! Estimate timestep based on Rayleigh-wave propagation speed
		static Real RayleighWaveTimeStep(const shared_ptr<Scene> rb=shared_ptr<Scene>());

		//! return 2d coordinates of a 3d point within plane defined by rotation axis and inclination of spiral, wrapped to the 0th period
		static boost::tuple<Real, Real, Real> spiralProject(const Vector3r& pt, Real dH_dTheta, int axis=2, Real periodStart=std::numeric_limits<Real>::quiet_NaN(), Real theta0=0);

		//! Calculate inscribed circle center of trianlge
		static Vector3r inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2);

		/// Get viscoelastic parameters kn,cn,ks,cs from analytical solution of
		/// a problem of interaction of pair disks with mass 1, collision
		/// time tc and restitution coefficients en,es.
	    static void getViscoelasticFromDisksInteraction(Real tc, Real en, Real es, shared_ptr<ViscElMat> b);

		//! Get unbalanced force of the whole simulation
		static Real unbalancedForce(bool useMaxForce=false, Scene* _rb=NULL);
		static Real kineticEnergy(Scene* _rb=NULL, Body::id_t* maxId=NULL);
		//! get total momentum of current simulation
		static Vector3r momentum();
		//! get total angular momentum of current simulation
		//static Vector3r angularMomentum(Vector3r origin = Vector3r::Zero());

		//static Vector3r totalForceInVolume(Real& avgIsoStiffness, Scene *_rb=NULL);


		//! create transientInteraction between 2 bodies, using existing Dispatcher in Omega
		static shared_ptr<Interaction> createExplicitInteraction(Body::id_t id1, Body::id_t id2, bool force);

		//! apply force on contact point on both bodies (reversed on body 2)
		static void applyForceAtContactPoint(const Vector2r& force, const Vector2r& contPt, Body::id_t id1, const Vector2r& pos1, Body::id_t id2, const Vector2r& pos2, Scene* scene);

		//! map scalar variable to 1d colorscale
		static Vector3r scalarOnColorScale(Real x, Real xmin=0., Real xmax=1.);



		//! Flip cell shear without affecting interactions; if flip is zeros, it will be computed such that abs of shear strain is minimal for each shear component
		//! Diagonal terms of flip are meaningless and ignored.
		//static Matrix3r flipCell(const Matrix3r& flip=Matrix3r::Zero());

		//! Class for storing stresses, affected on bodies, obtained from Interactions
		struct bodyState{
				Vector2r normStress, shearStress;
				bodyState (){
					normStress = Vector2r(0.0,0.0);
					shearStress = Vector2r(0.0,0.0);
				}
		};
		//! Function of getting stresses for each body
		static void getStressForEachBody(vector<Shop::bodyState>&);

		//! Define the exact average stress in each particle from contour integral ("LW" stands for Love-Weber, since this is what the contour integral gives).
		//static void getStressLWForEachBody(vector<Matrix3r>& bStresses);
		//static py::list getStressLWForEachBody();

		//! Function to compute overall ("macroscopic") stress.
		static Matrix3r getStress(Real volume=0);
		static Matrix3r getCapillaryStress(Real volume=0);
		static Matrix3r stressTensorOfPeriodicCell() { LOG_WARN("Shop::stressTensorOfPeriodicCelli is DEPRECATED: use getStress instead"); return Shop::getStress(); }
		//! Compute overall ("macroscopic") stress of periodic cell, returning 2 tensors
		//! (contribution of normal and shear forces)
		static py::tuple normalShearStressTensors(bool compressionPositive=false, bool splitNormalTensor=false, Real thresholdForce=NaN);

		//! Function to compute fabric tensor of periodic cell
		static void fabricTensor(Real& Fmean, Matrix3r& fabric, Matrix3r& fabricStrong, Matrix3r& fabricWeak, bool splitTensor=false, bool revertSign=false, Real thresholdForce=NaN);
		static py::tuple fabricTensor(bool splitTensor=false, bool revertSign=false, Real thresholdForce=NaN);

		//! Function to set translational and rotational velocities of all bodies to zero
		static void calm(const shared_ptr<Scene>& rb=shared_ptr<Scene>(), int mask=-1);

		//! Get a list of body-ids, which contacts the given body;
		static py::list getBodyIdsContacts(Body::id_t bodyID=-1);

		//! Set material and contact friction to the given value, non-dynamic bodies are not affected
		static void setContactFriction(Real angleRad);

		//! Homothetic change of sizes of disks and clumps
		static void growParticles(Real multiplier, bool updateMass, bool dynamicOnly);
		//! Change of size of a single disk or a clump
		// DEPREC, update wrt growParticles()
		static void growParticle(Body::id_t bodyID, Real multiplier, bool updateMass);

*/
		//! wrap floating number periodically to the given range
		static Real periodicWrap(Real x, Real x0, Real x1, long* period=NULL);
		//! Get unbalanced force of the whole simulation
		static Real unbalancedForce(bool useMaxForce=false, Scene* _rb=NULL);
		static py::tuple getStressAndTangent2D(Real z_dim = 1, bool symmetry=true);
		static py::tuple getStressTangentThermal2D(Real z_dim = 1, bool symmetry=true);
};
