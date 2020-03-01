// code not yet for use (6/12/2009); if long here uselessly, delete.
//

// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#include<sudodem/core/Serializable.hpp>
/*! Storage for general 3d view settings.

Is saved along with simulation and passed to every call to render(...).
Contains more or less what used to be inside OpenGLRenderer.

*/
class GLConfig: public Serializable{

	Vector3r lightPos,bgColor;
	Body::id_t currSel;
	bool dof,id,bbox,geom,wire,intrGeom,intrPhys;
	int mask;
	bool scaleDisplacements,scaleRotations;
	Vector3r displacementScale; Real rotationScale;
	vector<Se3r> clipPlaneSe3;
	vector<int> clipPlaneActive; // should be bool, but serialization doesn't handle vector<bool>

	// not saved
	Vector3r highlightEmission0;
	Vector3r highlightEmission1;

	// normalized saw signal with given periodicity, with values ∈ 〈0,1〉 */
	Real normSaw(Real t, Real period){ Real xi=(t-period*((int)(t/period)))/period; /* normalized value, (0-1〉 */ return (xi<.5?2*xi:2-2*xi); }
	Real normSquare(Real t, Real period){ Real xi=(t-period*((int)(t/period)))/period; /* normalized value, (0-1〉 */ return (xi<.5?0:1); }

	//! wrap number to interval x0…x1
	Real wrapCell(const Real x, const Real x0, const Real x1);
	//! wrap point to inside Scene's cell (identity if !Scene::isPeriodic)
	Vector3r wrapCellPt(const Vector3r& pt, Scene* rb);
	void drawPeriodicCell(Scene*);

	REGISTER_ATTRIBUTES(Serializable,(dof)(id)(bbox)(geom)(wire)(intrGeom)(intrPhys)(mask)(scaleDisplacements)(scaleRotations)(displacementScale)(rotationScale)(clipPlaneSe3)(clipPlaneActive));
};
