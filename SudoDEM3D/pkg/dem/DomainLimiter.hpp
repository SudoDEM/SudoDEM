
#include<sudodem/pkg/common/PeriodicEngines.hpp>
#include<sudodem/core/PartialEngine.hpp>
#include<sudodem/pkg/common/Sphere.hpp>

class DomainLimiter: public PeriodicEngine{
	public:
		virtual void action();
	SUDODEM_CLASS_BASE_DOC_ATTRS(DomainLimiter,PeriodicEngine,"Delete particles that are out of axis-aligned box given by *lo* and *hi*.",
		((Vector3r,lo,Vector3r(0,0,0),,"Lower corner of the domain."))
		((Vector3r,hi,Vector3r(0,0,0),,"Upper corner of the domain."))
		((long,nDeleted,0,Attr::readonly,"Cummulative number of particles deleted."))
		((Real,mDeleted,0,,"Mass of deleted particles."))
		((Real,vDeleted,0,,"Volume of deleted particles."))
		((int,mask,-1,,"If mask is defined, only particles with corresponding groupMask will be deleted."))
	);
};
REGISTER_SERIALIZABLE(DomainLimiter);

class LawTester: public PartialEngine{
	Body::id_t id1,id2; // shorthands for local use
	public:
		void init();
		virtual void action();
		void postLoad(LawTester&);
		void warnDeprec(const string& s1, const string& s2){ if(!warnedDeprecPtRot){ warnedDeprecPtRot=true; LOG_WARN("LawTester."<<s1<<" is deprecated, use LawTester."<<s2<<" instead.");} }
		Vector3r get_ptOurs(){ warnDeprec("ptOurs","uTest.head()"); return uTest.head<3>(); } Vector3r get_ptGeom(){ warnDeprec("ptGeom","uGeom.head()"); return uGeom.head<3>(); }
		Vector3r get_rotOurs(){ warnDeprec("rotOurs","uTest.tail()"); return uTest.tail<3>(); }  Vector3r get_rotGeom(){ warnDeprec("rotGeom","uGeom.tail()"); return uGeom.tail<3>(); }
	DECLARE_LOGGER;
	SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(LawTester,PartialEngine,"Prescribe and apply deformations of an interaction in terms of local normal and shear displacements and rotations (using either :yref:`disPpath<LawTester.disPath>` and :yref:`rotPath<LawTester.rotPath>` [or :yref:`path<LawTester.path>` in the future]). Supported :yref:`IGeom` types are :yref:`ScGeom`, :yref:`L3Geom` and :yref:`L6Geom`. \n\nSee :ysrc:`scripts/test/law-test.py` for an example.",
		((vector<Vector3r>,disPath,,Attr::triggerPostLoad,"Loading path, where each Vector3 contains desired normal displacement and two components of the shear displacement (in local coordinate system, which is being tracked automatically. If shorter than :yref:`rotPath<LawTester.rotPath>`, the last value is repeated."))
		((vector<Vector3r>,rotPath,,Attr::triggerPostLoad,"Rotational components of the loading path, where each item contains torsion and two bending rotations in local coordinates. If shorter than :yref:`path<LawTester.path>`, the last value is repeated."))
		((vector<string>,hooks,,,"Python commands to be run when the corresponding point in path is reached, before doing other things in that particular step. See also :yref:`doneHook<LawTester.doneHook>`. "))
		((Vector6r,uGeom,Vector6r::Zero(),,"Current generalized displacements (3 displacements, 3 rotations), as stored in the interation itself. They should corredpond to :yref:`uTest<LawTester.uTest>`, otherwise a bug is indicated."))
		((Vector6r,uTest,Vector6r::Zero(),,"Current generalized displacements (3 displacements, 3 rotations), as they should be according to this :yref:`LawTester`. Should correspond to :yref:`uGeom<LawTester.uGeom>`."))
		((Vector6r,uTestNext,Vector6r::Zero(),Attr::hidden,"The value of uTest in the next step; since uTest is computed before uGeom is updated (in the next time step), uTest and uGeom would be always shifted by one timestep."))
		((bool,warnedDeprecPtRot,false,Attr::hidden,"Flag to say that the user was already warned about using deprecated ptOurg/ptGeom/rotOurs/rotGeom."))
		((Vector3r,shearTot,Vector3r::Zero(),Attr::hidden,"Current displacement in global coordinates; only used internally with ScGeom."))
		((bool,displIsRel,true,,"Whether displacement values in *disPath* are normalized by reference contact length (r1+r2 for 2 spheres)."))
		//((Vector3i,forceControl,Vector3i::Zero(),,"Select which components of path (non-zero value) have force (stress) rather than displacement (strain) meaning."))
		((vector<int>,pathSteps,((void)"(constant step)",vector<int>(1,1)),Attr::triggerPostLoad,"Step number for corresponding values in :yref:`path<LawTester.path>`; if shorter than path, distance between last 2 values is used for the rest."))
		((vector<int>,_pathT,,(Attr::readonly|Attr::noSave),"Time value corresponding to points on path, computed from *pathSteps*. Length is always the same as path."))
		((vector<Vector6r>,_path,,(Attr::readonly|Attr::noSave),"Generalized displacement path values, computed from *disPath* and *rotPath* by appending zero initial value, and possibly repeating the last value to make lengths of *disPath* and *rotPath* match."))
		((shared_ptr<Interaction>,I,,(Attr::hidden),"Interaction object being tracked."))

		// axX, axY, axZ could be replaced by references to rows in trsf; perhaps do that in the future (in such a case, trsf would not be noSave)
		((Vector3r,axX,,Attr::hidden,"Local x-axis in global coordinates (normal of the contact) |yupdate|"))
		((Vector3r,axY,,Attr::hidden,"Local y-axis in global coordinates; perpendicular to axX; initialized arbitrarily, but tracked to be consistent. |yupdate|"))
		((Vector3r,axZ,,(Attr::hidden|Attr::noSave),"Local z-axis in global coordinates; computed from axX and axY. |yupdate|"))
		((Matrix3r,trsf,,(Attr::noSave|Attr::readonly),"Transformation matrix for the local coordinate system. |yupdate|"))
		((size_t,_interpPos,0,(Attr::readonly|Attr::hidden),"State cookie for the interpolation routine."))
		((Vector6r,uuPrev,Vector6r::Zero(),(Attr::readonly),"Generalized displacement values reached in the previous step, for knowing which increment to apply in the current step."))
		((int,step,1,,"Step number in which this engine is active; determines position in path, using pathSteps."))
		((string,doneHook,,,"Python command (as string) to run when end of the path is achieved. If empty, the engine will be set :yref:`dead<Engine.dead>`."))
		((Real,renderLength,0,,"Characteristic length for the purposes of rendering, set equal to the smaller radius."))
		((Real,refLength,0,(Attr::readonly),"Reference contact length, for rendering only."))
		((Vector3r,contPt,Vector3r::Zero(),Attr::hidden,"Contact point (for rendering only)"))
		((Real,idWeight,1,,"Float, usually ∈〈0,1〉, determining on how are displacements distributed between particles (0 for id1, 1 for id2); intermediate values will apply respective part to each of them. This parameter is ignored with 6-DoFs :yref:`IGeom`."))
		((Real,rotWeight,1,,"Float ∈〈0,1〉 determining whether shear displacement is applied as rotation or displacement on arc (0 is displacement-only, 1 is rotation-only). Not effective when mutual rotation is specified."))
		// reset force components along individual axes, instead of blocking DOFs which have no specific direction (for the force control)
		, /* deprec */ ((path,disPath,"LawTester.path will be used for generalized displacement (6-component) loading path in the future."))
		, /* init */
		, /* ctor */
		, /* py */ .add_property("ptOurs",&LawTester::get_ptOurs,"first 3 components of uTest |ydeprecated|") .add_property("ptGeom",&LawTester::get_ptGeom,"first 3 components of uGeom |ydeprecated|") .add_property("rotOurs",&LawTester::get_rotOurs,"last 3 components of uTest |ydeprecated|") .add_property("rotGeom",&LawTester::get_rotGeom,"last 3 components of uGeom |ydeprecated|")
	);
};
REGISTER_SERIALIZABLE(LawTester);

#ifdef SUDODEM_OPENGL
#include<sudodem/pkg/common/OpenGLRenderer.hpp>

class GlExtra_LawTester: public GlExtraDrawer{
	public:
	DECLARE_LOGGER;
	virtual void render();
	SUDODEM_CLASS_BASE_DOC_ATTRS(GlExtra_LawTester,GlExtraDrawer,"Find an instance of :yref:`LawTester` and show visually its data.",
		((shared_ptr<LawTester>,tester,,,"Associated :yref:`LawTester` object."))
	);
};
REGISTER_SERIALIZABLE(GlExtra_LawTester);

class GlExtra_OctreeCubes: public GlExtraDrawer{
	public:
	struct OctreeBox{ Vector3r center, extents; int fill; int level; };
	std::vector<OctreeBox> boxes;
	void postLoad(GlExtra_OctreeCubes&);
	virtual void render();
	SUDODEM_CLASS_BASE_DOC_ATTRS(GlExtra_OctreeCubes,GlExtraDrawer,"Render boxed read from file",
		((string,boxesFile,,Attr::triggerPostLoad,"File to read boxes from; ascii files with ``x0 y0 z0 x1 y1 z1 c`` records, where ``c`` is an integer specifying fill (0 for wire, 1 for filled)."))
		((Vector2i,fillRangeFill,Vector2i(2,2),,"Range of fill indices that will be filled."))
		((Vector2i,fillRangeDraw,Vector2i(-2,2),,"Range of fill indices that will be rendered."))
		((Vector2i,levelRangeDraw,Vector2i(-2,2),,"Range of levels that will be rendered."))
		((bool,noFillZero,true,,"Do not fill 0-fill boxed (those that are further subdivided)"))
	);
};
REGISTER_SERIALIZABLE(GlExtra_OctreeCubes);
#endif

