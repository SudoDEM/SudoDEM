// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

/*! Periodic cell parameters and routines. Usually instantiated as Scene::cell.

The Cell has current box configuration represented by the parallelepiped's base vectors (*hSize*). Lengths of the base vectors can be accessed via *size*.

* Matrix3r *trsf* is "deformation gradient tensor" F (http://en.wikipedia.org/wiki/Finite_strain_theory)
* Matrix3r *velGrad* is "velocity gradient tensor" (http://www.cs.otago.ac.nz/postgrads/alexis/FluidMech/node7.html)

The transformation is split between "normal" part and "rotation/shear" part for contact detection algorithms. The shearPt, unshearPt, getShearTrsf etc functions refer to both shear and rotation. This decomposition is frame-dependant and does not correspond to the rotation/stretch decomposition of mechanics (with stretch further decomposed into isotropic and deviatoric). Therefore, using the shearPt family in equations of mechanics is not recommended. Similarly, attributes assuming the existence of a "reference" state are considered deprecated (refSize, hSize0). It is better to not use them since there is no guarantee that the so-called "reference state" will be maintained consistently in the future.

*/

#pragma once
#include<sudodem/lib/base/Math.hpp>
#include<sudodem/lib/serialization/Serializable.hpp>
//#include<sudodem/lib/base/Math.hpp>

//#include<sudodem/pkg/dem/Shop.hpp>

class Cell: public Serializable{
	public:
	//! Get/set sizes of cell base vectors
	const Vector3r& getSize() const { return _size; }
	void setSize(const Vector3r& s){for (int k=0;k<3;k++) hSize.col(k)*=s[k]/hSize.col(k).norm(); refHSize=hSize;  postLoad(*this);}
	//! Return copy of the current size (used only by the python wrapper)
	Vector3r getSize_copy() const { return _size; }
	//! return vector of consines of skew angle in yz, xz, xy planes between respective transformed base vectors
	const Vector3r& getCos() const {return _cos;}
	//! transformation matrix applying pure shear&rotation (scaling removed)
	const Matrix3r& getShearTrsf() const { return _shearTrsf; }
	//! inverse of getShearTrsfMatrix().
	const Matrix3r& getUnshearTrsf() const {return _unshearTrsf;}
	//! transformation increment matrix applying arbitrary field (remove if not used in NewtonIntegrator! )
	// const Matrix3r& getTrsfInc() const { return _trsfInc; }

	/*! return pointer to column-major OpenGL 4x4 matrix applying pure shear. This matrix is suitable as argument for glMultMatrixd.

	Note: the order of OpenGL transoformations matters; for instance, if you draw sheared wire box of size *size*,
	centered at *center*, the order is:

		1. translation: glTranslatev(center);
		3. shearing: glMultMatrixd(scene->cell->getGlShearTrsfMatrix());
		2. scaling: glScalev(size);
		4. draw: glutWireCube(1);

	See also http://www.songho.ca/opengl/gl_transform.html#matrix .
	*/
	const double* getGlShearTrsfMatrix() const { return _glShearTrsfMatrix; }
	//! Whether any shear (non-diagonal) component of the strain matrix is nonzero.
	bool hasShear() const {return _hasShear; }

	// caches; private
	private:
		Matrix3r _invTrsf;
		Matrix3r _trsfInc;
		Matrix3r _vGradTimesPrevH;
		Vector3r _size, _cos;
		Vector3r _refSize;
		bool _hasShear;
		Matrix3r _shearTrsf, _unshearTrsf;
		double _glShearTrsfMatrix[16];
		void fillGlShearTrsfMatrix(double m[16]);
	public:

	DECLARE_LOGGER;

	//! "integrate" velGrad, update cached values used by public getter.
	void integrateAndUpdate(Real dt);
	/*! Return point inside periodic cell, even if shear is applied */
	Vector3r wrapShearedPt(const Vector3r& pt) const { return shearPt(wrapPt(unshearPt(pt))); }
	/*! Return point inside periodic cell, even if shear is applied; store cell coordinates in period. */
	Vector3r wrapShearedPt(const Vector3r& pt, Vector3i& period) const { return shearPt(wrapPt(unshearPt(pt),period)); }
	/*! Apply inverse shear on point; to put it inside (unsheared) periodic cell, apply wrapPt on the returned value. */
	Vector3r unshearPt(const Vector3r& pt) const { return _unshearTrsf*pt; }
	//! Apply shear on point.
	Vector3r shearPt(const Vector3r& pt) const { return _shearTrsf*pt; }
	/*! Wrap point to inside the periodic cell; don't compute number of periods wrapped */
	Vector3r wrapPt(const Vector3r& pt) const {
		Vector3r ret; for(int i=0;i<3;i++) ret[i]=wrapNum(pt[i],_size[i]); return ret;}
	/*! Wrap point to inside the periodic cell; period will contain by how many cells it was wrapped. */
	Vector3r wrapPt(const Vector3r& pt, Vector3i& period) const {
		Vector3r ret; for(int i=0; i<3; i++){ ret[i]=wrapNum(pt[i],_size[i],period[i]); } return ret;}
	/*! Wrap number to interval 0…sz */
	static Real wrapNum(const Real& x, const Real& sz) {
		Real norm=x/sz; return (norm-floor(norm))*sz;}
	/*! Wrap number to interval 0…sz; store how many intervals were wrapped in period */
	static Real wrapNum(const Real& x, const Real& sz, int& period) {
		Real norm=x/sz; period=(int)floor(norm); return (norm-period)*sz;}

	// relative position and velocity for interaction accross multiple cells
	Vector3r intrShiftPos(const Vector3i& cellDist) const { return hSize*cellDist.cast<Real>(); }
	Vector3r intrShiftVel(const Vector3i& cellDist) const { return _vGradTimesPrevH*cellDist.cast<Real>(); }
	// return body velocity while taking away mean field velocity (coming from velGrad) if the mean field velocity is applied on velocity
	Vector3r bodyFluctuationVel(const Vector3r& pos, const Vector3r& vel, const Matrix3r& prevVelGrad) const { return (vel-prevVelGrad*pos); }

	// get/set current shape; setting resets trsf to identity
	Matrix3r getHSize() const { return hSize; }
	void setHSize(const Matrix3r& m){ hSize=refHSize=m; postLoad(*this); }
	// set current transformation; has no influence on current configuration (hSize); sets display refHSize as side-effect
	Matrix3r getTrsf() const { return trsf; }
	void setTrsf(const Matrix3r& m){ trsf=m; postLoad(*this); }
	Matrix3r getVelGrad() const { return velGrad; }
	void setVelGrad(const Matrix3r& m){ nextVelGrad=m; velGradChanged=true;}
	//BEGIN Deprecated (see refSize property)
	// get undeformed shape
	Matrix3r getHSize0() const { return _invTrsf*hSize; }
	// edge lengths of the undeformed shape
	Vector3r getRefSize() const { Matrix3r h=getHSize0(); return Vector3r(h.col(0).norm(),h.col(1).norm(),h.col(2).norm()); }
	// temporary, will be removed in favor of more descriptive setBox(...)
	void setRefSize(const Vector3r& s){
		// if refSize is set to the current size and the cell is a box (used in older scripts), say it is not necessary
		Matrix3r hSizeEigen3=hSize.diagonal().asDiagonal();		//Eigen3 support
		if(s==_size && hSize==hSizeEigen3){ LOG_WARN("Setting O.cell.refSize=O.cell.size is useless, O.trsf=Matrix3.Identity is enough now."); }
		else {LOG_WARN("Setting Cell.refSize is deprecated, use Cell.setBox(...) instead.");}
		setBox(s); postLoad(*this);
	}
	//END Deprecated
	// set box shape of the cell
	void setBox(const Vector3r& size){ setHSize(size.asDiagonal()); trsf=Matrix3r::Identity(); postLoad(*this); }
	void setBox3(const Real& s0, const Real& s1, const Real& s2){ setBox(Vector3r(s0,s1,s2)); }


	// return current cell volume
	Real getVolume() const {return hSize.determinant();}
	void postLoad(Cell&){ integrateAndUpdate(0); }

	// to resolve overloads
	Vector3r wrapShearedPt_py(const Vector3r& pt) const { return wrapShearedPt(pt);}
	Vector3r wrapPt_py(const Vector3r& pt) const { return wrapPt(pt);}

	// strain measures
	Matrix3r getDefGrad() { return trsf; }
	Matrix3r getSmallStrain() { return .5*(trsf+trsf.transpose()) - Matrix3r::Identity(); }
	Matrix3r getRCauchyGreenDef() { return trsf.transpose()*trsf; }
	Matrix3r getLCauchyGreenDef() { return trsf*trsf.transpose(); }
	Matrix3r getLagrangianStrain() { return .5*(getRCauchyGreenDef()-Matrix3r::Identity()); }
	Matrix3r getEulerianAlmansiStrain() { return .5*(Matrix3r::Identity()-getLCauchyGreenDef().inverse()); }
	void computePolarDecOfDefGrad(Matrix3r& R, Matrix3r& U) { Matrix_computeUnitaryPositive(trsf,&R,&U); }
	boost::python::tuple getPolarDecOfDefGrad(){ Matrix3r R,U; computePolarDecOfDefGrad(R,U); return boost::python::make_tuple(R,U); }
	Matrix3r getRotation() { Matrix3r R,U; computePolarDecOfDefGrad(R,U); return R; }
	Matrix3r getLeftStretch() { Matrix3r R,U; computePolarDecOfDefGrad(R,U); return U; }
	Matrix3r getRightStretch() { Matrix3r R,U; computePolarDecOfDefGrad(R,U); return trsf*R.transpose(); }

	// stress measures
	//Matrix3r getStress() { return Shop::getStress(); }
	//Matrix3r getCauchyStress() { Matrix3r s=getStress(); return .5*(s+s.transpose()); }

	enum { HOMO_NONE=0, HOMO_POS=1, HOMO_VEL=2, HOMO_VEL_2ND=3 };
	SUDODEM_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(Cell,Serializable,"Parameters of periodic boundary conditions. Only applies if O.isPeriodic==True.",
		/* overridden below to be modified by getters/setters because of intended side-effects */
		((Matrix3r,trsf,Matrix3r::Identity(),,"[overridden]")) //"Current transformation matrix of the cell, which defines how far is the current cell geometry (:yref:`hSize<Cell.hSize>`) from the reference configuration. Changing trsf will not change :yref:`hSize<Cell.hSize>`, it serves only as accumulator for transformations applied via :yref:`velGrad<Cell.velGrad>`."))
		((Matrix3r,refHSize,Matrix3r::Identity(),,"Reference cell configuration, only used with :yref:`OpenGLRenderer.dispScale`. Updated automatically when :yref:`hSize<Cell.hSize>` or :yref:`trsf<Cell.trsf>` is assigned directly; also modified by :yref:`sudodem.utils.setRefSe3` (called e.g. by the ``Reference`` button in the UI)."))
		((Matrix3r,hSize,Matrix3r::Identity(),,"[overridden below]"))
		((Matrix3r,prevHSize,Matrix3r::Identity(),Attr::readonly,":yref:`hSize<Cell.hSize>` from the previous step, used in the definition of relative velocity across periods."))
		/* normal attributes */
		((Matrix3r,velGrad,Matrix3r::Zero(),,"[overridden below]"))
		((Matrix3r,nextVelGrad,Matrix3r::Zero(),Attr::readonly,"see :yref:`Cell.velGrad`."))
		((Matrix3r,prevVelGrad,Matrix3r::Zero(),Attr::readonly,"Velocity gradient in the previous step."))
		((bool,homoDeform,true,,"Deform (:yref:`velGrad<Cell.velGrad>`) the cell homothetically, by adjusting positions and velocities of bodies. The velocity change is obtained by deriving the expression v=∇v.x, where ∇v is the macroscopic velocity gradient, giving in an incremental form: Δv=Δ ∇v x + ∇v Δx. As a result, velocities are modified as soon as ``velGrad`` changes, according to the first term: Δv(t)=Δ ∇v x(t), while the 2nd term reflects a convective term: Δv'= ∇v v(t-dt/2)."))
		((bool,velGradChanged,false,Attr::readonly,"true when velGrad has been changed manually (see also :yref:`Cell.nextVelGrad`)")),
		/*deprec*/
		((Hsize,hSize,"conform to SudoDEM's names convention.")),
		/*init*/ ,
		/*ctor*/ _invTrsf=Matrix3r::Identity(); integrateAndUpdate(0),
		/*py*/
		// override some attributes above
		.add_property("hSize",&Cell::getHSize,&Cell::setHSize,"Base cell vectors (columns of the matrix), updated at every step from :yref:`velGrad<Cell.velGrad>` (:yref:`trsf<Cell.trsf>` accumulates applied :yref:`velGrad<Cell.velGrad>` transformations). Setting *hSize* during a simulation is not supported by most contact laws, it is only meant to be used at iteration 0 before any interactions have been created.")
		.add_property("size",&Cell::getSize_copy,&Cell::setSize,"Current size of the cell, i.e. lengths of the 3 cell lateral vectors contained in :yref:`Cell.hSize` columns. Updated automatically at every step. Assigning a value will change the lengths of base vectors (see :yref:`Cell.hSize`), keeping their orientations unchanged.")
		.add_property("refSize",&Cell::getRefSize,&Cell::setRefSize,"Reference size of the cell (lengths of initial cell vectors, i.e. column norms of :yref:`hSize<Cell.hSize>`).\n\n.. note::\n\t Modifying this value is deprecated, use :yref:`setBox<Cell.setBox>` instead.")
		// useful properties
		.add_property("trsf",&Cell::getTrsf,&Cell::setTrsf,"Current transformation matrix of the cell, obtained from time integration of :yref:`Cell.velGrad`.")
		.add_property("velGrad",&Cell::getVelGrad,&Cell::setVelGrad,"Velocity gradient of the transformation; used in :yref:`NewtonIntegrator`. Values of :yref:`velGrad<Cell.velGrad>` accumulate in :yref:`trsf<Cell.trsf>` at every step.\n\n NOTE: changing velGrad at the beginning of a simulation loop would lead to inacurate integration for one step, as it should normaly be changed after the contact laws (but before Newton). To avoid this problem, assignment is deferred automatically. The target value typed in terminal is actually stored in :yref:`Cell.nextVelGrad` and will be applied right in time by Newton integrator. \n\n.. note::\n\t Assigning individual components of velGrad is not possible (it will not return any error but it will have no effect). Instead, you can assign to :yref:`Cell.nextVelGrad`, as in O.cell.nextVelGrad[1,2]=1.")
		.def_readonly("size",&Cell::getSize_copy,"Current size of the cell, i.e. lengths of the 3 cell lateral vectors contained in :yref:`Cell.hSize` columns. Updated automatically at every step.")
		.add_property("volume",&Cell::getVolume,"Current volume of the cell.")
		// functions
		.def("setBox",&Cell::setBox,"Set :yref:`Cell` shape to be rectangular, with dimensions along axes specified by given argument. Shorthand for assigning diagonal matrix with respective entries to :yref:`hSize<Cell.hSize>`.")
		.def("setBox",&Cell::setBox3,"Set :yref:`Cell` shape to be rectangular, with dimensions along $x$, $y$, $z$ specified by arguments. Shorthand for assigning diagonal matrix with the respective entries to :yref:`hSize<Cell.hSize>`.")
		// debugging only
		.def("wrap",&Cell::wrapShearedPt_py,"Transform an arbitrary point into a point in the reference cell")
		.def("unshearPt",&Cell::unshearPt,"Apply inverse shear on the point (removes skew+rot of the cell)")
		.def("shearPt",&Cell::shearPt,"Apply shear (cell skew+rot) on the point")
		.def("wrapPt",&Cell::wrapPt_py,"Wrap point inside the reference cell, assuming the cell has no skew+rot.")
		.def("getDefGrad",&Cell::getDefGrad,"Returns deformation gradient tensor $\\mat{F}$ of the cell deformation (http://en.wikipedia.org/wiki/Finite_strain_theory)")
		.def("getSmallStrain",&Cell::getSmallStrain,"Returns small strain tensor $\\mat{\\varepsilon}=\\frac{1}{2}(\\mat{F}+\\mat{F}^T)-\\mat{I}$ of the cell (http://en.wikipedia.org/wiki/Finite_strain_theory)")
		.def("getRCauchyGreenDef",&Cell::getRCauchyGreenDef,"Returns right Cauchy-Green deformation tensor $\\mat{C}=\\mat{F}^T\\mat{F}$ of the cell (http://en.wikipedia.org/wiki/Finite_strain_theory)")
		.def("getLCauchyGreenDef",&Cell::getLCauchyGreenDef,"Returns left Cauchy-Green deformation tensor $\\mat{b}=\\mat{F}\\mat{F}^T$ of the cell (http://en.wikipedia.org/wiki/Finite_strain_theory)")
		.def("getLagrangianStrain",&Cell::getLagrangianStrain,"Returns Lagrangian strain tensor $\\mat{E}=\\frac{1}{2}(\\mat{C}-\\mat{I})=\\frac{1}{2}(\\mat{F}^T\\mat{F}-\\mat{I})=\\frac{1}{2}(\\mat{U}^2-\\mat{I})$ of the cell (http://en.wikipedia.org/wiki/Finite_strain_theory)")
		.def("getEulerianAlmansiStrain",&Cell::getEulerianAlmansiStrain,"Returns Eulerian-Almansi strain tensor $\\mat{e}=\\frac{1}{2}(\\mat{I}-\\mat{b}^{-1})=\\frac{1}{2}(\\mat{I}-(\\mat{F}\\mat{F}^T)^{-1})$ of the cell (http://en.wikipedia.org/wiki/Finite_strain_theory)")
		.def("getPolarDecOfDefGrad",&Cell::getPolarDecOfDefGrad,"Returns orthogonal matrix $\\mat{R}$ and symmetric positive semi-definite matrix $\\mat{U}$ as polar decomposition of deformation gradient $\\mat{F}$ of the cell ( $\\mat{F}=\\mat{RU}$ )")
		.def("getRotation",&Cell::getRotation,"Returns rotation of the cell (orthogonal matrix $\\mat{R}$ from polar decomposition $\\mat{F}=\\mat{RU}$ )")
		.def("getLeftStretch",&Cell::getLeftStretch,"Returns left (spatial) stretch tensor of the cell (matrix $\\mat{U}$ from polar decomposition $\\mat{F}=\\mat{RU}$ )")
		.def("getRightStretch",&Cell::getRightStretch,"Returns right (material) stretch tensor of the cell (matrix $\\mat{V}$ from polar decomposition $\\mat{F}=\\mat{RU}=\\mat{VR}\\ \\rightarrow\\ \\mat{V}=\\mat{FR}^T$ )")
		.def_readonly("shearTrsf",&Cell::_shearTrsf,"Current skew+rot transformation (no resize)")
		.def_readonly("unshearTrsf",&Cell::_unshearTrsf,"Inverse of the current skew+rot transformation (no resize)")
		.add_property("hSize0",&Cell::getHSize0,"Value of untransformed hSize, with respect to current :yref:`trsf<Cell::trsf>` (computed as :yref:`trsf<Cell::trsf>` ⁻¹ × :yref:`hSize<Cell::hSize>`.")
	);
};
REGISTER_SERIALIZABLE(Cell);
