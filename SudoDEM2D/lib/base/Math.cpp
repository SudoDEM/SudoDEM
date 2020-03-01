#include<sudodem/lib/base/Math.hpp>
template<> const Real Math<Real>::EPSILON = DBL_EPSILON;
template<> const Real Math<Real>::ZERO_TOLERANCE = 1e-20;
template<> const Real Math<Real>::MAX_REAL = DBL_MAX;
template<> const Real Math<Real>::PI = 4.0*atan(1.0);
template<> const Real Math<Real>::TWO_PI = 2.0*Math<Real>::PI;
template<> const Real Math<Real>::HALF_PI = 0.5*Math<Real>::PI;
template<> const Real Math<Real>::DEG_TO_RAD = Math<Real>::PI/180.0;
template<> const Real Math<Real>::RAD_TO_DEG = 180.0/Math<Real>::PI;

template<> int ZeroInitializer<int>(){ return (int)0; }
template<> Real ZeroInitializer<Real>(){ return (Real)0; }

#ifdef SUDODEM_MASK_ARBITRARY
bool operator==(const mask_t& g, int i) { return g == mask_t(i); }
bool operator==(int i, const mask_t& g) { return g == i; }
bool operator!=(const mask_t& g, int i) { return !(g == i); }
bool operator!=(int i, const mask_t& g) { return g != i; }
mask_t operator&(const mask_t& g, int i) { return g & mask_t(i); }
mask_t operator&(int i, const mask_t& g) { return g & i; }
mask_t operator|(const mask_t& g, int i) { return g | mask_t(i); }
mask_t operator|(int i, const mask_t& g) { return g | i; }
bool operator||(const mask_t& g, bool b) { return (g != 0) || b; }
bool operator||(bool b, const mask_t& g) { return g || b; }
bool operator&&(const mask_t& g, bool b) { return (g != 0) && b; }
bool operator&&(bool b, const mask_t& g) { return g && b; }
#endif

bool MatrixXr_pseudoInverse(const MatrixXr &a, MatrixXr &a_pinv, double epsilon){

	// see : http://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse#The_general_case_and_the_SVD_method
	if ( a.rows()<a.cols() ) return false;

	// SVD
	Eigen::JacobiSVD<MatrixXr> svdA;
	svdA.compute(a,Eigen::ComputeThinU|Eigen::ComputeThinV);
	MatrixXr vSingular = svdA.singularValues();

	// Build a diagonal matrix with the Inverted Singular values
	// The pseudo inverted singular matrix is easy to compute :
	// is formed by replacing every nonzero entry by its reciprocal (inversing).
	VectorXr vPseudoInvertedSingular(svdA.matrixV().cols(),1);
		
	for (int iRow =0; iRow<vSingular.rows(); iRow++){
		if(fabs(vSingular(iRow))<=epsilon) vPseudoInvertedSingular(iRow,0)=0.;
	   else vPseudoInvertedSingular(iRow,0)=1./vSingular(iRow);
	}

	// A little optimization here 
	MatrixXr mAdjointU = svdA.matrixU().adjoint().block(0,0,vSingular.rows(),svdA.matrixU().adjoint().cols());

	// Pseudo-Inversion : V * S * U'
	a_pinv = (svdA.matrixV() *  vPseudoInvertedSingular.asDiagonal()) * mAdjointU;

	return true;
}
