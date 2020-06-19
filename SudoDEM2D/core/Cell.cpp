
#include<sudodem/core/Cell.hpp>

CREATE_LOGGER(Cell);

void Cell::integrateAndUpdate(Real dt){
	//incremental displacement gradient
	_trsfInc=dt*velGrad;
	// total transformation; M = (Id+G).M = F.M
	trsf+=_trsfInc*trsf;
	_invTrsf=trsf.inverse();
	// hSize contains colums with updated base vectors
	prevHSize=hSize;
	_vGradTimesPrevH = velGrad*prevHSize;
	hSize+=_trsfInc*hSize;
	if(hSize.determinant()==0){ throw runtime_error("Cell is degenerate (zero volume)."); }
	// lengths of transformed cell vectors, skew cosines
	Matrix2r Hnorm; // normalized transformed base vectors
	for(int i=0; i<2; i++){
		Vector2r base(hSize.col(i));
		_size[i]=base.norm(); base/=_size[i]; //base is normalized now
		Hnorm(0,i)=base[0]; Hnorm(1,i)=base[1];};
	// skew cosines
		// sin between axes is cos of skew
		//_cos=(Hnorm.col(0).cross(Hnorm.col(1))).squaredNorm();
		//_cos=1 - pow(Hnorm.col(0).dot(Hnorm.col(1)),2);
		_cos=sqrt(1.0 - pow(Hnorm.col(0).dot(Hnorm.col(1)),2));//sin of the angle between the two basises
	// pure shear trsf: ones on diagonal
	_shearTrsf=Hnorm;
	// pure unshear transformation
	_unshearTrsf=_shearTrsf.inverse();
	// some parts can branch depending on presence/absence of shear
	_hasShear=(hSize(0,1)!=0 || hSize(1,0)!=0);
	// OpenGL shear matrix (used frequently)
	fillGlShearTrsfMatrix(_glShearTrsfMatrix);
}

void Cell::fillGlShearTrsfMatrix(double m[16]){
	m[0]=_shearTrsf(0,0); m[4]=_shearTrsf(0,1); m[8]=0; m[12]=0;
	m[1]=_shearTrsf(1,0); m[5]=_shearTrsf(1,1); m[9]=0; m[13]=0;
	m[2]=0;               m[6]=0;               m[10]=0;m[14]=0;
	m[3]=0;               m[7]=0;               m[11]=0;m[15]=1;
}
