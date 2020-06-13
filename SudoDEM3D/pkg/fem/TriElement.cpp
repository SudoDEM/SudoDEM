/*************************************************************************
*  Copyright (C) 2017 by Sway Zhao                                       *
*  zhswee@gmail.com                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include "TriElement.hpp"

#include<sudodem/pkg/common/InteractionLoop.hpp>

CREATE_LOGGER(TriElement);

TriElement::~TriElement()
{
}


void TriElement::postLoad(TriElement&)
{
	// if this fails, it means someone did vertices push_back, but they are resized to 3 at TriElement initialization already
	// in the future, a fixed-size array should be used instead of vector<Vector3r> for vertices
	// this is prevented by sudodem::serialization now IIRC
	if(vertices.size()!=3){ throw runtime_error(("TriElement must have exactly 3 vertices (not "+boost::lexical_cast<string>(vertices.size())+")").c_str()); }
	if(isnan(vertices[0][0])) return;  // not initialized, nothing to do
	Vector3r e[3] = {vertices[1]-vertices[0] ,vertices[2]-vertices[1] ,vertices[0]-vertices[2]};
	#define CHECK_EDGE(i) if(e[i].squaredNorm()==0){LOG_FATAL("TriElement has coincident vertices "<<i<<" ("<<vertices[i]<<") and "<<(i+1)%3<<" ("<<vertices[(i+1)%3]<<")!");}
		CHECK_EDGE(0); CHECK_EDGE(1);CHECK_EDGE(2);
	#undef CHECK_EDGE
	normal = e[0].cross(e[1]);
	area = .5*normal.norm();
	normal /= 2*area;
	for(int i=0; i<3; ++i){
		ne[i]=e[i].cross(normal); ne[i].normalize();
		vl[i]=vertices[i].norm();
		vu[i]=vertices[i]/vl[i];
	}
	Real p = e[0].norm()+e[1].norm()+e[2].norm();
	icr = e[0].norm()*ne[0].dot(e[2])/p;
}

SUDODEM_PLUGIN((TriElement));
/////////////
Vector3r TriElement::getNormal() const {
	return ((nodes[1]->pos-nodes[0]->pos).cross(nodes[2]->pos-nodes[0]->pos)).normalized();
}

Vector3r TriElement::getCentroid() const {
	return (1/3.)*(nodes[0]->pos+nodes[1]->pos+nodes[2]->pos);
}

Real TriElement::getArea() const {
	assert(numNodesOk());
	const Vector3r& A=nodes[0]->pos;
	const Vector3r& B=nodes[1]->pos;
	const Vector3r& C=nodes[2]->pos;
	return .5*((B-A).cross(C-A)).norm();
    //return triangleArea(nodes[0]->pos,nodes[1]->pos,nodes[2]->pos);
}

Real TriElement::getPerimeterSq() const {
	assert(numNodesOk());
	return (nodes[1]->pos-nodes[0]->pos).squaredNorm()+(nodes[2]->pos-nodes[1]->pos).squaredNorm()+(nodes[0]->pos-nodes[2]->pos).squaredNorm();

}

std::tuple<Vector3r,Vector3r,Vector3r> TriElement::getOuterVectors() const {
	assert(numNodesOk());
	// is not normalized
	Vector3r nn=(nodes[1]->pos-nodes[0]->pos).cross(nodes[2]->pos-nodes[0]->pos);
	return std::make_tuple((nodes[1]->pos-nodes[0]->pos).cross(nn),(nodes[2]->pos-nodes[1]->pos).cross(nn),(nodes[0]->pos-nodes[2]->pos).cross(nn));
}

bool TriElement::isPointInTriangle(const Vector3r& pt){
	//https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
	// u=P2−P1
    Vector3r u = nodes[1]->pos - nodes[0]->pos;
    // v=P3−P1
    Vector3r v = nodes[2]->pos - nodes[0]->pos;
    // n=u×v
    Vector3r n = u.cross(v);
    // w=P−P1
    Vector3r w = pt - nodes[0]->pos;
    // Barycentric coordinates of the projection P′of P onto T:
    // γ=[(u×w)⋅n]/n²
    float gamma = u.cross(w).dot(n) / n.dot(n);
	if (gamma > 1 || gamma < 0){//outside
		return false;
	}
    // β=[(w×v)⋅n]/n²
    float beta = w.cross(v).dot(n) / n.dot(n);
	if (beta > 1 || beta < 0){//outside
		return false;
	}
    float alpha = 1 - gamma - beta;
    // The point P′ lies inside T if:
    return ((0 <= gamma) && (gamma <= 1));
}

vector<Vector3r> TriElement::outerEdgeNormals() const{
	auto o=getOuterVectors();
	return {std::get<0>(o).normalized(),std::get<1>(o).normalized(),std::get<2>(o).normalized()};
}

#if 0
Vector3r TriElement::getNearestTrianglePt(const Vector3r& pt, const Vector3r[3] pts){
	const Vector3r& A(pts[0]); const Vector3r& B(pts[1]); const Vector3r& C(pts[2]);
	Vector3r n0=((B-A).cross(C-A));
	Vector3r n=n0.normalized();
	Vector3r outVec[3]={(B-A).cross(n0),(C-B).cross(n),(A-C).cross(n)};
	short w=0;
	for(int i:{0,1,2}) 
}
#endif

Vector3r TriElement::closestSegmentPt(const Vector3r& P, const Vector3r& A, const Vector3r& B, Real* normPos){
	Vector3r BA=B-A;
	if(unlikely(BA.squaredNorm()==0.)){ if(normPos) *normPos=0.; return A; }
	Real u=(P.dot(BA)-A.dot(BA))/(BA.squaredNorm());
	if(normPos) *normPos=u;
	return A+min(1.,max(0.,u))*BA;
}

Vector3r TriElement::closestSegmentPt(const Vector3r& P, const Vector3r& A, const Vector3r& B){
	Vector3r BA=B-A;
	Real u=(P.dot(BA)-A.dot(BA))/(BA.squaredNorm());
	return A+min(1.,max(0.,u))*BA;
}

Vector3r TriElement::getNearestPt(const Vector3r& pt) {
	// FIXME: actually no need to project, sign of the dot product will be the same regardless?!
	Vector3r fNormal=getNormal();
	Real planeDist=(pt-nodes[0]->pos).dot(fNormal);
	Vector3r fC=pt-planeDist*fNormal; // point projected to facet's plane
	Vector3r outVec[3];
	std::tie(outVec[0],outVec[1],outVec[2])=getOuterVectors();
	short w=0;
	for(int i:{0,1,2}) w|=(outVec[i].dot(fC-nodes[i]->pos)>0.?1:0)<<i;
	Vector3r contPt;
	switch(w){
		case 0: return fC; // ---: inside triangle
		case 1: return closestSegmentPt(fC,nodes[0]->pos,nodes[1]->pos); // +-- (n1)
		case 2: return closestSegmentPt(fC,nodes[1]->pos,nodes[2]->pos); // -+- (n2)
		case 4: return closestSegmentPt(fC,nodes[2]->pos,nodes[0]->pos); // --+ (n3)
		case 3: return nodes[1]->pos; // ++- (v1)
		case 5: return nodes[0]->pos; // +-+ (v0)
		case 6: return nodes[2]->pos; // -++ (v2)
		case 7: throw logic_error("TriElement::getNearestPt: Impossible sphere-facet intersection (all points are outside the edges). (please report bug)"); // +++ (impossible)
		default: throw logic_error("TriElement::getNearestPt: Nonsense intersection value. (please report bug)");
	}
}
//-----------------------------------------------------------------------------------------------------------
// http://en.wikipedia.org/wiki/Inertia_tensor_of_triangle
Matrix3r TriElement::triangleInertia(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2){
	Matrix3r V; V<<v0.transpose(),v1.transpose(),v2.transpose(); // rows!
	Real a=((v1-v0).cross(v2-v0)).norm(); // twice the triangle area
	Matrix3r S; S<<2,1,1, 1,2,1, 1,1,2; S*=(1/24.);
	Matrix3r C=a*V.transpose()*S*V;
	return Matrix3r::Identity()*C.trace()-C;
};
Real TriElement::triangleArea(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2){
	return (1/2.)*((v1-v0).cross(v2-v0)).norm();
}

Real TriElement::TetraVolume(const Vector3r& v){
    //given a point outside the triangle plane, v, which is selected as the base point
    //the volume can be obtained
    //(v1-v3).cross(v2-v3) is the normal (to the outside of the surface) by default.
    //the volume is positive when (v3-v) is (approx.) along the normal.
    //Volume with a sign is used to calculate the specimen volume after deformation.
    return (nodes[2]->pos-v).dot((nodes[0]->pos-nodes[2]->pos).cross(nodes[1]->pos-nodes[2]->pos))/6.;	
}

void TriElement::lumpMassInertia(Real density){
	
	if(!(halfThick>0)) return;
	for(int i:{0,1,2}){
	    // midpoint local coords
	    Vector3r vv[2]; for(int j:{1,2}) vv[i-1]=.5*nodes[i]->ori.conjugate()*(nodes[(i+j)%3]->pos-nodes[i]->pos);
	    Vector3r C = nodes[i]->ori.conjugate()*(getCentroid()-nodes[i]->pos);

	    //nodes[i]->inertia += density*(2*halfThick)*triangleInertia(Vector3r::Zero(),vv[0],vv[1]);//inertial is vector3r but the right hand is a matrix. FIXME when needed.
	    //nodes[i]->inertia += density*(2*halfThick)*triangleInertia(vv[0],vv[1],C);
	    nodes[i]->mass += density*(2*halfThick)*triangleArea(Vector3r::Zero(),vv[0],vv[1]);
	    nodes[i]->mass += density*(2*halfThick)*triangleArea(vv[0],vv[1],C);
    }
}

//-----------------------------------------------------------------------------------------------------------
void TriElement::stepUpdate(Real dt, bool rotIncr){
	if(!hasRefConf()) setRefConf();
	updateNode();
	computeNodalDisplacements(dt,rotIncr);
}

void TriElement::setRefConf(){
	// set initial node position and orientation
	pos = getCentroid();
	// for orientation, there is a few choices which one to pick
	enum {INI_CS_SIMPLE=0,INI_CS_NODE0_AT_X};
	const int ini_cs=INI_CS_NODE0_AT_X;
	//const int ini_cs=INI_CS_SIMPLE;
	switch(ini_cs){
		case INI_CS_SIMPLE:
			ori.setFromTwoVectors(Vector3r::UnitZ(),getNormal());
			break;
		case INI_CS_NODE0_AT_X:{
			Matrix3r T;
			T.col(0)=(nodes[0]->pos-pos).normalized();
			T.col(2)=this->getNormal();
			T.col(1)=T.col(2).cross(T.col(0));
			assert(T.col(0).dot(T.col(2))<1e-12);
			ori=Quaternionr(T);
			break;
		}
	};
	// reference nodal positions
	for(int i:{0,1,2}){
		Vector3r nl = glob2loc(nodes[i]->pos);
		assert(nl[2]<1e-6*(max(abs(nl[0]),abs(nl[1])))); // z-coord should be zero
		refPos.segment<2>(2*i)=nl.head<2>();
	}
	// reference nodal rotations
	refRot.resize(3);
	for(int i:{0,1,2}){
		// facet node orientation minus vertex node orientation, in local frame (read backwards)
		refRot[i]=nodes[i]->ori.conjugate()*ori;
		//LOG_WARN("refRot["<<i<<"]="<<AngleAxisr(refRot[i]).angle()<<"/"<<AngleAxisr(refRot[i]).axis().transpose());
	};
	// set displacements to zero
	uXy=phiXy=Vector6r::Zero();
	#ifdef MEMBRANE_DEBUG_ROT
		drill=Vector3r::Zero();
		currRot.resize(3);
	#endif
	// delete stiffness matrices to force their re-creating
	KKcst.resize(0,0); 
	KKdkt.resize(0,0);
};

void TriElement::updateNode(){
	assert(hasRefConf());
	pos = getCentroid();//duplicated!!
	// temporary orientation, just to be planar
	// see http://www.colorado.edu/engineering/cas/courses.d/NFEM.d/NFEM.AppC.d/NFEM.AppC.pdf
	Quaternionr ori0; ori0.setFromTwoVectors(Vector3r::UnitZ(),getNormal());
	Vector6r nxy0;
	for(int i:{0,1,2}){
        vertices[i] = nodes[i]->pos- pos;
		Vector3r xy0=ori0.conjugate()*vertices[i];
		#ifdef MEMBRANE_DEBUG_ROT
			if(xy0[2]>1e-5*(max(abs(xy0[0]),abs(xy0[1])))){
				LOG_ERROR("z-coordinate is not zero for node "<<i<<": ori0="<<AngleAxisr(ori0).axis()<<" !"<<AngleAxisr(ori0).angle()<<", xy0="<<xy0<<", position in global-oriented centroid-origin CS "<<(nodes[i]->pos-pos).transpose());
			}
		#else
			assert(xy0[2]<1e-5*(max(abs(xy0[0]),abs(xy0[1])))); // z-coord should be zero
		#endif
		nxy0.segment<2>(2*i)=xy0.head<2>();
	}
	// compute the best fit (C.8 of the paper)
	// split the fraction into numerator (y) and denominator (x) to get the correct quadrant
	Real tanTheta3_y=refPos[0]*nxy0[1]+refPos[2]*nxy0[3]+refPos[4]*nxy0[5]-refPos[1]*nxy0[0]-refPos[3]*nxy0[2]-refPos[5]*nxy0[4];
	Real tanTheta3_x=refPos.dot(nxy0);
	// rotation to be planar, plus rotation around plane normal to the element CR frame (read backwards)
	ori=ori0*AngleAxisr(atan2(tanTheta3_y,tanTheta3_x),Vector3r::UnitZ());

    //update some geometric measures of TriElement, and the following info would be stored for contact detection.
    Vector3r e[3] = {vertices[1]-vertices[0] ,vertices[2]-vertices[1] ,vertices[0]-vertices[2]};
	#define CHECK_EDGE(i) if(e[i].squaredNorm()==0){LOG_FATAL("Facet has coincident vertices "<<i<<" ("<<vertices[i]<<") and "<<(i+1)%3<<" ("<<vertices[(i+1)%3]<<")!");}
		CHECK_EDGE(0); CHECK_EDGE(1);CHECK_EDGE(2);
	#undef CHECK_EDGE
	normal = e[0].cross(e[1]);
	area = .5*normal.norm();
	normal /= 2*area;
	for(int i=0; i<3; ++i){
		ne[i]=e[i].cross(normal); ne[i].normalize();
		vl[i]=vertices[i].norm();
		vu[i]=vertices[i]/vl[i];
	}
	Real p = e[0].norm()+e[1].norm()+e[2].norm();
	icr = e[0].norm()*ne[0].dot(e[2])/p;
    
};

void TriElement::computeNodalDisplacements(Real dt, bool rotIncr){
	assert(hasRefConf());
	// supposes node is updated already
	for(int i:{0,1,2}){
		Vector3r xy = glob2loc(nodes[i]->pos);
		// relative tolerance of 1e-6 was too little, in some cases?!
		#ifdef MEMBRANE_DEBUG_ROT
			if(xy[2]>1e-5*(max(abs(xy[0]),abs(xy[1])))){
				LOG_ERROR("local z-coordinate is not zero for node "<<i<<": ori="<<AngleAxisr(ori).axis()<<" !"<<AngleAxisr(ori).angle()<<", xy="<<xy<<", position in global-oriented centroid-origin CS "<<(nodes[i]->pos-pos).transpose());
			}
		#else
			assert(xy[2]<1e-5*(max(abs(xy[0]),abs(xy[1]))));
		#endif
		// displacements
		uXy.segment<2>(2*i)=xy.head<2>()-refPos.segment<2>(2*i);
		// rotations
		if(rotIncr){
			// incremental
			Vector3r angVelL= glob2loc(nodes[i]->angVel); // angular velocity in element coords
			phiXy.segment<2>(2*i)-=dt*angVelL.head<2>();
			// throw std::runtime_error("Incremental rotation (In2_ElastMat_TriElement.rotIncr) is not yet implemented properly.");
		} else {
			// from total rotation difference
			AngleAxisr aa(refRot[i].conjugate()*(nodes[i]->ori.conjugate()*ori));
			/* aa sometimes gives angle close to 2π for rotations which are actually small and negative;
				I posted that question at http://forum.kde.org/viewtopic.php?f=74&t=110854 .
				Such a result is fixed by conditionally subtracting 2π:
			*/
			if(aa.angle()>M_PI) aa.angle()-=2*M_PI;
			Vector3r rot=Vector3r(aa.angle()*aa.axis()); // rotation vector in local coords
			// if(aa.angle()>3)
			/*if(rot.head<2>().squaredNorm()>3.1*3.1){ LOG_WARN("TriElement's in-plane rotation in a node is > 3.1 radians, expect unstability!");
                //output debug info
                AngleAxisr rr(refRot[i]);
				AngleAxisr cr(nodes[i]->ori.conjugate()*ori);
                std::cout<<"node "<<i<<" gid"<<nodes[i]->id<<"\n   refRot : "<<rr.axis()<<" !"<<rr.angle()<<"\n   currRot: "<<cr.axis()<<" !"<<cr.angle()<<"\n   diffRot: "<<aa.axis()<<" !"<<aa.angle()<<std::endl;
            }*/
			phiXy.segment<2>(2*i)=rot.head<2>(); // drilling rotation discarded
			#ifdef MEMBRANE_DEBUG_ROT
				AngleAxisr rr(refRot[i]);
				AngleAxisr cr(nodes[i]->ori.conjugate()*ori);
				LOG_TRACE("node "<<i<<"\n   refRot : "<<rr.axis()<<" !"<<rr.angle()<<"\n   currRot: "<<cr.axis()<<" !"<<cr.angle()<<"\n   diffRot: "<<aa.axis()<<" !"<<aa.angle());
				drill[i]=rot[2];
				currRot[i]=nodes[i]->ori.conjugate()*ori;
			#endif
		}
	}
};


void TriElement::ensureStiffnessMatrices(const Real& young, const Real& nu, const Real& thickness,bool bending, const Real& bendThickness){
	assert(hasRefConf());
	// do nothing if both matrices exist already
	if(KKcst.size()==36 && (!bending ||
		#ifdef MEMBRANE_CONDENSE_DKT		
			KKdkt.size()==54
		#else 
			KKdkt.size()==81
		#endif
	)) return; 
	// check thickness
	Real t=(isnan(thickness)?(2*halfThick):thickness);
	if(/*also covers NaN*/!(t>0)) throw std::runtime_error("TriElement::ensureStiffnessMatrices: Element thickness is not positive!");

	// plane stress stiffness matrix
	Matrix3r E;
	E<<1,nu,0, nu,1,0,0,0,(1-nu)/2;
	E*=young/(1-pow(nu,2));

	// strain-displacement matrix (CCT element)
	Real area = getArea();
	const Real& x1(refPos[0]); const Real& y1(refPos[1]); const Real& x2(refPos[2]); const Real& y2(refPos[3]); const Real& x3(refPos[4]); const Real& y3(refPos[5]);

	// Felippa: Introduction to FEM eq. (15.17), pg. 259
	Eigen::Matrix<Real,3,6> B;
	B<<
		(y2-y3),0      ,(y3-y1),0      ,(y1-y2),0      ,
		0      ,(x3-x2),0      ,(x1-x3),0      ,(x2-x1),
		(x3-x2),(y2-y3),(x1-x3),(y3-y1),(x2-x1),(y1-y2);
	B*=1/(2*area);
	KKcst.resize(6,6);
	KKcst=t*area*B.transpose()*E*B;

	// EBcst=E*B;

	if(!bending) return;

	// strain-displacement matrix (DKT element)

	Real dktT=(isnan(bendThickness)?t:bendThickness);
	assert(!isnan(dktT));

	// [Batoz,Bathe,Ho: 1980], Appendix A: Expansions for the DKT element
	// (using the same notation)
	Real x23=(x2-x3), x31(x3-x1), x12(x1-x2);
	Real y23=(y2-y3), y31(y3-y1), y12(y1-y2);
	Real l23_2=(pow(x23,2)+pow(y23,2)), l31_2=(pow(x31,2)+pow(y31,2)), l12_2=(pow(x12,2)+pow(y12,2));
	Real P4=-6*x23/l23_2, P5=-6*x31/l31_2, P6=-6*x12/l12_2;
	Real t4=-6*y23/l23_2, t5=-6*y31/l31_2, t6=-6*y12/l12_2;
	Real q4=3*x23*y23/l23_2, q5=3*x31*y31/l31_2, q6=3*x12*y12/l12_2;
	Real r4=3*pow(y23,2)/l23_2, r5=3*pow(y31,2)/l31_2, r6=3*pow(y12,2)/l12_2;

	// even though w1, w2, w3 are always zero due to our definition of triangle plane,
	// we need to get the corresponding nodal force, therefore cannot condense those DOFs away.

	auto Hx_xi=[&](const Real& xi, const Real& eta) -> Vector9r { Vector9r ret; ret<<
		P6*(1-2*xi)+(P5-P6)*eta,
		q6*(1-2*xi)-(q5+q6)*eta,
		-4+6*(xi+eta)+r6*(1-2*xi)-eta*(r5+r6),
		-P6*(1-2*xi)+eta*(P4+P6),
		q6*(1-2*xi)-eta*(q6-q4),
		-2+6*xi+r6*(1-2*xi)+eta*(r4-r6),
		-eta*(P5+P4),
		eta*(q4-q5),
		-eta*(r5-r4);
		return ret;
	};

	auto Hy_xi=[&](const Real& xi, const Real &eta) -> Vector9r { Vector9r ret; ret<<
		t6*(1-2*xi)+eta*(t5-t6),
		1+r6*(1-2*xi)-eta*(r5+r6),
		-q6*(1-2*xi)+eta*(q5+q6),
		-t6*(1-2*xi)+eta*(t4+t6),
		-1+r6*(1-2*xi)+eta*(r4-r6),
		-q6*(1-2*xi)-eta*(q4-q6),
		-eta*(t4+t5),
		eta*(r4-r5),
		-eta*(q4-q5);
		return ret;
	};

	auto Hx_eta=[&](const Real& xi, const Real& eta) -> Vector9r { Vector9r ret; ret<<
		-P5*(1-2*eta)-xi*(P6-P5),
		q5*(1-2*eta)-xi*(q5+q6),
		-4+6*(xi+eta)+r5*(1-2*eta)-xi*(r5+r6),
		xi*(P4+P6),
		xi*(q4-q6),
		-xi*(r6-r4),
		P5*(1-2*eta)-xi*(P4+P5),
		q5*(1-2*eta)+xi*(q4-q5),
		-2+6*eta+r5*(1-2*eta)+xi*(r4-r5);
		return ret;
	};

	auto Hy_eta=[&](const Real& xi, const Real& eta) -> Vector9r { Vector9r ret; ret<<
		-t5*(1-2*eta)-xi*(t6-t5),
		1+r5*(1-2*eta)-xi*(r5+r6),
		-q5*(1-2*eta)+xi*(q5+q6),
		xi*(t4+t6),
		xi*(r4-r6),
		-xi*(q4-q6),
		t5*(1-2*eta)-xi*(t4+t5),
		-1+r5*(1-2*eta)+xi*(r4-r5),
		-q5*(1-2*eta)-xi*(q4-q5);
		return ret;
	};

	// return B(xi,eta) [Batoz,Bathe,Ho,1980], (30)
	auto Bphi=[&](const Real& xi, const Real& eta) -> Eigen::Matrix<Real,3,9> {
		Eigen::Matrix<Real,3,9> ret;
		ret<<
			(y31*Hx_xi(xi,eta)+y12*Hx_eta(xi,eta)).transpose(),
			(-x31*Hy_xi(xi,eta)-x12*Hy_eta(xi,eta)).transpose(),
			(-x31*Hx_xi(xi,eta)-x12*Hx_eta(xi,eta)+y31*Hy_xi(xi,eta)+y12*Hy_eta(xi,eta)).transpose()
		;
		return ret/(2*area);
	};

	// bending elasticity matrix (the same as (t^3/12)*E, E being the plane stress matrix above)
	Matrix3r Db;
	Db<<1,nu,0, nu,1,0, 0,0,(1-nu)/2;
	Db*=young*pow(dktT,3)/(12*(1-pow(nu,2)));


	// assemble the matrix here
	// KKdkt0 is 9x9, then w_i dofs are condensed away, and KKdkt is only 9x6
	#ifdef MEMBRANE_CONDENSE_DKT
		MatrixXr KKdkt0(9,9); KKdkt0.setZero();
		// MatrixXr EBdkt0(9,6); EBdkt0.setZero();
	#else
		KKdkt.setZero(9,9);
		// EBdkt.setZero(9,6);
	#endif
	// gauss integration points and their weights
	Vector3r xxi(.5,.5,0), eeta(0,.5,.5);
	Real ww[]={1/3.,1/3.,1/3.};
	for(int i:{0,1,2}){
		for(int j:{0,1,2}){
			auto b=Bphi(xxi[i],eeta[j]); // evaluate B at the gauss point
			#ifdef MEMBRANE_CONDENSE_DKT
				KKdkt0
			#else
				KKdkt
			#endif
				+=(2*area*ww[j]*ww[i])*b.transpose()*Db*b;
			#if 0
				#ifdef MEMBRANE_CONDENSE_DKT
					EBdkt0
				#else
					EBdkt
				#endif
					+=(2*area*ww[j]*ww[i])*Db*b;
			#endif
		}
	}
	#ifdef MEMBRANE_CONDENSE_DKT
		// extract columns [_ 1 2 _ 4 5 _ 7 8]
		KKdkt.setZero(9,6);
		for(int i:{0,1,2}){
			KKdkt.col(2*i)=KKdkt0.col(3*i+1);
			KKdkt.col(2*i+1)=KKdkt0.col(3*i+2);
		}
	#endif
};

/*
void TriElement::addIntraStiffnesses(const shared_ptr<Node>& n,Vector3r& ktrans, Vector3r& krot) const {
	if(!hasRefConf()) return;
	int i=-1;
	for(int j=0; j<3; j++){ if(nodes[j].get()==n.get()){ i=j; break; }}
	if(i<0) throw std::logic_error("TriElement::addIntraStiffness:: node "+n->pyStr()+" not found within nodes of "+this->pyStr()+".");
	if(KKcst.size()==0) return;
	bool dkt=(KKdkt.size()>0);
	// local translational stiffness diagonal
	Vector3r ktl(KKcst(2*i,2*i),KKcst(2*i+1,2*i+1),dkt?abs(KKdkt(3*i)):0);
	ktrans+= ori*ktl;
	// local rotational stiffness, if needed
	if(dkt){
		#ifdef MEMBRANE_CONDENSE_DKT
			Vector3r krl(abs(KKdkt(3*i,2*i)),abs(KKdkt(3*i+1,2*i+1)),0);
		#else
			Vector3r krl(abs(KKdkt(3*i,3*i)),abs(KKdkt(3*i+1,3*i+1)),0);
		#endif
		krot+= ori*krl;
	}
}

*/

/////////////
SUDODEM_PLUGIN((Bo1_TriElement_Aabb));

void Bo1_TriElement_Aabb::go(	  const shared_ptr<Shape>& cm
				, shared_ptr<Bound>& bv
				, const Se3r& se3
				, const Body* b)
{
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
	TriElement* facet = static_cast<TriElement*>(cm.get());
	const Vector3r& O = facet->pos;//se3.position;
    //std::cout<<"trielement position:"<<facet->pos<<std::endl;
	Matrix3r facetAxisT=facet->ori.toRotationMatrix();
	const vector<Vector3r>& vertices=facet->vertices;
	/*if(!scene->isPeriodic){
		aabb->min=aabb->max = O + facetAxisT * vertices[0];
		for (int i=1;i<3;++i)
		{
			Vector3r v = O + facetAxisT * vertices[i];
			aabb->min = aabb->min.cwiseMin(v);
			aabb->max = aabb->max.cwiseMax(v);
		}
	} else {*/
		Real inf=std::numeric_limits<Real>::infinity();
		aabb->min=Vector3r(inf,inf,inf); aabb->max=Vector3r(-inf,-inf,-inf);
		for(int i=0; i<3; i++){
			//Vector3r v=scene->cell->unshearPt(/*O+*/facetAxisT*facet->nodes[i]->pos);
            Vector3r v=scene->cell->unshearPt(/*O+*/facet->nodes[i]->pos);
			aabb->min=aabb->min.cwiseMin(v);
			aabb->max=aabb->max.cwiseMax(v);
		}
	//}
}
//////////////

#ifdef SUDODEM_OPENGL

#include <sudodem/lib/opengl/OpenGLWrapper.hpp>
#include <sudodem/core/Scene.hpp>

SUDODEM_PLUGIN((Gl1_TriElement));


bool Gl1_TriElement::normals=false;
bool Gl1_TriElement::wire=false;

void Gl1_TriElement::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& ,bool wire2,const GLViewInfo&)
{
	TriElement* facet = static_cast<TriElement*>(cm.get());
	const vector<Vector3r>& vertices = facet->vertices;
    //const Vector3r& v1 = facet->node1->pos;
    //const Vector3r& v2 = facet->node2->pos;
    //const Vector3r& v3 = facet->node3->pos;
	const Vector3r* ne = facet->ne;
	const Real& icr = facet->icr;
    Vector3r normal = facet->normal; normal.normalize();

	if(/*cm->wire || */wire){
		// facet
		glBegin(GL_LINE_LOOP);
			glColor3v(normals ? Vector3r(1,0,0): cm->color);
		   glVertex3v(vertices[0]);
		   glVertex3v(vertices[1]);
		   glVertex3v(vertices[2]);
	    glEnd();
		if(!normals) return;
		// facet's normal
		glBegin(GL_LINES);
			glColor3(0.0,0.0,1.0);
			glVertex3(0.0,0.0,0.0);
			glVertex3v(Vector3r(normal*icr));
		glEnd();
		// normal of edges
		glColor3(0.0,0.0,1.0);
		glBegin(GL_LINES);
			glVertex3(0.0,0.0,0.0); glVertex3v(Vector3r(icr*ne[0]));
			glVertex3(0.0,0.0,0.0);	glVertex3v(Vector3r(icr*ne[1]));
			glVertex3(0.0,0.0,0.0);	glVertex3v(Vector3r(icr*ne[2]));
		glEnd();
	} else {
		glDisable(GL_CULL_FACE);
		//Vector3r normal=(facet->vertices[1]-facet->vertices[0]).cross(facet->vertices[2]-facet->vertices[1]); normal.normalize();
        //Vector3r normal=(v2-v1).cross(v3-v2); normal.normalize();
		glColor3v(cm->color);
		glBegin(GL_TRIANGLES);
			glNormal3v(normal); // this makes every triangle different WRT the light direction; important!
			glVertex3v(facet->vertices[0]);
			glVertex3v(facet->vertices[1]);
			glVertex3v(facet->vertices[2]);
		glEnd();
	}
}


#endif /* SUDODEM_OPENGL */
