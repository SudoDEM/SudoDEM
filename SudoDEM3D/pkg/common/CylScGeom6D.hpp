// 2011 © Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>
// 2012 © Kneib Francois <francois.kneib@irstea.fr>
#pragma once
#include<sudodem/pkg/dem/ScGeom.hpp>


class CylScGeom: public ScGeom {
public:
    /// Emulate a sphere whose position is the projection of sphere's center on cylinder sphere, and with motion linearly interpolated between nodes
    State fictiousState;
// 		shared_ptr<Interaction> duplicate;

    virtual ~CylScGeom () {};
    SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR(CylScGeom,ScGeom,"Geometry of a cylinder-sphere contact.",
                                   ((bool,onNode,false,,"contact on node?"))
                                   ((int,isDuplicate,0,,"this flag is turned true (1) automatically if the contact is shared between two chained cylinders. A duplicated interaction will be skipped once by the constitutive law, so that only one contact at a time is effective. If isDuplicate=2, it means one of the two duplicates has no longer geometric interaction, and should be erased by the constitutive laws."))
                                   ((int,trueInt,-1,,"Defines the body id of the cylinder where the contact is real, when :yref:`CylScGeom::isDuplicate`>0."))
                                   ((Vector3r,start,Vector3r::Zero(),,"position of 1st node |yupdate|"))
                                   ((Vector3r,end,Vector3r::Zero(),,"position of 2nd node |yupdate|"))
                                   ((Body::id_t,id3,0,,"id of next chained cylinder |yupdate|"))
                                   ((Real,relPos,0,,"position of the contact on the cylinder (0: node-, 1:node+) |yupdate|")),
                                   createIndex(); /*ctor*/
                                  );
    REGISTER_CLASS_INDEX(CylScGeom,ScGeom);
};
REGISTER_SERIALIZABLE(CylScGeom);


class CylScGeom6D: public ScGeom6D {
public:
    virtual ~CylScGeom6D() {};
    void precomputeRotations(const State& rbp1, const State& rbp2, bool isNew, bool creep=false) {
      initRotations(rbp1,rbp2);
    }
    void initRotations(const State& rbp1, const State& rbp2){
      initialOrientation1 = rbp1.ori;
      initialOrientation2 = rbp2.ori;
      twist=0;
      bending=Vector3r::Zero();
      twistCreep=Quaternionr(1.0,0.0,0.0,0.0);
    }
    State fictiousState;
    SUDODEM_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(CylScGeom6D,ScGeom6D,"Class representing :yref:`geometry<IGeom>` of two :yref:`bodies<Body>` in contact. The contact has 6 DOFs (normal, 2×shear, twist, 2xbending) and uses :yref:`ScGeom` incremental algorithm for updating shear.",
                                           ((bool,onNode,false,,"contact on node?"))
                                           ((int,isDuplicate,0,,"this flag is turned true (1) automatically if the contact is shared between two chained cylinders. A duplicated interaction will be skipped once by the constitutive law, so that only one contact at a time is effective. If isDuplicate=2, it means one of the two duplicates has no longer geometric interaction, and should be erased by the constitutive laws."))
                                           ((int,trueInt,-1,,"Defines the body id of the cylinder where the contact is real, when :yref:`CylScGeom::isDuplicate`>0."))
                                           ((Vector3r,start,Vector3r::Zero(),,"position of 1st node |yupdate|"))
                                           ((Vector3r,end,Vector3r::Zero(),,"position of 2nd node |yupdate|"))
                                           ((Body::id_t,id3,0,,"id of next chained cylinder |yupdate|"))
                                           ((Real,relPos,0,,"position of the contact on the cylinder (0: node-, 1:node+) |yupdate|")),
                                           /* extra initializers */,
                                           /* ctor */ createIndex();,
                                           /* py */
                                          );
    REGISTER_CLASS_INDEX(CylScGeom6D,ScGeom6D);
};
REGISTER_SERIALIZABLE(CylScGeom6D);

