//@2017 by Sway Zhao,  zhswee@gmail.com
//


#include<sudodem/pkg/fem/TriElement.hpp>
#include<sudodem/pkg/fem/Node.hpp>
#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Omega.hpp>
#include<sudodem/pkg/common/Sphere.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/lib/pyutil/doc_opts.hpp>
#include<cmath>

#include<numpy/ndarrayobject.h>
#include <boost/concept_check.hpp>

namespace py = boost::python;
//***********************************************************************************
//new membrane

shared_ptr<Body> NewMembrane(shared_ptr<Node> n1,shared_ptr<Node> n2,shared_ptr<Node> n3, shared_ptr<Material> mat){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<TriElement>(new TriElement());
	TriElement* A = static_cast<TriElement*>(body->shape.get());
    //A->nids = nids;
    A->node1 = n1;
    A->node2 = n2;
    A->node3 = n3;
    A->nodes.push_back(n1);
    A->nodes.push_back(n2);
    A->nodes.push_back(n3);
    A->pos = A->getCentroid();
    A->lumpMassInertia(mat->density);
    for(int i:{0,1,2}){ A->vertices[i] = A->nodes[i]->pos- A->pos;}
    body->state->pos = A->pos;
	body->state->blockedDOFs= Node::DOF_ALL;
	//body->state->mass=body->material->density*A->getVolume();
	//body->state->inertia =A->getInertia()*body->material->density;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(false);//we assume the mass and inertia of the TriElement is zero in DEM
	//

	return body;
}
//a cylindrical wall without ends
bool CylindricalWall(Vector3r bottomCenter, double radius, double height, int elm_res,shared_ptr<Material> mat){
    //bottomCenter: position of the bottom center of the cylindrical wall
    //radius: radius of the cylinder ends
    //height: height of the cyliner
    //elm_res: number of elements along the circumference
    double theta = 2.0*Mathr::PI/elm_res;
    //width of an element
    double elm_size_w = radius*theta;
    //guess the height of an element
    int elm_res_h = int(height/elm_size_w);
    double elm_size_h = height/elm_res_h;
    Vector3r p;
    double x,y,z;
    Scene* scene=Omega::instance().getScene().get();
    //create nodes
    std::cout<<"nodes creating ..."<<elm_res<<"x"<<elm_res_h+1<<std::endl;
    for(int j=0;j<elm_res_h+1;j++){//along height
        for(int i=0;i<elm_res;i++){//along circumference
            x = radius*cos(i*theta);
            y = radius*sin(i*theta);
            z = elm_size_h*j;
            p = Vector3r(x,y,z)+bottomCenter;
            //create a node
            shared_ptr<Node> node(new Node);
            node->pos = p;
            node->mass = 1000.0;
            node->inertia = Vector3r(1.,1.,1.);
            scene->nodes->insert(node);
            //std::cout<<"node id = "<<node->id<<std::endl;
        }
    }
    //create TriElements
    std::cout<<"TriElements creating ..."<<"2*"<<elm_res<<"x"<<elm_res_h+1<<std::endl;
    for(int j=0;j<elm_res_h;j++){//along height
        for(int i=0;i<elm_res-1;i++){//along circumference
            const shared_ptr<Node>& n1 = Node::byId(elm_res*j+i);
            const shared_ptr<Node>& n2 = Node::byId(elm_res*j+i+1);
            const shared_ptr<Node>& n3 = Node::byId(elm_res*(j+1)+i+1);
            const shared_ptr<Node>& n4 = Node::byId(elm_res*(j+1)+i);
            
            shared_ptr<Body> body1 = NewMembrane(n1,n2,n3,mat);
            shared_ptr<Body> body2 = NewMembrane(n1,n3,n4,mat);
            //std::cout<<i<<" "<<j<<std::endl;
            scene->bodies->insert(body1);
            scene->bodies->insert(body2);   
        }
        const shared_ptr<Node>& n1 = Node::byId(elm_res*j+(elm_res-1));
        const shared_ptr<Node>& n2 = Node::byId(elm_res*j);
        const shared_ptr<Node>& n3 = Node::byId(elm_res*(j+1));
        const shared_ptr<Node>& n4 = Node::byId(elm_res*(j+1)+(elm_res-1));
        
        shared_ptr<Body> body1 = NewMembrane(n1,n2,n3,mat);
        shared_ptr<Body> body2 = NewMembrane(n1,n3,n4,mat);
        scene->bodies->insert(body1);
        scene->bodies->insert(body2);  
        
    }
    return true;
}

//**********************************************************************************//


BOOST_PYTHON_MODULE(_fem_utils){
	// http://numpy.scipy.org/numpydoc/numpy-13.html mentions this must be done in module init, otherwise we will crash
	//import_array();

	SUDODEM_SET_DOCSTRING_OPTS;
    py::def("NewMembrane", NewMembrane, "creat a membrane element with three nodes.");
    py::def("CylindricalWall", CylindricalWall, "creat a membrane element with three nodes.");
}

