//@2016 Sway Zhao,  zhswee@gmail.com
//


#include"sudodem/pkg/dem/Superellipse.hpp"
#include"sudodem/pkg/common/Disk.hpp"

#include<sudodem/core/Scene.hpp>
#include<sudodem/core/Omega.hpp>
//#include<sudodem/pkg/common/Disk.hpp>
#include<sudodem/pkg/common/ElastMat.hpp>
#include<sudodem/lib/pyutil/doc_opts.hpp>
#include<cmath>

#include<numpy/ndarrayobject.h>
#include <boost/concept_check.hpp>

#include "svgobjects.hpp"//output svg images

namespace py = boost::python;
//***********************************************************************************

int Sign(double f){ if(f<0) return -1; if(f>0) return 1; return 0; }

//***********************************************************************************
shared_ptr<Body> NewSuperellipse(double x, double y, double ep, shared_ptr<Material> mat,bool rotate, bool isSphere, double z_dim =1.0){
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superellipse>(new Superellipse(x,y,ep));
	Superellipse* A = static_cast<Superellipse*>(body->shape.get());

  A->isSphere = isSphere;
	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getArea()*z_dim;//FIXME:with a thickness of unity
	body->state->inertia =A->getInertia()*body->material->density*z_dim;//FIXME:with a thickness of unity
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//
	if(rotate){
		double heading = double(rand())/RAND_MAX*2.0*Mathr::PI;//rotate around z
		Rotationr rot = Rotationr(heading);
	  A->setOrientation(A->getOrientation()*rot);
	 }
	body->state->ori = A->getOrientation();
	return body;
}

shared_ptr<Body> NewSuperellipse_rot(double x, double y, double ep, shared_ptr<Material> mat,double miniAngle, bool isSphere, double z_dim = 1.0){//with a specified quatanion
	shared_ptr<Body> body(new Body);
	body->material=mat;
	body->shape=shared_ptr<Superellipse>(new Superellipse(x,y,ep));
	Superellipse* A = static_cast<Superellipse*>(body->shape.get());
  A->isSphere = isSphere;

	body->state->pos= A->getPosition();
	body->state->mass=body->material->density*A->getArea()*z_dim;
	body->state->inertia =A->getInertia()*body->material->density*z_dim;
	body->shape->color = Vector3r(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->setAspherical(true);
	//
  Rotationr rot = Rotationr(miniAngle);
  A->setOrientation(rot);

	body->state->ori = A->getOrientation();
	return body;
}
//**********************************************************************************//
//output geometrical info of an assembly
void outputParticles(const char* filename){

	FILE * fin = fopen(filename,"w");//"a" - append;"r" -read; "w" -write


	//fprintf(fin,"volume\t%e\n",volume);
	fprintf(fin,"#rx\try\teps\tx\ty\trotation\n");

  //output superquadrics
  Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	//particlesVolume = 0;
	Rotationr rot;
	Vector2r pos,rxy;
	for ( ; bi!=biEnd; ++bi )
	{
		const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";

		if (b->shape->getClassName()=="Superellipse"){
			const shared_ptr<Superellipse>& A = SUDODEM_PTR_CAST<Superellipse> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());
			State* state=b->state.get();
      //Euler Angles from Quaternion
      rot = state->ori;
      pos = state->pos;
    	rxy = A->getrxy();
      //output
      fprintf(fin,"%e\t%e\t%e\t%e\t%e\t%e\n",rxy(0),rxy(1),A->geteps(),pos(0),pos(1),rot.angle());
		}

	}
	fclose(fin);
}

// -----------------------------------------------------------------------------
void write_radial_gradient_def( std::ostream & os, const std::string & name )
{
   //using namespace std ;
   os << "<radialGradient id='" << name << "' cx='50%' cy='50%' r='50%' fx='25%' fy='75%'>" << std::endl
      << "   <stop offset='0%'   style='stop-color:rgb(100%,100%,100%); stop-opacity:0.2' />" << std::endl
      << "   <stop offset='100%' style='stop-color:rgb(20%,20%,20%);    stop-opacity:0.2' />" << std::endl
      << "</radialGradient>" << std::endl ;
}

// -----------------------------------------------------------------------------

void write_radial_gradient_def_blue( std::ostream & os, const std::string & name )
{
   //using namespace std ;
   os << "<radialGradient id='" << name << "' cx='50%' cy='50%' r='50%' fx='50%' fy='50%'>" << std::endl
      << "   <stop offset='0%'   style='stop-color:rgb(100%,100%,100%); stop-opacity:0.2' />" << std::endl
      << "   <stop offset='100%' style='stop-color:rgb(0%,0%,50%);    stop-opacity:0.4' />" << std::endl
      << "</radialGradient>" << std::endl ;
}


ObjectsSet particleObjects(int slice, float line_width, float fill_opacity, bool draw_filled, bool draw_lines,bool color_po, Vector3r po_color, Vector3r solo_color){
	ObjectsSet objetos ;  // set of objects in the figure
	//add objects
	//output superquadrics
	Scene* scene=Omega::instance().getScene().get();
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	const bool isPeriodic = scene->isPeriodic;
	//particlesVolume = 0;
	Rotationr rot;
	Vector2r pos;
	double delta_theta = 2.0*Mathr::PI/slice;
	double line_width0 = -1.0;
	for ( ; bi!=biEnd; ++bi )
	{
		const shared_ptr<Body>& b = *bi;
		//std::cerr << "getClassName= " << b->shape->getClassName() << "\n";
		Vector3r color = (solo_color.norm()==0?b->shape->color:solo_color);
		double po_value = 0.0;//color for particle orientation
		std::shared_ptr<Polygon> poly = std::shared_ptr<Polygon>(new Polygon());

		State* state=b->state.get();
		//Euler Angles from Quaternion
		rot = state->ori;
		po_value = rot.smallestPositiveAngle();
		po_value = (po_value > Mathr::PI ? po_value-Mathr::PI : po_value);//normalized to [0,1) corresponding to angle in [0,pi).
		po_value = (po_value > 0.5*Mathr::PI ? abs(po_value-Mathr::PI) : po_value)/Mathr::PI*2.0;
		pos = state->pos;
		pos=(!isPeriodic ? pos : scene->cell->wrapShearedPt(pos));
		/*
		if(scene->isPeriodic){
			const Vector2r& cellSize(scene->cell->getSize());
			pos=scene->cell->unshearPt(pos); // remove the shear component
			pos += Vector2r(cellSize[0]*(pos[0]<0?1:0),cellSize[1]*(pos[1]<0?1:0));
			pos = scene->cell->shearPt(pos);
		}*/

		if (b->shape->getClassName()=="Superellipse"){
			const shared_ptr<Superellipse>& A = SUDODEM_PTR_CAST<Superellipse> (b->shape);
			//particlesVolume += A->getVolume();
			//Superquadrics* A = static_cast<Superquadrics*>(b->shape.get());

			if(line_width0 < 0){ line_width0 = min(A->getrxy()(0),A->getrxy()(1));}
			//output surface points
			for( unsigned i= 0 ; i < slice ; i++ )
			{
				float ang = i*delta_theta;
				//const shared_ptr<Superellipse>& A = SUDODEM_PTR_CAST<Superellipse> (b->shape);
				Vector2r p = pos + rot.toRotationMatrix()* A->getSurface(ang);
				//std::cout<<p<<std::endl;
				poly->points2D.push_back( p ) ;
			}
		}else if (b->shape->getClassName()=="Disk"){
			const shared_ptr<Disk>& A = SUDODEM_PTR_CAST<Disk> (b->shape);
			if(!A) break;
			double radius = A->radius;
			if(line_width0 < 0){line_width0 = radius;}
			//output surface points
			for( unsigned i= 0 ; i < slice ; i++ )
			{
				float ang = i*delta_theta;
		    double c = cos(double(ang));
		    double s = sin(double(ang));
			  Vector2r p = pos + radius*Vector2r(cos(ang),sin(ang));
				//std::cout<<p<<std::endl;
				poly->points2D.push_back( p ) ;
			}
		}else{
			std::cout<<"Please make sure that Superellipse or Disk are in the scene."<<std::endl;
		}
		poly->style.draw_filled = draw_filled ;
		poly->style.draw_lines  = draw_lines ;
		poly->style.close_lines = true ;
		poly->style.lines_width = line_width*line_width0 ;
		poly->style.lines_color = (color_po==true?po_value*po_color:color);
		poly->style.fill_color  = (color_po==true?po_value*po_color:color);
		poly->style.fill_opacity  = fill_opacity;
		objetos.add( poly );
	}

	return objetos;
}

ObjectsSet forceChainObjects(float line_width, float fill_opacity, Vector3r line_color){
	ObjectsSet objetos ;  // set of objects in the figure
	//add objects
	//output superquadrics
	Scene* scene=Omega::instance().getScene().get();
	const bool isPeriodic = scene->isPeriodic;
	double force=0.0;
	int force_num = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		//not considering clumped particles. CAUTION:SudoDEM does not prefer clumped particles.
		NormShearPhys* nsi=SUDODEM_CAST<NormShearPhys*> ( I->phys.get() );

		//Contact force
		double fn = (nsi->normalForce).norm();
		if (fn == 0)continue;

		std::shared_ptr<Polygon> poly = std::shared_ptr<Polygon>(new Polygon());
		//Vector2r ft = nsi->shearForce;
		force += fn ;//here we constructe a positive stress tensor
		force_num ++;
		Vector2r pos1 = Body::byId(I->getId1(),scene)->state->pos;
		Vector2r pos2 = Body::byId(I->getId2(),scene)->state->pos;
		Vector2r branch = pos2 - pos1;
		//shifting particle id2 only
		if(isPeriodic) branch += scene->cell->hSize*I->cellDist.cast<Real>();
		pos1=(!isPeriodic ? pos1 : scene->cell->wrapShearedPt(pos1));
		//if (branch.norm()>1e-2) std::cout<<"branch vector ..."<<std::endl;
		/*const Vector2r& cellSize(scene->cell->getSize());
		pos1=scene->cell->unshearPt(pos1); // remove the shear component
		pos1 += Vector2r(cellSize[0]*(pos1[0]<0?1:0),cellSize[1]*(pos1[1]<0?1:0));
		pos1 = scene->cell->shearPt(pos1);
		*/
		pos2 = branch + pos1;

		poly->points2D.push_back( pos1 ) ;
		poly->points2D.push_back( pos2 ) ;

		poly->style.draw_filled = false ;
		poly->style.draw_lines  = true ;
		poly->style.close_lines = false ;
		poly->style.lines_width = fn;//fnline_width*line_width0 ;
		poly->style.lines_color = line_color;;//line_value*line_color;
		//poly->style.fill_color  = line_value*line_color;
		poly->style.fill_opacity  = fill_opacity;
		objetos.add( poly );
	}
	//mean force
	force /=force_num;
	for(auto ob:objetos.objetos){
		Polygon *poly = dynamic_cast<Polygon*>(ob.get());
		double w = poly->style.lines_width * (line_width/force);

		if(w<0.1*line_width){
			w = 0.1*line_width;
		}else if(w>10.0*line_width){
			w = 10.0*line_width;
		}
		poly->style.lines_width = w;
		//std::cout<<poly->style.lines_width<<std::endl;
	}
	return objetos;
}

void drawSVG( const std::string & filename, float width_cm = 10.0, int slice = 20, float line_width = 0.1,
							float fill_opacity = 0.5, bool draw_filled = true, bool draw_lines = true,bool color_po = false,
							Vector3r po_color = Vector3r(1.0,0.0,0.0),Vector3r solo_color=Vector3r(0,0,0),float force_line_width=0.001, float force_fill_opacity=0,
							Vector3r force_line_color=Vector3r(1.0,1.0,1.0), bool draw_forcechain = false)
{
	 //width_cm: width (in centimeters) in the SVG header
	 //slice: to how many slices a Superellipse will be discretized
	 //line_width: the line width of the edge of a Superellipse
	 //fill_opacity: the opacity of the fill color in a Superellipse
	 //draw_lines: whether to draw the out profile of a Superellipse
	 //draw_filled: wheter to fill a Superellipse with a color
	 //color_po: whether to use a fill color to identify the orientation of a Superellipse
	 //po_color: the fill color to identify the orientation of a Superellipse
	 //solo_color: a solor color to fill a Superellipse. If its norm is zero, then use the color defined by po_color or shape->color
	 //force_line_width: the line width of the force chain with the average normal contact force
	 //force_line_color: the color of the force chain
	 //draw_forcechain: whether to draw the contact force chain that will be superposed with the particles

   std::vector< std::string > rad_fill_grad_names ; // name of radial fill gradients to output
   std::fstream fout( filename, std::ios_base::out) ;

   SVGContext ctx ;
   ctx.os = &fout ;
	 ObjectsSet objetos = particleObjects(slice,line_width,fill_opacity,draw_filled,draw_lines,color_po,po_color,solo_color);  // set of objects in the figure
	 ObjectsSet objetos2;//used for force chains
	 if(draw_forcechain){
		 objetos2 = forceChainObjects(force_line_width, force_fill_opacity, force_line_color);
		 /*for(auto ob:objetos2.objetos){
			 Polygon *poly = dynamic_cast<Polygon*>(ob);
			 objetos.add(poly);
		 }*/
	 }

	 //draw
   objetos.minmax();
   constexpr float fmrg = 0.01 ;
   //constexpr real fmrg = 0.01 ; // width del margen en X y en Y, expresado en porcentaje del width en X
   //std::cout << "objetos.max == " << objetos.max << ", objetos.min == " << objetos.min << std::endl ;

   const float ax       = objetos.max[0]-objetos.min[0] ;
   const Vector2r vmargen  = Vector2r( fmrg*ax, fmrg*ax );
   const Vector2r ptos_min = objetos.min,
              ptos_w   = objetos.max-objetos.min,
              box_min  = ptos_min-vmargen,
              box_w    = ptos_w+float(2.0)*vmargen ;

   //assert( 0 < box_w[0] && 0 < box_w[1] );
   float ratio = box_w[1]/box_w[0] ;

   const float wx = width_cm,      // width en centimetros
              wy = wx*ratio ;
	 //head
	 fout << "<!-- This file is written by SudoDEM2D. -->"<< std::endl;


   // cabecera svg
   fout << "<svg xmlns='http://www.w3.org/2000/svg' " << std::endl
        << "     width='" << wx << "cm' height='" << wy << "cm' " << std::endl
        << "     viewBox='" << box_min[0] << " " << box_min[1] << " " << box_w[0] << " " << box_w[1] << "'" << std::endl
        << ">" << std::endl ;


   //write_radial_gradient_def( fout, "hemisphereGradFill" );

   fout << "<defs>" << std::endl ;
   write_radial_gradient_def( fout, "hemisphereGradFill" );
   write_radial_gradient_def_blue( fout, "spherecapGradFill" );
   fout << "</defs>" << std::endl ;

   fout << "<g transform='translate(0.0 " << 2.0*box_min[1]+box_w[1] << ") scale(1.0 -1.0)'> <!-- transf global (inv y) -->"<< std::endl ;





   // objetos svg
	 if(draw_forcechain){
		 fout << "<g id='layer1' inkscape:label='Particles' inkscape:groupmode='layer'>" << std::endl ;
		 fout << "<g>" << std::endl ;
		 objetos.drawSVG( ctx );
		 fout << "</g>" << std::endl ;
		 fout << "</g>" << std::endl ;
		 fout << "<g id='layer2' inkscape:label='Force chains' inkscape:groupmode='layer'>" << std::endl ;
		 fout << "<g>" << std::endl ;
		 objetos2.drawSVG( ctx );
		 fout << "</g>" << std::endl ;
		 fout << "</g>" << std::endl ;
	 }else{
		 objetos.drawSVG( ctx );
	 }


   // pie svg
   fout << "</g>" << std::endl ;
   fout << "</svg>" << std::endl ;

   fout.close();
}



BOOST_PYTHON_MODULE(_superellipse_utils){
	// http://numpy.scipy.org/numpydoc/numpy-13.html mentions this must be done in module init, otherwise we will crash
	//import_array();

	SUDODEM_SET_DOCSTRING_OPTS;
	py::def("NewSuperellipse", NewSuperellipse,(py::args("z_dim")=1), "Generate a Superellipse.");
  py::def("NewSuperellipse_rot", NewSuperellipse_rot,(py::args("z_dim")=1), "Generate a Super-ellipsoid with parameters x, y, ep, .");
	py::def("outputParticles", outputParticles, "output particle info(rx,ry,eps,x,y,rotation angle) to a text file.");
	py::def("drawSVG", drawSVG, (py::args("width_cm")=10.0,py::args("slice")=20,py::args("line_width")=0.1,py::args("fill_opacity")=0.5,
															 py::args("draw_filled")=true,py::args("draw_lines")=true,py::args("color_po")=false,py::args("po_color")=Vector3r(1.0,0.0,0.0),
															 py::args("solo_color")=Vector3r(0.0,0.0,0.0),
															 py::args("force_line_width") = 0.001, py::args("force_fill_opacity")=0,py::args("force_line_color")=Vector3r(1.0,1.0,1.0),
															 py::args("draw_forcechain") = false),"draw Superellipses into a svg file.");
	//py::def("outputPOV",outputPOV,"output Superellipse into a POV file for post-processing using POV-ray.");
  //py::def("outputVTK",outputVTK,"output Superellipse into a VTK file for post-processing using Paraview.");
}
