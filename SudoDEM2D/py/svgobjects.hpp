// *********************************************************************
#ifndef SVGOBJECTS_HPP
#define SVGOBJECTS_HPP

#include <cassert>
#include <cmath>
#include <vector>
#include <string> // std::string, std::stoi
#include <iostream>
#include <fstream> // std::fstream
#include <sstream> // std::stringstream
#include <memory> //std::shared_ptr

#include <Eigen/Dense>

typedef Eigen::Matrix<double,2,1> Vector2r;
typedef Eigen::Matrix<double,3,1> Vector3r;

// *****************************************************************************
// SVGContext:
class SVGContext
{
  public:
  std::ostream * os ;
  SVGContext() ;
} ;

// *****************************************************************************
// clase para styles de polígonos y quizás otros tipos de objetos

class PathStyle
{
  public:
  PathStyle();
  void writeSVG( SVGContext & ctx ); // write style attrs to an svg file

  Vector3r  lines_color,      // lines color, when lines are drawn (draw_lines == true )
        fill_color ;      // fill color, when fill is drawn (draw_filled == true)
  bool  draw_lines,       // draw lines joining points (yes/no)
        close_lines,      // when draw_lines == true, join last point with first (yes/no)
        draw_filled ;     // fill the polygon (yes(no)
  float  lines_width ;
  float  fill_opacity ;    // fill opacity when draw_filled == true (0->transparent; 1->opaque)
  bool  dashed_lines  ;   // when draw_lines == true, draw dashed lines (yes/no)

  bool use_grad_fill ; // use gradient fill (when draw_filled == true)
  std::string grad_fill_name ; // name to use for gradient fill


} ;

// *****************************************************************************
// An abstract class for things which can be drawn to an SVG file and can be
// projected to 2d (must be projected before drawn)

class Object
{
  public:
  Object();
  virtual void drawSVG( SVGContext & ctx ) = 0 ;
  virtual void minmax( ) = 0 ;
  virtual ~Object() ;
  Vector2r min,max ;
} ;

// *****************************************************************************
// class Point
// A class for an isolated point, drawn with a given radius

class Point : public Object
{
  public:
  Point( Vector2r ppos2D, Vector3r color );
  virtual void drawSVG( SVGContext & ctx ) ;
  virtual void minmax() ;

  Vector3r color ;
  Vector2r pos2D ;
  float radius ; // 0.03 por defecto, se puede cambiar
};

// *****************************************************************************
// class Polygon
// Any Object which is described by a sequence of points
// it can be a polyline, a polygon, a filled polygon.

class Polygon : public Object // secuencia de points
{
   public:
   Polygon();
   virtual void minmax() ;
   virtual void drawSVG( SVGContext & ctx )  ;
   virtual ~Polygon() ;

   PathStyle         style ;
   std::vector<Vector2r> points2D ; // projected points
};
// *****************************************************************************
// class ObjectsSet
// A set of various objects

class ObjectsSet : public Object
{
   public:
   ObjectsSet();
   virtual void minmax() ;
   virtual void drawSVG( SVGContext & ctx )  ;
   virtual ~ObjectsSet() ;
   void add( const std::shared_ptr<Object>& pobj );  // add one object

   std::vector<std::shared_ptr<Object>> objetos ;
} ;


// *****************************************************************************
// class Segment
// A polygon with just two points.

class Segment : public Polygon
{
   public:
   Segment( const Vector2r & p0, const Vector2r & vd, const Vector3r & color, float width ) ;
   Segment( const Vector2r & p0, const Vector2r & p1  );
   Segment( const Point & p0, const Point & p1, float width );

};
// *****************************************************************************
// class Ellipse
// A planar ellipse in 3D

class Ellipse : public Polygon
{
   public:
   Ellipse( unsigned n, const Vector2r & center, const Vector2r & eje1, const Vector2r eje2  );
};
// *****************************************************************************
// class Figure
// A container for a set of objects which can be drawn to a SVG file.
// It can be written to a SVG file, the output includes the whole SVG elements
// (headers, footers, bounding box, etc....)
/*
class Figure
{
   public:

   Figure( ) ;
   void drawSVG( const std::string & nombre_arch ) ;

   ObjectsSet objetos ;  // set of objects in the figure
   float       width_cm ; // width (in centimeters) in the SVG header
   bool       flip_axes ; // true to flip axes (see Axes::Axes), false by default

   std::vector< std::string > rad_fill_grad_names ; // name of radial fill gradients to output
} ;
*/
#endif
