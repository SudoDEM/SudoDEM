// *********************************************************************
#include "svgobjects.hpp"

// Aux functions
// -----------------------------------------------------------------------------

void write_color( std::ostream & os, const Vector3r & col )
{

  //os << "rgb("<< col[0]*255.0 << "%,"<< col[1]*255.0 << "%,"<< col[2]*255.0 << "%)"  ;
  os << "rgb("<< col[0]*255.0 << ","<< col[1]*255.0 << ","<< col[2]*255.0 << ")"  ;
}

// -----------------------------------------------------------------------------

void write_stroke_color( std::ostream & os, const Vector3r & col )
{

  //os << "stroke='rgb("<< col[0]*255.0 << "%,"<< col[1]*255.0 << "%,"<< col[2]*255.0 << "%)' " << std::endl ;
  os << "stroke='rgb("<< col[0]*255.0 << ","<< col[1]*255.0 << ","<< col[2]*255.0 << ")' " << std::endl ;
}
// -----------------------------------------------------------------------------

void write_coord2( std::ostream & os, const Vector2r & p )
{
   //using namespace std ;
   if ( std::isnan(p[0]) && std::isnan(p[1]) )
      std::cout << "warning: not a number found in Vector2r. Not written to SVG output file." << std::endl ;
   else
       os << p[0] << " " << p[1] ;
}
// *****************************************************************************
// class SVGContext
// -----------------------------------------------------------------------------

SVGContext::SVGContext()
{
  os = nullptr ;
}

// *****************************************************************************
// class PathStyle
// -----------------------------------------------------------------------------

PathStyle::PathStyle()
{
  draw_lines   = true ;
  draw_filled  = true ;
  close_lines  = true ;
  lines_color  = Vector3r(0.0,0.0,0.0);
  fill_color   = Vector3r(1.0,0.5,0.5);
  fill_opacity = 0.5 ;
  lines_width  = 0.05 ;
  dashed_lines = false ;

  use_grad_fill = false ;
  grad_fill_name = "no name" ;
}
// -----------------------------------------------------------------------------

void PathStyle::writeSVG( SVGContext & ctx )
{
  //using namespace std ;
  std::ostream & os = *(ctx.os) ;

   os << "   style='" ;

   os << "fill:" ;
   if ( draw_filled )
   {
      if ( use_grad_fill )
      {
         if ( grad_fill_name == "no name" )
            std::cout << "WARNING: using 'no name' as grad fill name." << std::endl ;
         os << "url(#" << grad_fill_name << ")" ;
      }
      else
         write_color( os, fill_color ) ;
   }
   else
      os << "none" ;
   os << "; " ;

   if ( draw_filled )
      os << "fill-opacity:" << fill_opacity << "; " ;

   os << "stroke:" ;
   if ( draw_lines )
      write_color( os, lines_color ) ;
   else
      os << "none" ;
   os << "; " ;

   if ( draw_lines )
      os << "stroke-width:" << lines_width << "; " ;

   if ( draw_lines && dashed_lines )
      os << "stroke-dasharray:0.01,0.01; " ;

   os << "' " << std::endl ;

}

// *****************************************************************************
// class Object
// -----------------------------------------------------------------------------

Object::Object()
{
  min = max = Vector2r::Zero() ;
}
// -----------------------------------------------------------------------------

Object::~Object()
{

}

// *****************************************************************************
// class Puntos
// -----------------------------------------------------------------------------

Point::Point( Vector2r ppos2D, Vector3r pcolor )
{
  pos2D = ppos2D ;
  color = pcolor ;
  //min       = pos2D ;
  //max       = pos2D ;
  radius = 0.03 ;
}

void Point::minmax( )
{
  min       = pos2D ;
  max       = pos2D ;
}
// -----------------------------------------------------------------------------

void Point::drawSVG( SVGContext & ctx )
{
  assert( projected );
  //using namespace std ;


  PathStyle e ;
  e.draw_lines   = false ;
  e.draw_filled  = true ;
  e.fill_color   = color ;
  e.fill_opacity = 1.0;

  assert( ctx.os != nullptr );
  std::ostream & os = *(ctx.os) ;

  os << "<circle cx='" << pos2D[0] << "' cy='" << pos2D[1]
     <<       "' r='" << radius << "' " ;
  e.writeSVG( ctx );
  os << "/>" << std::endl ;

}
// *****************************************************************************
// class Polygon
// -----------------------------------------------------------------------------

Polygon::Polygon()
{

}
// -----------------------------------------------------------------------------

void Polygon::minmax( )
{
  for( const Vector2r & p2 : points2D )
  {

      min[0] = std::min( min[0], p2[0] );
      min[1] = std::min( min[1], p2[1] );
      max[0] = std::max( max[0], p2[0] );
      max[1] = std::max( max[1], p2[1] );
  }
}

Polygon::~Polygon()
{

}

// -----------------------------------------------------------------------------

void Polygon::drawSVG( SVGContext & ctx )
{
   //using namespace std ;
   assert( ctx.os != nullptr );


   if ( points2D.size() == 0 )
   {
      std::cout << "WARNING: attempting to draw an empty Polygon object" << std::endl ;
      return ;
   }

  //using namespace std ;
  std::ostream & os = *(ctx.os) ;

  os << "<path " << std::endl
     << "   stroke-linecap='round' stroke-linejoin='round'" << std::endl ;

  style.writeSVG( ctx );

  os << "   d=' M " ;
  write_coord2( os, points2D[0] );
  for( unsigned i = 1 ; i < points2D.size() ; i++ )
  {
      os << " L " ;
      write_coord2( os, points2D[i] );
      os << " " ;
      //cout << points_proy[i] << endl ;
  }
  if ( style.close_lines )
    os << " Z" ;

  os << "'/>" << std::endl;

}

// *****************************************************************************
// class ConjuntoPuntos
// -----------------------------------------------------------------------------

ObjectsSet::ObjectsSet()
{

}

ObjectsSet::~ObjectsSet()
{

}
// -----------------------------------------------------------------------------

void ObjectsSet::add( const std::shared_ptr<Object>& pobj )
{
  objetos.push_back( pobj );
}

// -----------------------------------------------------------------------------

void ObjectsSet::minmax()
{
   for( std::shared_ptr<Object> pobjeto : objetos )
   {
      assert( pobjeto != nullptr );
      pobjeto->minmax( );
      min[0] = std::min( min[0], pobjeto->min[0] );
      min[1] = std::min( min[1], pobjeto->min[1] );
      max[0] = std::max( max[0], pobjeto->max[0] );
      max[1] = std::max( max[1], pobjeto->max[1] );
  }

}

// -----------------------------------------------------------------------------

void ObjectsSet::drawSVG( SVGContext & ctx )
{
  for( std::shared_ptr<Object> pobjeto : objetos )
  {
    assert( pobjeto != nullptr );
    pobjeto->drawSVG( ctx );
  }
}

// *****************************************************************************
// class Segment
// -----------------------------------------------------------------------------

Segment::Segment( const Vector2r & p0, const Vector2r & vd, const Vector3r & color, float width )
{
  style.draw_lines   = true ;
  style.draw_filled  = false ;
  style.close_lines  = false ;
  style.lines_color  = color ;
  style.lines_width  = width ;

  points2D.push_back( p0 );
  points2D.push_back( p0+vd );
}

Segment::Segment( const Vector2r & p0, const Vector2r & p1 )
{
  style.draw_lines   = true ;
  style.draw_filled  = false ;
  style.close_lines  = false ;
  style.lines_color  = Vector3r(0.0,0.0,0.0) ;
  style.lines_width  = 0.007 ;

  points2D.push_back( p0 );
  points2D.push_back( p1 );
}

Segment::Segment( const Point & p0, const Point & p1, float width )
{
  style.draw_lines   = true ;
  style.draw_filled  = false ;
  style.close_lines  = false ;
  style.lines_color  = p0.color ;
  style.lines_width  = width ;

  points2D.push_back( p0.pos2D );
  points2D.push_back( p1.pos2D );
}

// *****************************************************************************
// class Ellipse
// an arbitrary ellipse, at any point, with any orientation

Ellipse::Ellipse( unsigned n, const Vector2r & center, const Vector2r & eje1, const Vector2r eje2 )
{
  for( unsigned i= 0 ; i < n ; i++ )
  {
    const float ang = float(i)*float(2.0)*float(M_PI)/float(n) ,
               c   = cos(double(ang)),
               s   = sin(double(ang));

    points2D.push_back( center + c*eje1 + s*eje2 ) ;
  }
  style.draw_filled = false ;
  style.draw_lines  = true ;
  style.close_lines = true ;
  style.lines_width = 0.006 ;
  style.lines_color = Vector3r( 0.0, 1.0, 0.0 );
}
