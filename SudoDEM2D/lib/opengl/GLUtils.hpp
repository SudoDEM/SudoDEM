// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
//
// header-only utility functions for GL (moved over from extra/Shop.cpp)
//
#pragma once

#include<sudodem/lib/opengl/OpenGLWrapper.hpp>
#include<sstream>
#include<iomanip>
#include<string>

struct GLUtils{
	// code copied from qglviewer
	struct QGLViewer{
		static void drawArrow(float length=1.0f, float radius=-1.0f, int nbSubdivisions=12);
		static void drawArrow(const Vector3r& from, const Vector3r& to, float radius=-1.0f, int nbSubdivisions=12);
	};
	// render wire of parallelepiped with sides given by vectors a,b,c; zero corner is at origin
	static void Parallelepiped(const Vector3r& a, const Vector3r& b, const Vector3r& c);
	static void Parallelogram(const Vector2r& a, const Vector2r& b);
	static void Square(const Real& length);
	static void GLDrawArrow(const Vector3r& from, const Vector3r& to, const Vector3r& color=Vector3r(1,1,1)){
		glEnable(GL_LIGHTING); glColor3v(color); QGLViewer::drawArrow(from,to);
	}
	static void GLDrawLine(const Vector3r& from, const Vector3r& to, const Vector3r& color=Vector3r(1,1,1)){
		glEnable(GL_LIGHTING); glColor3v(color);
		glBegin(GL_LINES); glVertex3v(from); glVertex3v(to); glEnd();
	}

	static void GLDrawNum(const Real& n, const Vector2r& pos, const Vector3r& color=Vector3r(1,1,1), unsigned precision=3){
		std::ostringstream oss; oss<<std::setprecision(precision)<< /* "w="<< */ (double)n;
		GLUtils::GLDrawText(oss.str(),pos,color);
	}

	static void GLDrawInt(long i, const Vector2r& pos, const Vector3r& color=Vector3r(1,1,1)){
		GLUtils::GLDrawText(boost::lexical_cast<std::string>(i),pos,color);
	}

	static void GLDrawText(const std::string& txt, const Vector2r& pos, const Vector3r& color=Vector3r(1,1,1)){
		glPushMatrix();
		glTranslatev(Vector3r(pos[0],pos[1],0));
		glColor3(color[0],color[1],color[2]);
		glRasterPos2i(0,0);
		for(unsigned int i=0;i<txt.length();i++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, txt[i]);
		glPopMatrix();
	}
};
