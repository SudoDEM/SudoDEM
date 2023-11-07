# enconding: utf-8
##########################################
#*************************************************************************
#*  Copyright (C) 2018 by Sway Zhao                                       *
#*  zhswee@gmail.com                                                      *
#*                                                                        *
#*  This program is free software; it is licensed under the terms of the  *
#*  GNU General Public License v2 or later. See file LICENSE for details. *
#*************************************************************************/
##########################################
"""
Module containing utinity functions for data visulization.
"""
## all exported names
__all__=['outputBoxPov','outPov']
import sudodem.qt
#from sudodem import *
from sudodem.wrapper import *
from sudodem import  _superquadrics_utils

#from sudodem._superquadrics_utils import *

import math,os
import numpy as np

#output *.POV of a cubic box for post-processing in Pov-ray
def outputBoxPov(filename,x,y,z,r=0.001,wallMask=[1,0,1,0,0,1]):#x=[x_min,x_max],y=[y_min,y_max],z=[z_min,z_max]
    fobj=open(filename,'w+')
    points = list()
    for k in z:
            for j in y:
                    for i in x:
                            points.append("<"+str(i)+","+str(j)+","+str(k)+">")
    print>>fobj, "// The radius of the cylinders to be used for each edge of the box"
    print>>fobj, "#declare r="+str(r)+";"
    print>>fobj, "#declare f1=finish{reflection 0.15 specular 0.3 ambient 0.42}"
    print>>fobj, "#declare t1=texture{pigment{rgbft <0.3,0.9,0.7,0,0.4>} finish{f1}}"
    print>>fobj, "#declare t2=texture{pigment{rgb <0.1,0.5,0.3>} finish{f1}}"
    #points include eight vertices of the bouding box

    print>>fobj, "// box"
    #faces
    print>>fobj, "union{"
    if wallMask[0] >0: print>>fobj, "//left wall\npolygon{4,\t\n"+points[0]+",\t\n"+points[4]+",\t\n"+points[6]+",\t\n"+points[2]+"\n}"
    if wallMask[1] >0: print>>fobj, "//right wall\npolygon{4,\t\n"+points[1]+",\t\n"+points[3]+",\t\n"+points[7]+",\t\n"+points[5]+"\n}"
    if wallMask[2] >0: print>>fobj, "//front wall\npolygon{4,\t\n"+points[0]+",\t\n"+points[1]+",\t\n"+points[5]+",\t\n"+points[4]+"\n}"
    if wallMask[3] >0: print>>fobj, "//back wall\npolygon{4,\t\n"+points[3]+",\t\n"+points[2]+",\t\n"+points[6]+",\t\n"+points[7]+"\n}"
    if wallMask[4] >0: print>>fobj, "//top wall\npolygon{4,\t\n"+points[4]+",\t\n"+points[5]+",\t\n"+points[7]+",\t\n"+points[6]+"\n}"
    if wallMask[5] >0: print>>fobj, "//bottom wall\npolygon{4,\t\n"+points[0]+",\t\n"+points[2]+",\t\n"+points[3]+",\t\n"+points[1]+"\n}"
    print>>fobj, "texture{t1}\n}"
    print>>fobj, "union{"
    #vertices
    for v in points:
        print>>fobj,"sphere{"+v+",r}"
    #edges at the bottom
    pts = points[:4]
    tmp =pts[2]
    pts[2] = pts[3]
    pts[3] = tmp
    for i in range(4):
            print>>fobj, "cylinder{"+pts[i-1]+","+pts[i]+",r}"
    #edges at the top
    pts = points[4:]
    tmp =pts[2]
    pts[2] = pts[3]
    pts[3] = tmp
    for i in range(4):
            print>>fobj, "cylinder{"+pts[i-1]+","+pts[i]+",r}"
    #edges at the side
    for i in range(4):
            print>>fobj, "cylinder{"+points[i]+","+points[i+4]+",r}"

    #print>>fobj, "    pigment{rgb <0.7,0.95,1>} finish{specular 0.5 ambient 0.42}"
    print>>fobj, "texture{t2}\n}"
    fobj.close()

#Poly-EllipsoidPov
def outPov(fname,transmit=0,phong=1.0,singleColor=True,ids=[],scale=1.0, withbox = True):
    views=sudodem.qt._GLViewer.views()
    if len(views) == 0:
        raise (RuntimeError,"Open a view by clicking the show3D button.")
    view = views[0]
    fout = open("./"+fname+".pov",'w')
    loc = view.eyePosition*scale
    ps = view.lookAt*scale
    #view.axes=True
    #print view.upVector.dot(view.viewDir)
    #print view.viewDir
    #print box
    if(withbox):
        x = [O.bodies[0].state.pos[0],O.bodies[1].state.pos[0]]
        y = [O.bodies[2].state.pos[1],O.bodies[3].state.pos[1]]
        z = [O.bodies[4].state.pos[2],O.bodies[5].state.pos[2]]
        outputBoxPov("./"+fname+"box.inc",x,y,z,r=0.001,wallMask=[1,1,1,1,1,1])
    print>>fout,"//using the command: povray **.pov +A0.01 +W1600 +H1200"
    print>>fout,"#include \"colors.inc\""
    print>>fout,"#declare singleColor = ", "true;" if singleColor else "false;"
    print>>fout,"#declare Random_r = seed (1432);"
    print>>fout,"#declare Random_g = seed (7242);"
    print>>fout,"#declare Random_b = seed (9912);"
    print>>fout,"#declare para_trans = ", transmit,";"
    print>>fout,"#declare para_phong = ", phong,";"
    print>>fout,"camera {"
    print>>fout,"location <",loc[0],",",loc[1],",",loc[2],">"
    print>>fout,"sky z"
    print>>fout,"right -x*image_width/image_height"
    print>>fout,"look_at <",ps[0],",",ps[1],",",ps[2],">"
    print>>fout,"}"
    #print>>fout,"#include \"axes.inc\""
    print>>fout,"#include \""+fname+".inc\""
    if(withbox):
        print>>fout,"#include \""+fname+"box.inc\""
    print>>fout,"light_source{<4,8,5>*10 color rgb <1,1,1>}"
    print>>fout,"light_source{<12,-6>*10 color rgb <1,1,1>}"
    print>>fout,"light_source { <0, 2, 10> White }"
    print>>fout,"background{rgb 1}"
    print>>fout,"plane { z, -5 pigment { checker Green White }}"
    _superquadrics_utils.PolySuperellipsoidPOV("./"+fname+".inc",ids=ids,scale=scale)
    fout.close()
