# enconding: utf-8
##########################################
#*************************************************************************
#*  Copyright (C) 2016 by Sway Zhao                                       *
#*  zhswee@gmail.com                                                      *
#*                                                                        *
#*  This program is free software; it is licensed under the terms of the  *
#*  GNU General Public License v2 or later. See file LICENSE for details. *
#*************************************************************************/
##########################################
"""
Module containing utinity functions for data post-processing.
"""
## all exported names
__all__=['outputBoxPov','save_particleinfo','save_contactinfo','save_wallinfo','save_modelinfo','exportForceChains']

from sudodem import *
from sudodem.wrapper import *
from sudodem import  _superellipse_utils

from sudodem._superellipse_utils import *

import math,os
import numpy as np
##########################################
#some auxiliary functions
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
            print>>fobj,"disk{"+v+",r}"
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
        
def readPos(filename):     
    if os.path.isfile(filename):##not exist ,false
        fobj=open(filename)
        
        while 1:
            line=fobj.readline()
            if not line:break
            s=line.split()
            yield ([float(s[0]),float(s[1]),float(s[2])],[float(s[3]),float(s[4]),float(s[5]),float(s[6])])
        fobj.close()    
def CN(i):
         return sum([1 for m in i.intrs() if not m.phys.normalForce.norm() ==0]) 
      
def density(z=0.1):
    #calculate the number of particles below a specified height
    count = 0
    v = 0.0
    avg_cn = 0.0
    for i in O.bodies:
        if isinstance(i.shape,Superellipse):
            pos = i.state.pos
            if pos[2]<= z and pos[2]>=0.0:
                if pos[0]<= 0.35 and pos[0]>= 0.0:
                    if pos[1]<= 0.35 and pos[1]>= 0.0:
                        #yes, including this particle
                        count += 1
                        v += i.shape.getVolume()
                        avg_cn += CN(i)
                        #print vv
                        #v += vv
    #print count
    v_box = 0.35*0.35*z
    #print v
    #print v_box
    #print z
    return v/v_box,avg_cn/count
def frange(x, y, jump):
  while x < y:
    yield x
    x += jump       
def avg_density(zmin,zmax,num):
    step = (zmax-zmin)/(1.0*num)
    #print step
    zz = frange(zmin,zmax,step)
    dens = 0.0
    cn = 0.0
    for i in zz:
        data=density(z=i)
        dens += data[0]
        cn += data[1]
    return dens/num,cn/num
    
    
def CoordinateNum(BoxSize=[0.35,0.35,0.2],circle=False):
    #para:BoxSize,xlength,ylength,height are x,y,z length of the container respectively.
    #out of date para:r, the search range of a particle
    #para:meaBox,measure box which specifies which particle should be taken into account.
    #para:all,if all particles are used to calculate.
    #caution:To get the right result, one should execute 'O.run(1)' before this function.Now I also don't know the reason.
    bodies=[]
   
    #get the polyhedrons, then put them into bodies
    if circle==True:
        for i in O.bodies:
            if isinstance(i.shape,Superellipse):
                rr=math.sqrt((i.state.pos[0]-0.1)**2.0+(i.state.pos[1]-0.1)**2.0)
                if rr<=0.1:
                    bodies.append(i)
    else:
        for i in O.bodies:
            if isinstance(i.shape,Superellipse):
                bodies.append(i)
    #process the boides
    CNlist=[]   #store CN
    for i in bodies:
        count=0
        #check others around the particle
        intrs=i.intrs() #intrs around this particle
        for j in intrs:
            if j.id1==i.id:
                id2=j.id2
            else:
                id2=j.id1

            #is it polyhedron?
            if isinstance(O.bodies[id2].shape,Facet):#facet?
                count+=1
            else:
                obj1=i
                obj2=O.bodies[id2]

                #if sum([(i.state.pos[k]-j.state.pos[k])**2.0 for k in range(3)])<=r**2.0:#check particle within the range of radius r
                    #need to chech if there is an interaction between the two paticles
            if check(obj1.shape,obj2.shape,obj1.state,obj2.state) :
                    #if True,then count it.
                    count+=1
        CNlist.append(count)
    #process CN,maybe need to do a statistics
    FC=dict()   #store frenquency count
    for i in CNlist:
        if i in FC.keys():
            FC[i]=FC[i]+1
        else:
            FC[i]=1
    return FC
def FC_CN(z=0.2):
    #process the boides
    CNlist=[]   #store CN
    for i in O.bodies:
        if isinstance(i.shape,Superellipse):
            cn1 = CN(i)
            if i.state.pos[2] <= z:
                if cn1 >= 1:
                    CNlist.append(cn1)
    #process CN,maybe need to do a statistics
    FC=dict()   #store frenquency count
    for i in range(30):#
        FC[i]=0
    for i in CNlist:
        if i in FC.keys():
            FC[i]=FC[i]+1
        else:
            FC[i]=1
    return FC
def coordNumDistr(maxCN=30):
    "coordination number distribution"
    #process the boides
    CNlist=[]   #store CN
    for i in O.bodies:
        if not isinstance(i.shape,Wall):
            cn1 = CN(i)
            if cn1 >= 1:
                CNlist.append(cn1)
    #process CN,maybe need to do a statistics
    FC=dict()   #store frenquency count
    tot_p = sum(CNlist)
    for i in range(maxCN):#
        FC[i]=0
    for i in CNlist:
        if i in FC.keys():
            FC[i]=FC[i]+1.0/tot_p
        else:
            FC[i]=1/tot_p
    return FC

def exportCoordNumDistr(filename,maxCN=30):
    "export coordination number to a file"
    f = open(filename,'w')
    FC = coordNumDistr(maxCN=maxCN)
    for i in range(maxCN):
        print>>f,i,FC[i]
    f.close()
      
###############################################################################
def GetTop1(r=0.01):
    zmax = 0.
    for i in O.bodies:
        if isinstance(i.shape,Superellipse):
            if i.state.pos[2]>zmax:
                zmax = i.state.pos[2]
    if zmax==0:
            return 0.5*r
    return zmax+r
def GetTop(r=0.01):
    zmax = 0.
    for i in O.bodies:
        if isinstance(i.shape,Superellipse):
            cn = sum([1 for m in i.intrs() if not m.phys.normalForce.norm() ==0])
            if cn == 0:
                i.state.vel=(0,0,-1)
            if i.state.pos[2]>zmax:
                zmax = i.state.pos[2]
    if zmax==0:
            return 0.5*r
    return zmax+r
def vectordot(v):
    temp=np.zeros([3,3])
    modv=math.sqrt(sum([i**2. for i in v]))
    if modv!=0:
        v=[i/modv for i in v]
    for i in range(3):
        for j in range(3):
            temp[i,j]=v[i]*v[j]
    return temp
################outfut anisotropy lambda_d from the original data


def anisotropyfromRAW(vectorslist):
    subtensor=list()
    count=len(vectorslist)
    subtensor = [vectordot(v/v.norm()) for v in vectorslist]
    fabtensor=np.zeros([3,3])
    fabtensor=sum(subtensor)/count
    b,c=np.linalg.eig(fabtensor)
    return math.sqrt((b[0]-b[1])**2.+(b[0]-b[2])**2.+(b[1]-b[2])**2.)/math.sqrt(2.)  
def output_fabrics():
    #normal contacts
    vlist_nc = list()
    #branch vectors
    vlist_bv = list()
    for i in O.interactions:
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <6:#wall
            continue
        if fn.norm()>0.0:
            vlist_nc.append(fn)
            p = O.bodies[id1].state.pos - O.bodies[id2].state.pos
            vlist_bv.append(p)
    #calculate the anisotropy
    ani_nc = anisotropyfromRAW(vlist_nc)
    ani_bv = anisotropyfromRAW(vlist_bv)
    return (ani_nc,ani_bv)
def output_CNrotation():
    cns = 0.0
    rot = 0.0
    for i in O.bodies:
        if i.id>5:#not wall
            cn = CN(i)
            if cn >=2:
                cns += cn
            rot += i.state.angVel.norm()
    p_num = float(len(O.bodies))
    cns_avg = cns/p_num
    rot_avg = rot/p_num
    return (cns_avg,rot_avg)
def outorientation(a,prefix,startstep,pressure,sourcepath,destinationpath,sp=False):
    for i in a:
        filename=str(pressure)+'-'+str(i*2)
        filename=destinationpath+'/'+filename
        O.reset()
        O.load(sourcepath+'/'+prefix+str(startstep+i*1600)+'.xml.bz2')
        #vtkExporter = VTKExporter(filename)
        obj=open(filename+'.dat','w')
        if sp:
            ids=12
        else:
            ids=8
        #vtkExporter.exportInteractions(what=[('fn','i.phys.normalForce.norm()')],idskip=ids)
        for i in O.bodies:
            if i.id>=ids:
                vv=i.state.se3[1].toRotationMatrix()*i.shape.v[0]
                print>>obj,vv[0],vv[1],vv[2]
        obj.close()


def output_FC(directory = "m0.67"):
    filename=directory +'data-z0.5.dat'
    f = open(filename,'w')
    ff = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5]
    cns = list()
    for i in ff:
            O.reset()
            O.load(directory+"/"+"a"+str(i)+directory+"/finalpacking.xml.bz2")
            C=FC_CN(z=0.5)
            cns.append(C)
            print i," finished!"
            
    for C in cns:
            out = ' '        
            for i in range(30):
                    out += str(C[i])+' '
            print>>f,out 
    f.close()        
def notBoundaryP(i):
    for itr in i.intrs():
        if not itr.phys.normalForce.norm() == 0:
            #exclude the interaction with a wall
            id = min(itr.id1,itr.id2)
            if id < 6:#wall
                return False
    return True
           
def output_ori(directory = "m0.67",ff = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5]):
    cns = list()
    for ii in ff:
        fname = directory+"/"+"a"+str(ii)+directory
        f = open(fname+'/orientation.dat','w')    
        O.reset()
        O.load(fname +"/finalpacking.xml.bz2")
        ore_v = Vector3(1,0,0)#the maximum semi-axis is along (1,0,0) when generating a particle
        if ii > 1.0:
            ore_v = Vector3(0,0,1)
        for i in O.bodies:
            if isinstance(i.shape,Superellipse) and notBoundaryP(i):
                ori = i.state.se3[1].toRotationMatrix()*ore_v
                print>>f,ori[0],ori[1],ori[2]
                print>>f,-ori[0],-ori[1],-ori[2]
        print ii," finished!"
        f.close()
#######################################################################
def save_particleinfo(filename):
    "output particle info to a file"
    f = open(filename,'w')
    print>>f,"###Info of each particle:id radius x y z spin"
    for i in O.bodies:
        if isinstance(i.shape,Superellipse):#focus on superball 
            p = i.state.pos
            r = i.shape.rx
            print>>f,i.id,r,p[0],p[1],p[2],i.state.angVel.norm()
        if isinstance(i.shape,Disk):#focus on superball 
            p = i.state.pos
            r = i.shape.radius
            print>>f,i.id,r,p[0],p[1],p[2],i.state.angVel.norm()
    f.close()

def save_contactinfo(filename):
    "output contact info to a file"
    f = open(filename,'w')
    print>>f,"###c_type obj1.id obj2.id nforce(x,y,z) sforce(x,y,z)"
    print>>f,"### Info of contacts including particle-particle (named 1) and particle-wall (named 2)"
    w_contact = list()#store contacts on walls
    for itr in O.interactions:
        fn = itr.phys.normalForce
        if fn.norm() > 0.0:#repulsive force considered
            #this is a real contact
            fs = itr.phys.shearForce
            id1 = itr.id1#id1 is less than id2 by default
            id2 = itr.id2
            out = ''
            if isinstance(O.bodies[id1].shape,Wall):#wall?
                out = '2 '+str(id1)+' '+str(id2)+' '+str(fn[0])+' '+str(fn[1])+' '+str(fn[2])+' '+str(fs[0])+' '+str(fs[1])+' '+str(fs[2])
                w_contact.append(out)            
            else:#particle-particle contact            
                out = '1 '+str(id1)+' '+str(id2)+' '+str(fn[0])+' '+str(fn[1])+' '+str(fn[2])+' '+str(fs[0])+' '+str(fs[1])+' '+str(fs[2])
                print>>f, out
    for i in w_contact:
        print>>f,i
    f.close()
def save_wallinfo(filename):
    "output wall positions to a file"
    f = open(filename,'w')
    for i in range(6):
        p =O.bodies[i].state.pos        
        print>>f,p[0],p[1],p[2]
    f.close()

def save_modelinfo(path,step=0):
    "output all model info to three individual files corresponding to info of particles, contacts and walls."
    #create the path if it does not exists
    path = os.path.dirname(path+'/')
    if not os.path.exists(path):#not exists
        os.makedirs(path)
    save_particleinfo(path+'/particleinfo_'+str(step)+'.txt')
    save_contactinfo(path+'/contactinfo_'+str(step)+'.txt')
    save_wallinfo(path+'/wallinfo_'+str(step)+'.txt')
    print "Model info output finished!"

def save_strain(filename,step):
    f = open(filename,'a')
    triax=O.engines[3]
    print>>f,step,(1.0-triax.height/triax.height0)*100.0
    f.close()
def save_init(path,step):
    O.reset()
    #O.load(path+'/final_con.xml.bz2')
    path = path+'/shear'
    O.load(path+'/shear'+str(step)+'.xml.bz2')
    path = path+'/outdata'
    #save_strain(path+'/strainstep.txt',step)
    save_particleinfo(path+'/ballinfo_'+str(step)+'.txt')
    save_contactinfo(path+'/contactinfo_'+str(step)+'.txt')
    save_wallinfo(path+'/wallinfo_'+str(step)+'.txt')

################################################################
#  Contact force chains
################################################################
def PointonPlan(para,point):
    #para:[A,B,C,D], parameters of a plane equation
    #point:[x0,y0,z0], the coordinate of a point outside the given plane
    t = (sum([para[i]*point[i] for i in range(3)]) - para[3])/(sum([para[i]**2. for i in range(3)])) 
    return [point[i]-para[i]*t for i in range(3)]

def ForceChains(inputfile,ct_filename,w_filename, comment="comment"):
    '''
    inputfile: ballinfo
    ct_filename:contactinfo
    w_filename:wallinfo
    '''
    wall_file = open(w_filename, 'r')
    lines = wall_file.readlines()
    wall_file.close()
    w_pos = [i[:-2].split(' ') for i in lines]
    wall_planes=[
            [1,0,0,float(w_pos[0][0])],#wall 3,x
            [1,0,0,float(w_pos[1][0])],#wall 4,x
            [0,1,0,float(w_pos[2][1])],#wall 5,y
            [0,1,0,float(w_pos[3][1])], #wall 6,y
            [0,0,1,float(w_pos[4][2])],#wall 1,z
                [0,0,1,float(w_pos[5][2])]#wall 2,z
            ]
    f = open(inputfile,'r')
    # output file
    fName = ct_filename +'.vtp'
    outContactFile = open(fName, 'w')
    lines = f.readlines()[1:]########## the range may need to be changed
    f.close()
    #open the original file including the info of all balls
    nBodies = 0
    radius = dict()
    position =dict()
    for l in lines:
        l1 = l.lstrip()
        l1 = l1[:-2].split(' ')
        #print l1
        b_id = int(l1[0])
        position[b_id] = [float(l1[2]),float(l1[3]),float(l1[4])]
        radius[b_id] = float(l1[1])
        nBodies += 1
    # output file
    contact_file = open(ct_filename, 'r')
    lines = contact_file.readlines()[2:]
    contact_file.close()
    ##########################################################################
    ###contacts
    nIntrs = len(lines)    
    # head
    outContactFile.write("<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>\n<PolyData>\n")
    outContactFile.write("<Piece NumberOfPoints='%s' NumberOfVerts='0' NumberOfLines='%s' NumberOfStrips='0' NumberOfPolys='0'>\n"%(str(2*nIntrs),str(nIntrs)))
    # write coords of intrs bodies (also taking into account possible periodicity
    outContactFile.write("<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n")
    for l in lines:
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        #print l1[:3],'sss',l1[3:]
        [contact_type,id1,id2] = [int(i) for i in l1[:3]]
        #[fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        #[fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        if contact_type < 2:#1: ball-ball contact
            #find positions of two touching balls
            pos = position[id1]
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
            pos = position[id2]
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
        else:#2 ball_wall contact
            
            pos = PointonPlan(wall_planes[id1],position[id2])    ###???
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2])) 
            pos = position[id2]
            outContactFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
        #print contact_type,id1,id2,fn_x,fn_y,fn_z,fs_x,fs_y,fs_z
    
    outContactFile.write("</DataArray>\n</Points>\n<Lines>\n<DataArray type='Int32' Name='connectivity' format='ascii'>\n")
        
    ss=''
    for con in range(2*nIntrs):
        ss+=' '+str(con)
    outContactFile.write(ss+'\n')
    outContactFile.write("</DataArray>\n<DataArray type='Int32' Name='offsets' format='ascii'>\n")
    ss=''
    for con in range(nIntrs):
        ss+=' '+str(con*2+2)
    outContactFile.write(ss)
    outContactFile.write("\n</DataArray>\n</Lines>\n")
    ##
    name = 'Fn'
    outContactFile.write("<PointData Scalars='%s'>\n<DataArray type='Float32' Name='%s' format='ascii'>\n"%(name,name))
    for l in lines:
        l=l.lstrip()
        l1 = l[:-2].split(' ')
        #[contact_type,id1,id2] = [int(i) for i in l1[:3]]
        [fn_x,fn_y,fn_z,fs_x,fs_y,fs_z] = [float(i) for i in l1[3:]]
        fn=(fn_x**2.+fn_y**2.+fn_z**2.)**0.5
        outContactFile.write("%g %g\n"%(fn,fn))
    outContactFile.write("</DataArray>\n</PointData>")
    outContactFile.write("\n</Piece>\n</PolyData>\n</VTKFile>")
    outContactFile.close()


def exportDisk(inputfile, comment="comment"):
    "export disks to a VTK file."
    f = open(inputfile,'r')
    outDiskFile = open(inputfile+'.vtk', 'w')
    lines = f.readlines()[1:]########## the range may need to be changed
    f.close()
    #open the original file including the info of all balls
    nBodies = 0
    radius = dict()
    position =dict()
    for l in lines:
        l1 = l.lstrip()
        l1 = l1[:-2].split(' ')
        #print l1
        b_id = int(l1[0])
        position[b_id] = [float(l1[2]),float(l1[3]),float(l1[4])]
        radius[b_id] = float(l1[1])
        nBodies += 1
    # head
    outDiskFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n"%(comment,nBodies))
    # write position of disks
    for key in position.keys():
        pos = position[key]
        outDiskFile.write("%g %g %g\n"%(pos[0],pos[1],pos[2]))
    # write radius
    outDiskFile.write("\nPOINT_DATA %d\nSCALARS radius double 1\nLOOKUP_TABLE default\n"%(nBodies))
    for key in position.keys():
        outDiskFile.write("%g\n"%(radius[key]))

    outDiskFile.close()


def exportForceChains(path,step=0):
    "export force chains with path of input files and step number."
    path = os.path.dirname(path+'/')
    if not os.path.exists(path):#not exists
        print "The typed path dose not exist!"
        return False
    particlefile = os.path.join(path,'particleinfo_'+str(step)+'.txt')
    contactfile = os.path.join(path,'contactinfo_'+str(step)+'.txt')
    wallfile = os.path.join(path,'wallinfo_'+str(step)+'.txt')
    #check if files exist
    if not (os.path.isfile(particlefile) and os.path.isfile(contactfile) and os.path.isfile(wallfile)):
        print "Input files are missing! Please check the path of particleinfo*.txt, contactinfo*.txt and wallinfo*.txt is set correctly."
        return False
    ForceChains(particlefile,contactfile,wallfile, comment="comment")
    return True
################################################################
##Zhao, Shiwei, and Xiaowen Zhou. 2017. "Effects of Particle Asphericity on the Macro-
##and Micro-Mechanical Behaviors of Granular Assemblies." Granular Matter 19 (2):38.
##Zhao, Shiwei, T. Matthew Evans, and Xiaowen Zhou. 2018. "Three-Dimensional Voronoi 
##Analysis of Monodisperse Ellipsoids during Triaxial Shear." Powder Technology 323:323-36.


class VTK3Dhist():
    "A class for 3d histogram exporting to a vtk file for visulization with Paraview."
    def __init__(self,dimension = 8):
        #self.file_in,self.file_outpath,self.file_outname = file_in, file_outpath,file_outname
        self.dimension = dimension
        self.dimension2 = self.dimension*self.dimension
        self.elements_num = 6*self.dimension2
        self.delta = 2.0/self.dimension
        self.xx = 0
        self.yy = 0
        self.zz = 0
    #mapping from cube to disk
    #A direct mapping from cube to (its inscribed) disk: 
    #http://catlikecoding.com/unity/tutorials/cube-disk/
    def cube2sp(self,x,y,z):
        m,n = x.shape
        xx = np.zeros([m,n])
        yy = np.zeros([m,n])
        zz = np.zeros([m,n])
        for i in range(m):
            for j in range(n):
                x1 = x[i,j]
                y1 = y[i,j]
                z1 = z[i,j]
                #print 1.-x1**2/2-z1**2/2+x1**2.*z1**2/3
                xx[i,j] = x1*np.sqrt(1.-y1**2/2-z1**2/2+y1**2.*z1**2/3)
                yy[i,j] = y1*np.sqrt(1.-x1**2/2-z1**2/2+x1**2.*z1**2/3)
                zz[i,j] = z1*np.sqrt(1.-y1**2/2-x1**2/2+y1**2.*x1**2/3)
        #return xx,yy,zz
        self.xx = xx
        self.yy = yy
        self.zz = zz

    def generate_cube(self):
        " "    
        #we construct a 3-d array for storing info from six faces of the cube.
        #x=0,x=1,y=0,y=1,z=0,z=1
        #cubeelements = zeros(6,self.dimension,self.dimension);%
        #cubenodes = zeros(6,self.dimension+1,self.dimension+1);
        #nodes_num = 6*(self.dimension-1)**2+12*(self.dimension-1)+8;
        #elements_num = 6*self.dimension2
        #update the private member variables
        self.dimension2 = self.dimension*self.dimension
        self.elements_num = 6*self.dimension2
        self.delta = 2.0/self.dimension

        elements = np.zeros([12,self.elements_num])
        #
        face_elements = np.zeros([2,self.dimension2])
        
        for i in range(self.dimension):
            for j in range(self.dimension):
                face_elements[:,i*self.dimension+j]=np.array([i*self.delta,j*self.delta]) 

        one = np.ones(self.dimension2);

        face1 = face_elements[0,:];
        face2 = face_elements[1,:];
        #num = self.dimension*self.dimension
        #print face1,face2,one
        temp = np.array([one,face1,face2])
        #print 'temp = ', temp
        indices = [0,0,1,1,2,2]
        for i in range(6):
            I = np.eye(3)
            index =  indices[i] 
            tmp =  np.copy(I[:,index])#Caution: The problem is that Numpy basic slicing does not create copies of the actual data
            #print "tmp==",tmp
            I[:, index]=I[:,0]
            I[:,0]=tmp
            #print 'tmp2=',tmp
            #print "i====",i+1,index
            #print I
            #print "I temp=",I.dot(np.array([one,face1+self.delta,face2]))
            I[index,0]=np.mod(i,2)*2
            elements[0:3,i*self.dimension2:(i+1)*self.dimension2]=I.dot(temp)
            elements[3:6,i*self.dimension2:(i+1)*self.dimension2]=I.dot(np.array([one,face1+self.delta,face2]))
            elements[6:9,i*self.dimension2:(i+1)*self.dimension2]=I.dot(np.array([one,face1+self.delta,face2+self.delta]))
            elements[9: ,i*self.dimension2:(i+1)*self.dimension2]=I.dot(np.array([one,face1,face2+self.delta]))
    
        elements = elements-np.ones([12,self.elements_num])
        cube_x = elements[[0,3,6,9],:]
        cube_y = elements[[1,4,7,10],:]
        cube_z = elements[[2,5,8,11],:]
    
        # =========================================================================
        # Mapping nodes from cube to disk.
        # =========================================================================
        #print cube_x
        #print cube_y
        self.cube2sp(cube_x,cube_y,cube_z)
        #self.writeVTK('cube.vtk',self.xx,self.yy,self.zz)
    def transform(self,x,y,z):
        #print x.shape
        NumData=len(x)
        #print NumData
        #transform the coordinates (x,y,z) on a disk to a cubical surface (cube_x,cube_y,cube_z)
        cube_x=np.zeros(NumData)#coordinate x in the cubical system
        cube_y=np.zeros(NumData)# y
        cube_z=np.zeros(NumData)# z
        for i in range(NumData):
            #find the maximum of |x[i]|,|y[i]| and |z[i]|
            abc = [abs(x[i]),abs(y[i]),abs(z[i])]
            index = [0,1,2]#the first one corresponds to the maximum of a,b,and c
            if abc[0] < abc[1]:#swap the indices of a and b
                index[0],index[1] = index[1],index[0]
                if abc[1] < abc[2]:
                    index[0],index[2] = index[2],index[0]
            elif abc[0] < abc[2]:
                index[0],index[2] = index[2],index[0]
            xyz = [x[i],y[i],z[i]]
            xyz_t = [xyz[j] for j in index]
            xyz2 = [2.0*j**2 for j in xyz_t]
            tmp = [0.0,0.0,0.0]
            xyz_r = [0.,0.,0.]
            #print xyz2
            delta = -np.sqrt((-xyz2[1]+xyz2[2]-3)**2-12*xyz2[1])
            tmp[0] = 1.0
            tmp[1] = np.sqrt((delta+xyz2[1]-xyz2[2]+3)/2)
            tmp[2] = np.sqrt((delta-xyz2[1]+xyz2[2]+3)/2)
            for j in range(3):
                xyz_r[index[j]] = tmp[j]*np.sign(xyz[index[j]])
            #print xyz_r,i
            cube_x[i],cube_y[i],cube_z[i] = xyz_r
        return cube_x,cube_y,cube_z

    def local_index(self,x,delta):#find the element index that x falls in
        return np.fix((x+1)/self.delta)+1 #index of the element array
    def find_j_ab(self,n,x,y,z):
        a1 = (n+0.5*(x+1))*self.dimension**2
        a2 = self.dimension*(self.local_index(z,self.delta)-1)+self.local_index(y,self.delta)
        return int(a1+a2)
          
    def find_j(self,x,y,z):
        #x_ind = self.local_index(x,self.delta)
        #y_ind = self.local_index(y,self.delta)
        #z_ind = self.local_index(z,self.delta)
        n = 0
        if abs(y) == 1:
            x,y = y,x
            n = 2
        elif abs(z) == 1:
            x,z = z,x
            n = 4
        return self.find_j_ab(n,x,y,z),self.find_j_ab(n,-x,-y,-z)
        ##done
        if (0):
            if abs(x) == 1:#face 1 when x=-1, and face 2 when x=1
                #y_ind = self.local_index(y,self.delta)
                #z_ind = self.local_index(z,self.delta)
                #y_ind1 = self.local_index(-y,self.delta)
                #z_ind1 = self.local_index(-z,self.delta)
                #a = int(0.5*(x+1)*self.dimension**2+self.dimension*(z_ind-1)+y_ind)
                #b = int(0.5*(1-x)*self.dimension**2+self.dimension*(z_ind1-1)+y_ind1)
                #a1 = 0.5*(x+1)*self.dimension**2
                #a2 = self.dimension*(self.local_index(z,self.delta)-1)+self.local_index(y,self.delta)
                a = self.find_j_ab(0,x,y,z)
                #b1 = 0.5*(1-x)*self.dimension**2
                #b2 = self.dimension*(self.local_index(-z,self.delta)-1)+self.local_index(-y,self.delta)
                b = self.find_j_ab(0,-x,-y,-z)
                return a,b

            if abs(y) == 1:#face 3 when y=-1, and face 4 when y=1
                #x_ind = self.local_index(x,self.delta)
                #z_ind = self.local_index(z,self.delta)
                #x_ind1 = self.local_index(-x,self.delta)
                #z_ind1 = self.local_index(-z,self.delta)
                #return int((2+0.5*(y+1))*self.dimension**2+self.dimension*(z_ind-1)+x_ind),int((2+0.5*(1-y))*self.dimension**2+self.dimension*(z_ind1-1)+x_ind1)
                a = self.find_j_ab(2,y,x,z)
                b = self.find_j_ab(2,-y,-x,-z)
                return a,b

            if abs(z) == 1:#face 5 when z = -1, and face 6 when z=1
                #y_ind = self.local_index(y,self.delta)
                #x_ind = self.local_index(x,self.delta)
                #y_ind1 = self.local_index(-y,self.delta)
                #x_ind1 = self.local_index(-x,self.delta)
                #return int((4+0.5*(z+1))*self.dimension**2+self.dimension*(x_ind-1)+y_ind),int((4+0.5*(1-z))*self.dimension**2+self.dimension*(x_ind1-1)+y_ind1)
                a = self.find_j_ab(4,z,y,x)
                b = self.find_j_ab(4,-z,-y,-x)
                return a,b



    def writeVTK(self,filename,xx,yy,zz):
        print("writing Polydata into a VTK file")
        f_out = open(filename,'w')
    
        #writing the file head
        print>>f_out,"# vtk DataFile Version 2.0"
        print>>f_out,"Sway data processing"
        print>>f_out,"ASCII"
        print>>f_out,"DATASET POLYDATA"

        m,n = xx.shape
        print>>f_out,"POINTS   ", n*4," float"  ###index is from 0 in VTK
        #output points
        for i in range(n):#n elements or polygons
            for j in range(m):
                print>>f_out, xx[j,i],yy[j,i],zz[j,i]
        #output polygons
        print>>f_out,"POLYGONS ", n,n*5  #polygon_num polygon_num*5
        for i in range(n):
            print>>f_out,"4 ",i*4,i*4+1,i*4+2,i*4+3
        f_out.close()
    def drawBar(self,p1,p2,p3,p4,p5):
        #p1~p4: vertices on the top face
        #p5: the origin
        #convert to string
        p1_s = ' ' + str(p1)
        p2_s = ' ' + str(p2)
        p3_s = ' ' + str(p3)
        p4_s = ' ' + str(p4)
        p5_s = ' ' + str(p5)

        line = ""
        #top face
        line += "4 "+ p1_s + p2_s + p3_s + p4_s
        #side faces
        line += " 3 "+ p5_s + p2_s + p1_s  #5 2 1
        line += " 3 "+ p5_s + p1_s + p4_s  #5 1 4
        line += " 3 "+ p5_s + p4_s + p3_s  #5 4 3
        line += " 3 "+ p5_s + p3_s + p2_s  #5 3 2
        return line


    def write3DhistogramVTK(self,filename,coeff):
        print("writing Polydata into a VTK file")
        f_out = open(filename,'w')
    
        #writing the file head
        print>>f_out,"# vtk DataFile Version 2.0"
        print>>f_out,"Sway data processing"
        print>>f_out,"ASCII"
        print>>f_out,"DATASET POLYDATA"

        m,n = self.xx.shape
        print>>f_out,"POINTS   ", n*4 + 1," float"  ###index is from 0 in VTK; here we include the origin
        #output points
        print>>f_out, 0,0,0  #output the origin
        for i in range(n):#n elements or polygons
            for j in range(m):
                print>>f_out, self.xx[j,i]*coeff[i],self.yy[j,i]*coeff[i],self.zz[j,i]*coeff[i]
        #output polygons
        print>>f_out,"POLYGONS ", n*5,n*21  #polygon_num polygon_num*5
        for i in range(n):
            #print>>f_out,"4 ",i*4,i*4+1,i*4+2,i*4+3
            print>>f_out,self.drawBar(i*4+1,i*4+2,i*4+3,i*4+4,0)
        #output magnitude
        print>>f_out,"CELL_DATA ", n*5  #
        print>>f_out,"SCALARS cell_scalars float 1"
        print>>f_out,"LOOKUP_TABLE default"
        #print values of all cells (polygons)
        for i in range(n):
            print>>f_out,coeff[i],coeff[i],coeff[i],coeff[i],coeff[i]
        f_out.close()

    def VTK3DhistogramCN(self,file_input,file_output,shift):#for orientation of Voronoi cell
        #shift = 3; % 0 v_max, shortest axis direction; 3 v_min, longest axis direction
        self.dimension = 10;
        self.generate_cube()#we construct a 3-d array for storing info from six faces of the cube.
        # Reading contacts' orientations.
        orient = np.loadtxt(file_input,delimiter=' ', skiprows=1)
        #v_max and v_min are imported
        x=orient[:,1+shift];
        y=orient[:,2+shift];
        z=orient[:,3+shift];
        #print x.shape
        NumData=len(x)
        cube_x,cube_y,cube_z = self.transform(x,y,z)
     
        # Cheking to see an orientation is inside which element in cubical system.
        
        coeff = np.zeros(self.elements_num)
        for i in range(NumData):
            j = self.find_j(cube_x[i],cube_y[i],cube_z[i],dimension)
            #print j
            coeff[j-1] = coeff[j-1]+1;     # Number of normal forces pointing within an element.


        coeff= (coeff/NumData)/(4.0*np.pi/self.elements_num)
        if shift == 0:#v_max
            suffix = 'shortest.vtk'
        elif shift == 3:#v_min
            suffix = 'longest.vtk'
        else:
            suffix ='.vtk'
        self.write3DhistogramVTK(file_output+suffix,coeff)



    def VTK3Dhistogram(self,file_input,file_outputpath,file_outputname):
        "file_input: path of the input file.\nfile_outpath: directory of the output file.\nfile_outputname: name of the output file."
        shift = 0
        #shift = 3; % 0 v_max, shortest axis direction; 3 v_min, longest axis direction
        self.generate_cube()
        # =========================================================================
        # Reading contacts' orientations.
        data = np.loadtxt(file_input,delimiter=' ', skiprows=2) #No. c_type object1.id object2.id nforce(x,y,z) sforce(x,y,z)
        #find the range of data
        ids = data[:,0]
        m = int(2*len(ids)-np.sum(ids)) #total number of lines with particle-particle contacts
        #data1 = data[:m,3+shift:6+shift]
        fn_v = data[:m,3:6]
        ft_v = data[:m,6:9]
        fn = np.sqrt(fn_v[:,0]**2.0+fn_v[:,1]**2.0+fn_v[:,2]**2.0)
        ft = np.sqrt(ft_v[:,0]**2.0+ft_v[:,1]**2.0+ft_v[:,2]**2.0)
        avg_fn1 = np.average(fn)
    
        #v_max and v_min are imported
        #contact normal: [x,y,z]
        x=fn_v[:,0]/fn;
        y=fn_v[:,1]/fn;
        z=fn_v[:,2]/fn;
        for i in range(len(fn)):
            if fn[i]>10.0*avg_fn1:#the value larger than 10<Fn> has a pretty low probability of less than 1e-5
                fn[i]=10.0*avg_fn1
        #print "size of data: ",orient.shape
        #concatenate vectors 
        #x = np.concatenate((x,-x))
        #y = np.concatenate((y,-y))
        #z = np.concatenate((z,-z))
    
        ##########
        NumData = len(x)
        cube_x,cube_y,cube_z = self.transform(x,y,z)
        # Checking to see an orientation is inside which element in cubical system.
    
        cn_coeff = np.zeros(self.elements_num)
        fn_coeff = np.zeros(self.elements_num)
        ft_coeff = np.zeros(self.elements_num)
        #print elements_num
        for i in range(NumData):
            #print mappx[i],mappy[i],mappz[i]
            j1,j2 = self.find_j(cube_x[i],cube_y[i],cube_z[i])
            #contact normal
            cn_coeff[j1-1] += 1     # Number of normal forces pointing within an element.
            cn_coeff[j2-1] += 1
            #normal contact force
            fn_coeff[j1-1] += fn[i]     # normal forces pointing within an element.
            fn_coeff[j2-1] += fn[i]
            #tangential contact force
            ft_coeff[j1-1] += ft[i]     # Number of normal forces pointing within an element.
            ft_coeff[j2-1] += ft[i]

        cn_coeff1= (cn_coeff/NumData/2.0)/(4.0*np.pi/self.elements_num)
    
        fn_coeff = fn_coeff/cn_coeff
        avg_fn = sum(fn_coeff)/self.elements_num
        fn_coeff /=avg_fn
        ft_coeff = ft_coeff/cn_coeff/avg_fn
    
        fname = file_outputpath+'/cn'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,cn_coeff1)
        fname = file_outputpath+'/fn'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,fn_coeff)
        fname = file_outputpath+'/ft'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,ft_coeff)


    def VTK3Dhistogram2(self,file_input,file_outputpath,file_outputname):
        shift = 0
        #shift = 3; % 0 v_max, shortest axis direction; 3 v_min, longest axis direction
        self.generate_cube()
        # =========================================================================
        # Reading contacts' orientations.
        data = np.loadtxt(file_input,delimiter=' ', skiprows=2) #No. c_type object1.id object2.id nforce(x,y,z) sforce(x,y,z)
        #find the range of data
        ids = data[:,0]
        m = 2*len(ids)-np.sum(ids) #total number of lines with particle-particle contacts
        #data1 = data[:m,3+shift:6+shift]
        fn_v = data[:m,3:6]
        ft_v = data[:m,6:9]
        fn1 = np.sqrt(fn_v[:,0]**2.0+fn_v[:,1]**2.0+fn_v[:,2]**2.0)
        ft1 = np.sqrt(ft_v[:,0]**2.0+ft_v[:,1]**2.0+ft_v[:,2]**2.0)
        #v_max and v_min are imported
        #contact normal: [x,y,z]
        x = list()
        y = list()
        z = list()
        fn = list()
        ft = list()
        for i in range(len(ft1)):
            #print ft_v[i,0],ft[i]
            if ft1[i] >0:
                x.append(ft_v[i,0]/ft1[i])
                y.append(ft_v[i,1]/ft1[i])
                z.append(ft_v[i,2]/ft1[i])
                fn.append(fn1[i])
                ft.append(ft1[i])
                #print "ft=0"
            #print ft_v[i,0]/ft[i]
        #x=ft_v[:,0]/ft;
        #y=ft_v[:,1]/ft;
        #z=ft_v[:,2]/ft;
        #print "size of data: ",orient.shape
        #concatenate vectors 
        #x = np.concatenate((x,-x))
        #y = np.concatenate((y,-y))
        #z = np.concatenate((z,-z))
    
        #print x.shape
        NumData=len(x)
        #print NumData
        cube_x,cube_y,cube_z = self.transform(x,y,z)
        # Cheking to see an orientation is inside which element in cubical system.
    
        cn_coeff = np.zeros(elements_num)
        fn_coeff = np.zeros(elements_num)
        ft_coeff = np.zeros(elements_num)
        #print elements_num
        for i in range(NumData):
            j1,j2 = self.find_j(cube_x[i],cube_y[i],cube_z[i])
            #contact normal
            cn_coeff[j1-1] += 1     # Number of normal forces pointing within an element.
            cn_coeff[j2-1] += 1
            #normal contact force
            fn_coeff[j1-1] += fn[i]     # normal forces pointing within an element.
            fn_coeff[j2-1] += fn[i]
            #tangential contact force
            ft_coeff[j1-1] += ft[i]     # Number of normal forces pointing within an element.
            ft_coeff[j2-1] += ft[i]

        cn_coeff1= (cn_coeff/NumData/2.0)/(4.0*np.pi/elements_num)
        avg_fn = np.average(fn)
        fn_coeff = fn_coeff/cn_coeff/avg_fn
        ft_coeff = ft_coeff/cn_coeff/avg_fn
    
        #fname = file_outputpath+'/cn'+file_outputname+'.vtk'
        #self.write3DhistogramVTK(fname,xx,yy,zz,cn_coeff1)
        #fname = file_outputpath+'/fn'+file_outputname+'.vtk'
        #self.write3DhistogramVTK(fname,xx,yy,zz,fn_coeff)
        fname = file_outputpath+'/ft'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,ft_coeff)
#############################excute####################
########the following is for the asphericity project
#eta = '1.0'
#numlist = [1,2,12,130]
#eta = '1.05'
#numlist = [1,2,12,130]
#eta = '1.2'
#numlist = [1,12,128]
#eta = '1.4'
#numlist = [1,14,127]
#eta = '1.6'
#numlist = [1,13,126]
#for num in numlist:
#    file_in = eta+'/shear/outdata/contactinfo_'+str(num)+'.txt'
    #file_out = '3DhistOutput'+'/'+eta+'-'+str(num)
#    file_outpath,file_outname = 'histogramtest/',eta+'-'+str(num)
#    VTK3Dhistogram(file_in,file_outpath,file_outname)
#file_in = "contactinfo_12.txt"
#file_outpath = "./"
#file_outname = "3dhistg"
#hist = VTK3Dhist()
#hist.VTK3Dhistogram(file_in,file_outpath,file_outname)
##########the following is for the aspectratio project
#d = [10,12,15,17,20]
#num = [0,130]
#shift = [0,3]
#for i in d:
#    for j in num:
#        file_in = str(i)+'/'+str(j)+'reduced.txt'
#        file_out = '3DhistOutput'+'/'+str(i)+'-'+str(j)
#        for k in shift:
#            print i,j,k
#            VTK3DhistogramCN(file_in,file_out,k)


################################################################
def out_anisotropy_cn(path,steps):
    f = open(path+'/anisoCn.txt','w')
    for m in range(steps):
        O.reset()
        print path,m
        #O.load(path+'/final_con.xml.bz2')
        O.load(path+'/shear/shear'+str(m+1)+'.xml.bz2')
        #path = path+'/outdata'
        #find CompressionEngine
        triax = O.engines[3]
        for e in O.engines:
            if isinstance(e,CompressionEngine):
                triax = e
                break
        zstrain = (1.0-triax.height/triax.height0)*100.0
        ani_nc,ani_bv = output_fabrics()
        cns_avg,rot_avg = output_CNrotation()
        print>>f,m+1,zstrain,ani_nc,ani_bv,cns_avg,rot_avg
    f.close()
def out_porosity(path='.'):
    f = open(path+'/porosity.txt','w')
    steps = [130,130,128,127,126]
    mm= ['1.0','1.05','1.2','1.4','1.6']
    for i in range(len(steps)):
        m = mm[i]
        O.reset()
        print path,m
        #O.load(path+'/'+m+'/final_con.xml.bz2')
        O.load(path+'/'+m+'/shear/shear'+str(steps[i])+'.xml.bz2')
        #path = path+'/outdata'
        #find CompressionEngine
        triax = O.engines[3]
        for e in O.engines:
            if isinstance(e,CompressionEngine):
                triax = e
                break
        zstrain = (1.0-triax.height/triax.height0)*100.0
        v=0.0
        for i in O.bodies:
            if isinstance(i.shape,Superellipse):
                v += i.shape.getVolume()
        box_v = triax.height*triax.width*triax.depth
        print>>f,1.0-v/box_v
    f.close()   
#######################################################
########## new anisotropy: stress-force-fabric
#######################################################
def tensor_deviator(tensor):
    I = np.zeros((3,3))
    I[0,0]=I[1,1]=I[2,2]=1.0
    I1 = np.trace(tensor)/3.0
    return tensor - I1*I
def tensor_dot_tensor(a,b):#return a scalar
    return np.einsum('ij,ij',a,b)
def vector_dyadic(v1,v2):
    tensor = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            tensor[i,j] = v1[i]*v2[j]
    return tensor    
def new_anisotropy(ani_tensor,deviator_stress):
    sign = np.sign(tensor_dot_tensor(ani_tensor,deviator_stress))
    return sign*math.sqrt(1.5*tensor_dot_tensor(ani_tensor,ani_tensor))
def calc_anisotropy():
    #get contact force and branch vectors
    #normal contacts
    vlist_fn = list()
    vlist_fs = list()
    #branch vectors
    vlist_bv = list()
    avg_fn = 0.0
    vlist_fsfn = list()
    c_num_strong = 0.0
    for i in O.interactions:
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <6:#wall
            continue
        if fn.norm()>0.0:
            vlist_fn.append(fn)
            fs = i.phys.shearForce
            vlist_fsfn.append(fs.norm()/fn.norm())
            avg_fn += fn.norm()
            vlist_fs.append(fs)
            p = O.bodies[id1].state.pos - O.bodies[id2].state.pos
            vlist_bv.append(p)
    c_num=len(vlist_fn)
    avg_fn /=c_num
    for f in vlist_fn:
        if f.norm() >= avg_fn:#strong contact
            c_num_strong += 1
    c_num_sliding = sum([1 for eta in vlist_fsfn if eta >= 0.499999999])
    pro_sliding = float(c_num_sliding)/c_num #proportion of sliding contacts
    c_num_weak = c_num - c_num_strong
    pro_weak = float(c_num_weak)/c_num #proportion of weak contacts
    #get stress tensor
    tensor_stress = np.zeros([3,3])
    for c in range(c_num):
        c_f = vlist_fn[c] + vlist_fs[c]
        c_bv = -vlist_bv[c]
        
        v_dya1 = vector_dyadic(c_f,c_bv)
       
        tensor_stress += v_dya1
        
    tensor_stress = tensor_stress/c_num
    p = np.einsum('ii',tensor_stress)/3.0
    q1 = tensor_deviator(tensor_stress)
    q = math.sqrt(1.5*tensor_dot_tensor(q1,q1))   
    #fabric tensor of contact normals
    tensor_cn = np.zeros([3,3])
    tensor_cn_strong = np.zeros([3,3])
    tensor_cn_weak = np.zeros([3,3])
    for v in vlist_fn:
        modv=v.norm()
        if modv!=0:
            v=[i/modv for i in v]
        dya1 = vector_dyadic(v,v)
        tensor_cn += dya1
        if modv >= avg_fn:#strong contact
            tensor_cn_strong += dya1
        else:#weak contact
            tensor_cn_weak += dya1
    
    tensor_cn = tensor_cn/c_num
    tensor_cn_strong /= c_num_strong
    tensor_cn_weak /= c_num_weak
    #deviator of cn
    ani_cn = tensor_deviator(tensor_cn)*15.0/2.0
    ani_cn_strong = tensor_deviator(tensor_cn_strong)*15.0/2.0
    ani_cn_weak = tensor_deviator(tensor_cn_weak)*15.0/2.0
    ################
    #anisotropy of contact force
    #fabric tensor of normal contact forces
    tensor_fn = np.zeros([3,3])
    tensor_fs = np.zeros([3,3])
    tensor_fn_strong = np.zeros([3,3])
    tensor_fs_strong = np.zeros([3,3])
    tensor_fn_weak = np.zeros([3,3])
    tensor_fs_weak = np.zeros([3,3])
    for c in range(c_num):
        c_fn = vlist_fn[c]
        c_fs = vlist_fs[c]
        fn1 = c_fn.norm()
        fs1 = c_fs.norm()
        
            
        c_fn=c_fn/fn1
        c_fs=c_fs/fs1
        v_dya1 = vector_dyadic(c_fn,c_fn)
        v_dya2 = vector_dyadic(c_fs,c_fn)
        
        tensor_fn += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        tensor_fs += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        if fn1 >= avg_fn:#strong contact
            tensor_fn_strong += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn_strong,v_dya1))
            tensor_fs_strong += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn_strong,v_dya1))
        else:
            tensor_fn_weak += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn_weak,v_dya1))
            tensor_fs_weak += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn_weak,v_dya1))
    tensor_fn = tensor_fn/c_num
    tensor_fs = tensor_fs/c_num
    tensor_fn_strong /= c_num_strong
    tensor_fs_strong /= c_num_strong
    tensor_fn_weak /= c_num_weak
    tensor_fs_weak /= c_num_weak
    #deviator of fn
    f0_bar = np.einsum('ii',tensor_fn)
    ani_fn = tensor_deviator(tensor_fn)*15.0/2.0/f0_bar
    ani_fs = tensor_deviator(tensor_fs)*15.0/3.0/f0_bar
    f0_bar_strong = np.einsum('ii',tensor_fn_strong)
    ani_fn_strong = tensor_deviator(tensor_fn_strong)*15.0/2.0/f0_bar_strong
    ani_fs_strong = tensor_deviator(tensor_fs_strong)*15.0/3.0/f0_bar_strong
    f0_bar_weak = np.einsum('ii',tensor_fn_weak)
    ani_fn_weak = tensor_deviator(tensor_fn_weak)*15.0/2.0/f0_bar_weak
    ani_fs_weak = tensor_deviator(tensor_fs_weak)*15.0/3.0/f0_bar_weak
    #b,c=np.linalg.eig(fabtensor)
    ##################
    ###anisotropy of branch vectors
    #fabric tensor of branch vectors
    tensor_bvn = np.zeros([3,3])
    tensor_bvs = np.zeros([3,3])
    tensor_bvn_strong = np.zeros([3,3])
    tensor_bvs_strong = np.zeros([3,3])
    tensor_bvn_weak = np.zeros([3,3])
    tensor_bvs_weak = np.zeros([3,3])
    for c in range(c_num):
        c_fn = vlist_fn[c]
        c_bv = vlist_bv[c]
        fn1 = c_fn.norm()
        c_bv_n = c_bv.dot(c_fn)*c_fn/fn1/fn1###project bv on cn
        c_bv_s = c_bv - c_bv_n
        
        bvn1 = c_bv_n.norm()
        bvs1 = c_bv_s.norm()
        c_fn=c_fn/fn1
        c_bv_n=c_bv_n/bvn1
        c_bv_s=c_bv_s/bvs1
        v_dya1 = vector_dyadic(c_fn,c_fn)
        v_dya2 = vector_dyadic(c_bv_s,c_fn)
        tensor_bvn += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        tensor_bvs += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        if fn1 >= avg_fn:#strong contact
            tensor_bvn_strong += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn_strong,v_dya1))
            tensor_bvs_strong += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn_strong,v_dya1))
        else:
            tensor_bvn_weak += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn_weak,v_dya1))
            tensor_bvs_weak += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn_weak,v_dya1))
    tensor_bvn = tensor_bvn/c_num
    tensor_bvs = tensor_bvs/c_num
    tensor_bvn_strong /= c_num_strong
    tensor_bvs_strong /= c_num_strong
    tensor_bvn_weak /= c_num_weak
    tensor_bvs_weak /= c_num_weak
    #deviator of fn
    b0_bar = np.einsum('ii',tensor_bvn)
    ani_bvn = tensor_deviator(tensor_bvn)*15.0/2.0/b0_bar
    ani_bvs = tensor_deviator(tensor_bvs)*15.0/3.0/b0_bar  
    b0_bar_strong = np.einsum('ii',tensor_bvn_strong)
    ani_bvn_strong = tensor_deviator(tensor_bvn_strong)*15.0/2.0/b0_bar_strong
    ani_bvs_strong = tensor_deviator(tensor_bvs_strong)*15.0/3.0/b0_bar_strong
    b0_bar_weak = np.einsum('ii',tensor_bvn_weak)
    ani_bvn_weak = tensor_deviator(tensor_bvn_weak)*15.0/2.0/b0_bar_weak
    ani_bvs_weak = tensor_deviator(tensor_bvs_weak)*15.0/3.0/b0_bar_weak
    #########output the deviatoric invariants
    base_tensor = q1
    ani_cn1 = new_anisotropy(ani_cn,base_tensor)  
    ani_bvn1 = new_anisotropy(ani_bvn,base_tensor) 
    ani_bvs1 = new_anisotropy(ani_bvs,base_tensor) 
    ani_fn1 = new_anisotropy(ani_fn,base_tensor) 
    ani_fs1 = new_anisotropy(ani_fs,base_tensor)
    ani_cn1_strong = new_anisotropy(ani_cn_strong,base_tensor)  
    ani_bvn1_strong = new_anisotropy(ani_bvn_strong,base_tensor) 
    ani_bvs1_strong = new_anisotropy(ani_bvs_strong,base_tensor) 
    ani_fn1_strong = new_anisotropy(ani_fn_strong,base_tensor) 
    ani_fs1_strong = new_anisotropy(ani_fs_strong,base_tensor)
    ani_cn1_weak = new_anisotropy(ani_cn_weak,base_tensor)  
    ani_bvn1_weak = new_anisotropy(ani_bvn_weak,base_tensor) 
    ani_bvs1_weak = new_anisotropy(ani_bvs_weak,base_tensor) 
    ani_fn1_weak = new_anisotropy(ani_fn_weak,base_tensor) 
    ani_fs1_weak = new_anisotropy(ani_fs_weak,base_tensor)
    sigma0 = q/p
    sigma =  0.4*(ani_cn1+ani_bvn1+1.5*ani_bvs1+ani_fn1+1.5*ani_fs1)  
    ######
    return (ani_cn1,ani_bvn1,ani_bvs1,ani_fn1,ani_fs1,sigma,sigma0,avg_fn,pro_sliding,pro_weak,ani_cn1_strong,ani_bvn1_strong,ani_bvs1_strong,ani_fn1_strong,ani_fs1_strong,ani_cn1_weak,ani_bvn1_weak,ani_bvs1_weak,ani_fn1_weak,ani_fs1_weak)
def calc_anisotropy2():
    #get contact force and branch vectors
    #normal contacts
    vlist_fn = list()
    vlist_fs = list()
    #branch vectors
    vlist_bv = list()
    avg_fn = 0.0
    vlist_fsfn = list()
    c_num_strong = 0.0
    for i in O.interactions:
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <6:#wall
            continue
        if fn.norm()>0.0:
            vlist_fn.append(fn)
            fs = i.phys.shearForce
            vlist_fsfn.append(fs.norm()/fn.norm())
            avg_fn += fn.norm()
            vlist_fs.append(fs)
            p = O.bodies[id1].state.pos - O.bodies[id2].state.pos
            vlist_bv.append(p)
    c_num=len(vlist_fn)
    avg_fn /=c_num
    for f in vlist_fn:
        if f.norm() >= avg_fn:#strong contact
            c_num_strong += 1
    c_num_sliding = sum([1 for eta in vlist_fsfn if eta >= 0.499999999])
    pro_sliding = float(c_num_sliding)/c_num #proportion of sliding contacts
    c_num_weak = c_num - c_num_strong
    pro_weak = float(c_num_weak)/c_num #proportion of weak contacts
    #get stress tensor
    tensor_stress = np.zeros([3,3])
    for c in range(c_num):
        c_f = vlist_fn[c] + vlist_fs[c]
        c_bv = -vlist_bv[c]
        
        v_dya1 = vector_dyadic(c_f,c_bv)
       
        tensor_stress += v_dya1
        
    tensor_stress = tensor_stress/c_num
    p = np.einsum('ii',tensor_stress)/3.0
    q1 = tensor_deviator(tensor_stress)
    q = math.sqrt(1.5*tensor_dot_tensor(q1,q1))   
    #fabric tensor of contact normals
    tensor_cn = np.zeros([3,3])
    tensor_cn_strong = np.zeros([3,3])
    tensor_cn_weak = np.zeros([3,3])
    for v in vlist_fn:
        modv=v.norm()
        if modv!=0:
            v=[i/modv for i in v]
        dya1 = vector_dyadic(v,v)
        tensor_cn += dya1
        if modv >= avg_fn:#strong contact
            tensor_cn_strong += dya1
        else:#weak contact
            tensor_cn_weak += dya1
    
    tensor_cn = tensor_cn/c_num
    tensor_cn_strong /= c_num_strong
    tensor_cn_weak /= c_num_weak
    #deviator of cn
    ani_cn = tensor_deviator(tensor_cn)*15.0/2.0
    ani_cn_strong = tensor_deviator(tensor_cn_strong)*15.0/2.0
    ani_cn_weak = tensor_deviator(tensor_cn_weak)*15.0/2.0
    ################
    #anisotropy of contact force
    #fabric tensor of normal contact forces
    tensor_fn = np.zeros([3,3])
    tensor_fs = np.zeros([3,3])
    tensor_fn_strong = np.zeros([3,3])
    tensor_fs_strong = np.zeros([3,3])
    tensor_fn_weak = np.zeros([3,3])
    tensor_fs_weak = np.zeros([3,3])
    for c in range(c_num):
        c_fn = vlist_fn[c]
        c_fs = vlist_fs[c]
        fn1 = c_fn.norm()
        fs1 = c_fs.norm()
        
            
        c_fn=c_fn/fn1
        c_fs=c_fs/fs1
        v_dya1 = vector_dyadic(c_fn,c_fn)
        v_dya2 = vector_dyadic(c_fs,c_fn)
        
        tensor_fn += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        tensor_fs += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        if fn1 >= avg_fn:#strong contact
            tensor_fn_strong += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            tensor_fs_strong += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        else:
            tensor_fn_weak += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            tensor_fs_weak += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
    tensor_fn = tensor_fn/c_num
    tensor_fs = tensor_fs/c_num
    tensor_fn_strong /= c_num_strong
    tensor_fs_strong /= c_num_strong
    tensor_fn_weak /= c_num_weak
    tensor_fs_weak /= c_num_weak
    #deviator of fn
    f0_bar = np.einsum('ii',tensor_fn)
    ani_fn = tensor_deviator(tensor_fn)*15.0/2.0/f0_bar
    ani_fs = tensor_deviator(tensor_fs)*15.0/3.0/f0_bar
    f0_bar_strong = np.einsum('ii',tensor_fn_strong)
    ani_fn_strong = tensor_deviator(tensor_fn_strong)*15.0/2.0/f0_bar_strong
    ani_fs_strong = tensor_deviator(tensor_fs_strong)*15.0/3.0/f0_bar_strong
    f0_bar_weak = np.einsum('ii',tensor_fn_weak)
    ani_fn_weak = tensor_deviator(tensor_fn_weak)*15.0/2.0/f0_bar_weak
    ani_fs_weak = tensor_deviator(tensor_fs_weak)*15.0/3.0/f0_bar_weak
    #b,c=np.linalg.eig(fabtensor)
    ##################
    ###anisotropy of branch vectors
    #fabric tensor of branch vectors
    tensor_bvn = np.zeros([3,3])
    tensor_bvs = np.zeros([3,3])
    tensor_bvn_strong = np.zeros([3,3])
    tensor_bvs_strong = np.zeros([3,3])
    tensor_bvn_weak = np.zeros([3,3])
    tensor_bvs_weak = np.zeros([3,3])
    for c in range(c_num):
        c_fn = vlist_fn[c]
        c_bv = vlist_bv[c]
        fn1 = c_fn.norm()
        c_bv_n = c_bv.dot(c_fn)*c_fn/fn1/fn1###project bv on cn
        c_bv_s = c_bv - c_bv_n
        
        bvn1 = c_bv_n.norm()
        bvs1 = c_bv_s.norm()
        c_fn=c_fn/fn1
        c_bv_n=c_bv_n/bvn1
        c_bv_s=c_bv_s/bvs1
        v_dya1 = vector_dyadic(c_fn,c_fn)
        v_dya2 = vector_dyadic(c_bv_s,c_fn)
        tensor_bvn += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        tensor_bvs += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        if fn1 >= avg_fn:#strong contact
            tensor_bvn_strong += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            tensor_bvs_strong += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
        else:
            tensor_bvn_weak += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            tensor_bvs_weak += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
    tensor_bvn = tensor_bvn/c_num
    tensor_bvs = tensor_bvs/c_num
    tensor_bvn_strong /= c_num_strong
    tensor_bvs_strong /= c_num_strong
    tensor_bvn_weak /= c_num_weak
    tensor_bvs_weak /= c_num_weak
    #deviator of fn
    b0_bar = np.einsum('ii',tensor_bvn)
    ani_bvn = tensor_deviator(tensor_bvn)*15.0/2.0/b0_bar
    ani_bvs = tensor_deviator(tensor_bvs)*15.0/3.0/b0_bar  
    b0_bar_strong = np.einsum('ii',tensor_bvn_strong)
    ani_bvn_strong = tensor_deviator(tensor_bvn_strong)*15.0/2.0/b0_bar_strong
    ani_bvs_strong = tensor_deviator(tensor_bvs_strong)*15.0/3.0/b0_bar_strong
    b0_bar_weak = np.einsum('ii',tensor_bvn_weak)
    ani_bvn_weak = tensor_deviator(tensor_bvn_weak)*15.0/2.0/b0_bar_weak
    ani_bvs_weak = tensor_deviator(tensor_bvs_weak)*15.0/3.0/b0_bar_weak
    #########output the deviatoric invariants
    base_tensor = q1
    ani_cn1 = new_anisotropy(ani_cn,base_tensor)  
    ani_bvn1 = new_anisotropy(ani_bvn,base_tensor) 
    ani_bvs1 = new_anisotropy(ani_bvs,base_tensor) 
    ani_fn1 = new_anisotropy(ani_fn,base_tensor) 
    ani_fs1 = new_anisotropy(ani_fs,base_tensor)
    ani_cn1_strong = new_anisotropy(ani_cn_strong,base_tensor)  
    ani_bvn1_strong = new_anisotropy(ani_bvn_strong,base_tensor) 
    ani_bvs1_strong = new_anisotropy(ani_bvs_strong,base_tensor) 
    ani_fn1_strong = new_anisotropy(ani_fn_strong,base_tensor) 
    ani_fs1_strong = new_anisotropy(ani_fs_strong,base_tensor)
    ani_cn1_weak = new_anisotropy(ani_cn_weak,base_tensor)  
    ani_bvn1_weak = new_anisotropy(ani_bvn_weak,base_tensor) 
    ani_bvs1_weak = new_anisotropy(ani_bvs_weak,base_tensor) 
    ani_fn1_weak = new_anisotropy(ani_fn_weak,base_tensor) 
    ani_fs1_weak = new_anisotropy(ani_fs_weak,base_tensor)
    sigma0 = q/p
    sigma =  0.4*(ani_cn1+ani_bvn1+1.5*ani_bvs1+ani_fn1+1.5*ani_fs1)  
    ######
    return (ani_cn1,ani_bvn1,ani_bvs1,ani_fn1,ani_fs1,sigma,sigma0,avg_fn,pro_sliding,pro_weak,ani_cn1_strong,ani_bvn1_strong,ani_bvs1_strong,ani_fn1_strong,ani_fs1_strong,ani_cn1_weak,ani_bvn1_weak,ani_bvs1_weak,ani_fn1_weak,ani_fs1_weak)

###eta network
def aniso_etanetwork(base_tensor,vlist_fn,vlist_fs,vlist_bv,avg_fn):
    #fabric tensor of contact normals
    tensor_cn = np.zeros([3,3])
    tensor_cn_weak = np.zeros([3,3])
    c_num = len(vlist_fn)
    c_num_weak = 0
    for v in vlist_fn:
        modv=v.norm()
        if modv!=0:
            v=[i/modv for i in v]
        dya1 = vector_dyadic(v,v)
        tensor_cn += dya1
        if modv <= avg_fn:#weak contact, eta network
            tensor_cn_weak += dya1
            c_num_weak += 1
    
    tensor_cn = tensor_cn/c_num
    tensor_cn_weak /= c_num_weak
    tensor_cn = tensor_cn_weak
    #print "tot weak",c_num,c_num_weak
    #deviator of cn
    ani_cn = tensor_deviator(tensor_cn)*15.0/2.0
    ani_cn_weak = tensor_deviator(tensor_cn_weak)*15.0/2.0
    ################
    #anisotropy of contact force
    #fabric tensor of normal contact forces
    tensor_fn_weak = np.zeros([3,3])
    tensor_fs_weak = np.zeros([3,3])
    #fabric tensor of branch vectors
    tensor_bvn_weak = np.zeros([3,3])
    tensor_bvs_weak = np.zeros([3,3])
    c_num_sliding = 0
    for c in range(c_num):
        c_fn0 = vlist_fn[c]
        c_fs0 = vlist_fs[c]
        fn1 = c_fn0.norm()
        fs1 = c_fs0.norm()
        if fn1 <= avg_fn:#weak contact, eta network
            #sliding contacts
            if fs1/fn1 > 0.49999:#sliding contact
                c_num_sliding += 1        
            #####contact force fabric
            c_fn=c_fn0/fn1
            c_fs=c_fs0/fs1
            v_dya1 = vector_dyadic(c_fn,c_fn)
            v_dya2 = vector_dyadic(c_fs,c_fn)

            tensor_fn_weak += fn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            tensor_fs_weak += fs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            #branch vector fabric
            c_bv = vlist_bv[c]
            c_bv_n = c_bv.dot(c_fn0)*c_fn0/fn1/fn1###project bv on cn
            c_bv_s = c_bv - c_bv_n
        
            bvn1 = c_bv_n.norm()
            bvs1 = c_bv_s.norm()
            
            c_bv_n=c_bv_n/bvn1
            c_bv_s=c_bv_s/bvs1
            v_dya1 = vector_dyadic(c_fn,c_fn)
            v_dya2 = vector_dyadic(c_bv_s,c_fn)

            tensor_bvn_weak += bvn1*v_dya1/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))
            tensor_bvs_weak += bvs1*v_dya2/(1.0 + tensor_dot_tensor(ani_cn,v_dya1))

    tensor_fn_weak /= c_num_weak
    tensor_fs_weak /= c_num_weak
    tensor_bvn_weak /= c_num_weak
    tensor_bvs_weak /= c_num_weak
    #deviator of fn
    f0_bar_weak = np.einsum('ii',tensor_fn_weak)
    ani_fn_weak = tensor_deviator(tensor_fn_weak)*15.0/2.0/f0_bar_weak
    ani_fs_weak = tensor_deviator(tensor_fs_weak)*15.0/3.0/f0_bar_weak

    b0_bar_weak = np.einsum('ii',tensor_bvn_weak)
    ani_bvn_weak = tensor_deviator(tensor_bvn_weak)*15.0/2.0/b0_bar_weak
    ani_bvs_weak = tensor_deviator(tensor_bvs_weak)*15.0/3.0/b0_bar_weak
    #########output the deviatoric invariants
    ani_cn1_weak = new_anisotropy(ani_cn_weak,base_tensor)  
    ani_bvn1_weak = new_anisotropy(ani_bvn_weak,base_tensor) 
    ani_bvs1_weak = new_anisotropy(ani_bvs_weak,base_tensor) 
    ani_fn1_weak = new_anisotropy(ani_fn_weak,base_tensor) 
    ani_fs1_weak = new_anisotropy(ani_fs_weak,base_tensor)
    ######
    proportion_sliding = float(c_num_sliding)/c_num
    return     ani_cn1_weak,ani_fn1_weak,ani_fs1_weak,ani_bvn1_weak,ani_bvs1_weak,proportion_sliding
def calc_aniso_etanetwork(outfile,eta=[0.025,0.05,0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0]):
    #get contact force and branch vectors
    #normal contacts
    vlist_fn = list()
    vlist_fs = list()
    #branch vectors
    vlist_bv = list()
    avg_fn = 0.0
    vlist_fsfn = list()
    c_num_strong = 0.0
    for i in O.interactions:
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <6:#wall
            continue
        if fn.norm()>0.0:
            vlist_fn.append(fn)
            fs = i.phys.shearForce
            vlist_fsfn.append(fs.norm()/fn.norm())
            avg_fn += fn.norm()
            vlist_fs.append(fs)
            p = O.bodies[id1].state.pos - O.bodies[id2].state.pos
            vlist_bv.append(p)
    c_num=len(vlist_fn)
    avg_fn /=c_num
    #get stress tensor
    tensor_stress = np.zeros([3,3])
    for c in range(c_num):
        c_f = vlist_fn[c] + vlist_fs[c]
        c_bv = -vlist_bv[c]
        
        v_dya1 = vector_dyadic(c_f,c_bv)
       
        tensor_stress += v_dya1
        
    tensor_stress = tensor_stress/c_num
    p = np.einsum('ii',tensor_stress)/3.0
    q1 = tensor_deviator(tensor_stress)
    q = math.sqrt(1.5*tensor_dot_tensor(q1,q1))   
    #################################
    ###eta network
    f_out = open(outfile,'w')
    print>>f_out,"eta ani_cn,ani_fn,ani_fs,ani_bvn,ani_bvs,ani_tot" ##file header
    for eta_i in eta:
        print eta_i
        fn_threshold = eta_i*avg_fn
        ani_cn,ani_fn,ani_fs,ani_bvn,ani_bvs,p_sliding = aniso_etanetwork(q1,vlist_fn,vlist_fs,vlist_bv,fn_threshold)
        ani_tot = 0.4*(ani_cn+ani_bvn+1.5*ani_bvs+ani_fn+1.5*ani_fs)  
        print>>f_out,eta_i,ani_cn,ani_fn,ani_fs,ani_bvn,ani_bvs,ani_tot,p_sliding
    f_out.close()
    #return

def calc_stressratio_etanetwork(outfile,eta=[0.025,0.05,0.1,0.25,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0]):
    #get contact force and branch vectors
    #normal contacts
    vlist_fn = list()
    vlist_fs = list()
    #branch vectors
    vlist_bv = list()
    avg_fn = 0.0
    avg_f = 0.0
    vlist_fsfn = list()
    c_num_strong = 0.0
    vlist_f = list()
    for i in O.interactions:
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <6:#wall
            continue
        if fn.norm()>0.0:
            vlist_fn.append(fn)
            fs = i.phys.shearForce
            vlist_f.append(fn+fs)
            avg_f += (fn+fs).norm()
            vlist_fsfn.append(fs.norm()/fn.norm())
            avg_fn += fn.norm()
            vlist_fs.append(fs)
            p = O.bodies[id1].state.pos - O.bodies[id2].state.pos
            vlist_bv.append(p)
    c_num=len(vlist_fn)
    avg_fn /=c_num
    avg_f /=c_num
    #################################
    #get volume of the assembly
    triax = O.engines[3]
    for e in O.engines:
        if isinstance(e,CompressionEngine):
            triax = e
            break
    Vol = triax.height*triax.width*triax.depth
    
    #get the total stress tensor
    tensor_stress = np.zeros([3,3])
    for c in range(c_num):
        c_f = vlist_fn[c] + vlist_fs[c]
        c_bv = -vlist_bv[c]
        
        v_dya1 = vector_dyadic(c_f,c_bv)
       
        tensor_stress += v_dya1
        
    tensor_stress = tensor_stress/Vol
    p = np.einsum('ii',tensor_stress)/3.0  

    ###eta network
    f_out = open(outfile,'w')
    print>>f_out,"eta stressratio slidingcontacts" ##file header
    
    for eta_i in eta:
        fn_threshold = eta_i*avg_fn
        #get stress tensor
        tensor_stress = np.zeros([3,3])
        c_num_sliding = 0
        for c in range(c_num):
            c_fn0 = vlist_fn[c]
            c_fs0 = vlist_fs[c]
            
            c_f = vlist_fn[c] + vlist_fs[c]
            fn1 = c_fn0.norm()
            fs1 = c_fs0.norm()
            if fn1 <= fn_threshold:
                #sliding contacts
                if fs1/fn1 > 0.49999:#sliding contact
                    c_num_sliding += 1    
                #stress tensor    
                c_bv = -vlist_bv[c]
                v_dya1 = vector_dyadic(c_f,c_bv)
                tensor_stress += v_dya1
        
        tensor_stress = tensor_stress/Vol
        #p = np.einsum('ii',tensor_stress)/3.0
        q1 = tensor_deviator(tensor_stress)
        q = math.sqrt(1.5*tensor_dot_tensor(q1,q1))   
        prop_sliding = float(c_num_sliding)/c_num
        print eta_i
        print>>f_out,eta_i,q/p,prop_sliding
    f_out.close()
    #return


def calc_sliding_etanetwork(outfile,eta=[0.025,0.05,0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0]):
    #get contact force and branch vectors
    #normal contacts
    vlist_fn = list()
    vlist_fs = list()
    #branch vectors
    vlist_bv = list()
    avg_fn = 0.0
    avg_f = 0.0
    vlist_fsfn = list()
    c_num_strong = 0.0
    vlist_f = list()
    for i in O.interactions:
        fn = i.phys.normalForce
        id1 = i.id1
        id2 = i.id2
        if id1 <6:#wall
            continue
        if fn.norm()>0.0:
            vlist_fn.append(fn)
            fs = i.phys.shearForce
            vlist_f.append(fn+fs)
            avg_f += (fn+fs).norm()
            vlist_fsfn.append(fs.norm()/fn.norm())
            avg_fn += fn.norm()
            vlist_fs.append(fs)
    c_num=len(vlist_fn)
    avg_fn /=c_num
    avg_f /=c_num
    #################################

    ###eta network
    f_out = open(outfile,'w')
    print>>f_out,"eta stressratio slidingcontacts" ##file header
    
    for eta_i in eta:
        fn_threshold = eta_i*avg_fn
        c_num_sliding = 0
        for c in range(c_num):
            c_fn0 = vlist_fn[c]
            c_fs0 = vlist_fs[c]
            
            c_f = vlist_fn[c] + vlist_fs[c]
            fn1 = c_fn0.norm()
            fs1 = c_fs0.norm()
            if fn1 <= fn_threshold:
                #sliding contacts
                if fs1/fn1 > 0.49999:#sliding contact
                    c_num_sliding += 1    
        prop_sliding = float(c_num_sliding)/c_num
        print eta_i
        print>>f_out,eta_i,prop_sliding
    f_out.close()

#friction mobilization, calc average friction mobilization index <I_m>
def calc_avgFricMobilIndex(path,steps):
    filename = 'fricMobiIndex.txt'
    f = open(path+'/'+filename,'w')
    print>>f,'m,zstrain,avg_Im'
    f.close()
    for m in range(steps):
        O.reset()
        print path,m
        #O.load(path+'/final_con.xml.bz2')
        O.load(path+'/shear/shear'+str(m+1)+'.xml.bz2')
        #path = path+'/outdata'
        #find CompressionEngine
        triax = O.engines[3]
        for e in O.engines:
            if isinstance(e,CompressionEngine):
                triax = e
                break
        zstrain = log(triax.height0/triax.height)*100.0
        ####calc the average Im
        c_num = 0.0
        Im_tot = 0.0
        for i in O.interactions:
            fn = i.phys.normalForce
            id1 = i.id1
            id2 = i.id2
            if id1 <6:#wall
                continue
            if fn.norm()>0.0:
                fs = i.phys.shearForce
                Im = fs.norm()/fn.norm()*2.0  #friction coefficient is 0.5
                Im_tot += Im
                c_num += 1
            
        avg_Im = Im_tot/c_num
        ####
        f = open(path+'/'+filename,'a')
        print>>f,m+1,zstrain,avg_Im
    f.close()




def new_out_anisotropy_cn(path,steps):
    filename = 'aniso-latest.txt'
    f = open(path+'/'+filename,'w')
    print>>f,'m,zstrain,ani_cn,ani_bvn,ani_bvs,ani_fn,ani_fs,sigma,sigma0,avg_fn,pro_sliding,pro_weak,ani_cn_s,ani_bvn_s,ani_bvs_s,ani_fn_s,ani_fs_s,ani_cn_w,ani_bvn_w,ani_bvs_w,ani_fn_w,ani_fs_w,cns_avg,rot_avg'
    f.close()
    for m in range(steps):
        O.reset()
        print path,m
        #O.load(path+'/final_con.xml.bz2')
        O.load(path+'/shear/shear'+str(m+1)+'.xml.bz2')
        #path = path+'/outdata'
        #find CompressionEngine
        triax = O.engines[3]
        for e in O.engines:
            if isinstance(e,CompressionEngine):
                triax = e
                break
        zstrain = log(triax.height0/triax.height)*100.0
        #ani_cn,ani_bvn,ani_bvs,ani_fn,ani_fs,sigma,sigma0,avg_fn,pro_sliding,pro_weak,ani_cn_s,ani_bvn_s,ani_bvs_s,ani_fn_s,ani_fs_s,ani_cn_w,ani_bvn_w,ani_bvs_w,ani_fn_w,ani_fs_w = calc_anisotropy()
        #invoke calc_anisotropy2()
        ani_cn,ani_bvn,ani_bvs,ani_fn,ani_fs,sigma,sigma0,avg_fn,pro_sliding,pro_weak,ani_cn_s,ani_bvn_s,ani_bvs_s,ani_fn_s,ani_fs_s,ani_cn_w,ani_bvn_w,ani_bvs_w,ani_fn_w,ani_fs_w = calc_anisotropy2()
        cns_avg,rot_avg = output_CNrotation()
        f = open(path+'/'+filename,'a')
        print>>f,m+1,zstrain,ani_cn,ani_bvn,ani_bvs,ani_fn,ani_fs,sigma,sigma0,avg_fn,pro_sliding,pro_weak,ani_cn_s,ani_bvn_s,ani_bvs_s,ani_fn_s,ani_fs_s,ani_cn_w,ani_bvn_w,ani_bvs_w,ani_fn_w,ani_fs_w,cns_avg,rot_avg
        f.close()

def new_out_anisotropy_eta(m,step):
    file_out = 'Anisotropy_eta'+'/'+m+'-'+str(step)+'aniso.dat'
    O.reset()
    #print path,m
    #O.load(path+'/final_con.xml.bz2')
    O.load(m+'/shear/shear'+str(step)+'.xml.bz2')
    calc_aniso_etanetwork(file_out)

def new_out_stressratio_eta(m,step,eta=[]):
    file_out = 'Anisotropy_eta'+'/'+m+'-'+str(step)+'stressratiopeak.dat'
    O.reset()
    #print path,m
    #O.load(path+'/final_con.xml.bz2')
    O.load(m+'/shear/shear'+str(step)+'.xml.bz2')
    calc_stressratio_etanetwork(file_out,eta=eta)

def new_out_sliding_eta(m,step,eta=[]):
    file_out = 'Anisotropy_eta'+'/'+m+'-'+str(step)+'slidingcontact.dat'
    O.reset()
    #print path,m
    #O.load(path+'/final_con.xml.bz2')
    O.load(m+'/shear/shear'+str(step)+'.xml.bz2')
    calc_sliding_etanetwork(file_out,eta=eta)
####################################################### 
#####################################################
#########
#########output data
#for i in ['1.0','1.2','1.4','1.6']:
#i='1.6'
#steps_num = 128#1.2
#steps_num = 127#1.4
#steps_num = 126#1.6
#steps_num = 130#1.0,1.05
#calc_avgFricMobilIndex(i,steps_num)
#for m in range(129):
#    print i,m
#    save_init(i,m+1)
#new_out_anisotropy_cn(i,steps_num)    
#save_init('1.0',130)
#save_init('1.2',128)
#save_init('1.4',127)
#save_init('1.6',126)
#save_init('1.05',1)
#out_porosity(path='.')


#eta = '1.0'
#numlist = [1,12,130]
#m = '1.05'
#numlist = [1,12,130]
#eta = '1.2'
#numlist = [1,12,128]
#eta = '1.4'
#numlist = [1,14,127]
#eta = '1.6'
#numlist = [1,14,126]
#m = '1.05'
#step = 2
#mm = ['1.0','1.05','1.2','1.4','1.6']
#ss = [130,130,128,127,126]
#ss = [12,12,12,14,13]
#mm=['1.6']
#ss=[13]
#etas = np.arange(0.025,10.0,0.025)
#ss = [2,2,2,2,2]
#for i in range(1):
    #m = mm[i]
    #step = ss[i]
    #new_out_anisotropy_eta(m,step)
    #new_out_stressratio_eta(m,step,eta=etas)
    #new_out_sliding_eta(m,step,eta=etas)
