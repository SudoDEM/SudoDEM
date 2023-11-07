# 2013 Jan Elias, http://www.fce.vutbr.cz/STM/elias.j/, elias.j@fce.vutbr.cz
# https://www.vutbr.cz/www_base/gigadisk.php?i=95194aa9a

"""
Auxiliary functions for polyhedra
"""


import math,random,doctest,geom,numpy
from sudodem import Vector3
from sudodem.wrapper import *
#from miniEigen import *
try: # use psyco if available
    import psyco
    psyco.full()
except ImportError: pass


# c++ implementations for performance reasons
from sudodem._polyhedra_utils import *

#**********************************************************************************
def randomColor(seed=None):
    random.seed(seed);
    #Return random Vector3 with each component in interval 0...1 (uniform distribution)
    return Vector3(random.random(),random.random(),random.random())

#**********************************************************************************
#create polyhedra, one can specify vertices directly, or leave it empty for random shape
def polyhedra(material,size=Vector3(1,1,1),seed=None,v=[],mask=1,fixed=False, color=[-1,-1,-1]):
    """create polyhedra, one can specify vertices directly, or leave it empty for random shape.

    :param Material material: material of new body
    :param Vector3 size: size of new body (see Polyhedra docs)
    :param float seed: seed for random operations
    :param [Vector3] v: list of body vertices (see Polyhedra docs)
    """
    b=Body()
    random.seed(seed);
    b.aspherical = True
    if len(v)>0:
        b.shape = Polyhedra(v=v)
    else:
        b.shape = Polyhedra(size = size, seed=random.randint(0,1E6))
    if color[0] == -1:
        b.shape.color = randomColor(seed=random.randint(0,1E6))
    else:
        b.shape.color = color
    b.mat = material
    b.state.mass = b.mat.density*b.shape.GetVolume()
    b.state.inertia = b.shape.GetInertia()*b.mat.density
    b.state.ori = b.shape.GetOri()
    b.state.pos = b.shape.GetCentroid()
    b.mask=mask
    if fixed:
        b.state.blockedDOFs = 'xyzXYZ'
    return b

#**********************************************************************************
#creates polyhedra having N vertices and resembling disk
def polyhedralBall(radius, N, material, center,mask=1):
    """creates polyhedra having N vertices and resembling disk

    :param float radius: ball radius
    :param int N: number of vertices
    :param Material material: material of new body
    :param Vector3 center: center of the new body
    """
    pts = []

    inc = math.pi * (3. - math.sqrt(5.))
    off = 2. / float(N)
    for k in range(0, N):
        y = k * off - 1. + (off / 2.)
        r = math.sqrt(1. - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r*radius, y*radius, math.sin(phi)*r*radius])

    ball = polyhedra(material,v=pts)
    ball.state.pos = center
    return ball

#**********************************************************************************
def polyhedraTruncIcosaHed(radius, material, centre,mask=1):
    pts = []

    p = (1.+math.sqrt(5.))/2.
    f = radius/math.sqrt(9.*p + 1.)
    A = [[0.,1.,3.*p],[2.,1.+2.*p,p],[1.,2.+p,2.*p]]
    for a in A:
        a = [a[0]*f,a[1]*f,a[2]*f]
        B = [a,[a[1],a[2],a[0]],[a[2],a[0],a[1]]]
        for b in B:
            pts.append(b)
            if not b[0] == 0:
                pts.append([-b[0], b[1], b[2]])
                if not b[1] == 0:
                    pts.append([-b[0],-b[1], b[2]])
                    if not b[2] == 0: pts.append([-b[0],-b[1],-b[2]])
                if not b[2] == 0:
                    pts.append([-b[0], b[1],-b[2]])
            if not b[1] == 0:
                pts.append([ b[0],-b[1], b[2]])
                if not b[2] == 0:
                    pts.append([ b[0],-b[1],-b[2]])
            if not b[2] == 0: pts.append([ b[0], b[1],-b[2]])
    ball = polyhedra(material,v=pts)
    ball.state.pos = centre
    return ball

#**********************************************************************************
def polyhedraSnubCube(radius, material, centre, mask=1):
    pts = []

    f = radius/1.3437133737446
    c1 = 0.337754
    c2 = 1.14261
    c3 = 0.621226
    A = [[c2,c1,c3],[c1,c3,c2],[c3,c2,c1],[-c1,-c2,-c3],[-c2,-c3,-c1],[-c3,-c1,-c2]]
    for a in A:
        a = [a[0]*f,a[1]*f,a[2]*f]
        pts.append([-a[0],-a[1], a[2]])
        pts.append([ a[0],-a[1],-a[2]])
        pts.append([-a[0], a[1],-a[2]])
        pts.append([ a[0], a[1], a[2]])
    ball = polyhedra(material,v=pts)
    ball.state.pos = centre
    return ball
#**********************************************************************************
#fill box [mincoord, maxcoord] by non-overlaping polyhedrons with random geometry and sizes within the range (uniformly distributed)
def fillBox(mincoord, maxcoord,material,sizemin=[1,1,1],sizemax=[1,1,1],ratio=[0,0,0],seed=None,mask=1):
    """fill box [mincoord, maxcoord] by non-overlaping polyhedrons with random geometry and sizes within the range (uniformly distributed)
    :param Vector3 mincoord: first corner
    :param Vector3 maxcoord: second corner
    :param Vector3 sizemin: minimal size of bodies
    :param Vector3 sizemax: maximal size of bodies
    :param Vector3 ratio: scaling ratio
    :param float seed: random seed
    """
    random.seed(seed);
    v = fillBox_cpp(mincoord, maxcoord, sizemin,sizemax,  ratio, random.randint(0,1E6), material)
    #lastnan = -1
    #for i in range(0,len(v)):
    #    if(math.isnan(v[i][0])):
    #        O.bodies.append(polyhedra(material,seed=random.randint(0,1E6),v=v[lastnan+1:i],mask=1,fixed=False))
    #        lastnan = i

#**********************************************************************************
#fill box [mincoord, maxcoord] by non-overlaping polyhedrons with random geometry and sizes within the range (uniformly distributed)
def fillBoxByBalls(mincoord, maxcoord,material,sizemin=[1,1,1],sizemax=[1,1,1],ratio=[0,0,0],seed=None,mask=1,numpoints=60):
    random.seed(seed);
    v = fillBoxByBalls_cpp(mincoord, maxcoord, sizemin,sizemax,  ratio, random.randint(0,1E6), material,numpoints)
    #lastnan = -1
    #for i in range(0,len(v)):
    #    if(math.isnan(v[i][0])):
    #        O.bodies.append(polyhedra(material,seed=random.randint(0,1E6),v=v[lastnan+1:i],mask=1,fixed=False))
    #        lastnan = i


