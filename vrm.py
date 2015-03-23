bl_info = {
    "name": "VRM",
    "description": "Vector Rendering Method script",
    "author": "Antonio Ospite",
    "version": (0, 3),
    "blender": (2, 73, 0),
    "wiki_url": "http://vrm.ao2.it",
    "category": "Render"
}

# ---------------------------------------------------------------------
#    Copyright (c) 2006, 2007, 2008, 2009, 2012 Antonio Ospite
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# ---------------------------------------------------------------------
#
# Additional credits:
#   Thanks to Emilio Aguirre for S2flender from which I took inspirations :)
#   Thanks to Nikola Radovanovic, the author of the original VRM script,
#       the code you read here has been rewritten _almost_ entirely
#       from scratch but Nikola gave me the idea, so I thank him publicly.
#   Thanks to Folkert de Vries for Freestyle SVG Exporter
#
# ---------------------------------------------------------------------
#
# Things TODO for a next release:
#   - Shadeless shader
#   - FIX the issue with negative scales in object tranformations!
#   - Use a better depth sorting algorithm
#   - Review how selections are made (this script uses selection states of
#     primitives to represent visibility infos)
#   - Use a data structure other than Mesh to represent the 2D image?
#     Think to a way to merge (adjacent) polygons that have the same color.
#     Or a way to use paths for silhouettes and contours.
#   - Consider SMIL for animation handling instead of ECMA Script? (Firefox do
#     not support SMIL for animations)
#   - Switch to the Mesh structure, should be considerably faster
#     (partially done, but with Mesh we cannot sort faces, yet)
#   - Implement Edge Styles (silhouettes, contours, etc.) (partially done).
#   - Implement Shading Styles? (partially done, to make more flexible).
#   - Add Vector Writers other than SVG.
#   - set the background color!
#   - Check memory use!!
#
# ---------------------------------------------------------------------

import bpy
from bpy import context
from bpy.props import *
from bpy.types import Operator, AddonPreferences
#import Blender
#from Blender import Scene, Object, Mesh, NMesh, Material, Lamp, Camera, Window
#from Blender.Mathutils import *
import mathutils
from mathutils import *
from math import *
import sys
import time
from functools import cmp_to_key
import bmesh
from freestyle.types import (
        StrokeShader,
        Interface0DIterator,
        Operators,
        )
import parameter_editor

try:
    set()
except NameError:
    from sets import Set as set


def uniq(alist):
    tmpdict = dict()
    return [tmpdict.setdefault(e, e) for e in alist if e not in tmpdict]
    # in python > 2.4 we ca use the following
    #return [ u for u in alist if u not in locals()['_[1]'] ]


# Constants
EPS = 10e-5

# We use a global progress Indicator Object
progress = None


# Config class for global settings

def GetSupportedFormat(scene, context):
  global SWFSupported, PDFSupported
  supportedFormats = [("SVG", "SVG", "")]
  if SWFSupported:
    supportedFormats.append(("SWF", "SWF", ""))
  if PDFSupported:
    supportedFormats.append(("PDF", "PDF", ""))
  return supportedFormats

def update_outputFORMAT(self, context):
  vrm = bpy.context.scene.VRM
  if vrm.outputPATH!="":
    if vrm.outputFORMAT=="SVG" and not vrm.outputPATH.upper().endswith(".SVG"):
      vrm.outputPATH = vrm.outputPATH[:-4]+".svg"
    elif vrm.outputFORMAT=="SWF" and not vrm.outputPATH.upper().endswith(".SWF"):
      vrm.outputPATH = vrm.outputPATH[:-4]+".swf"

class VRM(bpy.types.PropertyGroup):
    bl_idname = "RENDER_PT_VRM"
    
    polygonsSHOW = BoolProperty(default=True,name='Filled Polygons',description='Render filled polygons')
    polygonsSHADING = EnumProperty(items = [('FLAT', 'FLAT', ''), ('TOON', 'TOON', '')],default="FLAT", name='Shading Style',description='Choose the shading style')
    polygonsHSR = EnumProperty(items = [('PAINTER', 'PAINTER', ''), ('NEWELL', 'NEWELL', '')],default="PAINTER")
    # Hidden to the user for now
    polygonsEXPANSION_TRICK = BoolProperty(default=True,)

    polygonsTOON_LEVELS = 2

    freestyleSHOW = BoolProperty(default=False,name='Show FreeStyle\'s Line',description='Render FreeStyle\'s Line')
    edgesSHOW = BoolProperty(default=False,name='Show Edges',description='Render polygon edges')
    edgesSHOW_HIDDEN = BoolProperty(default=False,name='Show Hidden Edges',description='Render hidden edges as dashed lines')
    edgesSTYLE = EnumProperty(items = [('MESH', 'MESH', ''), ('SILHOUETTE', 'SILHOUETTE', '')],default="MESH",name='Style',description='Choose the edge style')
    edgesWIDTH = FloatProperty(default=2,min=0,subtype='PIXEL', name="",description='Change Edge Width')
    edgesCOLOR = IntVectorProperty(subtype="COLOR", min=0, max=255, name="",description='Change Edge Color')

    outputFORMAT = EnumProperty(items=GetSupportedFormat, name='Select the output Format:',description='Choose the Output Format', update=update_outputFORMAT)
    outputANIMATION = BoolProperty(default=False,name='Animation',description='Toggle rendering of animations')
    outputJOIN_OBJECTS = BoolProperty(default=True,name='Join objects',description='Join objects in the rendered file')
    outputPATH = StringProperty(subtype='FILE_PATH',name='path')

class VRMRenderOp(bpy.types.Operator):
    bl_idname = "vrm.render_op"
    bl_label = "Render"

    def execute(self, context):
      vectorize(bpy.path.abspath(bpy.context.scene.VRM.outputPATH))
      return {'FINISHED'}

class VRMPanel(bpy.types.Panel):
    #bl_idname = __name__
    bl_idname = "RENDER_PT_VRMPanel"
    bl_space_type = 'PROPERTIES'
    bl_label = "Vector Rendering Method"
    bl_region_type = 'WINDOW'
    bl_context = "render"
    
    def draw(self, context):
        layout = self.layout
        vrm = context.scene.VRM
        
        box = layout.box()
        box.label("Select the output format")
        box.prop(vrm, "outputPATH")
        box.prop(vrm, "outputFORMAT")
        box.prop(vrm, "outputANIMATION")
        box.prop(vrm, "outputJOIN_OBJECTS")
        
        layout.operator(VRMRenderOp.bl_idname)
        
        box = layout.box()
        box.label("Rendering Style")
        box.prop(vrm, "polygonsSHOW")
        box.prop(vrm, "polygonsSHADING")
        #layout.prop(self, "polygonsHSR")
        #layout.prop(self, "polygonsEXPANSION_TRICK")
        #layout.prop(self, "polygonsTOON_LEVELS")
        box.prop(vrm, "freestyleSHOW")
        #box.prop(vrm, "edgesSHOW")
        #box.prop(vrm, "edgesSTYLE")
        #row = box.row()
        #row.prop(vrm, "edgesWIDTH")
        #row.prop(vrm, "edgesCOLOR")
        #box.prop(vrm, "edgesSHOW_HIDDEN")


# Utility functions
print_debug = False


def dumpfaces(flist, filename):
    """Dump a single face to a file.
    """
    if not print_debug:
        return

    class tmpmesh:
        pass

    m = tmpmesh()
    m.faces = flist

    writerobj = SVGVectorWriter(filename)

    writerobj.open()
    writerobj._printPolygons(m)

    writerobj.close()


def debug(msg):
    if print_debug:
        sys.stderr.write(msg)


def EQ(v1, v2):
    return (abs(v1[0] - v2[0]) < EPS and
            abs(v1[1] - v2[1]) < EPS)
by_furthest_z = (lambda f1, f2:
    cmp(max([v.co[2] for v in f1]), max([v.co[2] for v in f2]) + EPS)
    )


def sign(x):

    if x < -EPS:
    #if x < 0:
        return -1
    elif x > EPS:
    #elif x > 0:
        return 1
    else:
        return 0


# ---------------------------------------------------------------------
#
## HSR Utility class
#
# ---------------------------------------------------------------------

EPS = 10e-5
INF = 10e5


class HSR:
    """A utility class for HSR processing.
    """

    def is_nonplanar_quad(face):
        """Determine if a quad is non-planar.

        From: http://mathworld.wolfram.com/Coplanar.html

        Geometric objects lying in a common plane are said to be coplanar.
        Three noncollinear points determine a plane and so are trivially
        coplanar. Four points are coplanar iff the volume of the tetrahedron
        defined by them is 0,

            | x_1 y_1 z_1 1 |
            | x_2 y_2 z_2 1 |
            | x_3 y_3 z_3 1 |
            | x_4 y_4 z_4 1 | == 0

        Coplanarity is equivalent to the statement that the pair of lines
        determined by the four points are not skew, and can be equivalently
        stated in vector form as (x_3-x_1).[(x_2-x_1)x(x_4-x_3)]==0.

        An arbitrary number of n points x_1, ..., x_n can be tested for
        coplanarity by finding the point-plane distances of the points
        x_4, ..., x_n from the plane determined by (x_1,x_2,x_3)
        and checking if they are all zero.
        If so, the points are all coplanar.

        We here check only for 4-point complanarity.
        """
        n = len(face)

        # assert(n>4)
        if n < 3 or n > 4:
            print("ERROR a mesh in Blender can't have more than 4 vertices or less than 3")
            raise AssertionError

        elif n == 3:
            # three points must be complanar
            return False
        else:  # n == 4
            x1 = Vector(face[0].co)
            x2 = Vector(face[1].co)
            x3 = Vector(face[2].co)
            x4 = Vector(face[3].co)

            v = (x3 - x1) * CrossVecs((x2 - x1), (x4 - x3))
            if v != 0:
                return True

        return False

    is_nonplanar_quad = staticmethod(is_nonplanar_quad)

    def pointInPolygon(poly, v):
        return False

    pointInPolygon = staticmethod(pointInPolygon)

    def edgeIntersection(s1, s2, do_perturbate=False):

        (x1, y1) = s1[0].co[0], s1[0].co[1]
        (x2, y2) = s1[1].co[0], s1[1].co[1]

        (x3, y3) = s2[0].co[0], s2[0].co[1]
        (x4, y4) = s2[1].co[0], s2[1].co[1]

        #z1 = s1[0].co[2]
        #z2 = s1[1].co[2]
        #z3 = s2[0].co[2]
        #z4 = s2[1].co[2]

        # calculate delta values (vector components)
        dx1 = x2 - x1
        dx2 = x4 - x3
        dy1 = y2 - y1
        dy2 = y4 - y3

        #dz1 = z2 - z1
        #dz2 = z4 - z3

        C = dy2 * dx1 - dx2 * dy1  # cross product
        if C == 0:  # parallel
            return None

        dx3 = x1 - x3  # combined origin offset vector
        dy3 = y1 - y3

        a1 = (dy3 * dx2 - dx3 * dy2) / C
        a2 = (dy3 * dx1 - dx3 * dy1) / C

        # check for degeneracies
        #print_debug("\n")
        #print_debug(str(a1)+"\n")
        #print_debug(str(a2)+"\n\n")

        if (a1 == 0 or a1 == 1 or a2 == 0 or a2 == 1):
            # Intersection on boundaries, we consider the point external?
            return None

        elif (a1 > 0.0 and a1 < 1.0 and a2 > 0.0 and a2 < 1.0):  # lines cross
            x = x1 + a1 * dx1
            y = y1 + a1 * dy1

            #z = z1 + a1 * dz1
            z = 0
            return (NMesh.Vert(x, y, z), a1, a2)

        else:
            # lines have intersections but not those segments
            return None

    edgeIntersection = staticmethod(edgeIntersection)

    def isVertInside(self, v):
        winding_number = 0
        coincidence = False

        # Create point at infinity
        point_at_infinity = NMesh.Vert(-INF, v.co[1], -INF)

        for i in range(len(self.v)):
            s1 = (point_at_infinity, v)
            s2 = (self.v[i - 1], self.v[i])

            if EQ(v.co, s2[0].co) or EQ(v.co, s2[1].co):
                coincidence = True

            if HSR.edgeIntersection(s1, s2, do_perturbate=False):
                winding_number += 1

        # Check even or odd
        if (winding_number % 2) == 0:
            return False
        else:
            if coincidence:
                return False
            return True

    isVertInside = staticmethod(isVertInside)

    def det(a, b, c):
        return ((b[0] - a[0]) * (c[1] - a[1]) -
                (b[1] - a[1]) * (c[0] - a[0]))

    det = staticmethod(det)

    def pointInPolygon(q, P):
        is_in = False

        point_at_infinity = NMesh.Vert(-INF, q.co[1], -INF)

        det = HSR.det

        for i in range(len(P.v)):
            p0 = P.v[i - 1]
            p1 = P.v[i]
            if (det(q.co, point_at_infinity.co, p0.co) < 0) != (det(q.co, point_at_infinity.co, p1.co) < 0):
                if det(p0.co, p1.co, q.co) == 0:
                    #print "On Boundary"
                    return False
                elif (det(p0.co, p1.co, q.co) < 0) != (det(p0.co, p1.co, point_at_infinity.co) < 0):
                    is_in = not is_in

        return is_in

    pointInPolygon = staticmethod(pointInPolygon)

    def projectionsOverlap(f1, f2):
        """ If you have nonconvex, but still simple polygons, an acceptable method
        is to iterate over all vertices and perform the Point-in-polygon test[1].
        The advantage of this method is that you can compute the exact
        intersection point and collision normal that you will need to simulate
        collision. When you have the point that lies inside the other polygon, you
        just iterate over all edges of the second polygon again and look for edge
        intersections. Note that this method detects collsion when it already
        happens. This algorithm is fast enough to perform it hundreds of times per
        sec.  """

        for i in range(len(f1.v)):

            # If a point of f1 in inside f2, there is an overlap!
            v1 = f1.v[i]
            #if HSR.isVertInside(f2, v1):
            if HSR.pointInPolygon(v1, f2):
                return True

            # If not the polygon can be ovelap as well, so we check for
            # intersection between an edge of f1 and all the edges of f2

            v0 = f1.v[i - 1]

            for j in range(len(f2.v)):
                v2 = f2.v[j - 1]
                v3 = f2.v[j]

                e1 = v0, v1
                e2 = v2, v3

                intrs = HSR.edgeIntersection(e1, e2)
                if intrs:
                    #print_debug(str(v0.co) + " " + str(v1.co) + " " +
                    #        str(v2.co) + " " + str(v3.co) )
                    #print_debug("\nIntersection\n")

                    return True

        return False

    projectionsOverlap = staticmethod(projectionsOverlap)

    def midpoint(p1, p2):
        """Return the midpoint of two vertices.
        """
        m = MidpointVecs(Vector(p1), Vector(p2))
        mv = NMesh.Vert(m[0], m[1], m[2])

        return mv

    midpoint = staticmethod(midpoint)

    def facesplit(P, Q, facelist, nmesh):
        """Split P or Q according to the strategy illustrated in the Newell's
        paper.
        """

        by_furthest_z = (lambda f1, f2:
                cmp(max([v.co[2] for v in f1]), max([v.co[2] for v in f2]) + EPS)
                )

        # Choose if split P on Q plane or vice-versa

        n = 0
        for Pi in P:
            d = HSR.Distance(Vector(Pi), Q)
            if d <= EPS:
                n += 1
        pIntersectQ = (n != len(P))

        n = 0
        for Qi in Q:
            d = HSR.Distance(Vector(Qi), P)
            if d >= -EPS:
                n += 1
        qIntersectP = (n != len(Q))

        newfaces = []

        # 1. If parts of P lie in both half-spaces of Q
        # then splice P in two with the plane of Q
        if pIntersectQ:
            #print "We split P"
            f = P
            plane = Q

            newfaces = HSR.splitOn(plane, f)

        # 2. Else if parts of Q lie in both half-space of P
        # then splice Q in two with the plane of P
        if qIntersectP and newfaces == None:
            #print "We split Q"
            f = Q
            plane = P

            newfaces = HSR.splitOn(plane, f)
            #print "After"

        # 3. Else slice P in half through the mid-point of
        # the longest pair of opposite sides
        if newfaces == None:

            print("We ignore P...")
            facelist.remove(P)
            return facelist

            #f = P

            #if len(P)==3:
            #    v1 = midpoint(f[0], f[1])
            #    v2 = midpoint(f[1], f[2])
            #if len(P)==4:
            #    v1 = midpoint(f[0], f[1])
            #    v2 = midpoint(f[2], f[3])
            #vec3 = (Vector(v2)+10*Vector(f.normal))
            #
            #v3 = NMesh.Vert(vec3[0], vec3[1], vec3[2])

            #plane = NMesh.Face([v1, v2, v3])
            #
            #newfaces = splitOn(plane, f)

        if newfaces == None:
            print("Big FAT problem, we weren't able to split POLYGONS!")
            raise AssertionError

        #print newfaces
        if newfaces:
            #for v in f:
            #    if v not in plane and v in nmesh.verts:
            #        nmesh.verts.remove(v)
            for nf in newfaces:

                nf.mat = f.mat
                nf.sel = f.sel
                nf.col = [f.col[0]] * len(nf.v)

                nf.smooth = 0

                for v in nf:
                    nmesh.verts.append(v)
                # insert pieces in the list
                facelist.append(nf)

            facelist.remove(f)

        # and resort the faces
        facelist.sort(by_furthest_z)
        facelist.sort(lambda f1, f2: cmp(f1.smooth, f2.smooth))
        facelist.reverse()

        #print [ f.smooth for f in facelist ]

        return facelist

    facesplit = staticmethod(facesplit)

    def isOnSegment(v1, v2, p, extremes_internal=False):
        """Check if point p is in segment v1v2.
        """

        l1 = (v1 - p).length
        l2 = (v2 - p).length

        # Should we consider extreme points as internal ?
        # The test:
        # if p == v1 or p == v2:
        if l1 < EPS or l2 < EPS:
            return extremes_internal

        l = (v1 - v2).length

        # if the sum of l1 and l2 is circa l, then the point is on segment,
        if abs(l - (l1 + l2)) < EPS:
            return True
        else:
            return False

    isOnSegment = staticmethod(isOnSegment)

    def Distance(point, face):
        """ Calculate the distance between a point and a face.

        An alternative but more expensive method can be:

            ip = Intersect(Vector(face[0]), Vector(face[1]), Vector(face[2]),
                    Vector(face.no), Vector(point), 0)

            d = Vector(ip - point).length

        See: http://mathworld.wolfram.com/Point-PlaneDistance.html
        """

        p = Vector(point)
        plNormal = Vector(face.no)
        plVert0 = Vector(face.v[0])

        d = (plVert0 * plNormal) - (p * plNormal)

        #d = plNormal * (plVert0 - p)

        #print "\nd: %.10f - sel: %d, %s\n" % (d, face.sel, str(point))

        return d

    Distance = staticmethod(Distance)

    def makeFaces(vl):
        #
        # make one or two new faces based on a list of vertex-indices
        #
        newfaces = []

        if len(vl) <= 4:
            nf = NMesh.Face()

            for v in vl:
                nf.v.append(v)

            newfaces.append(nf)

        else:
            nf = NMesh.Face()

            nf.v.append(vl[0])
            nf.v.append(vl[1])
            nf.v.append(vl[2])
            nf.v.append(vl[3])
            newfaces.append(nf)

            nf = NMesh.Face()
            nf.v.append(vl[3])
            nf.v.append(vl[4])
            nf.v.append(vl[0])
            newfaces.append(nf)

        return newfaces

    makeFaces = staticmethod(makeFaces)

    def splitOn(Q, P, return_positive_faces=True, return_negative_faces=True):
        """Split P using the plane of Q.
        Logic taken from the knife.py python script
        """

        # Check if P and Q are parallel
        u = CrossVecs(Vector(Q.no), Vector(P.no))
        ax = abs(u[0])
        ay = abs(u[1])
        az = abs(u[2])

        if (ax + ay + az) < EPS:
            print("PARALLEL planes!!")
            return

        # The final aim is to find the intersection line between P
        # and the plane of Q, and split P along this line

        nP = len(P.v)

        # Calculate point-plane Distance between vertices of P and plane Q
        d = []
        for i in range(0, nP):
            d.append(HSR.Distance(P.v[i], Q))

        newVertList = []

        posVertList = []
        negVertList = []
        for i in range(nP):
            d0 = d[i - 1]
            V0 = P.v[i - 1]

            d1 = d[i]
            V1 = P.v[i]

            #print "d0:", d0, "d1:", d1

            # if the vertex lies in the cutplane
            if abs(d1) < EPS:
                #print "d1 On cutplane"
                posVertList.append(V1)
                negVertList.append(V1)
            else:
                # if the previous vertex lies in cutplane
                if abs(d0) < EPS:
                    #print "d0 on Cutplane"
                    if d1 > 0:
                        #print "d1 on positive Halfspace"
                        posVertList.append(V1)
                    else:
                        #print "d1 on negative Halfspace"
                        negVertList.append(V1)
                else:
                    # if they are on the same side of the plane
                    if (d1 * d0) > 0:
                        #print "On the same half-space"
                        if d1 > 0:
                            #print "d1 on positive Halfspace"
                            posVertList.append(V1)
                        else:
                            #print "d1 on negative Halfspace"
                            negVertList.append(V1)

                    # the vertices are not on the same side of the plane, so we have an intersection
                    else:
                        #print "Intersection"

                        e = Vector(V0), Vector(V1)
                        tri = Vector(Q[0]), Vector(Q[1]), Vector(Q[2])

                        inters = Intersect(tri[0], tri[1], tri[2], e[1] - e[0], e[0], 0)
                        if inters == None:
                            print("Split Break")
                            break

                        #print "Intersection", inters

                        nv = NMesh.Vert(inters[0], inters[1], inters[2])
                        newVertList.append(nv)

                        posVertList.append(nv)
                        negVertList.append(nv)

                        if d1 > 0:
                            posVertList.append(V1)
                        else:
                            negVertList.append(V1)

        # uniq for python > 2.4
        #posVertList = [ u for u in posVertList if u not in locals()['_[1]'] ]
        #negVertList = [ u for u in negVertList if u not in locals()['_[1]'] ]

        # a more portable way
        posVertList = uniq(posVertList)
        negVertList = uniq(negVertList)

        # If vertex are all on the same half-space, return
        #if len(posVertList) < 3:
        #    print "Problem, we created a face with less that 3 vertices??"
        #    posVertList = []
        #if len(negVertList) < 3:
        #    print "Problem, we created a face with less that 3 vertices??"
        #    negVertList = []

        if len(posVertList) < 3 or len(negVertList) < 3:
            #print "RETURN NONE, SURE???"
            return None

        if not return_positive_faces:
            posVertList = []
        if not return_negative_faces:
            negVertList = []

        newfaces = HSR.addNewFaces(posVertList, negVertList)

        return newfaces

    splitOn = staticmethod(splitOn)

    def addNewFaces(posVertList, negVertList):
        # Create new faces resulting from the split
        outfaces = []
        if len(posVertList) or len(negVertList):

            #newfaces = [posVertList] + [negVertList]
            newfaces = ([[NMesh.Vert(v[0], v[1], v[2]) for v in posVertList]] +
                    [[NMesh.Vert(v[0], v[1], v[2]) for v in negVertList]])

            for nf in newfaces:
                if nf and len(nf) > 2:
                    outfaces += HSR.makeFaces(nf)

        return outfaces

    addNewFaces = staticmethod(addNewFaces)


# ---------------------------------------------------------------------
#
## Mesh Utility class
#
# ---------------------------------------------------------------------

class MeshUtils:

    def buildEdgeFaceUsersCache(me):
        '''
        Takes a mesh and returns a list aligned with the meshes edges.
        Each item is a list of the faces that use the edge
        would be the equiv for having ed.face_users as a property

        Taken from .blender/scripts/bpymodules/BPyMesh.py,
        thanks to ideasman_42.
        '''

        def sorted_edge_indicies(ed):
            i1 = ed.vertices[0]
            i2 = ed.vertices[1]
            if i1 > i2:
                i1, i2 = i2, i1
            return i1, i2

        me.update(calc_edges=False, calc_tessface=True)
        face_edges_dict = dict([(sorted_edge_indicies(ed), (ed.index, [])) for ed in me.edges])
        for f in me.polygons:#.tessfaces:
            fvi = f.vertices  # face vert idx's
            for i in range(len(f.vertices)):
                i1 = fvi[i]
                i2 = fvi[i - 1]

                if i1 > i2:
                    i1, i2 = i2, i1

                face_edges_dict[i1, i2][1].append(f)

        face_edges = [None] * len(me.edges)
        for ed_index, ed_faces in face_edges_dict.values():
            face_edges[ed_index] = ed_faces

        return face_edges

    def isMeshEdge(adjacent_faces):
        """Mesh edge rule.

        A mesh edge is visible if _at_least_one_ of its adjacent faces is selected.
        Note: if the edge has no adjacent faces we want to show it as well,
        useful for "edge only" portion of objects.
        """

        if len(adjacent_faces) == 0:
            return True

        selected_faces = [f for f in adjacent_faces if f.select]

        if len(selected_faces) != 0:
            return True
        else:
            return False

    def isSilhouetteEdge(adjacent_faces):
        """Silhuette selection rule.

        An edge is a silhuette edge if it is shared by two faces with
        different selection status or if it is a boundary edge of a selected
        face.
        """

        if ((len(adjacent_faces) == 1 and adjacent_faces[0].select == True) or
            (len(adjacent_faces) == 2 and
                adjacent_faces[0].select != adjacent_faces[1].select)
            ):
            return True
        else:
            return False

    buildEdgeFaceUsersCache = staticmethod(buildEdgeFaceUsersCache)
    isMeshEdge = staticmethod(isMeshEdge)
    isSilhouetteEdge = staticmethod(isSilhouetteEdge)


# ---------------------------------------------------------------------
#
## Shading Utility class
#
# ---------------------------------------------------------------------

class ShadingUtils:

    shademap = None

    @staticmethod
    def toonShadingMapSetup():
        vrm = bpy.context.scene.VRM
        levels = vrm.polygonsTOON_LEVELS

        texels = 2 * levels - 1
        ShadingUtils.shademap = [0.0] + [(i) / float(texels - 1) for i in range(1, texels - 1)] + [1.0]
        print("toonMap: "+str(ShadingUtils.shademap))

        return

    @staticmethod
    def toonShading(u):

        if not ShadingUtils.shademap:
            ShadingUtils.toonShadingMapSetup()
        
        shademap = ShadingUtils.shademap

        v = 1.0
        for i in range(0, len(shademap) - 1):
            pivot = (shademap[i] + shademap[i + 1]) / 2.0
            j = int(u > pivot)

            v = shademap[i + j]

            if v < shademap[i + 1]:
                return v

        return v

    #toonShadingMapSetup = staticmethod(toonShadingMapSetup)
    #toonShading = staticmethod(toonShading)


# ---------------------------------------------------------------------
#
## Projections classes
#
# ---------------------------------------------------------------------

class Projector:
    """Calculate the projection of an object given the camera.

    A projector is useful to so some per-object transformation to obtain the
    projection of an object given the camera.

    The main method is #doProjection# see the method description for the
    parameter list.
    """

    def __init__(self, cameraObj, canvasRatio):
        """Calculate the projection matrix.

        The projection matrix depends, in this case, on the camera settings.
        TAKE CARE: This projector expects vertices in World Coordinates!
        """

        camera = cameraObj.data

        aspect = float(canvasRatio[0]) / float(canvasRatio[1])
        near = camera.clip_start
        far = camera.clip_end

        scale = float(camera.draw_size)

        fovy = atan(0.5 / aspect / (camera.lens / 32))
        fovy = fovy * 360.0 / pi

        camPersp = 'PERSP'
        camOrtho = 'ORTHO'

        # What projection do we want?
        if camera.type == camPersp:
            mP = self._calcPerspectiveMatrix(fovy, aspect, near, far)
        elif camera.type == camOrtho:
            mP = self._calcOrthoMatrix(fovy, aspect, near, far, scale)

        # View transformation
        #cam = Matrix(camera.getInverseMatrix())
        #cam.transpose()
        cam = mathutils.Matrix(cameraObj.matrix_world.copy()) #http://blenderartists.org/forum/showthread.php?199557-getInverseMatrix()-equivalence-in-2-5
        cam.invert()
        cam.transpose()
        cam.transpose()

        mP = mP * cam

        self.projectionMatrix = mP

    ##
    # Public methods
    #

    def doProjection(self, v):
        """Project the point on the view plane.

        Given a vertex calculate the projection using the current projection
        matrix.
        """

        # Note that we have to work on the vertex using homogeneous coordinates
        # From blender 2.42+ we don't need to resize the vector to be 4d
        # when applying a 4x4 matrix, but we do that anyway since we need the
        # 4th coordinate later
        p = self.projectionMatrix * Vector(v).to_4d()

        # Perspective division
        if p[3] != 0:
            p[0] = p[0] / p[3]
            p[1] = p[1] / p[3]
            p[2] = p[2] / p[3]

        # restore the size
        p[3] = 1.0
        p.resize_3d()

        return p

    ##
    # Private methods
    #

    def _calcPerspectiveMatrix(self, fovy, aspect, near, far):
        """Return a perspective projection matrix.
        """

        top = near * tan(fovy * pi / 360.0)
        bottom = -top
        left = bottom * aspect
        right = top * aspect
        x = (2.0 * near) / (right - left)
        y = (2.0 * near) / (top - bottom)
        a = (right + left) / (right - left)
        b = (top + bottom) / (top - bottom)
        c = - ((far + near) / (far - near))
        d = - ((2 * far * near) / (far - near))

        m = mathutils.Matrix((
                (x,   0.0,    a,    0.0),
                (0.0,   y,    b,    0.0),
                (0.0, 0.0,    c,      d),
                (0.0, 0.0, -1.0,    0.0)))

        return m

    def _calcOrthoMatrix(self, fovy, aspect, near, far, scale):
        """Return an orthogonal projection matrix.
        """

        # The 11 in the formula was found emiprically
        top = near * tan(fovy * pi / 360.0) * (scale * 11)
        bottom = -top
        left = bottom * aspect
        right = top * aspect
        rl = right - left
        tb = top - bottom
        fn = near - far
        tx = -((right + left) / rl)
        ty = -((top + bottom) / tb)
        tz = ((far + near) / fn)

        m = Matrix(
                [2.0 / rl, 0.0,      0.0,       tx],
                [0.0,      2.0 / tb, 0.0,       ty],
                [0.0,      0.0,      2.0 / fn,  tz],
                [0.0,      0.0,      0.0,      1.0])

        return m


# ---------------------------------------------------------------------
#
## Progress Indicator
#
# ---------------------------------------------------------------------

class Progress:
    """A model for a progress indicator.

    Do the progress calculation calculation and
    the view independent stuff of a progress indicator.
    """
    def __init__(self, steps=0):
        self.name = ""
        self.steps = steps
        self.completed = 0
        self.progress = 0

    def setSteps(self, steps):
        """Set the number of steps of the activity wich we want to track.
        """
        self.steps = steps

    def getSteps(self):
        return self.steps

    def setName(self, name):
        """Set the name of the activity wich we want to track.
        """
        self.name = name

    def getName(self):
        return self.name

    def getProgress(self):
        return self.progress

    def reset(self):
        self.completed = 0
        self.progress = 0

    def update(self):
        """Update the model, call this method when one step is completed.
        """
        if self.progress == 100:
            return False

        self.completed += 1
        self.progress = (float(self.completed) / float(self.steps)) * 100
        self.progress = int(self.progress)

        return True


class ProgressIndicator:
    """An abstraction of a View for the Progress Model
    """
    def __init__(self):

        # Use a refresh rate so we do not show the progress at
        # every update, but every 'self.refresh_rate' times.
        self.refresh_rate = 10
        self.shows_counter = 0

        self.quiet = False

        self.progressModel = None

    def setQuiet(self, value):
        self.quiet = value

    def setActivity(self, name, steps):
        """Initialize the Model.

        In a future version (with subactivities-progress support) this method
        could only set the current activity.
        """
        self.progressModel = Progress()
        self.progressModel.setName(name)
        self.progressModel.setSteps(steps)

    def getActivity(self):
        return self.progressModel

    def update(self):
        """Update the model and show the actual progress.
        """
        assert(self.progressModel)

        if self.progressModel.update():
            if self.quiet:
                return

            self.show(self.progressModel.getProgress(),
                    self.progressModel.getName())

        # We return always True here so we can call the update() method also
        # from lambda funcs (putting the call in logical AND with other ops)
        return True

    def show(self, progress, name=""):
        self.shows_counter = (self.shows_counter + 1) % self.refresh_rate
        if self.shows_counter != 0:
            return

        if progress == 100:
            self.shows_counter = -1


class ConsoleProgressIndicator(ProgressIndicator):
    """Show a progress bar on stderr, a la wget.
    """
    def __init__(self):
        ProgressIndicator.__init__(self)

        self.swirl_chars = ["-", "\\", "|", "/"]
        self.swirl_count = -1

    def show(self, progress, name):
        ProgressIndicator.show(self, progress, name)

        bar_length = 70
        bar_progress = int((progress / 100.0) * bar_length)
        bar = ("=" * bar_progress).ljust(bar_length)

        self.swirl_count = (self.swirl_count + 1) % len(self.swirl_chars)
        swirl_char = self.swirl_chars[self.swirl_count]

        progress_bar = "%s |%s| %c %3d%%" % (name, bar, swirl_char, progress)

        sys.stderr.write(progress_bar + "\r")
        if progress == 100:
            sys.stderr.write("\n")


class GraphicalProgressIndicator(ProgressIndicator):
    """Interface to the Blender.Window.DrawProgressBar() method.
    """
    def __init__(self):
        ProgressIndicator.__init__(self)

        #self.swirl_chars = ["-", "\\", "|", "/"]
        # We have to use letters with the same width, for now!
        # Blender progress bar considers the font widths when
        # calculating the progress bar width.
        self.swirl_chars = ["\\", "/"]
        self.swirl_count = -1

    def show(self, progress, name):
        ProgressIndicator.show(self, progress)

        self.swirl_count = (self.swirl_count + 1) % len(self.swirl_chars)
        swirl_char = self.swirl_chars[self.swirl_count]

        progress_text = "%s - %c %3d%%" % (name, swirl_char, progress)

        # Finally draw  the Progress Bar
        Window.WaitCursor(1)  # Maybe we can move that call in the constructor?
        Window.DrawProgressBar(progress / 100.0, progress_text)

        if progress == 100:
            Window.DrawProgressBar(1, progress_text)
            Window.WaitCursor(0)

# ---------------------------------------------------------------------
#
## 2D Object representation class
#
# ---------------------------------------------------------------------

# TODO: a class to represent the needed properties of a 2D vector image
# For now just using a [N]Mesh structure.

# ---------------------------------------------------------------------
#
## Vector Drawing Classes
#
# ---------------------------------------------------------------------

## A generic Writer

class VectorWriter:
    """
    A class for printing output in a vectorial format.

    Given a 2D representation of the 3D scene the class is responsible to
    write it is a vector format.

    Every subclasses of VectorWriter must have at last the following public
    methods:
        - open(self)
        - close(self)
        - printCanvas(self, scene,
            doPrintPolygons=True, doPrintEdges=False, showHiddenEdges=False):
    """

    def __init__(self, fileName):
        """Set the output file name and other properties"""

        try:
            global writerSETTING
            writerSETTING
        except:
            writerSETTING = True

        self.outputFileName = fileName

        render = bpy.context.scene.render
        self.canvasSize = (int(render.resolution_x * render.resolution_percentage / 100), int(render.resolution_y * render.resolution_percentage / 100))

        self.fps = render.fps

        self.startFrame = 1
        self.endFrame = 1
        self.animation = False

    ##
    # Public Methods
    #

    def open(self, startFrame=1, endFrame=1):
        if startFrame != endFrame:
            self.startFrame = startFrame
            self.endFrame = endFrame
            self.animation = True

        print("Outputting to: ", self.outputFileName)

        return

    def close(self):
        return

    def printCanvas(self, scene, doPrintPolygons=True, doPrintEdges=False,
            showHiddenEdges=False):
        """This is the interface for the needed printing routine.
        """
        return


## SVG Writer

class SVGVectorWriter(VectorWriter):
    """A concrete class for writing SVG output.
    """

    def __init__(self, fileName):
        """Simply call the parent Contructor.
        """
        VectorWriter.__init__(self, fileName)

        self.file = None

    ##
    # Public Methods
    #

    def open(self, startFrame=1, endFrame=1):
        """Do some initialization operations.
        """
        VectorWriter.open(self, startFrame, endFrame)

        self.file = open(self.outputFileName, "w")

        self._printHeader()

    def close(self):
        """Do some finalization operation.
        """
        self._printFooter()

        if self.file:
            self.file.close()

        # remember to call the close method of the parent as last
        VectorWriter.close(self)
        print("close...")

    def printCanvas(self, scene, doPrintPolygons=True, doPrintEdges=False,
            showHiddenEdges=False):
        """Convert the scene representation to SVG.
        """

        Objects = scene.objects

        render = scene.render
        framenumber = scene.frame_current

        if self.animation:
            framestyle = "display:none"
        else:
            framestyle = "display:block"

        # Assign an id to this group so we can set properties on it using DOM
        self.file.write("<g id=\"frame%d\" style=\"%s\">\n" %
                (framenumber, framestyle))

        for obj in Objects:

            if obj.type != 'MESH':
                continue

            self.file.write("<g id=\"%s\">\n" % obj.name)

            mesh = obj.data

            if doPrintPolygons:
                self._printPolygons(mesh)

            if False:#doPrintEdges:
                self._printEdges(mesh, showHiddenEdges)

            self.file.write("</g>\n")

    def finalFrame(self):
        self.file.write("</g>\n")

    ##
    # Private Methods
    #

    def _calcCanvasCoord(self, v):
        """Convert vertex in scene coordinates to canvas coordinates.
        """

        pt = Vector([0, 0, 0])

        mW = float(self.canvasSize[0]) / 2.0
        mH = float(self.canvasSize[1]) / 2.0

        # rescale to canvas size
        pt[0] = v.co[0] * mW + mW
        pt[1] = v.co[1] * mH + mH
        pt[2] = v.co[2]

        # For now we want (0,0) in the top-left corner of the canvas.
        # Mirror and translate along y
        pt[1] *= -1
        pt[1] += self.canvasSize[1]

        return pt

    def _printHeader(self):
        """Print SVG header."""

        self.file.write("<?xml version=\"1.0\"?>\n")
        self.file.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"\n")
        self.file.write("\t\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n")
        self.file.write("<svg version=\"1.0\"\n")
        self.file.write("\txmlns=\"http://www.w3.org/2000/svg\"\n")
        self.file.write("\twidth=\"%d\" height=\"%d\">\n\n" %
                self.canvasSize)

        if self.animation:
            delay = 1000 / self.fps

            self.file.write("""\n<script type="text/javascript"><![CDATA[
            globalStartFrame=%d;
            globalEndFrame=%d;

            timerID = setInterval("NextFrame()", %d);
            globalFrameCounter=%d;
            \n""" % (self.startFrame, self.endFrame, delay, self.startFrame))

            self.file.write("""\n
            function NextFrame()
            {
              currentElement  = document.getElementById('frame'+globalFrameCounter)
              previousElement = document.getElementById('frame'+(globalFrameCounter-1))

              if (!currentElement)
              {
                return;
              }

              if (globalFrameCounter > globalEndFrame)
              {
                clearInterval(timerID)
              }
              else
              {
                if(previousElement)
                {
                    previousElement.style.display="none";
                }
                currentElement.style.display="block";
                globalFrameCounter++;
              }
            }
            \n]]></script>\n
            \n""")

    def _printFooter(self):
        """Print the SVG footer."""

        self.file.write("\n</svg>\n")

    def _printPolygons(self, mesh):
        """Print the selected (visible) polygons.
        """
        
        print("_printPolygons: ",mesh)
        #mesh.update(calc_edges=False, calc_tessface=True)

        if len(mesh.polygons) == 0:
            return

        self.file.write("<g>\n")
        backface = mesh.polygon_layers_int["backface"].data
        rlayer = mesh.polygon_layers_int["rlayer"].data
        glayer = mesh.polygon_layers_int["glayer"].data
        blayer = mesh.polygon_layers_int["blayer"].data
        alayer = mesh.polygon_layers_int["alayer"].data

        for face in mesh.polygons:
            print("_printPolygons: ",mesh.polygons)
            if not backface[face.index].value:
                continue
            
            print("  selected!")

            self.file.write("<path d=\"")

            #p = self._calcCanvasCoord(face.verts[0])
            p = self._calcCanvasCoord(mesh.vertices[face.vertices[0]])
            self.file.write("M %g,%g L " % (p[0], p[1]))

            for idx in face.vertices[1:]:
                p = self._calcCanvasCoord(mesh.vertices[idx])
                self.file.write("%g,%g " % (p[0], p[1]))

            # get rid of the last blank space, just cosmetics here.
            #self.file.seek(-1, 1)
            self.file.write(" z\"\n")

            # take as face color the first vertex color
            #if mesh.vertex_colors.active.data[face.vertices[0]]:
            #    fcol = mesh.vertex_colors.active.data[face.vertices[0]].color
            #    print("out color: ",fcol)
            #    if len(fcol)==3:
            #      color = [int(fcol[0]*255), int(fcol[1]*255), int(fcol[2]*255), 255]
            #    else:
            #      color = [int(fcol[0]*255), int(fcol[1]*255), int(fcol[2]*255), int(fcol[3]*255)]
            #else:
            #    color = [255, 255, 255, 255]
            color = [rlayer[face.index].value, glayer[face.index].value, blayer[face.index].value, alayer[face.index].value]
            
            print("  color: ",color)

            # Convert the color to the #RRGGBB form
            str_col = "#%02X%02X%02X" % (color[0], color[1], color[2])

            # Handle transparent polygons
            opacity_string = ""
            if color[3] != 255:
                opacity = float(color[3]) / 255.0
                opacity_string = " fill-opacity: %g; stroke-opacity: %g; opacity: 1;" % (opacity, opacity)
                #opacity_string = "opacity: %g;" % (opacity)

            self.file.write("\tstyle=\"fill:" + str_col + ";")
            self.file.write(opacity_string)

            # use the stroke property to alleviate the "adjacent edges" problem,
            # we simulate polygon expansion using borders,
            # see http://www.antigrain.com/svg/index.html for more info
            stroke_width = 1.0

            vrm = bpy.context.scene.VRM
            # EXPANSION TRICK is not that useful where there is transparency
            if vrm.polygonsEXPANSION_TRICK and color[3] == 255:
                # str_col = "#000000" # For debug
                self.file.write(" stroke:%s;\n" % str_col)
                self.file.write(" stroke-width:" + str(stroke_width) + ";\n")
                self.file.write(" stroke-linecap:round;stroke-linejoin:round")

            self.file.write("\"/>\n")

        self.file.write("</g>\n")
        
    def _printFreeStyleStroke(self, stroke, thickness, caps, r, g, b, a):
      self.file.write("<g>\n <path fill=\"none\" ")
      self.file.write("stroke-width=\"" + str(thickness) + "\" ")
      #self.file.write("stroke=\"rgba(" + str(int(r*255)) + "," + str(int(g*255)) + "," + str(int(b*255)) + "," + str(a) + ")\" ")
      self.file.write("stroke=\"rgb(" + str(int(r*255)) + "," + str(int(g*255)) + "," + str(int(b*255)) + ")\" ")
      self.file.write("d=\"")
      bFirst = True
      for v in stroke:
        x, y = v.point
        y = self.canvasSize[1]-y
        if bFirst:
          bFirst = False
          self.file.write("M " + str(x) + " " + str(y) + " ")
        else:
          self.file.write("L " + str(x) + " " + str(y) + " ")
      self.file.write("\"/>\n</g>\n")

    def _printEdges(self, mesh, showHiddenEdges=False):
        """Print the wireframe using mesh edges.
        """

        vrm = bpy.context.scene.VRM
        stroke_width = vrm.edgesWIDTH
        stroke_col = vrm.edgesCOLOR
        
        self.file.write("<g>\n")

        for e in mesh.edges:

            hidden_stroke_style = ""

            if e.select == 0:
                if showHiddenEdges == False:
                    continue
                else:
                    hidden_stroke_style = ";\n stroke-dasharray:3, 3"

            p1 = self._calcCanvasCoord(mesh.vertices[e.vertices[0]])
            p2 = self._calcCanvasCoord(mesh.vertices[e.vertices[1]])

            self.file.write("<line x1=\"%g\" y1=\"%g\" x2=\"%g\" y2=\"%g\"\n"
                    % (p1[0], p1[1], p2[0], p2[1]))
            self.file.write(" style=\"stroke:rgb(" + str(stroke_col[0]) + "," + str(stroke_col[1]) + "," + str(stroke_col[2]) + ");")
            self.file.write(" stroke-width:" + str(stroke_width) + ";\n")
            self.file.write(" stroke-linecap:round;stroke-linejoin:round")
            self.file.write(hidden_stroke_style)
            self.file.write("\"/>\n")

        self.file.write("</g>\n")


## SWF Writer

try:
    from ming import *
    SWFSupported = True
except:
    SWFSupported = False


class SWFVectorWriter(VectorWriter):
    """A concrete class for writing SWF output.
    """

    def __init__(self, fileName):
        """Simply call the parent Contructor.
        """
        Ming_useSWFVersion(8)
        VectorWriter.__init__(self, fileName)

        self.movie = None
        self.sprite = None
        self.s = None
        self.scale = 10

    ##
    # Public Methods
    #

    def open(self, startFrame=1, endFrame=1):
        """Do some initialization operations.
        """
        VectorWriter.open(self, startFrame, endFrame)
        #print("self.movie = SWFMovie()")
        self.movie = SWFMovie()
        self.movie.setBackground(255, 255, 255);
        #print("self.movie.setDimension("+str(self.canvasSize[0])+", "+str(self.canvasSize[1])+")")
        sm = max(self.canvasSize[0], self.canvasSize[1])
        if sm!=0:
          self.scale = 16000.0/sm
          Ming_setScale(self.scale)
        #Ming_setScale(10)
        self.movie.setDimension(self.canvasSize[0], self.canvasSize[1])
        if self.animation:
            self.movie.setRate(self.fps)
            numframes = endFrame - startFrame + 1
            self.movie.setFrames(numframes)

    def close(self):
        """Do some finalization operation.
        """
        #print("self.movie.save(self.outputFileName)")
        self.movie.save(self.outputFileName)

        # remember to call the close method of the parent
        VectorWriter.close(self)

    def printCanvas(self, scene, doPrintPolygons=True, doPrintEdges=False,
            showHiddenEdges=False):
        """Convert the scene representation to SVG.
        """
        render = scene.render
        framenumber = scene.frame_current

        Objects = scene.objects

        if self.sprite:
            #print("self.movie.remove(self.sprite)")
            self.movie.remove(self.sprite)

        #print("sprite = SWFSprite()")
        sprite = SWFSprite()
        self.s = sprite

        for obj in Objects:

            if(obj.type != 'MESH'):
                continue

            mesh = obj.data

            if doPrintPolygons:
                self._printPolygons(mesh, sprite)

            if False:#doPrintEdges:
                self._printEdges(mesh, sprite, showHiddenEdges)

    def finalFrame(self):
        #print("sprite.nextFrame()")
        self.s.nextFrame()
        #print("i = self.movie.add(sprite)")
        i = self.movie.add(self.s)
        # Remove the instance the next time
        #print("self.sprite = i")
        self.sprite = i
        if self.animation:
            self.movie.nextFrame()

    ##
    # Private Methods
    #

    def _calcCanvasCoord(self, v):
        """Convert vertex in scene coordinates to canvas coordinates.
        """

        pt = Vector([0, 0, 0])

        mW = float(self.canvasSize[0]) / 2.0
        mH = float(self.canvasSize[1]) / 2.0

        # rescale to canvas size
        pt[0] = v.co[0] * mW + mW
        pt[1] = v.co[1] * mH + mH
        pt[2] = v.co[2]

        # For now we want (0,0) in the top-left corner of the canvas.
        # Mirror and translate along y
        pt[1] *= -1
        pt[1] += self.canvasSize[1]

        return pt

    def _printPolygons(self, mesh, sprite):
        """Print the selected (visible) polygons.
        """

        if len(mesh.polygons) == 0:
            return
        
        backface = mesh.polygon_layers_int["backface"].data
        rlayer = mesh.polygon_layers_int["rlayer"].data
        glayer = mesh.polygon_layers_int["glayer"].data
        blayer = mesh.polygon_layers_int["blayer"].data
        alayer = mesh.polygon_layers_int["alayer"].data

        for face in mesh.polygons:
            if not backface[face.index].value:
                continue

            color = [rlayer[face.index].value, glayer[face.index].value, blayer[face.index].value, alayer[face.index].value]
            
            #print("s = SWFShape()")
            s = SWFShape()
            #print("f = s.addFill("+str(color[0])+", "+str(color[1])+", "+str(color[2])+", "+str(color[3])+")")
            f = s.addFill(color[0], color[1], color[2], color[3])
            #print("s.setRightFill(f)")
            s.setLeftFill(f)

            # The starting point of the shape
            p0 = self._calcCanvasCoord(mesh.vertices[face.vertices[0]])
            #print("s.movePenTo("+str(p0[0])+", "+str(p0[1])+")")
            s.movePenTo(p0[0], p0[1])

            for idx in face.vertices[1:]:
                p = self._calcCanvasCoord(mesh.vertices[idx])
                #print("s.drawLineTo("+str(p[0])+", "+str(p[1])+")")
                s.drawLineTo(p[0], p[1])

            # Closing the shape
            #print("s.drawLineTo("+str(p0[0])+", "+str(p0[1])+")")
            s.drawLineTo(p0[0], p0[1])

            #print("s.end()")
            s.end()
            #print("sprite.add(s)")
            sprite.add(s)
    
    def _printFreeStyleStroke(self, stroke, thickness, caps, r, g, b, a):
      s = SWFShape()
      if caps=='ROUND':
        c = mingc.SWF_LINESTYLE_FLAG_ENDCAP_ROUND
      elif caps=='SQUARE':
        c = mingc.SWF_LINESTYLE_FLAG_ENDCAP_SQUARE
      else:
        c = mingc.SWF_LINESTYLE_FLAG_ENDCAP_NONE
      
      ti = int(thickness)
      if ti<1:
        ti = 1
      Ming_setScale(1)
      s.setLine2(int(thickness), int(r*255), int(g*255), int(b*255), int(a*255), c, 1)
      Ming_setScale(self.scale)
      bFirst = True
      
      for v in stroke:
        x, y = v.point
        y = self.canvasSize[1]-y
        
        if bFirst:
          s.movePenTo(x, y)
          bFirst = False
        else:
          s.drawLineTo(x, y)
      s.end()
      self.s.add(s)

    def _printEdges(self, mesh, sprite, showHiddenEdges=False):
        """Print the wireframe using mesh edges.
        """

        vrm = bpy.context.scene.VRM
        stroke_width = vrm.edgesWIDTH
        stroke_col = vrm.edgesCOLOR

        s = SWFShape()

        for e in mesh.edges:

            # Next, we set the line width and color for our shape.
            #print(stroke_width)
            s.setLine(int(stroke_width), stroke_col[0], stroke_col[1], stroke_col[2], 255)
            #s.setLine(2, 0,0,0, 255)

            if e.select == 0:
                if showHiddenEdges == False:
                    continue
                else:
                    # SWF does not support dashed lines natively, so -for now-
                    # draw hidden lines thinner and half-trasparent
                    s.setLine(int(stroke_width / 2), stroke_col[0], stroke_col[1], stroke_col[2], 128)
                    #s.setLine(2, 0,0,0, 128)

            p1 = self._calcCanvasCoord(mesh.vertices[e.vertices[0]])
            p2 = self._calcCanvasCoord(mesh.vertices[e.vertices[1]])

            s.movePenTo(p1[0], p1[1])
            s.drawLineTo(p2[0], p2[1])

        s.end()
        sprite.add(s)


## PDF Writer

try:
    from reportlab.pdfgen import canvas
    PDFSupported = True
except:
    PDFSupported = False


class PDFVectorWriter(VectorWriter):
    """A concrete class for writing PDF output.
    """

    def __init__(self, fileName):
        """Simply call the parent Contructor.
        """
        VectorWriter.__init__(self, fileName)

        self.canvas = None

    ##
    # Public Methods
    #

    def open(self, startFrame=1, endFrame=1):
        """Do some initialization operations.
        """
        VectorWriter.open(self, startFrame, endFrame)
        size = (self.canvasSize[0], self.canvasSize[1])
        self.canvas = canvas.Canvas(self.outputFileName, pagesize=size, bottomup=0)

    def close(self):
        """Do some finalization operation.
        """
        self.canvas.save()

        # remember to call the close method of the parent
        VectorWriter.close(self)

    def printCanvas(self, scene, doPrintPolygons=True, doPrintEdges=False,
            showHiddenEdges=False):
        """Convert the scene representation to SVG.
        """
        context = scene.getRenderingContext()
        framenumber = context.currentFrame()

        Objects = scene.objects

        for obj in Objects:

            if(obj.type != 'MESH'):
                continue

            mesh = obj.data

            if doPrintPolygons:
                self._printPolygons(mesh)

            if doPrintEdges:
                self._printEdges(mesh, showHiddenEdges)

        self.canvas.showPage()

    ##
    # Private Methods
    #

    def _calcCanvasCoord(self, v):
        """Convert vertex in scene coordinates to canvas coordinates.
        """

        pt = Vector([0, 0, 0])

        mW = float(self.canvasSize[0]) / 2.0
        mH = float(self.canvasSize[1]) / 2.0

        # rescale to canvas size
        pt[0] = v.co[0] * mW + mW
        pt[1] = v.co[1] * mH + mH
        pt[2] = v.co[2]

        # For now we want (0,0) in the top-left corner of the canvas.
        # Mirror and translate along y
        pt[1] *= -1
        pt[1] += self.canvasSize[1]

        return pt

    def _printPolygons(self, mesh):
        """Print the selected (visible) polygons.
        """

        if len(mesh.faces) == 0:
            return

        for face in mesh.faces:
            if not face.sel:
                continue

            if face.col:
                fcol = face.col[0]
                color = [fcol.r / 255.0, fcol.g / 255.0, fcol.b / 255.0,
                        fcol.a / 255.0]
            else:
                color = [1, 1, 1, 1]

            self.canvas.setFillColorRGB(color[0], color[1], color[2])
            # For debug
            self.canvas.setStrokeColorRGB(0, 0, 0)

            path = self.canvas.beginPath()

            # The starting point of the path
            p0 = self._calcCanvasCoord(face.verts[0])
            path.moveTo(p0[0], p0[1])

            for v in face.verts[1:]:
                p = self._calcCanvasCoord(v)
                path.lineTo(p[0], p[1])

            # Closing the shape
            path.close()

            self.canvas.drawPath(path, stroke=0, fill=1)

    def _printEdges(self, mesh, showHiddenEdges=False):
        """Print the wireframe using mesh edges.
        """

        stroke_width = config.edges['WIDTH']
        stroke_col = config.edges['COLOR']

        self.canvas.setLineCap(1)
        self.canvas.setLineJoin(1)
        self.canvas.setLineWidth(stroke_width)
        self.canvas.setStrokeColorRGB(stroke_col[0] / 255.0, stroke_col[1] / 255.0,
            stroke_col[2] / 255)

        for e in mesh.edges:

            self.canvas.setLineWidth(stroke_width)

            if e.sel == 0:
                if showHiddenEdges == False:
                    continue
                else:
                    # PDF does not support dashed lines natively, so -for now-
                    # draw hidden lines thinner
                    self.canvas.setLineWidth(stroke_width / 2.0)

            p1 = self._calcCanvasCoord(e.v1)
            p2 = self._calcCanvasCoord(e.v2)

            self.canvas.line(p1[0], p1[1], p2[0], p2[1])

def cmp(a, b):
  if a==b:
    return 0
  elif a < b:
    return -1
  else:
    return 1


# ---------------------------------------------------------------------
#
## Rendering Classes
#
# ---------------------------------------------------------------------

# A dictionary to collect different shading style methods
shadingStyles = dict()
shadingStyles['FLAT'] = None
shadingStyles['TOON'] = None

# A dictionary to collect different edge style methods
edgeStyles = dict()
edgeStyles['MESH'] = MeshUtils.isMeshEdge
edgeStyles['SILHOUETTE'] = MeshUtils.isSilhouetteEdge

# A dictionary to collect the supported output formats
outputWriters = dict()
outputWriters['SVG'] = SVGVectorWriter
if SWFSupported:
    outputWriters['SWF'] = SWFVectorWriter
if PDFSupported:
    outputWriters['PDF'] = PDFVectorWriter

gOutput = None


class Renderer:
    """Render a scene viewed from the active camera.

    This class is responsible of the rendering process, transformation and
    projection of the objects in the scene are invoked by the renderer.

    The rendering is done using the active camera for the current scene.
    """

    def __init__(self):
        """Make the rendering process only for the current scene by default.

        We will work on a copy of the scene, to be sure that the current scene do
        not get modified in any way.
        """

        # Render the current Scene, this should be a READ-ONLY property
        self._SCENE = bpy.context.scene

        # Use the aspect ratio of the scene rendering context
        render = self._SCENE.render

        aspect_ratio = float(render.resolution_x) / float(render.resolution_y)
        self.canvasRatio = (float(render.pixel_aspect_x) * aspect_ratio,
                            float(render.pixel_aspect_y)
                            )
        print(render.resolution_x)
        print(render.resolution_y)
        print(render.pixel_aspect_x)
        print(render.pixel_aspect_y)
        print(aspect_ratio)

        # Render from the currently active camera
        #self.cameraObj = self._SCENE.objects.camera

        self.lights = []

    ##
    # Public Methods
    #

    def doRendering(self, outputWriter, animation=False):
        global gOutput
        """Render picture or animation and write it out.

        The parameters are:
            - a Vector writer object that will be used to output the result.
            - a flag to tell if we want to render an animation or only the
              current frame.
        """

        render = self._SCENE.render
        origCurrentFrame = self._SCENE.frame_current

        # Handle the animation case
        if not animation:
            startFrame = origCurrentFrame
            endFrame = startFrame
            outputWriter.open()
        else:
            startFrame = self._SCENE.frame_start
            endFrame = self._SCENE.frame_end
            outputWriter.open(startFrame, endFrame)
        
        delTmpScene = True
        gOutput = outputWriter

        # Do the rendering process frame by frame
        print("Start Rendering of %d frames" % (endFrame - startFrame + 1))
        for f in range(startFrame, endFrame + 1):
            print("\n\nFrame: %d" % f)

            # FIXME To get the correct camera position we have to use +1 here.
            # Is there a bug somewhere in the Scene module?
            self._SCENE.frame_current = f + 1

            # Use some temporary workspace, a full copy of the scene
            #inputScene = self._SCENE.copy()
            
            bpy.ops.scene.new(type='FULL_COPY')
            newScene = bpy.data.scenes[-1]
            bpy.context.screen.scene = newScene
            inputScene = newScene
            self.cameraObj = inputScene.camera

            # To get the objects at this frame remove the +1 ...
            ctx = inputScene.render
            inputScene.frame_current = f
            

            # Get a projector for this camera.
            # NOTE: the projector wants object in world coordinates,
            # so we should remember to apply modelview transformations
            # _before_ we do projection transformations.
            self.proj = Projector(self.cameraObj, self.canvasRatio)
            
            try:
                renderedScene = self.doRenderScene(inputScene)
            except:
                print("There was an error! Aborting.")
                import traceback
                print(traceback.print_exc())
                gOutput = None

                bpy.context.screen.scene = self._SCENE #self._SCENE.makeCurrent()
                if delTmpScene:
                  bpy.data.scenes.remove(inputScene) #Scene.Unlink(renderedScene)
                  #del inputScene
                return
            
            vrm = bpy.context.scene.VRM
            outputWriter.printCanvas(renderedScene,
                    doPrintPolygons=vrm.polygonsSHOW,
                    doPrintEdges=vrm.edgesSHOW,
                    showHiddenEdges=vrm.edgesSHOW_HIDDEN)

            # delete the rendered scene
            bpy.context.screen.scene = self._SCENE #self._SCENE.makeCurrent()
            if delTmpScene:
              print(inputScene)
              bpy.data.scenes.remove(inputScene) #Scene.Unlink(renderedScene)
              #del inputScene
            
            bpy.ops.render.render(scene=self._SCENE.name)
            
            outputWriter.finalFrame()

        gOutput = None
        outputWriter.close()
        print("Done!")
        self._SCENE.frame_set(origCurrentFrame)

    def doRenderScene(self, workScene):
        """Control the rendering process.

        Here we control the entire rendering process invoking the operation
        needed to transform and project the 3D scene in two dimensions.
        """

        print("doRenderScene")
        # global processing of the scene
        #self._splitFaces(workScene)

        self._filterHiddenObjects(workScene)

        self._buildLightSetup(workScene)

        self._doSceneClipping(workScene)

        self._doConvertGeometricObjsToMesh(workScene)
        
        vrm = bpy.context.scene.VRM
        if vrm.outputJOIN_OBJECTS:
            self._joinMeshObjectsInScene(workScene)

        self._doSceneDepthSorting(workScene)

        # Per object activities

        Objects = workScene.objects

        print("Total Objects: %d" % len(Objects))
        for i, obj in enumerate(Objects):
            print("\n\n-------")
            print("Rendering Object: %d" % i)

            if obj.type != 'MESH':
                print("Only Mesh supported! - Skipping type:", obj.type)
                continue

            print("Rendering: ", obj.name)

            mesh = obj.data
            #mesh.update(calc_edges=False, calc_tessface=True)

            self._doModelingTransformation(mesh, obj.matrix_world)

            self._doBackFaceCulling(mesh)

            # When doing HSR with NEWELL we may want to flip all normals
            # toward the viewer
            if vrm.polygonsHSR == "NEWELL":
                for f in mesh.faces:
                    f.sel = 1 - f.sel
                mesh.flipNormals()
                for f in mesh.faces:
                    f.sel = 1

            self._doLighting(mesh)

            # Do "projection" now so we perform further processing
            # in Normalized View Coordinates
            self._doProjection(mesh, self.proj)

            self._doViewFrustumClipping(mesh, obj, workScene)

            self._doHiddenSurfaceRemoval(mesh, obj, workScene)

            self._doEdgesStyle(mesh, edgeStyles[vrm.edgesSTYLE])

            # Update the object data, important! :)
            #mesh.update(calc_edges=False, calc_tessface=True)

        return workScene

    ##
    # Private Methods
    #

    # Utility methods
    
    def _splitFaces(self, workScene):
      for o in workScene.objects:
       if o.type == 'MESH':
          mesh = o.data
          bpy.context.scene.objects.active = o
          bpy.ops.object.select_all(action='DESELECT')
          o.select = True
          bpy.context.scene.objects.active = o
          #bpy.ops.object.mode_set({"scene": workScene}, mode='EDIT')
          bm = bmesh.new()
          bm.from_mesh(mesh)
          
          bmesh.ops.triangulate(bm, faces=bm.faces)
          #bmesh.ops.split_edges(bm, edges=bm.edges, use_verts=True)
          #bmesh.ops.split(bm, geom=bm.faces, dest=bm)
          bmesh.ops.split_edges(bm, edges=bm.edges[:])
          
          #bmesh.update_edit_mesh(mesh, True)
          
          bm.to_mesh(mesh)
          bm.free()
          #bpy.ops.object.mode_set({"scene": workScene}, mode='OBJECT')
          mesh.update(calc_edges=False, calc_tessface=True)

    def _getObjPosition(self, obj):
        """Return the obj position in World coordinates.
        """
        return obj.matrix_world.to_translation()

    def _cameraViewVector(self):
        """Get the View Direction form the camera matrix.
        """
        print(self.cameraObj.matrix_world)
        v = Vector(self.cameraObj.matrix_world.transposed()[2].copy())
        v.resize_3d()
        return v

    # Faces methods

    def _isFaceVisible(self, face, mesh):
        """Determine if a face of an object is visible from the current camera.

        The view vector is calculated from the camera location and one of the
        vertices of the face (expressed in World coordinates, after applying
        modelview transformations).

        After those transformations we determine if a face is visible by
        computing the angle between the face normal and the view vector, this
        angle has to be between -90 and 90 degrees for the face to be visible.
        This corresponds somehow to the dot product between the two, if it
        results > 0 then the face is visible.

        There is no need to normalize those vectors since we are only interested in
        the sign of the cross product and not in the product value.

        NOTE: here we assume the face vertices are in WorldCoordinates, so
        please transform the object _before_ doing the test.
        """

        normal = Vector(face.normal)
        camPos = self._getObjPosition(self.cameraObj)
        view_vect = None

        # View Vector in orthographics projections is the view Direction of
        # the camera
        if self.cameraObj.data.type == 'ORTHO':
            view_vect = self._cameraViewVector()

        # View vector in perspective projections can be considered as
        # the difference between the camera position and one point of
        # the face, we choose the farthest point from the camera.
        if self.cameraObj.data.type == 'PERSP':
            vv = max([((camPos - Vector(mesh.vertices[idx].co)).length, (camPos - Vector(mesh.vertices[idx].co))) for idx in face.vertices])
            view_vect = vv[1]

        # if d > 0 the face is visible from the camera
        d = view_vect * normal

        if d > 0:
            return True
        else:
            return False

    # Scene methods

    def _filterHiddenObjects(self, scene):
        """Discard object that are on hidden layers in the scene.
        """

        for o in scene.objects:
            if o.hide_render:
                scene.objects.unlink(o)

        scene.update()

    def _buildLightSetup(self, scene):
        # Get the list of lighting sources
        obj_lst = scene.objects
        self.lights = [o for o in obj_lst if o.type == 'LAMP']

        # When there are no lights we use a default lighting source
        # that have the same position of the camera
        if len(self.lights) == 0:
            l = Lamp.New('Lamp')
            lobj = Object.New('Lamp')
            lobj.loc = self.cameraObj.loc
            lobj.link(l)
            self.lights.append(lobj)

    def _doSceneClipping(self, scene):
        """Clip whole objects against the View Frustum.

        For now clip away only objects according to their center position.
        """

        cam_pos = self._getObjPosition(self.cameraObj)
        view_vect = self._cameraViewVector()
        
        near = self.cameraObj.data.clip_start
        far = self.cameraObj.data.clip_end

        aspect = float(self.canvasRatio[0]) / float(self.canvasRatio[1])
        fovy = atan(0.5 / aspect / (self.cameraObj.data.lens / 32))
        fovy = fovy * 360.0 / pi

        Objects = scene.objects

        for o in Objects:
            if o.type != 'MESH':
                continue

            """
            obj_vect = Vector(cam_pos) - self._getObjPosition(o)

            d = obj_vect*view_vect
            theta = AngleBetweenVecs(obj_vect, view_vect)

            # if the object is outside the view frustum, clip it away
            if (d < near) or (d > far) or (theta > fovy):
                scene.objects.unlink(o)
            """

            # Use the object bounding box
            # (whose points are already in WorldSpace Coordinate)
            
            bbox = []
            for [x,y,z] in o.bound_box:
              bbox.append(Vector((x,y,z)))

            bb = [o.matrix_world.copy() * v.copy() for v in bbox]

            points_outside = 0
            for p in bb:
                p_vect = Vector(cam_pos) - Vector(p)

                d = p_vect * view_vect
                theta = p_vect.angle(view_vect)

                # Is this point outside the view frustum?
                if (d < near) or (d > far) or (theta > fovy):
                    points_outside += 1

            # If the bb is all outside the view frustum we clip the whole
            # object away
            if points_outside == len(bb):
                scene.objects.unlink(o)

    def _doConvertGeometricObjsToMesh(self, scene):
        """Convert all "geometric" objects to mesh ones.
        """
        geometricObjTypes = ['MESH', 'SURFACE', 'CURVE', 'TEXT']
        #geometricObjTypes = ['Mesh', 'Surf', 'Curve']

        Objects = scene.objects
        
        objList = [o for o in Objects if o.type in geometricObjTypes]
        for obj in objList:
            old_obj = obj
            obj = self._convertToRawMeshObj(obj)
            scene.objects.link(obj)
            scene.objects.unlink(old_obj)
            print(old_obj.type)

            # XXX Workaround for Text and Curve which have some normals
            # inverted when they are converted to Mesh, REMOVE that when
            # blender will fix that!!
            if old_obj.type in ['CURVE', 'TEXT']:
                me = obj.getData(mesh=1)

                for f in me.faces:
                    f.sel = 1
                for v in me.verts:
                    v.sel = 1

                me.remDoubles(0)
                me.triangleToQuad()
                me.recalcNormals()
                me.update()

    def _doSceneDepthSorting(self, scene):
        """Sort objects in the scene.

        The object sorting is done accordingly to the object centers.
        """

        c = self._getObjPosition(self.cameraObj)

        by_obj_center_pos = (lambda o1, o2:
                (o1.type == 'MESH' and o2.type == 'MESH') and
                cmp((self._getObjPosition(o1) - Vector(c)).length,
                    (self._getObjPosition(o2) - Vector(c)).length)
            )

        # Implement sorting by bounding box, the object with the bb
        # nearest to the camera should be drawn as last.
        by_nearest_bbox_point = (lambda o1, o2:
                (o1.type == 'MESH' and o2.type == 'MESH') and
                cmp(min([(Vector(p) - Vector(c)).length for p in [o1.matrix_world.copy() * v.copy() for v in [Vector(o1.bound_box[0:3]), Vector(o1.bound_box[3:6]), Vector(o1.bound_box[6:9]), Vector(o1.bound_box[9:12]), Vector(o1.bound_box[12:15]), Vector(o1.bound_box[15:18], Vector(o1.bound_box[18:21]), Vector(o1.bound_box[21:24]))]]]),
                    min([(Vector(p) - Vector(c)).length for p in [o2.matrix_world.copy() * v.copy() for v in [Vector(o2.bound_box[0:3]), Vector(o2.bound_box[3:6]), Vector(o2.bound_box[6:9]), Vector(o2.bound_box[9:12]), Vector(o2.bound_box[12:15]), Vector(o2.bound_box[15:18], Vector(o2.bound_box[18:21]), Vector(o2.bound_box[21:24]))]]])
                )
            )

        Objects = list(scene.objects)

        #Objects.sort(by_obj_center_pos)
        Objects.sort(key=cmp_to_key(by_nearest_bbox_point))

        # update the scene
        for o in Objects:
            scene.objects.unlink(o)
            scene.objects.link(o)

    def _joinMeshObjectsInScene(self, scene):
        """Merge all the Mesh Objects in a scene into a single Mesh Object.
        """

        bpy.ops.object.select_all(action='DESELECT')
        
        oList = [o for o in scene.objects if o.type == 'MESH']

        # FIXME: Object.join() do not work if the list contains 1 object
        if len(oList) == 1:
            return
        
        mesh = bpy.data.meshes.new('BigOneMesh')
        bigObj = bpy.data.objects.new('BigOne', mesh)
        scene.objects.link(bigObj)
        scene.objects.active = bigObj
        bigObj.select = True
        
        bpy.ops.object.select_by_type(extend=False, type='MESH')
        #print (' bpy.context.active_object =', bpy.context.active_object)
        #print (' bpy.context.selected_objects =',bpy.context.selected_objects)
        bpy.ops.object.join()
        

        for o in bpy.context.selected_objects:
          if not o.name.startswith('BigOne'):
            scene.objects.unlink(o)

        scene.update()

    # Per object/mesh methods

    def _convertToRawMeshObj(self, object):
        """Convert geometry based object to a mesh object.
        """
        #me = bpy.data.meshes.new('RawMesh_' + object.name)
        me = object.to_mesh(bpy.context.scene, False, 'PREVIEW')

        newObject = bpy.data.objects.new('RawMesh_' + object.name, me)

        # If the object has no materials set a default material
        if not me.materials:
            tmpMat = bpy.data.materials.new('Tmp')
            me.materials.append(tmpMat)
            #for f in me.faces: f.mat = 0

        newObject.matrix_world = object.matrix_world

        return newObject

    def _doModelingTransformation(self, mesh, matrix):
        """Transform object coordinates to world coordinates.

        This step is done simply applying to the object its tranformation
        matrix and recalculating its normals.
        """
        # XXX FIXME: blender do not transform normals in the right way when
        # there are negative scale values
        if matrix[0][0] < 0 or matrix[1][1] < 0 or matrix[2][2] < 0:
            print("WARNING: Negative scales, expect incorrect results!")

        mesh.transform(matrix, True)

    def _doBackFaceCulling(self, mesh):
        """Simple Backface Culling routine.

        At this level we simply do a visibility test face by face and then
        select the vertices belonging to visible faces.
        """
        #mesh.update(calc_edges=False, calc_tessface=True)
        mesh.polygon_layers_int.new("backface")
        backface = mesh.polygon_layers_int["backface"].data

        # Select all vertices, so edges can be displayed even if there are no
        # faces
        for v in mesh.vertices:
            v.select = 1

        bpy.context.tool_settings.mesh_select_mode=[False,False,True] #face
        # Loop on faces
        for f in mesh.polygons:
            backface[f.index].value = 0
            if self._isFaceVisible(f, mesh):
                backface[f.index].value = 1

    def _doLighting(self, mesh):
        """Apply an Illumination and shading model to the object.

        The model used is the Phong one, it may be inefficient,
        but I'm just learning about rendering and starting from Phong seemed
        the most natural way.
        """

        vrm = bpy.context.scene.VRM
        #mesh.update(calc_edges=False, calc_tessface=True)
        
        # If the mesh has vertex colors already, use them,
        # otherwise turn them on and do some calculations
        if len(mesh.vertex_colors):
            return
        mesh.vertex_colors.new()

        materials = mesh.materials

        camPos = self._getObjPosition(self.cameraObj)
        mesh.polygon_layers_int.new("rlayer")
        mesh.polygon_layers_int.new("glayer")
        mesh.polygon_layers_int.new("blayer")
        mesh.polygon_layers_int.new("alayer")
        backface = mesh.polygon_layers_int["backface"].data
        rlayer = mesh.polygon_layers_int["rlayer"].data
        glayer = mesh.polygon_layers_int["glayer"].data
        blayer = mesh.polygon_layers_int["blayer"].data
        alayer = mesh.polygon_layers_int["alayer"].data
        
        Iamb = Vector(bpy.context.scene.world.ambient_color)

        # We do per-face color calculation (FLAT Shading), we can easily turn
        # to a per-vertex calculation if we want to implement some shading
        # technique. For an example see:
        # http://www.miralab.unige.ch/papers/368.pdf
        for f in mesh.polygons:
            if not backface[f.index].value:
                continue

            mat = None
            if materials:
                mat = materials[f.material_index]

            # A new default material
            if mat == None:
                mat = bpy.data.materials.new('defMat')

            # Check if it is a shadeless material
            elif mat.use_shadeless:
                I = mat.diffuse_color
                # Convert to a value between 0 and 255
                tmp_col = [c for c in I]
                
                rlayer[f.index].value = tmp_col[0]*255
                glayer[f.index].value = tmp_col[1]*255
                blayer[f.index].value = tmp_col[2]*255
                alayer[f.index].value = 255

                continue

            # do vertex color calculation

            TotDiffSpec = Vector([0.0, 0.0, 0.0])

            for l in self.lights:
                light_obj = l
                light_pos = self._getObjPosition(l)
                light = light_obj.data

                L = Vector(light_pos).normalized()
                #print("L: ",L)

                V = (Vector(camPos) - Vector(f.center)).normalized()
                #print("V: ",V)

                N = Vector(f.normal).normalized()
                #print("N: ",N)

                if vrm.polygonsSHADING == 'TOON':
                    NL = ShadingUtils.toonShading(N * L)
                else:
                    NL = (N * L)
                #print("NL: ",NL)

                # Should we use NL instead of (N*L) here?
                R = 2 * (N * L) * N - L
                #print("R: ",R)

                Ip = light.energy
                #print("Ip: ",Ip)

                # Diffuse co-efficient
                kd = mat.diffuse_intensity * Vector(mat.diffuse_color)
                #print("kd: ",kd)
                for i in [0, 1, 2]:
                    kd[i] *= light.color[i]
                #print("kd: ",kd)

                Idiff = Ip * kd * max(0, NL)
                #print("Idiff: ",Idiff)

                # Specular component
                ks = mat.specular_intensity * Vector(mat.specular_color)
                #print("ks: ",ks)
                ns = mat.specular_hardness
                #print("ns: ",ns)
                Ispec = Ip * ks * pow(max(0, (V * R)), ns)
                #print("Ispec: ",Ispec)

                TotDiffSpec += (Idiff + Ispec)
                #print("TotDiffSpec: ",TotDiffSpec)

            #print("")
            # Ambient component
            #print("Iamb: ",Iamb)
            ka = mat.ambient
            #print("ka: ",ka)

            # Emissive component (convert to a triplet)
            ki = Vector([mat.emit] * 3)
            #print("ki: ",ki)

            #I = ki + Iamb + (Idiff + Ispec)
            I = ki + (ka * Iamb) + TotDiffSpec
            #print("I: ",I)

            # Set Alpha component
            I = list(I)
            I.append(mat.alpha)

            # Clamp I values between 0 and 1
            I = [min(c, 1) for c in I]
            I = [max(0, c) for c in I]
            #print("I: ",I)

            # Convert to a value between 0 and 255
            tmp_col = [c for c in I]
            #print("vert color: ",tmp_col)
            
            rlayer[f.index].value = tmp_col[0]*255
            glayer[f.index].value = tmp_col[1]*255
            blayer[f.index].value = tmp_col[2]*255
            
            alayer[f.index].value = tmp_col[3]*255

            #for idx in f.vertices:
            #mesh.tessface_vertex_colors.active.data[f.index].color1 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #mesh.tessface_vertex_colors.active.data[f.index].color2 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #mesh.tessface_vertex_colors.active.data[f.index].color3 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #mesh.tessface_vertex_colors.active.data[f.index].color4 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #vcolTmp = mesh.tessface_vertex_colors[f.index].data
            #vcolTmp.color1 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #vcolTmp.color2 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #vcolTmp.color3 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #vcolTmp.color4 = (tmp_col[0], tmp_col[1], tmp_col[2])
            #vcolTmp.color4.a = tmp_col[3]
            
            #for vidx in f.vertices:
            #  mesh.vertex_colors.active.data[vidx].color = (tmp_col[0], tmp_col[1], tmp_col[2])
            
            #self.printVertexColors(mesh)
            
            #mesh.update(calc_edges=False, calc_tessface=True)
            #print("        : ",mesh.tessface_vertex_colors.active.data[f.index].color1)
    
    def printVertexColors(self, mesh):
      print("printVertexColors: ")
      for f in mesh.polygons:
        print("  face: ",f)
        for vidx in f.vertices:
          print("    ",mesh.vertex_colors.active.data[vidx].color)

    def _doProjection(self, mesh, projector):
        """Apply Viewing and Projection tranformations.
        """

        for v in mesh.vertices:
            p = projector.doProjection(v.co[:])
            v.co[0] = p[0]
            v.co[1] = p[1]
            v.co[2] = p[2]

        #mesh.recalcNormals()
        #mesh.update()

        # We could reeset Camera matrix, since now
        # we are in Normalized Viewing Coordinates,
        # but doung that would affect World Coordinate
        # processing for other objects

        #self.cameraObj.data.type = 1
        #self.cameraObj.data.scale = 2.0
        #m = Matrix().identity()
        #self.cameraObj.setMatrix(m)

    def _doViewFrustumClipping(self, mesh, ob, workScene):
        """Clip faces against the View Frustum.
        """

        # The Canonical View Volume, 8 vertices, and 6 faces,
        # We consider its face normals pointing outside

        v1 = [1, 1, -1]
        v2 = [1, -1, -1]
        v3 = [-1, -1, -1]
        v4 = [-1, 1, -1]
        v5 = [1, 1, 1]
        v6 = [1, -1, 1]
        v7 = [-1, -1, 1]
        v8 = [-1, 1, 1]

        cvv = []
        f1 = [v1, v4, v3, v2]
        cvv.append(f1)
        f2 = [v5, v6, v7, v8]
        cvv.append(f2)
        f3 = [v1, v2, v6, v5]
        cvv.append(f3)
        f4 = [v2, v3, v7, v6]
        cvv.append(f4)
        f5 = [v3, v4, v8, v7]
        cvv.append(f5)
        f6 = [v4, v1, v5, v8]
        cvv.append(f6)
        
        #mode = bpy.context.tool_settings.mesh_select_mode
        #bpy.context.tool_settings.mesh_select_mode = [False, False, True]

        nmesh = bpy.data.meshes[mesh.name]
        #nmesh.update(calc_edges=False, calc_tessface=True)
        #bpy.context.screen.scene = ob.users_scene[0]
        bpy.context.scene.objects.active = ob
        bpy.ops.object.select_all(action='DESELECT')
        ob.select = True
        bpy.context.scene.objects.active = ob
        #bpy.ops.object.mode_set(mode='OBJECT')
        #bpy.ops.object.mode_set({"scene": workScene}, mode='EDIT')
        #bpy.ops.mesh.select_all(action='DESELECT')
        
        bm = bmesh.new()
        bm.from_mesh(nmesh)
        #clippedfaces = nmesh.tessfaces[:]
        clippedfaces = bm.faces[:]
        facelist = clippedfaces[:]

        for clipface in cvv:

            clippedfaces = []

            for f in facelist:

                #newfaces = HSR.splitOn(clipface, f, return_positive_faces=False)
                newfaces = None

                if not newfaces:
                    # Check if the face is all outside the view frustum
                    # TODO: Do this test before, it is more efficient
                    points_outside = 0
                    for ve in f.verts:
                        #v = nmesh.vertices[idx].co
                        v = ve.co
                        if abs(v[0]) > (1 - EPS) or abs(v[1]) > (1 - EPS) or abs(v[2]) > (1 - EPS):
                            points_outside += 1

                    if points_outside != len(f.verts):
                        clippedfaces.append(f)
                    else:
                        #clippedfaces.select = True
                        #bmesh.ops.delete(bm, f, DEL_ONLYFACES)
                        #print(f)
                        bmesh.ops.delete(bm, geom=[f], context=2)
                else:
                    for nf in newfaces:
                        for v in nf:
                            nmesh.verts.append(v)

                        nf.mat = f.mat
                        nf.sel = f.sel
                        nf.col = [f.col[0]] * len(nf.v)

                        clippedfaces.append(nf)
            facelist = clippedfaces[:]
            #bpy.ops.mesh.delete(type='ONLY_FACE')

        #nmesh.tessfaces = facelist
        #nmesh.from_pydata(nmesh.vertices[:],nmesh.edges[:], facelist)
        #nmesh.tessfaces.foreach_set(facelist)
        #bmesh.update_edit_mesh(nmesh, True)
        bm.to_mesh(nmesh)
        bm.free()
        nmesh.update(calc_edges=False, calc_tessface=True)
        #bpy.ops.object.mode_set({"scene": workScene}, mode='OBJECT')
        
        #bpy.context.tool_settings.mesh_select_mode = mode

    # HSR routines
    def __simpleDepthSort(self, mesh, ob, workScene):
        """Sort faces by the furthest vertex.

        This simple mesthod is known also as the painter algorithm, and it
        solves HSR correctly only for convex meshes.
        """

        #global progress
        #mesh.update(calc_edges=False, calc_tessface=True)
        bpy.context.scene.objects.active = ob
        bpy.ops.object.select_all(action='DESELECT')
        ob.select = True
        bpy.context.scene.objects.active = ob
        #bpy.ops.object.mode_set({"scene": workScene}, mode='EDIT')
        bm = bmesh.new()
        bm.from_mesh(mesh)

        # The sorting requires circa n*log(n) steps
        n = len(mesh.polygons)
        #progress.setActivity("HSR: Painter", n * log(n))

        by_furthest_z = (lambda f1, f2: #progress.update() and
                cmp(max([v.co[2] for v in f1]), max([v.co[2] for v in f2]) + EPS)
                )

        # FIXME: using NMesh to sort faces. We should avoid that!
        #nmesh = NMesh.GetRaw(mesh.name)

        # remember that _higher_ z values mean further points
        #nmesh.faces.sort(by_furthest_z)
        #nmesh.faces.reverse()
        bm.faces.sort(key=lambda f: max([v.co[2] for v in f.verts]), reverse=True)
        bm.faces.index_update()

        #bmesh.update_edit_mesh(mesh, True)
        bm.to_mesh(mesh)
        bm.free()
        #mesh.update(calc_edges=False, calc_tessface=True)
        #bpy.ops.object.mode_set({"scene": workScene}, mode='OBJECT')
        #nmesh.update()

    def __newellDepthSort(self, mesh):
        """Newell's depth sorting.

        """

        #global progress

        # Find non planar quads and convert them to triangle
        #for f in mesh.faces:
        #    f.sel = 0
        #    if is_nonplanar_quad(f.v):
        #        print "NON QUAD??"
        #        f.sel = 1

        # Now reselect all faces
        for f in mesh.faces:
            f.sel = 1
        mesh.quadToTriangle()

        # FIXME: using NMesh to sort faces. We should avoid that!
        nmesh = NMesh.GetRaw(mesh.name)

        # remember that _higher_ z values mean further points
        nmesh.faces.sort(by_furthest_z)
        nmesh.faces.reverse()

        # Begin depth sort tests

        # use the smooth flag to set marked faces
        for f in nmesh.faces:
            f.smooth = 0

        facelist = nmesh.faces[:]
        maplist = []

        # The steps are _at_least_ equal to len(facelist), we do not count the
        # feces coming out from splitting!!
        progress.setActivity("HSR: Newell", len(facelist))
        #progress.setQuiet(True)

        while len(facelist):
            debug("\n----------------------\n")
            debug("len(facelits): %d\n" % len(facelist))
            P = facelist[0]

            pSign = sign(P.normal[2])

            # We can discard faces parallel to the view vector
            #if P.normal[2] == 0:
            #    facelist.remove(P)
            #    continue

            split_done = 0
            face_marked = 0

            for Q in facelist[1:]:

                debug("P.smooth: " + str(P.smooth) + "\n")
                debug("Q.smooth: " + str(Q.smooth) + "\n")
                debug("\n")

                qSign = sign(Q.normal[2])
                # TODO: check also if Q is parallel??

                # Test 0: We need to test only those Qs whose furthest vertex
                # is closer to the observer than the closest vertex of P.

                zP = [v.co[2] for v in P.v]
                zQ = [v.co[2] for v in Q.v]
                notZOverlap = min(zP) > max(zQ) + EPS

                if notZOverlap:
                    debug("\nTest 0\n")
                    debug("NOT Z OVERLAP!\n")
                    if Q.smooth == 0:
                        # If Q is not marked then we can safely print P
                        break
                    else:
                        debug("met a marked face\n")
                        continue

                # Test 1: X extent overlapping
                xP = [v.co[0] for v in P.v]
                xQ = [v.co[0] for v in Q.v]
                #notXOverlap = (max(xP) <= min(xQ)) or (max(xQ) <= min(xP))
                notXOverlap = min(xQ) >= (max(xP) - EPS) or min(xP) >= (max(xQ) - EPS)

                if notXOverlap:
                    debug("\nTest 1\n")
                    debug("NOT X OVERLAP!\n")
                    continue

                # Test 2: Y extent Overlapping
                yP = [v.co[1] for v in P.v]
                yQ = [v.co[1] for v in Q.v]
                #notYOverlap = max(yP) <= min(yQ) or max(yQ) <= min(yP)
                notYOverlap = min(yQ) >= (max(yP) - EPS) or min(yP) >= (max(yQ) - EPS)

                if notYOverlap:
                    debug("\nTest 2\n")
                    debug("NOT Y OVERLAP!\n")
                    continue

                # Test 3: P vertices are all behind the plane of Q
                n = 0
                for Pi in P:
                    d = qSign * HSR.Distance(Vector(Pi), Q)
                    if d <= EPS:
                        n += 1
                pVerticesBehindPlaneQ = (n == len(P))

                if pVerticesBehindPlaneQ:
                    debug("\nTest 3\n")
                    debug("P BEHIND Q!\n")
                    continue

                # Test 4: Q vertices in front of the plane of P
                n = 0
                for Qi in Q:
                    d = pSign * HSR.Distance(Vector(Qi), P)
                    if d >= -EPS:
                        n += 1
                qVerticesInFrontPlaneP = (n == len(Q))

                if qVerticesInFrontPlaneP:
                    debug("\nTest 4\n")
                    debug("Q IN FRONT OF P!\n")
                    continue

                # Test 5: Check if projections of polygons effectively overlap,
                # in previous tests we checked only bounding boxes.

                #if not projectionsOverlap(P, Q):
                if not (HSR.projectionsOverlap(P, Q) or HSR.projectionsOverlap(Q, P)):
                    debug("\nTest 5\n")
                    debug("Projections do not overlap!\n")
                    continue

                # We still can't say if P obscures Q.

                # But if Q is marked we do a face-split trying to resolve a
                # difficulty (maybe a visibility cycle).
                if Q.smooth == 1:
                    # Split P or Q
                    debug("Possibly a cycle detected!\n")
                    debug("Split here!!\n")

                    facelist = HSR.facesplit(P, Q, facelist, nmesh)
                    split_done = 1
                    break

                # The question now is: Does Q obscure P?

                # Test 3bis: Q vertices are all behind the plane of P
                n = 0
                for Qi in Q:
                    d = pSign * HSR.Distance(Vector(Qi), P)
                    if d <= EPS:
                        n += 1
                qVerticesBehindPlaneP = (n == len(Q))

                if qVerticesBehindPlaneP:
                    debug("\nTest 3bis\n")
                    debug("Q BEHIND P!\n")

                # Test 4bis: P vertices in front of the plane of Q
                n = 0
                for Pi in P:
                    d = qSign * HSR.Distance(Vector(Pi), Q)
                    if d >= -EPS:
                        n += 1
                pVerticesInFrontPlaneQ = (n == len(P))

                if pVerticesInFrontPlaneQ:
                    debug("\nTest 4bis\n")
                    debug("P IN FRONT OF Q!\n")

                # We don't even know if Q does obscure P, so they should
                # intersect each other, split one of them in two parts.
                if not qVerticesBehindPlaneP and not pVerticesInFrontPlaneQ:
                    debug("\nSimple Intersection?\n")
                    debug("Test 3bis or 4bis failed\n")
                    debug("Split here!!2\n")

                    facelist = HSR.facesplit(P, Q, facelist, nmesh)
                    split_done = 1
                    break

                facelist.remove(Q)
                facelist.insert(0, Q)
                Q.smooth = 1
                face_marked = 1
                debug("Q marked!\n")
                break

            # Write P!
            if split_done == 0 and face_marked == 0:
                facelist.remove(P)
                maplist.append(P)
                dumpfaces(maplist, "dump" + str(len(maplist)).zfill(4) + ".svg")

                progress.update()

            if len(facelist) == 870:
                dumpfaces([P, Q], "loopdebug.svg")

            #if facelist == None:
            #    maplist = [P, Q]
            #    print [v.co for v in P]
            #    print [v.co for v in Q]
            #    break

            # end of while len(facelist)

        nmesh.faces = maplist
        #for f in nmesh.faces:
        #    f.sel = 1

        nmesh.update()

    def _doHiddenSurfaceRemoval(self, mesh, ob, workScene):
        """Do HSR for the given mesh.
        """
        vrm = bpy.context.scene.VRM
        #mesh.update(calc_edges=False, calc_tessface=True)
        if len(mesh.polygons) == 0:
            return

        if vrm.polygonsHSR == 'PAINTER':
            print("\nUsing the Painter algorithm for HSR.")
            self.__simpleDepthSort(mesh, ob, workScene)

        elif vrm.polygonsHSR == 'NEWELL':
            print("\nUsing the Newell's algorithm for HSR.")
            self.__newellDepthSort(mesh)

    def _doEdgesStyle(self, mesh, edgestyleSelect):
        """Process Mesh Edges accroding to a given selection style.

        Examples of algorithms:

        Contours:
            given an edge if its adjacent faces have the same normal (that is
            they are complanar), than deselect it.

        Silhouettes:
            given an edge if one its adjacent faces is frontfacing and the
            other is backfacing, than select it, else deselect.
        """

        bpy.context.tool_settings.mesh_select_mode = (False, True, False) #Mesh.Mode(Mesh.SelectModes['EDGE'])

        edge_cache = MeshUtils.buildEdgeFaceUsersCache(mesh)

        for i, edge_faces in enumerate(edge_cache):
            mesh.edges[i].select = 0
            if edgestyleSelect(edge_faces):
                mesh.edges[i].select = 1

        """
        for e in mesh.edges:

            e.sel = 0
            if edgestyleSelect(e, mesh):
                e.sel = 1
        """
        #



# A wrapper function for the vectorizing process
def vectorize(filename):
    """The vectorizing process is as follows:

     - Instanciate the writer and the renderer
     - Render!
     """

    if filename == "":
        print("\nERROR: invalid file name!")
        return

    #from Blender import Window
    #editmode = Window.EditMode()
    #if editmode:
    #    Window.EditMode(0)
    vrm = bpy.context.scene.VRM
    
    if vrm.edgesSHOW:
      vrm.edgesSHOW = False

    actualWriter = outputWriters[vrm.outputFORMAT]
    writer = actualWriter(filename)

    renderer = Renderer()
    renderer.doRendering(writer, vrm.outputANIMATION)

    #if editmode:
    #    Window.EditMode(1)


class VRMPathShader(StrokeShader):
  def __init__(self, _thickness, _caps, _r, _g, _b, _a):
    StrokeShader.__init__(self)
    self.thickness = _thickness
    self.caps = _caps
    self.r = _r
    self.g = _g
    self.b = _b
    self.a = _a

  @classmethod
  def from_lineset(cls, lineset):
    linestyle = lineset.linestyle
    return cls(linestyle.thickness, linestyle.caps, linestyle.color[0], linestyle.color[1], linestyle.color[2], linestyle.alpha)
  
  def shade(self, stroke):
    global gOutput
    #print("VRMPathShader::shade")
    gOutput._printFreeStyleStroke(stroke, self.thickness, self.caps, self.r, self.g, self.b, self.a)


class ParameterEditorCallback(object):
    """Object to store callbacks for the Parameter Editor in"""
    def lineset_pre(self, scene, layer, lineset):
        raise NotImplementedError()

    def modifier_post(self, scene, layer, lineset):
        raise NotImplementedError()

    def lineset_post(self, scene, layer, lineset):
        raise NotImplementedError()

class VRMPathShaderCallback(ParameterEditorCallback):
    @classmethod
    def modifier_post(cls, scene, layer, lineset):
        global gOutput
        #print("VRMPathShaderCallback::modifier_post")
        
        vrm = bpy.context.scene.VRM
        if gOutput is None:
            return []
        if not (scene.render.use_freestyle and vrm.freestyleSHOW):
            return []

        cls.shader = VRMPathShader.from_lineset(lineset)
        return [cls.shader]


def register():
  bpy.utils.register_class(VRMRenderOp)
  bpy.utils.register_class(VRM)
  bpy.utils.register_class(VRMPanel)
  bpy.types.Scene.VRM = PointerProperty(type=VRM)
  
  parameter_editor.callbacks_modifiers_post.append(VRMPathShaderCallback.modifier_post)

def unregister():
  bpy.utils.unregister_class(VRMRenderOp)
  bpy.utils.unregister_class(VRMPanel)
  bpy.utils.unregister_class(VRM)
  del bpy.types.Scene.VRM
  
  parameter_editor.callbacks_modifiers_post.remove(VRMPathShaderCallback.modifier_post)

# Here the main
if __name__ == "__main__":
    register()