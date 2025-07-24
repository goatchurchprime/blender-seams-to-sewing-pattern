import bpy
from collections import defaultdict
from bpy.types import Operator
import bmesh
import mathutils
import math
from bpy.props import (
    BoolProperty,
    FloatProperty,
    IntProperty,
    EnumProperty,
)

# needs to select the area to flatten (triangle select and L)
# then go back to object mode, select this option, turn off generate flat mesh
# run, then edit mode, attributes panel, see CCC, and viewport shading color by attribute CCC

import sys
sys.path.append("/nix/store/izs69w0zy2wimfkw6yfrrkazra631lid-freecad-1.0.0/lib")
sys.path.append("/nix/store/4chbw98xmm4bglfl8f8fmij17ajvf4yi-python3.11-numpy-2.2.4/lib/python3.11/site-packages/")
import numpy, flatmesh, math

class Freecad_flatten_component(Operator):
    bl_idname = "object.freecad_flatten_component"
    bl_label = "Freecad flatten"
    bl_description = (
        "Flattens a seamed component using FreeCAD flattener"
    )
    bl_options = {'REGISTER', 'UNDO'}

    apply_modifiers: BoolProperty(
        name="Apply modifiers",
        description="Applies all modifiers before operating.",
        default=True,
    )
    distort_colors: BoolProperty(
        name="Distort colors",
        description="Color triangles by distortion.",
        default=True,
    )
    make_flatmesh: BoolProperty(
        name="Make flat mesh",
        default=True,
    )
    stretch_range: FloatProperty(
        name="Stretch range",
        min=0.000001,
        default = 0.1 
    )

    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self, width=250)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Unfolds this mesh by cutting along seams.", icon='INFO')
        layout.separator()
        layout.row()

        layout.row()
        row = layout.row()
        row.prop(self, "apply_modifiers")
        row = layout.row()
        row.prop(self, "distort_colors")
        row = layout.row()
        row.prop(self, "make_flatmesh")
        row = layout.row()
        row.prop(self, "stretch_range")


    def execute(self, context):
        if self.apply_modifiers:
            bpy.ops.object.convert(target='MESH')
            obj = bpy.context.active_object
        bpy.ops.object.mode_set(mode='EDIT')
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)

        fg = { f for f in bm.faces if f.select }
        if len(fg) == 0:
            print("** no faces selected, do this and hit L for the component up to the seam before converting to Object mode")
        bmverts = bm.verts[:]
        g = sorted(set().union(*((v.index  for v in f.verts) for f in fg)))
        verts = [ tuple(bmverts[i].co)  for i in g ]
        gmap = dict(zip(g, range(len(g))))

        print("faces ", len(bm.faces), len(verts), len(fg))
        seam_edges = list(filter(lambda e: e.seam, bm.edges))
        facverts = facevertstodoubleup(seam_edges, fg)
        tris = trilist(gmap, verts, facverts, fg)

        flattener = flatmesh.FaceUnwrapper(numpy.array(verts), numpy.array(tris))
        flattener.findFlatNodes(10, 0.95)
        fpts = [ (ze[0], ze[1], 0)  for ze in flattener.ze_nodes ]
        if math.isnan(fpts[0][0]) or math.isnan(fpts[0][1]) or math.isnan(fpts[0][2]):
            print("fluffnumber flattening failed", fluffnumber)

        stretchrange =  self.stretch_range
        print("stretchrange ", stretchrange)
        if self.distort_colors:
            if len(bm.loops.layers.color) > 0:
                color_layer = bm.loops.layers.color[obj.data.vertex_colors.active_index]
            else:
                color_layer = bm.loops.layers.color.new("CCC")
            fvals = triangledistortions(numpy.array(verts), numpy.array(fpts), numpy.array(tris))
            print("min max distort", min(fvals), max(fvals))
            for face, d in zip(fg, fvals):
                l = min(1.0, abs(d)/stretchrange)
                c = col1*(1-l) + (col0 if d < 0 else col2)*l
                for loop in face.loops:
                    loop[color_layer] = c
            #bm.to_mesh(obj.data)

        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        if self.make_flatmesh:
            mesh = bpy.data.meshes.new("flatmesh")
            mesh.from_pydata(fpts, [], tris)
            mobj = bpy.data.objects.new(mesh.name, mesh)
            mcoll = bpy.data.collections.new(mesh.name)
            bpy.context.scene.collection.children.link(mcoll)
            mcoll.objects.link(mobj)
            bpy.context.view_layer.objects.active = mobj

        return {'FINISHED'}

def trianglearea(p0, p1, p2):
    a = p1 - p0
    b = p2 - p0
    c = numpy.array((a[1]*b[2] - b[1]*a[2], -a[0]*b[2] + b[0]*a[2], a[0]*b[1] - b[0]*a[1]))
    return math.sqrt(sum(c*c))/2

col0 = numpy.array((1.0,0.0,0.0,1.0))
col1 = numpy.array((0.9,0.9,0.9,1.0))
col2 = numpy.array((0.0,0.0,1.0,1.0))


def triangledistortions(verts, fpts, tris):
    fvals = [ ]
    for tri in tris:
        av = trianglearea(verts[tri[0]], verts[tri[1]], verts[tri[2]])
        fv = trianglearea(fpts[tri[0]], fpts[tri[1]], fpts[tri[2]])
        fvals.append(fv/av - 1.0)
    return fvals

"""
import bpy, bmesh

mesh = bpy.context.object.data
bm = bmesh.new()
bm.from_mesh(mesh)
color_layer = bm.loops.layers.color.new("CCC")
for vert in bm.verts:
    for loop in vert.link_loops:
        loop[color_layer] = tuple(vert.co)+(1,)
        
bm.to_mesh(mesh)


mesh = bpy.context.object.data
bm = bmesh.new()
bm.from_mesh(mesh)
color_layer = bm.loops.layers.color[mesh.vertex_colors.active_index]
for vert in bm.verts[:10]:
    for loop in vert.link_loops:
        loop[color_layer] = (0,0,1,1)

for face in bm.faces[:10]:
    for loop in face.loops:
        loop[color_layer] = (0,1,0,1)
        
bm.to_mesh(mesh)

colv0 = -0.1
colv2 = 0.1
colv1 = (colv0 + colv2)/2
col0 = P3(1,0,0)
col1 = P3(0.9,0.9,0.9)
col2 = P3(0,0,1)
baspectratio = True
def convcl(v, v0, v1, cl0, cl1):
    lam = max(0, min(1, (v - v0)/(v1 - v0)))
    c = cl0*(255*(1-lam)) + cl1*(255*lam)
    c = cl0*(255*(1-lam)) + cl1*(255*lam)
    return (int(c[0])<<24) + (int(c[1])<<16) + (int(c[2])<<8) + 255
def convcol(v):
    return convcl(v, colv0, colv1, col0, col1) if v < colv1 else convcl(v, colv1, colv2, col1, col2)
def colorbydistortion(mtriangulation, mflattened):
    facets1 = mtriangulation.Mesh.Facets
    facets2 = mflattened.Mesh.Facets
    assert(len(facets1) == len(facets2))
    FaceColors = [ ]
    cs = [ ]
    for f1, f2 in zip(facets1, facets2):
        if baspectratio:
            c = f1.AspectRatio/f2.AspectRatio - 1.0
        else:
            c = f1.Area/f2.Area - 1.0
        cs.append(c)
        FaceColors.append(convcol(c))


"""

def trilist(gmap, verts, facverts, fg):
    for i in range(len(facverts)):
        gverts = { }
        for v in facverts[i][1]:
            vi = gmap[v.index]
            gverts[vi] = len(verts)
            verts.append(verts[vi])   # or add in tinydisplacement here, to protect remeshing closing it up if we do this
        facverts[i][1] = gverts
    
    tris = [ ]
    for f in fg:
        tri = [ ]
        for v in f.verts:
            k = gmap[v.index]
            for lfaces, gverts in facverts:
                if f in lfaces and k in gverts:
                    k = gverts[k]
            tri.append(k)
        tris.append(tri)
    return tris

def facevertstodoubleup(seam_edges, fg):
    res = [ ]
    internal_seams = set(e  for e in seam_edges  if sum(int(f in fg)  for f in e.link_faces) == 2)
    print([[e.verts[0].index, e.verts[1].index]  for e in internal_seams])
    while internal_seams:
        lfaces = [ ]
        lverts = set()
        e0 = internal_seams.pop()
        f0 = e0.link_faces[0]
        lfaces.append(f0)
        print("beginning face ", f0, e0)
        for e, f, v in [ (e0, f0, e0.verts[0]), (e0, f0, e0.verts[1]) ]:
            es = e
            for i in range(1000):
                print(i, e.seam, e, f, v)
                e = next(ef  for ef in f.edges  if ef != e and v in ef.verts)
                if e == es:
                    print("wrap endpoint")
                    lfaces.pop()
                    break  # wrap end point
                if e.seam:
                    lverts.add(v)
                    if not (e in internal_seams):
                        print("hit boundary")
                        break # hit boundary
                    es = e
                    internal_seams.remove(es)
                    v = next(vf  for vf in e.verts  if vf != v)
                else:
                    f = next(ff  for ff in e.link_faces  if ff != f)
                    lfaces.append(f)
        print([v.index  for v in lverts])
        for f in lfaces:
            print("  ", [v.index for v in f.verts])
        res.append([lfaces, lverts])
    return res
