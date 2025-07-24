import bpy
from collections import defaultdict
from bpy.types import Operator
import bmesh
import mathutils
import math
from bpy.props import (
    BoolProperty,
    IntProperty,
    EnumProperty,
)


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
        layout.row()

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
        print(tris)
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        flattener = flatmesh.FaceUnwrapper(numpy.array(verts), numpy.array(tris))
        flattener.findFlatNodes(10, 0.95)
        fpts = [ (ze[0], ze[1], 0)  for ze in flattener.ze_nodes ]
        if math.isnan(fpts[0][0]) or math.isnan(fpts[0][1]) or math.isnan(fpts[0][2]):
            print("fluffnumber flattening failed", fluffnumber)

        mesh = bpy.data.meshes.new("flatmesh")
        #fpts = [ (0.,0.,0.),(1.,0.,0.),(0.,1.,0.) ]
        #tris = [ (0,1,2) ]
        mesh.from_pydata(fpts, [], tris)
        mobj = bpy.data.objects.new(mesh.name, mesh)
        mcoll = bpy.data.collections.new(mesh.name)
        bpy.context.scene.collection.children.link(mcoll)
        mcoll.objects.link(mobj)

        bpy.context.view_layer.objects.active = mobj

        return {'FINISHED'}


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
