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

# Now no longer necessary using the blender made with nix develop in magic blender file
#import sys
#sys.path.append("/nix/store/izs69w0zy2wimfkw6yfrrkazra631lid-freecad-1.0.0/lib")
#sys.path.append("/nix/store/4chbw98xmm4bglfl8f8fmij17ajvf4yi-python3.11-numpy-2.2.4/lib/python3.11/site-packages/")

import numpy
import flatmesh
import math

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
        default = 0.01 
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
            print("Do we need to copy before applying modifiers?")
            obj = bpy.context.active_object
            with bpy.context.temp_override(object=obj):
                mod_names = [mod.name for mod in obj.modifiers]
                for mod_name in mod_names:
                    bpy.ops.object.modifier_apply(modifier=mod_name)

        bpy.ops.object.mode_set(mode='EDIT')
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)

        bmverts = bm.verts[:]
        faces = set(bm.faces[:])

        fg = { f for f in bm.faces if f.select }
        face_groups = []
        if not fg:
            print("** no faces selected, separating all groups by seams")
            faces = set(bm.faces[:])
            seam_edges = list(filter(lambda e: e.seam, bm.edges))
            while faces:
                bpy.ops.mesh.select_all(action='DESELECT')
                face = faces.pop()
                face.select = True
                bpy.ops.mesh.select_linked()
                selected_faces = {f for f in faces if f.select}
                selected_faces.add(face) # this or bm.faces above?
                face_groups.append(selected_faces)
                faces -= selected_faces
        else:
            face_groups.append(fg)
        
        flatteneddata = [ ]
        for i, fg in enumerate(face_groups):
            g = sorted(set().union(*((v.index  for v in f.verts) for f in fg)))
            verts = [ tuple(bmverts[j].co)  for j in g ]
            gmap = dict(zip(g, range(len(g))))

            print(i, "faces ", len(fg), len(verts))
            seam_edges = list(filter(lambda e: e.seam, bm.edges))
            facverts = facevertstodoubleup(seam_edges, fg)
            tris = trilist(gmap, verts, facverts, fg)

            npverts = numpy.array(verts)
            nptris = numpy.array(tris)
            flattener = flatmesh.FaceUnwrapper(npverts, nptris)
            flattener.findFlatNodes(10, 0.95)
            fpts = [ (ze[0], ze[1], i*0.01)  for ze in flattener.ze_nodes ]
            npfpts = numpy.array(fpts)
            if math.isnan(fpts[0][0]) or math.isnan(fpts[0][1]) or math.isnan(fpts[0][2]):
                print("fluffnumber flattening failed")
                continue
            flatteneddata.append([fg, npverts, npfpts, nptris])
            
        if self.distort_colors:
            if len(bm.loops.layers.color) > 0:
                color_layer = bm.loops.layers.color[obj.data.vertex_colors.active_index]
            else:
                color_layer = bm.loops.layers.color.new("CCC")
            stretchrange = self.stretch_range
            triangdistortcols = [ ]
            for i in range(len(flatteneddata)):
                triangdistortcols.append([])
                fg, npverts, npfpts, nptris = flatteneddata[i]
                fvals = triangledistortions(npverts, npfpts, nptris)
                print("min max distort", min(fvals), max(fvals))
                for face, d in zip(fg, fvals):
                    l = min(1.0, abs(d)/stretchrange)
                    c = col1*(1-l) + (col0 if d < 0 else col2)*l
                    triangdistortcols[-1].append(c)
                    for loop in face.loops:
                        loop[color_layer] = c
                #bm.to_mesh(obj.data)

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.context.view_layer.objects.active = None
            
        if self.make_flatmesh:
            fi = bpy.data.collections.find("flatmeshes")
            if fi == -1:
                mcoll = bpy.data.collections.new("flatmeshes")
                bpy.context.scene.collection.children.link(mcoll)
            else:
                mcoll = bpy.data.collections[fi]
                [ mcoll.objects.unlink(o)  for o in mcoll.objects ]
            lmcoll = bpy.context.view_layer.layer_collection.children[mcoll.name]
            lmcoll.exclude = False
            for i in range(len(flatteneddata)):
                fg, npverts, npfpts, nptris = flatteneddata[i]
                mesh = bpy.data.meshes.new("flatmesh%d" % len(fg))
                mesh.from_pydata(list(npfpts), [], list(nptris))
                mobj = bpy.data.objects.new(mesh.name, mesh)
                mcoll.objects.link(mobj)

                print("mmm ", mobj, list(mcoll.objects))
                bpy.context.view_layer.objects.active = mobj
                #mobj.select_set(True)
                bpy.ops.object.mode_set(mode='EDIT')
                if self.distort_colors:
                    print("ggg", bpy.context.edit_object)
                    fbm = bmesh.from_edit_mesh(bpy.context.edit_object.data)
                    fcolor_layer = fbm.loops.layers.color.new("CCD")
                    print("fcolor_layer", fcolor_layer)
                    for face, c in zip(fbm.faces, triangdistortcols[i]):
                        for loop in face.loops:
                            loop[fcolor_layer] = c
                bpy.ops.object.mode_set(mode='OBJECT')
                #mobj.select_set(False)
                bpy.context.view_layer.objects.active = None

        return {'FINISHED'}


def trianglearea(p0, p1, p2):
    a = p1 - p0
    b = p2 - p0
    c = numpy.array((a[1]*b[2] - b[1]*a[2], -a[0]*b[2] + b[0]*a[2], a[0]*b[1] - b[0]*a[1]))
    return math.sqrt(sum(c*c))/2

col0 = numpy.array((0.0,0.0,1.0,1.0))
col2 = numpy.array((1.0,0.0,0.0,1.0))
col1 = numpy.array((0.9,0.9,0.9,1.0))

def triangledistortions(verts, fpts, tris):
    fvals = [ ]
    for tri in tris:
        av = trianglearea(verts[tri[0]], verts[tri[1]], verts[tri[2]])
        fv = trianglearea(fpts[tri[0]], fpts[tri[1]], fpts[tri[2]])
        fvals.append(fv/av - 1.0)
    return fvals

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
    #print("internal seams vertex pairs", [[e.verts[0].index, e.verts[1].index]  for e in internal_seams])
    while internal_seams:
        lfaces = [ ]
        lverts = set()
        e0 = internal_seams.pop()
        f0 = e0.link_faces[0]
        lfaces.append(f0)
        #print("beginning face ", f0, e0)
        for e, f, v in [ (e0, f0, e0.verts[0]), (e0, f0, e0.verts[1]) ]:
            es = e
            for i in range(1000):
                #print(i, e.seam, e, f, v)
                e = next(ef  for ef in f.edges  if ef != e and v in ef.verts)
                if e == es:
                    #print("wrap endpoint")
                    lfaces.pop()
                    break  # wrap end point
                if e.seam:
                    lverts.add(v)
                    if not (e in internal_seams):
                        #print("hit boundary")
                        break # hit boundary
                    es = e
                    internal_seams.remove(es)
                    v = next(vf  for vf in e.verts  if vf != v)
                else:
                    try:
                        f = next(ff  for ff in e.link_faces  if ff != f)
                        lfaces.append(f)
                    except StopIteration:
                        pass
        #print([v.index  for v in lverts])
        #for f in lfaces:
        #    print("  ", [v.index for v in f.verts])
        res.append([lfaces, lverts])
    return res
