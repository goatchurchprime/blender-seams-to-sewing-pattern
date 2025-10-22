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
    distort_colors_pertriangle: BoolProperty(
        name="Distort colors per triangle",
        description="Color triangles by distortion.",
        default=True,
    )
    distort_colors_peredge: BoolProperty(
        name="Distort colors per edge",
        description="Color triangles corner/sides distortion.",
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
        row.prop(self, "distort_colors_pertriangle")
        row = layout.row()
        row.prop(self, "distort_colors_peredge")
        row = layout.row()
        row.prop(self, "make_flatmesh")
        row = layout.row()
        row.prop(self, "stretch_range")


    def execute(self, context):
        if self.apply_modifiers:
            print("Do we need to copy before applying modifiers?")
            print(" otherwise we need to Ctrl-Z several times to return to original state")
            obj = bpy.context.active_object
            with bpy.context.temp_override(object=obj):
                mod_names = [ mod.name  for mod in obj.modifiers  if mod.show_viewport ]
                for mod_name in mod_names:
                    bpy.ops.object.modifier_apply(modifier=mod_name)


        # Separate out the faces by connectivity
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_mode(type='FACE')
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)
        #bm = bmesh.from_edit_mesh(bpy.context.edit_object.data)

        face_groups = []
        fg = { f for f in bm.faces if f.select }
        if not fg:
            faces = set(bm.faces[:])
            print("** no faces selected, separating all groups by seams nfaces=", len(faces))
            while faces:
                bpy.ops.mesh.select_all(action='DESELECT')
                face = faces.pop()
                face.select = True
                bpy.ops.mesh.select_linked()
                selected_faces = {f for f in faces if f.select}
                selected_faces.add(face) # this or bm.faces above?
                print("facegroup", len(face_groups), " nfaces=", len(selected_faces))
                face_groups.append(selected_faces)
                faces -= selected_faces
        else:
            face_groups.append(fg)
        
        # run the flattener on each face group
        flatteneddata = [ ]
        bmverts = bm.verts[:]
        seam_edges = list(filter(lambda e: e.seam, bm.edges))  # from which internal seams are derived
        for i, fg in enumerate(face_groups):
            g = sorted(set().union(*((v.index  for v in f.verts) for f in fg)))
            verts = [ tuple(bmverts[j].co)  for j in g ]
            gmap = dict(zip(g, range(len(g))))

            print(i, "faces ", len(fg), len(verts))
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
            flatteneddata.append([ list(fg), npverts, npfpts, nptris ])
            
        # colour the facegroups
        if self.distort_colors_pertriangle:
            if len(bm.loops.layers.color) > 0:
                color_layer = bm.loops.layers.color[obj.data.vertex_colors.active_index]
            else:
                color_layer = bm.loops.layers.color.new("CCC")
            obj.data.vertex_colors.active_index = obj.data.vertex_colors.find("CCC")
            stretchrangeedge = self.stretch_range
            stretchrangetriang = self.stretch_range*0.25

            for i in range(len(flatteneddata)):
                fg, npverts, npfpts, nptris = flatteneddata[i]
                for j in range(len(fg)):
                    face = fg[j]
                    tri = nptris[j]
                    fv = trianglearea(npfpts[tri[0]], npfpts[tri[1]], npfpts[tri[2]])
                    av = trianglearea(npverts[tri[0]], npverts[tri[1]], npverts[tri[2]])
                    ft = fv/av - 1.0
                    c = stretchtocolour(ft, stretchrangetriang)
                    for k, loop in enumerate(face.loops):
                        loop[color_layer] = c

        # generate the flattened meshes from the flattened data
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
                cfpts, cpolys, cpolysvcolours = trianglestocolourednpolys(npfpts, nptris, npverts)

                mesh = bpy.data.meshes.new("flatmesh%d" % len(fg))
                #mesh.from_pydata(list(npfpts), [], list(nptris))
                mesh.from_pydata(list(cfpts), [], list(cpolys))
                mobj = bpy.data.objects.new(mesh.name, mesh)
                mcoll.objects.link(mobj)

                # apply flattened colours to the new mesh
                if self.distort_colors_peredge:
                    bpy.context.view_layer.objects.active = mobj
                    bpy.ops.object.mode_set(mode='EDIT')
                    print("ggg", bpy.context.edit_object)
                    fbm = bmesh.from_edit_mesh(bpy.context.edit_object.data)
                    fcolor_layer = fbm.loops.layers.color.new("CCD")
                    mobj.data.vertex_colors.active_index = mobj.data.vertex_colors.find("CCD")
                    print("fcolor_layer", fcolor_layer)
                    for j, face in enumerate(fbm.faces):
                        for k, loop in enumerate(face.loops):
                            loop[fcolor_layer] = stretchtocolour(cpolysvcolours[j][k], stretchrangeedge)
                    bpy.ops.object.mode_set(mode='OBJECT')
                    bpy.context.view_layer.objects.active = None

        return {'FINISHED'}


def trianglearea(p0, p1, p2):
    a = p1 - p0
    b = p2 - p0
    c = numpy.array((a[1]*b[2] - b[1]*a[2], -a[0]*b[2] + b[0]*a[2], a[0]*b[1] - b[0]*a[1]))
    return math.sqrt(sum(c*c))/2

def fedgelen(p0, p1):
    a = p1 - p0
    return math.sqrt(sum(a*a))/2

def stretchtocolour(d, dstretchrange):
    l = min(1.0, abs(d)/dstretchrange)
    return col1*(1-l) + (col0 if d < 0 else col2)*l

col1 = numpy.array((0.9,0.9,0.9,1.0))
col2 = numpy.array((0.0,0.0,1.0,1.0))
col0 = numpy.array((1.0,0.0,0.0,1.0))

def edgesplitverts(i0, i1, cfpts, edgesplitmap):
    edgemidpoints = edgesplitmap.pop((i1, i0), None)
    if edgemidpoints is None:
        midpoint0 = cfpts[i0]*0.875 + cfpts[i1]*0.125
        midpoint1 = cfpts[i0]*0.125 + cfpts[i1]*0.875
        edgemidpoints = [ len(cfpts), len(cfpts)+1 ]
        cfpts.append(midpoint0)
        cfpts.append(midpoint1)
        edgesplitmap[(i0, i1)] = edgemidpoints
    else:
        edgemidpoints.reverse()
    return edgemidpoints
    
def trianglestocolourednpolys(fpts, tris, pts):
    cfpts = list(fpts)
    cpolys = [ ]
    cpolysvcolours = [ ]
    edgesplitmap = { }
    for tri in tris:
        cpolyvvalues = [ ]
        cpolyvcolours = [ ]
        fe01 = fedgelen(fpts[tri[0]], fpts[tri[1]])/fedgelen(pts[tri[0]], pts[tri[1]]) - 1.0
        fe12 = fedgelen(fpts[tri[1]], fpts[tri[2]])/fedgelen(pts[tri[1]], pts[tri[2]]) - 1.0
        fe20 = fedgelen(fpts[tri[2]], fpts[tri[0]])/fedgelen(pts[tri[2]], pts[tri[0]]) - 1.0

        cpoly = [ tri[0] ]
        cpolyvcolours = [ 0.0*(fe01 + fe20 - fe12) ]
        cpoly.extend(edgesplitverts(tri[0], tri[1], cfpts, edgesplitmap))
        cpolyvcolours.extend([ fe01 ]*(len(cpoly) - len(cpolyvcolours)))
        cpoly.append(tri[1])
        cpolyvcolours.append(0.0*(fe12 + fe01 - fe20))
        cpoly.extend(edgesplitverts(tri[1], tri[2], cfpts, edgesplitmap))
        cpolyvcolours.extend([ fe12 ]*(len(cpoly) - len(cpolyvcolours)))
        cpoly.append(tri[2])
        cpolyvcolours.append(0.0*(fe20 + fe12 - fe01))
        cpoly.extend(edgesplitverts(tri[2], tri[0], cfpts, edgesplitmap))
        cpolyvcolours.extend([ fe20 ]*(len(cpoly) - len(cpolyvcolours)))

        cpolys.append(cpoly)
        cpolysvcolours.append(cpolyvcolours)
    
    return cfpts, cpolys, cpolysvcolours

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
                if e.seam or len(e.link_faces) == 1:
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
