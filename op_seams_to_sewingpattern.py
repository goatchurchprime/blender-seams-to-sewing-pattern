import bpy
from collections import defaultdict
from bpy.types import Operator
import bmesh
import mathutils
import math
from bpy.props import (
    BoolProperty,
    IntProperty,
    FloatProperty,
    EnumProperty,
)

if bpy.app.version >= (3, 0, 0):
    from . import function_wrapper_3_0 as function_wrapper
elif bpy.app.version >= (2, 90, 0):
    from . import function_wrapper_2_9 as function_wrapper
else:
    from . import function_wrapper_2_8 as function_wrapper


"""
import bmesh
me = bpy.context.edit_object.data
bm = bmesh.from_edit_mesh(me)
bpy.ops.uv.unwrap(method='CONFORMAL', margin=0.02)
Anvbyss
art by ashley
Fid_
"""

import numpy
import flatmesh

def freecadflatten(bm, me):
    uv_layer = bm.loops.layers.uv.active
    bpy.ops.uv.select(deselect_all=True)
    bpy.ops.mesh.select_mode(type="FACE")
    bpy.ops.mesh.select_all(action='SELECT')
    faces = set(bm.faces[:])
    faceGroups = [ ]
    while faces:
        bpy.ops.uv.select(deselect_all=True)
        eselface = [l.face  for f in faces  for l in f.loops  if l[uv_layer].select]  # should be empty from function above!
        if not eselface:
            face = faces.pop()
            face.loops[0][uv_layer].select = True
        else:
            face = eselface[0]
        bpy.ops.uv.select_linked()
        selected_faces = set(l.face  for f in faces  for l in f.loops  if l[uv_layer].select)
        selected_faces.add(face)
        faceGroups.append(selected_faces)
        print("sss ", len(selected_faces), len(faces))
        faces -= selected_faces
    
    for fg in faceGroups:
        flattenfacegroup(fg, uv_layer)
    bpy.ops.mesh.select_all(action='DESELECT')
    
# Now in each facegroup the loops are the unique face-vertex pairs
# We need to find which ones are the same by vertex index when the edge 
# between them is not a seam

def flattenfacegroup(fg, uv_layer):
    print("fcf fg", len(fg))
    # each contained non-seam edge equates 2 corresponding loops (face-vertex pairs)
    jedges = set()
    for f in fg:
        for e in f.edges:
            if len(e.link_faces) == 2:
                if e.link_faces[0] in fg and e.link_faces[1] in fg:
                    if not e.seam:
                        jedges.add(e)

    # derive the equivalence classes from the equivalence relation
    allloops = set().union(*(f.loops for f in fg))
    loopconns = dict((l, (l,))  for l in allloops)
    for e in jedges:
        tloops = list(e.link_faces[0].loops) + list(e.link_faces[1].loops)
        tloops.sort(key=lambda l: l.vert.index)
        assert (len(tloops) == 6)
        for i in range(5):
            if tloops[i].vert == tloops[i+1].vert:
                l1, l2 = tloops[i], tloops[i+1]
                if l1 not in loopconns[l2]:
                    comb = loopconns[l1] + loopconns[l2]
                    for l in comb:
                        loopconns[l] = comb

    # generate the points from equivalence classes and map the triangles 
    looppoints = list(set(loopconns.values()))
    npverts = numpy.array([ looppoint[0].vert.co  for looppoint in looppoints ])
    connlookup = dict([ (looppoint, i)  for i, looppoint in enumerate(looppoints) ])
    nptris = numpy.array( [ [ connlookup[loopconns[l]] for l in f.loops ] for f in fg ] )

    # flatten and map back to the uv values
    flattener = flatmesh.FaceUnwrapper(npverts, nptris)
    flattener.findFlatNodes(10, 0.95)
    fpts = [ mathutils.Vector((ze[0], ze[1]))  for ze in flattener.ze_nodes ]
    for f in fg:
        for l in f.loops:
            l[uv_layer].uv = fpts[connlookup[loopconns[l]]]*0.5 + mathutils.Vector((0.5,0.5))


class Seams_To_SewingPattern(Operator):
    bl_idname = "object.seams_to_sewingpattern"
    bl_label = "Seams to Sewing Pattern"
    bl_description = (
        "Converts a manifold mesh with seams into a swewing pattern for cloth"
        " simulation"
    )
    bl_options = {'REGISTER', 'UNDO'}

    do_unwrap: EnumProperty(
        name="Unwrap",
        description=(
            "Perform an unwrap before unfolding. Identical to UV > Unwrap"
        ),
        items=(
            ('ANGLE_BASED', "Angle based", ""),
            ('CONFORMAL', "Conformal", ""),
            ('KEEP', "Keep existing (advanced)", ""),
        ),
        default='ANGLE_BASED',
    )
    keep_original: BoolProperty(
        name="Work on duplicate",
        description=(
            "Creates a duplicate of the selected object and operates on that"
            " instead. This keeps your original object intact."
        ),
        default=True,
    )
    
    bevel_skim: BoolProperty(
        name="BevelSkim",
        description="Bevel seams to double and create edge line gap between",
        default=True,
    )
    bevel_offset: FloatProperty(
        name="Bevel offset",
        description="For the gaps applied to the seams before projection",
        default=0.0002,
    )
    dont_flatten: BoolProperty(
        name="Dontflatten",
        description="Don't actually flatten the panels of the new mesh",
        default=False,
    )
    dont_project: BoolProperty(
        name="Dontproject",
        description="Don't project the flattened UVs back to the 3D surfaces",
        default=False,
    )

    freecad_flattener: BoolProperty(
        name="FreecadFlattener",
        description="Flatten using the FreeCAD flattener",
        default=True,
    )

    use_remesh_seams: BoolProperty(
        name="Remesh seams",
        description="Subdivide along the seams",
        default=True,
    )
    use_remesh: BoolProperty(
        name="Remesh",
        description="Use Boundary Aligned Remesh to remesh",
        default=True,
    )
    apply_modifiers: BoolProperty(
        name="Apply modifiers",
        description="Applies all modifiers before operating.",
        default=True,
    )
    target_tris: IntProperty(
        name="Target number of triangles",
        description="Actual number of triangle migh be a bit off",
        default=5000,
    )

    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self, width=250)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(
            text="Unfolds this mesh by cutting along seams.", icon='INFO'
        )
        layout.separator()
        layout.row()
        layout.row()
        row = layout.row()
        row.prop(self, "do_unwrap")
        if(self.do_unwrap == 'KEEP'):
            row = layout.row()
            row.alignment = 'EXPAND'
            row.label(text="Ensure your seams match your UV's!", icon='EDGESEL')

        layout.row()
        row = layout.row()
        row.prop(self, "keep_original")
        row = layout.row()
        row.prop(self, "apply_modifiers")
        row = layout.row()
        row.prop(self, "bevel_skim")
        row = layout.row()
        row.prop(self, "bevel_offset")
        row = layout.row()
        row.prop(self, "dont_flatten")
        row = layout.row()
        row.prop(self, "dont_project")
        row = layout.row()
        row.prop(self, "freecad_flattener")
        row = layout.row()
        row.prop(self, "use_remesh_seams")
        row = layout.row()
        row.prop(self, "use_remesh")
        row = layout.row()
        row.prop(self, "target_tris")
        row.enabled = self.use_remesh
        layout.row()

    def execute(self, context):
        if self.keep_original:
            # Duplicate selection to keep original.
            src_obj = bpy.context.active_object
            obj = src_obj.copy()
            obj.data = src_obj.data.copy()
            obj.animation_data_clear()
            bpy.context.collection.objects.link(obj)

            obj.select_set(True)
            src_obj.select_set(False)
            bpy.context.view_layer.objects.active = obj
            bpy.ops.mesh.customdata_custom_splitnormals_clear()
            bpy.context.object.data.materials.clear()  # try this next

        if self.apply_modifiers:
            if self.keep_original:
                with bpy.context.temp_override(object=obj):
                    mod_names = [mod.name for mod in obj.modifiers]
                    for mod_name in mod_names:
                        bpy.ops.object.modifier_apply(modifier=mod_name)
            else:
                print("no modifiers applied if not keeping original")

        wm = bpy.context.window_manager
        bpy.ops.object.mode_set(mode='EDIT')

        obj = bpy.context.edit_object
        me = obj.data

        bpy.ops.mesh.select_mode(type="EDGE")

        bpy.ops.mesh.select_all(action='SELECT')
        if (self.do_unwrap != 'KEEP'):
            bpy.ops.uv.unwrap(method=self.do_unwrap, margin=0.02)
            
            
        bpy.ops.mesh.select_all(action='DESELECT')

        bm = bmesh.from_edit_mesh(me)

        obj["S2S_InitialVolume"] = bm.calc_volume()

        function_wrapper.do_update_edit_mesh(me)

        if not self.dont_flatten and self.freecad_flattener:
            freecadflatten(bm, me)
            #function_wrapper.do_update_edit_mesh(me)


        # Calculate edge length based on a surface of equilateral triangles.
        if (self.use_remesh_seams):
            current_area = sum(f.calc_area() for f in bm.faces)
            target_triangle_count = self.target_tris
            area_per_triangle = current_area / target_triangle_count

            max_edge_length = math.sqrt(area_per_triangle/(math.sqrt(3)/4))

            # A bias to compensate for stretching.
            self.ensure_edgelength(max_edge_length * 0.8, bm, wm)


        warn_any_seam = False
        for e in bm.edges:
            if e.seam:
                e.select = True
                warn_any_seam = True

        if not warn_any_seam:
            self.report(
                {'ERROR'},
                (
                    'There are no seams in this mesh. Please add seams where'
                    ' you want to cut the model.'
                )
            )
            return {'CANCELLED'}

        # the bevel thing makes those empty quads so corresponding vertexes can be connected between the faces
        if self.bevel_skim:
            function_wrapper.do_bevel(self.bevel_offset)
            # fix fanning seams
            degenerate_edges = list()
            for f in list(filter(lambda f: (f.select), bm.faces)):
                is_degenerate = False
                for v in f.verts:
                    vert_degenerate = True
                    for e in v.link_edges:
                        if e.seam:
                            vert_degenerate = False
                    if vert_degenerate:
                        is_degenerate = True

                for e in f.edges:
                    if e.is_boundary:
                        is_degenerate = False

                if is_degenerate:
                    for e in f.edges:
                        degenerate_edges.append(e)

            print("degenerate edges to collapse", degenerate_edges)
            bmesh.ops.collapse(bm, edges=list(set(degenerate_edges)), uvs=True)


            print("quitting early")
            return {'FINISHED'}

            bpy.ops.mesh.delete(type='ONLY_FACE')


        bpy.ops.mesh.select_mode(type="FACE")
        faceGroups = []

        # isolate all face islands, and UV unwrap each island

        faces = set(bm.faces[:])
        wm.progress_begin(0, 99)
        progress_max = len(faces)
        progress = 0
        while faces:
            bpy.ops.mesh.select_all(action='DESELECT')
            face = faces.pop()
            face.select = True
            bpy.ops.mesh.select_linked()
            selected_faces = {f for f in faces if f.select}
            selected_faces.add(face)  # this or bm.faces above?
            faceGroups.append(selected_faces)
            faces -= selected_faces

            progress += len(selected_faces)
            wm.progress_update((progress / progress_max))

        print("We have found ", len(faceGroups), " facegroups.")

        uv_layer = bm.loops.layers.uv.active

        progress = 0
        area_before = 0
        area_after = 0

        function_wrapper.do_update_edit_mesh(me)
        for g in faceGroups:
            if self.dont_flatten or self.dont_project:
                break
            progress += 1
            wm.progress_update((progress / len(faceGroups)))
            bpy.ops.mesh.select_mode(type='FACE')
            bpy.ops.mesh.select_all(action='DESELECT')
            average_position = mathutils.Vector((0, 0, 0))
            facenum = 0

            # calculate the area, average position

            for f in g:
                f.select = True
                area_before += f.calc_area()
                average_position += f.calc_center_median()
                facenum += 1

            average_position /= facenum

            average_tangent = mathutils.Vector((0, 0, 0))
            average_bitangent = mathutils.Vector((0, 0, 0))

            # calculate a rough tangent and a bitangent

            average_uv_position = mathutils.Vector((0, 0))
            uv_position_samples = 0

            for face in g:
                for loop in face.loops:
                    uv = loop[uv_layer].uv
                    uv_position_samples += 1
                    average_uv_position += uv
                    delta = loop.vert.co - average_position
                    average_tangent += delta * (uv.x - 0.5)
                    average_bitangent += delta * (uv.y - 0.5)

            # reorient the tangent and bitangent

            average_uv_position /= uv_position_samples
            average_tangent = average_tangent.normalized()
            average_bitangent = average_bitangent.normalized()
            average_normal = average_tangent.cross(
                average_bitangent
            ).normalized()
            halfvector = average_bitangent + average_tangent
            halfvector /= 2
            halfvector.normalize()
            # straighten out half vector
            halfvector = average_normal.cross(halfvector)
            halfvector = average_normal.cross(halfvector)
            cw = mathutils.Matrix.Rotation(math.radians(45.0), 4, average_normal)
            ccw = mathutils.Matrix.Rotation(math.radians(-45.0), 4, average_normal)

            average_tangent = mathutils.Vector(halfvector)
            average_tangent.rotate(ccw)

            average_bitangent = mathutils.Vector(halfvector)
            average_bitangent.rotate(cw)
            print("facegroup UV projected tangents ", average_tangent, average_bitangent, average_tangent.length, average_bitangent.length)

            # offset each face island by their UV value, using the tangent and
            # bitangent to recreate the flat shape defined by the UVs in space

            for face in g:
                for loop in face.loops:
                    uv = loop[uv_layer].uv
                    vert = loop.vert
                    pos = mathutils.Vector((0, 0, 0))
                    pos += average_position
                    pos += average_tangent * -(uv.x - average_uv_position.x)
                    pos += average_bitangent * -(uv.y - average_uv_position.y)
                    # arbitrary - should probably depend on object scale?
                    pos += average_normal * 0.3
                    
                    # commenting out this avoids projecting it out so we can find what the processed subdivided surface looks like 
                    vert.co = pos
                    vert.normal = average_normal

            function_wrapper.do_update_edit_mesh(me)
            area_after += sum(f.calc_area() for f in g)

        # done
        if not (self.dont_flatten or self.dont_project):
            area_ratio = math.sqrt(area_before / area_after)
            bpy.ops.mesh.select_all(action='SELECT')
            previous_pivot = bpy.context.scene.tool_settings.transform_pivot_point
            bpy.context.scene.tool_settings.transform_pivot_point = ('INDIVIDUAL_ORIGINS')
            bpy.ops.transform.resize(value=(area_ratio, area_ratio, area_ratio))
            bpy.context.scene.tool_settings.transform_pivot_point = previous_pivot
            obj["S2S_UVtoWORLDscale"] = area_ratio

        function_wrapper.do_update_edit_mesh(me)
        bpy.ops.mesh.select_all(action='SELECT')

        bpy.ops.mesh.remove_doubles(threshold=0.0004, use_unselected=False)

        if (self.use_remesh):
            bpy.ops.mesh.dissolve_limited(angle_limit=0.01)
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            bpy.ops.remesh.boundary_aligned_remesh(edge_length=max_edge_length, iterations=10, reproject=False)

        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        wm.progress_end()

        # fix 2.9 wm.progress problem
        bpy.context.window.cursor_set('NONE')
        bpy.context.window.cursor_set('DEFAULT')

        return {'FINISHED'}

    def ensure_edgelength(self, max_length, mesh, wm):
        seam_edges = list(filter(lambda e: e.seam, mesh.edges))
        edge_groups = defaultdict(list)
        for e in seam_edges:
            edge_groups[math.floor(e.calc_length() / max_length)].append(e)

        wm.progress_begin(0, 99)

        # A little weird, but by grouping the edges by number of required cuts,
        # subdivide_edges() can work a lot more effecient

        # this seems okay
        for progress, k in enumerate(sorted(edge_groups.keys(), reverse=True)):
            if k:
                eg = edge_groups[k]
                #print("subdivide edge", k, len(eg))
                wm.progress_update((progress / len(edge_groups)))
                bmesh.ops.subdivide_edges(mesh, edges=eg, cuts=k)
            else:
                print("hi k0", k)
        bmesh.ops.triangulate(mesh, faces=mesh.faces, quad_method='BEAUTY', ngon_method='BEAUTY')
        # done
