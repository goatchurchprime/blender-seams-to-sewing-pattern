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

if bpy.app.version >= (3, 0, 0):
    from . import function_wrapper_3_0 as function_wrapper
elif bpy.app.version >= (2, 90, 0):
    from . import function_wrapper_2_9 as function_wrapper
else:
    from . import function_wrapper_2_8 as function_wrapper


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
    dont_flatten: BoolProperty(
        name="Dontflatten",
        description="Don't actually flatten the panels of the new mesh",
        default=False,
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
        default=25000,
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
        row.prop(self, "dont_flatten")
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
            bpy.ops.object.convert(target='MESH')
            obj = bpy.context.active_object

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

        if self.bevel_skim:
            function_wrapper.do_bevel()
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

            bmesh.ops.collapse(bm, edges=degenerate_edges, uvs=True)

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
            if self.dont_flatten:
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

            # offset each face island by their UV value, using the tangent and
            # bitangent

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
        if not self.dont_flatten:
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

        return{'FINISHED'}

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
