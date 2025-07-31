import bpy
from os.path import basename
from xml.sax.saxutils import escape
from bpy.props import (
    StringProperty,
    BoolProperty,
    EnumProperty,
    IntVectorProperty,
    FloatProperty,
)
import bmesh
import mathutils
import random

class Export_Flatpattern(bpy.types.Operator):
    """Export Sewingpattern to .SVG file format. This should be called after the Seams to Sewing Pattern operator"""

    bl_idname = "object.export_flatpattern"
    bl_label = "Export Sewing Pattern"
    bl_options = {'REGISTER', 'UNDO'}

    filepath: StringProperty(
        subtype='FILE_PATH',
    )
    alignment_markers: EnumProperty(
        items=(
            ('OFF', "Off",
             "No alignment markers"),
            ('SEAM', "Marked as seam",
             "Use sewing edges manually marked as seam"),
            ('AUTO', "Autodetect + seam",
             "Finds sewing edges of corners automatically and marks them as seam"),
        ),
        name="Alignment markers",
        description="Exports matching colored lines on the borders of sewing patterns to assist with alignment",
        default='AUTO',
    )
    file_format: EnumProperty(
        items=(
            ('SVG', "Scalable Vector Graphic (.svg)",
             "Export the sewing pattern to a .SVG file"),
        ),
        name="Format",
        description="File format to export the UV layout to",
        default='SVG',
    )

    @classmethod
    def poll(cls, context):
        return True

    def invoke(self, context, event):
        #stuff to check / set before goes here :)
        self.filepath = self.get_default_file_name(context) + "." + self.file_format.lower()
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def get_default_file_name(self, context):
        return context.active_object.name

    def check(self, context):
        if any(self.filepath.endswith(ext) for ext in (".png", ".eps", ".svg")):
            self.filepath = self.filepath[:-4]

        ext = "." + self.file_format.lower()
        self.filepath = bpy.path.ensure_ext(self.filepath, ext)
        return True

    def execute(self, context):
        i = bpy.data.collections.find("flatmeshes")
        if i == -1:
            print("no flatmeshes collection")
            return {'CANCELLED'}
        mcoll = bpy.data.collections[i]

        filepath = self.filepath
        filepath = bpy.path.ensure_ext(filepath, "." + self.file_format.lower())
        self.export(filepath, mcoll)
        return {'FINISHED'}
    
    def export(self, filepath, mcoll):
        document_scale = 1000.0 #millimeter
        document_width = 1000.0
        fout = open(filepath, "w")
        fout.write('<svg xmlns="http://www.w3.org/2000/svg"\n viewBox="0 0 ' + str(document_scale) + ' ' + str(document_scale) +'" ')
        fout.write('width="' + str(document_width) + 'mm" height="' + str(document_width) + 'mm">\n')
        fout.write('<defs><style>.seam{stroke: #000; stroke-width:0.5px; fill:white} .sewinguide{stroke-width:0.5px;}</style></defs>\n')

        for obj in mcoll.objects:
            bm = bmesh.new()
            bm.from_mesh(obj.data)

            bls = { }
            for f in bm.faces:
                for l in f.loops:
                    if l.edge.is_boundary:
                        bls[l.vert] = l.edge.verts[0] if l.vert != l.edge.verts[0] else l.edge.verts[1]
            vbound = list(bls.popitem())
            while vbound[-1] in bls:
                vbound.append(bls.pop(vbound[-1]))

            fout.write('<g><path class="seam" d="M')
            for v in vbound:
                fout.write('%.1f,%.1f ' % (v.co[0]*1000, (-v.co[1])*1000))
            fout.write('"/></g>\n')
        fout.write('</svg>\n')
            
        bpy.ops.object.mode_set(mode='OBJECT')
        
