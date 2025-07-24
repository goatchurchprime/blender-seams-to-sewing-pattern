import bpy
from os.path import basename
import bmesh
from bpy.props import (StringProperty,BoolProperty,EnumProperty,IntVectorProperty,FloatProperty)
import mathutils
import random, struct, json
from functools import reduce

tinydisplacement = mathutils.Vector((0,0,0))
tinydisplacement = mathutils.Vector((0,0.00001,0))

class Export_SeamsSTL(bpy.types.Operator):
    """Export Seams encoded into .STL file format."""

    bl_idname = "object.export_seamsstl"
    bl_label = "Export Seams STL/JSON"
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
            ('STL', "Stereolithography (.stl)",
             "Export the seams 3D to a .STL file"),
            ('JSON', "JSON (.json)",
             "Export parts at seams .json file"),
        ),
        name="Format",
        description="File format to export the UV layout to",
        default='JSON',
    )

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.data.uv_layers

    def invoke(self, context, event):
        #stuff to check / set before goes here :)
        self.filepath = self.get_default_file_name(context) + "." + self.file_format.lower()
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def get_default_file_name(self, context):
        return context.active_object.name

    def check(self, context):
        if any(self.filepath.endswith(ext) for ext in (".stl",)):
            self.filepath = self.filepath[:-4]

        ext = "." + self.file_format.lower()
        self.filepath = bpy.path.ensure_ext(self.filepath, ext)
        return True

    def execute(self, context):
        obj = context.active_object
        is_editmode = (obj.mode == 'EDIT')
        if is_editmode:
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        filepath = self.filepath
        filepath = bpy.path.ensure_ext(filepath, "." + self.file_format.lower())
        
        self.export(filepath)

        if is_editmode:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)

        return {'FINISHED'}
    
    def export(self, filepath):
        #get loops:
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_mode(type="FACE")
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)

        print("self.file_format", self.file_format)
        if self.file_format == "STL":
            fout = open(filepath, "wb")
            header = bytearray("Furcut seams in STL".encode())
            while len(header) < 80:
                header.append(0x20)
            fout.write(header)
            fout.write(struct.pack("i", len(bm.faces)))
            for f in bm.faces:
                assert (len(f.verts) == 3), "Mesh is not triangulated!"
                ls = f.loops[:]
                n = (int(ls[0].edge.seam), int(ls[1].edge.seam), int(ls[2].edge.seam))
                fout.write(struct.pack("fff", *n))
                for l in f.loops:
                    fout.write(struct.pack("fff", *(l.vert.co*1000)))
                fout.write(struct.pack("h", 1)) # material
            fout.close()

        elif self.file_format == "JSON":
            fout = open(filepath, "w")
            face_groups = []
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

            bmverts = bm.verts[:]
            jdata = [ ]
            for fluffnumber, fg in enumerate(face_groups):
                g = set()
                for f in fg:
                    g.update(v.index  for v in f.verts)
                g = sorted(g)  # now a list
                verts = [ tuple(bmverts[i].co)  for i in g ]
                gmap = dict(zip(g, range(len(g))))

                facverts = facevertstodoubleup(seam_edges, fg)
                for i in range(len(facverts)):
                    gverts = { }
                    for v in facverts[i][1]:
                        vi = gmap[v.index]
                        gverts[vi] = len(verts)
                        verts.append(verts[vi] + tinydisplacement)
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

                jdata.append({"verts":verts, "tris":tris})
            json.dump(jdata, fout)
            fout.close()
            
        else:
            print("unknown file format", self.file_format)

        bpy.ops.object.mode_set(mode='OBJECT')

        
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
