import bpy
import bmesh

def do_bevel(offset=0.0002):
    bpy.ops.mesh.bevel(vertex_only=False, offset=offset)

def do_update_edit_mesh(m):
    bmesh.update_edit_mesh(m, False)
