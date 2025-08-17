bl_info = {
        'name': 'Seams to Sewing Pattern',
        'author': 'Thomas KKKKKole',
        'version': (1, 0),
        'blender': (2, 80, 0),
        'category': 'Cloth',
        'description': 'Converts a mesh with seams into sewing patterns with sewing edges',
        'location': 'Object > Seams to Sewing Pattern > ...',
        'wiki_url': 'https://blenderartists.org/t/1248713'}

if "bpy" in locals():
    import importlib
    importlib.reload(op_seams_to_sewingpattern)
    importlib.reload(op_export_sewingpattern)
    importlib.reload(op_export_flatpattern)
    importlib.reload(op_export_stl_seams)
    importlib.reload(op_quick_clothsim)
    importlib.reload(op_boundary_alinged_remesh)
    importlib.reload(op_clean_up_edges)
    importlib.reload(op_freecad_flatten_component)
else:
    from . import op_seams_to_sewingpattern
    from . import op_export_sewingpattern
    from . import op_export_flatpattern
    from . import op_export_stl_seams
    from . import op_quick_clothsim
    from . import op_boundary_alinged_remesh
    from . import op_clean_up_edges
    from . import op_freecad_flatten_component

print("RELOADING SEAMS TO SEWINGPATTERN ADDON*** YYYYYYYY ***")

import bpy
from bpy.types import Menu

def clean_up_func(self, context):
    self.layout.separator()
    self.layout.operator("mesh.clean_up_knife_cut")

def menu_func(self, context):
    lay_out = self.layout
    lay_out.operator_context = 'INVOKE_REGION_WIN'

    lay_out.separator()
    lay_out.menu("VIEW3D_MT_object_seams_to_sewing_pattern_menu",
                text="Seams to Sewing Pattern")
    
class VIEW3D_MT_object_seams_to_sewing_pattern_menu(Menu):
    bl_idname = "VIEW3D_MT_object_seams_to_sewing_pattern_menu"
    bl_label = "Seams to Sewing Pattern"

    def draw(self, context):
        layout = self.layout
        layout.operator("object.freecad_flatten_component", text="Freecad flatten component", icon="OUTLINER_DATA_SURFACE")
        layout.separator()
        layout.operator("object.export_flatpattern", text="Export flat Pattern", icon="OUTLINER_DATA_SURFACE")
        layout.separator()
        layout.operator("object.seams_to_sewingpattern", text="Seams to Sewing Pattern", icon="OUTLINER_DATA_SURFACE")
        layout.separator()
        layout.operator("object.export_sewingpattern", text="Export Sewing Pattern (.svg)", icon="EXPORT")
        layout.separator()
        layout.operator("object.export_seamsstl", text="Export JSON with seams", icon="EXPORT")
        layout.separator()
        layout.operator("object.quick_clothsim", text="Quick Clothsim", icon="MOD_CLOTH")

# Register
classes = [
    VIEW3D_MT_object_seams_to_sewing_pattern_menu,
    op_seams_to_sewingpattern.Seams_To_SewingPattern,
    op_export_sewingpattern.Export_Sewingpattern,
    op_export_flatpattern.Export_Flatpattern,
    op_export_stl_seams.Export_SeamsSTL, 
    op_quick_clothsim.QuickClothsim,
    op_boundary_alinged_remesh.Remesher,
    op_clean_up_edges.CleanUpEdges,
    op_freecad_flatten_component.Freecad_flatten_component
    ]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    # Adds submenu in View3D > Seams to Sewing Pattern
    bpy.types.VIEW3D_MT_object.append(menu_func)
    bpy.types.VIEW3D_MT_edit_mesh_clean.append(clean_up_func)
    bpy.types.VIEW3D_MT_edit_mesh_edges.append(clean_up_func)
    bpy.types.VIEW3D_MT_edit_mesh_context_menu.append(clean_up_func)



def unregister():
    
    bpy.types.VIEW3D_MT_edit_mesh_clean.remove(clean_up_func)
    bpy.types.VIEW3D_MT_edit_mesh_edges.remove(clean_up_func)
    bpy.types.VIEW3D_MT_edit_mesh_context_menu.append(clean_up_func)
    bpy.types.VIEW3D_MT_object.remove(menu_func)
    # Removes submenu
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

if __name__ == "__main__":
    register()
