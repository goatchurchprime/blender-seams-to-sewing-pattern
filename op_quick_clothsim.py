import bpy
from bpy.types import Operator
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
)


class ObjectModeOperator:
    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT'

class QuickClothsim(ObjectModeOperator, Operator):
    """Add clothsim setup to the selected objects"""
    bl_idname = "object.quick_clothsim"
    bl_label = "Quick Clothsim"
    bl_options = {'REGISTER', 'UNDO'}

    use_sewing: BoolProperty(
        name="Sewing",
        description="Enable Shape > Sewing",
        default=True,
    )
    use_gravity: BoolProperty(
        name="Gravity",
        description="Use world gravity",
        default=False,
    )
    pressure_style: EnumProperty(
        name="Pressure",
        items=(
            ('OFF', "Off", "No pressure"),
            ('MEDIUM', "Medium", "Medium pressure"),
            ('HIGH', "High", "High pressure")
        ),
        default='MEDIUM',
    )
    air_visc: BoolProperty(
        name="Increase Air Viscosity",
        description="Increas Physical Properties > Air Viscosity. Use when the objects collapses in on itself during simulation",
        default=True,
    )

    def execute(self, context):
        objects = bpy.context.selected_objects
        if objects is not None :
            for obj in objects:
                    cloth_mod = obj.modifiers.new(name = 'Cloth', type = 'CLOTH')
                    cloth_mod = obj.modifiers.get('Cloth')
                    #pressure
                    if (self.pressure_style != 'OFF'):
                        cloth_mod.settings.use_pressure = True
                    if (self.pressure_style == 'MEDIUM'):
                        cloth_mod.settings.uniform_pressure_force = 10
                    if (self.pressure_style == 'HIGH'):
                        cloth_mod.settings.uniform_pressure_force = 50
                        
                    #sewing
                    cloth_mod.settings.use_sewing_springs = self.use_sewing
                    if (self.pressure_style == 'MEDIUM'):
                        cloth_mod.settings.sewing_force_max = 5
                    if (self.pressure_style == 'HIGH'):
                        cloth_mod.settings.sewing_force_max = 15
                    
                    #viscosity
                    if (self.air_visc == True):
                        cloth_mod.settings.air_damping = 10
                    if (self.use_gravity == False):
                        cloth_mod.settings.effector_weights.gravity = 0
        return {'FINISHED'}
