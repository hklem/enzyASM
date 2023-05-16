from pymol import cmd
import os

PDB_code = '3zr4'

pdb_file = '/{}.pdb'.format(PDB_code)

cat_center_file = '/catalytic_center.pdb'
active_site_core_file = '/active_site_radius_5.pdb'
active_site_truncated_file = '/truncated_active_site_radius_5_protonated.pdb'

cwd_path = os.getcwd()

path_to_originalPDB = cwd_path+pdb_file
path_to_cat_center = cwd_path+cat_center_file
path_to_core = cwd_path+active_site_core_file
path_to_model = cwd_path+active_site_truncated_file

cmd.load(path_to_originalPDB,'Original_PDB')
cmd.color('gray70','Original_PDB')
cmd.hide('sticks','Original_PDB')
cmd.orient()
cmd.scene('Whole protein', 'store')

cmd.load(path_to_cat_center,'Catalytic_Center')
cmd.show_as('sticks','Catalytic_Center')
cmd.color('cyan','Catalytic_Center')
cmd.hide("everything", "Original_PDB")
cmd.zoom("visible")
cmd.orient("visible")
cmd.scene('Catalytic center', 'store')

cmd.load(path_to_core,'Active_Site')
cmd.hide("everything", "Catalytic_Center")
arg1 = 'Active_Site'
cmd.show("surface",arg1)
cmd.set("surface_color","white") 
cmd.set("cartoon_transparency",0.5,arg1)
cmd.show("cartoon", "Original_PDB")
cmd.show("sticks", arg1)
cmd.show("spheres", arg1)
cmd.color("gray85","elem C and "+arg1)
cmd.color("gray98","elem H and "+arg1)
cmd.color("slate","elem N and "+arg1)
cmd.set("stick_radius",0.07, arg1)
cmd.set("sphere_scale",0.18, arg1)
cmd.set("sphere_scale",0.13, arg1+" and elem H")
cmd.set("dash_gap",0.01, arg1)
cmd.set("dash_radius",0.07, arg1)
cmd.set("stick_color","black", arg1)
cmd.set("dash_gap",0.01)
cmd.set("dash_radius",0.035)
cmd.hide("nonbonded", arg1)
cmd.hide("lines", arg1)
cmd.zoom("visible")
cmd.scene('Active site', 'store')

cmd.load(path_to_model,'Truncated_Model')
cmd.hide("everything", "Active_Site")
cmd.hide("everything", "Original_PDB")
arg1 = 'Truncated_Model'
cmd.show("sticks", arg1)
cmd.show("spheres", arg1)
cmd.color("gray85","elem C and "+arg1)
cmd.color("gray98","elem H and "+arg1)
cmd.color("slate","elem N and "+arg1)
cmd.set("stick_radius",0.07, arg1)
cmd.set("sphere_scale",0.18, arg1)
cmd.set("sphere_scale",0.13, arg1+" and elem H")
cmd.set("dash_gap",0.01, arg1)
cmd.set("dash_radius",0.07, arg1)
cmd.set("stick_color","black", arg1)
cmd.set("dash_gap",0.01)
cmd.set("dash_radius",0.035)
cmd.hide("nonbonded", arg1)
cmd.hide("lines", arg1)
cmd.zoom("visible")
cmd.scene('Truncated model', 'store')

cmd.show("cartoon", "Original_PDB")
cmd.show("surface","Active_Site")
cmd.set("cartoon_transparency",0.7,"Original_PDB")
cmd.zoom("all")
cmd.scene('Overview', 'store')

