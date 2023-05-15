################################################
# Code written by Heidi Klem while at 
# Colorado State University as a graduate student
# in the Paton and McCullagh groups. 
#################################################

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel as ob

def define_catalytic_center(pdb_file, catalytic_center_resname=None, \
catalytic_center_resnum=None, catalytic_center_chain='*'):
    protein_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False,proximityBonding=False)
    catalytic_center_mol = Chem.RWMol(protein_mol)
    count = 0
    if catalytic_center_resname==None and catalytic_center_resnum==None:
        #catalytic_center_resname='CCC'
        resnumbers = []
        for atom in reversed(protein_mol.GetAtoms()):
            if res_info(atom,'chain') == '*':
                if res_info(atom,'number') in resnumbers:
                    continue
                else:
                    resnumbers.append(res_info(atom,'number'))
                    continue
            else:
                catalytic_center_mol.RemoveAtom(atom.GetIdx())
        if len(resnumbers) == 0:
            print("WARNING: No atoms found matching catalytic center definition.")
        if len(resnumbers) > 1:
            print("WARNING: There were {} residues with distinct residue numbers ({}) containing atom(s) corresponding to the specified chain {}. Double check this is what was intended, otherwise please specify residue number in addition to the residue name to ensure uniqueness.".format(len(resnumbers),resnumbers,catalytic_center_chain))


    if catalytic_center_resname!=None and catalytic_center_resnum==None:
        print("Catalytic center defined by residue name {} in the PDB file".\
        format(catalytic_center_resname))
        resnumbers = []
        for atom in reversed(protein_mol.GetAtoms()):
            if res_info(atom,'name') == catalytic_center_resname:
                if res_info(atom,'number') in resnumbers:
                    continue
                else:
                    resnumbers.append(res_info(atom,'number'))
                    continue
            else:
                catalytic_center_mol.RemoveAtom(atom.GetIdx())
        if len(resnumbers) == 0:
            print("WARNING: No atoms found matching catalytic center definition.")
        if len(resnumbers) > 1:
            print("WARNING: There were {} residues with distinct residue numbers ({}) containing atom(s) corresponding to the specified residue name {}. Double check this is what was intended, otherwise please specify residue number in addition to the residue name to ensure uniqueness.".format(len(resnumbers),resnumbers,catalytic_center_resname))

    if catalytic_center_resname==None and catalytic_center_resnum!=None:
        print("Catalytic center defined by residue number {} in the PDB file".\
        format(catalytic_center_resnum))
        resnames = []
        for atom in reversed(protein_mol.GetAtoms()):
            if res_info(atom,'number') == catalytic_center_resnum:
                if res_info(atom,'name') in resnames:
                    continue
                else:
                    resnames.append(res_info(atom,'name'))
                    continue
            else:
                catalytic_center_mol.RemoveAtom(atom.GetIdx())
        if len(resnames) == 0:
            print("WARNING: No atoms found matching catalytic center"+\
            " definition.")
        if len(resnames) > 1:
            print("WARNING: There were {}".format(len(resnames)) +" residues"+\
            " with distinct residue names ({})".format(resnames) +\
            " containing atom(s) corresponding to the residue number {}."\
            .format(catalytic_center_resnum))
            print("Double check this is what was intended, otherwise please" +\
            " specify residue name in addition to the residue number to" +\
            " ensure uniqueness.")

            #print("WARNING: There were {} residues with distinct residue names"+\
            #" ({}) containing atom(s) corresponding to the residue number {}."+\
            #" Double check this is what was intended, otherwise please specify"+\
            #" residue name in addition to the resisude number to ensure"+\
            #" uniqueness.".format(len(resnames),resnames,catalytic_center_resnum))

    if catalytic_center_resname!=None and catalytic_center_resnum!=None:
        print("Catalytic center defined by residue number {}".\
        format(catalytic_center_resnum) + " and residue name {}".\
        format(catalytic_center_resname) + " in the PDB file.")

        for atom in reversed(protein_mol.GetAtoms()):
            if res_info(atom,'number') == catalytic_center_resnum and res_info\
            (atom,'name') == catalytic_center_resname:

                count += 1
                continue
            else:
                catalytic_center_mol.RemoveAtom(atom.GetIdx())
        if count == 0:
            print("WARNING: No atoms found matching catalytic center definition.")
    print('Catalytic center contains {} atoms.'.format(catalytic_center_mol.\
    GetNumAtoms()))

    return catalytic_center_mol

def residue_shell(base_mol,center_mol,radius,centroid=False,include_residues=[]):
    #keep_residues = []
    res_name, res_number, res_chain, res_atom, atom_type = [], [], [], [], []
    N_termini_interactions = []
    side_chain_interactions = []
    C_termini_interactions = []
    current_res = 'nothing'
    backbone_atoms=[' O  ',' C  ',' N  ',' CA ']
    
    if centroid == True:
        centroid_coords = np.asarray(Chem.rdMolTransforms.\
            ComputeCentroid((center_mol.GetConformer())))
        distances = [np.linalg.norm(np.asarray(center_mol.GetConformer().\
            GetAtomPosition(atom.GetIdx()))-centroid_coords) for atom in \
            center_mol.GetAtoms()]
        radius_buffer = np.max(distances)
        
        for atom1 in base_mol.GetAtoms():
            if str(res_info(atom1,'name'))+str(res_info(atom1,'number')) in include_residues:
                atomic_distance = radius     
            else:
                coords1 = np.asarray(base_mol.GetConformer().GetAtomPosition(atom1.GetIdx()))
                atomic_distance = np.linalg.norm(coords1-centroid_coords)
            if atomic_distance < radius+radius_buffer:
                res_name.append(res_info(atom1,'name'))
                res_number.append(res_info(atom1,'number'))
                res_chain.append(res_info(atom1,'chain'))
                res_atom.append(res_info(atom1,'atom_name'))
                if res_info(atom1,'atom_name') in backbone_atoms:
                    atom_type.append('Backbone')
                else: 
                    atom_type.append('Sidechain')
                    
    already_added = []                
    if centroid == False:
        for i,atom1 in enumerate(base_mol.GetAtoms()):
            if str(res_info(atom1,'name'))+str(res_info(atom1,'number')) in include_residues:
                keep_atom = True
            else:
                keep_atom = False
                coords1 = np.asarray(base_mol.GetConformer().GetAtomPosition(atom1.GetIdx()))
            for j,atom2 in enumerate(center_mol.GetAtoms()):
                if keep_atom == True:
                    atomic_distance = radius-1 
                if keep_atom == False:
                    coords2 = np.asarray(center_mol.GetConformer().GetAtomPosition(atom2.GetIdx()))
                    atomic_distance = np.linalg.norm(coords1-coords2)    
                if atomic_distance < radius:
                    if str(res_info(atom1,'atom_name'))+str(res_info(atom1,'number')) in already_added:
                        continue
                    res_name.append(res_info(atom1,'name'))
                    res_number.append(res_info(atom1,'number'))
                    res_chain.append(res_info(atom1,'chain'))
                    res_atom.append(res_info(atom1,'atom_name'))
                    already_added.append(str(res_info(atom1,'atom_name'))+str(res_info(atom1,'number')))
                    if res_info(atom1,'atom_name') in backbone_atoms:
                        atom_type.append('Backbone')
                    else: 
                        atom_type.append('Sidechain')   
                                         
    new_mol = Chem.RWMol(base_mol)
    
    res_dict = {'Residue Name':res_name,'Residue Number':res_number,'Residue Atom':res_atom,'Atom Type':atom_type,'Residue Chain':res_chain}
    
    for atom in reversed(base_mol.GetAtoms()):
        if res_info(atom,'number') in res_number:
            if res_info(atom,'chain') == res_chain[res_number.index(res_info(atom,'number'))]:
                continue
            else:
                new_mol.RemoveAtom(atom.GetIdx())
        else:
            new_mol.RemoveAtom(atom.GetIdx())

    return new_mol, res_dict

def truncate(base_mol, skip_residues=['HOH'], skip_resnumbers=[], remove_resnumbers=[], remove_atoms=[' O  ',' C  ',' N  '], constrain_atoms=[' CA ']):

    new_mol = Chem.RWMol(base_mol)
    proline_count = 0
    CA_list = []
    
    
    
    for atom in reversed(base_mol.GetAtoms()):

        if res_info(atom,'name') in skip_residues:
            continue
        if res_info(atom,'number') in skip_resnumbers:
            continue
        if res_info(atom,'number') in remove_resnumbers:
            new_mol.RemoveAtom(atom.GetIdx())
            continue
        if res_info(atom,'atom_name') in remove_atoms:
            if res_info(atom,'name') == 'PRO' and \
                    res_info(atom,'atom_name') == ' N  ':
                proline_count += 1
                continue
            new_mol.RemoveAtom(atom.GetIdx())

    for atom in new_mol.GetAtoms():
        if res_info(atom,'atom_name') in constrain_atoms:
            CA_list.append(atom.GetIdx()+1)
            
    if proline_count > 0:
        print('WARNING: active site model contains {} proline(s), which cannot be methyl capped. N atom was kept and one artificial H was added.'.format(proline_count))
    return new_mol, CA_list


def obabel_protonate(active_site_pdb_file):

    ob_conv = ob.OBConversion()
    ob_conv.SetInFormat('pdb')
    ob_mol = ob.OBMol()
    ob_conv.ReadFile(ob_mol, str(active_site_pdb_file))

    ob_mol.AddHydrogens(False, True, 7.0) # polaronly, correctForPH, pH
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('pdb', 'pdb')

    his_ids = []

    active_site_mol = Chem.MolFromPDBFile(str(active_site_pdb_file),removeHs=False,sanitize=False,proximityBonding=False)
    for atom in active_site_mol.GetAtoms():
        if res_info(atom,'name') == 'HIS':
            his_ids.append(res_info(atom,'number'))

    if len(his_ids) > 0:
        print('WARNING: There are {} HIS residues, user should check protonation state and adjust pdb output file accordingly. Default is epsilon_N.'.format(len(np.unique(his_ids))))

    file_name_base = str(active_site_pdb_file).split('.pdb')[0]
    conv.WriteFile(ob_mol, str(file_name_base)+'_protonated'+'.pdb')


def res_info(atom,info):
    ''' Options for info are name, number, or chain.
    '''
    if str(info) == 'name':
        return atom.GetPDBResidueInfo().GetResidueName()
    if info == 'atom_name':
        return atom.GetPDBResidueInfo().GetName()
    if info == 'number':
        return atom.GetPDBResidueInfo().GetResidueNumber()
    if info == 'chain':
        return atom.GetPDBResidueInfo().GetChainId()


#def evaluate_interactions(res_dict):
