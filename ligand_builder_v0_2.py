# This module is created to collect different methods for ligand
# generation. *Under development*
#
#
# The data structure of the excel file should be:
#   Data in Sheet1
#   Chosen names of ligands in column A
#   Smiles stored in column B, append name (for now)
#   n lables in column A should equal n smiles in column B
#
# ----------------------------------------------------------------------
# ligand_builder v0.2, beta pre-version 1, Mary Pitman, Mobley Lab UCI
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# ----------------------------------------------------------------------

__doc__="""
        Modules designed to generate ligands topologies and setup files
        for pmxworkflow via OpenBabel or Openeye. Reads data from excel
        input file.
        """
__all__ = ['collect_excel_data', 'smi_to_sdf']
__version__ = '0.2'
__author__ = 'M. Pitman'


import os
import sys
import subprocess as sub

import openpyxl
from openforcefield.topology import Molecule, Topology

# Edit for conditional import.
from openeye import oechem
from openeye.oechem import *

# ----------------------------------------------------------------------

def collect_excel_data(input_excel_file: str) -> str:
    """ reads in 'my_data.xlsx' """
    try:
        global MY_DATA, SHEET
        MY_DATA = openpyxl.load_workbook(input_excel_file)
        SHEET = MY_DATA['Sheet1']
    except:
        print("The argument for read_excel_file is 'my_data.xlsx' .")
        sys.exit()
        
    # Test that the shape of the lists are equal.
    if len(SHEET['A']) != len(SHEET['B']):
        print("Unequal number of ligand names and smiles *.xlsx. Assign a name to each smiles")
        sys.exit()
        
    # Retrieve the SMILES from the my_data.xlsx.
    global SMILES_COLLECTOR
    SMILES_COLLECTOR = ""
    for cell in SHEET['B']:
        SMILES_COLLECTOR += str('-:"') + str(cell.value) + str('" ')
          
    # Retrieve the names of each ligands from the my_data.xlsx.
    global NAME_COLLECTOR
    NAME_COLLECTOR = ""
    for cell in SHEET['A']:
        NAME_COLLECTOR += str(cell.value) + str(' ')
   
# Not totally sure the technical output type here. Need to double check
def smi_to_sdf(toolkit_to_use: str, input_excel_file: str) -> str:
    """ calculate a .sdf for each ligand using openeye or openbabel
    example: smi_to_sdf('openbabel', 'my_data.xlsx')
    """
    collect_excel_data(input_excel_file)
    
    # Generate the .sdf container using obabel.
    if toolkit_to_use is str("openbabel"):
        import openbabel
        for cell in SHEET['B']:
            smiles_string = str('-:"') + str(cell.value) + str('" ')
            cell_shift = cell.row - 1
            lig_name = SHEET['A'][cell_shift].value
            conversion = str("obabel {} -osdf >  {}.sdf --gen3d"
                            .format(smiles_string, lig_name))
            sub.call(conversion, shell=True)
     

    # Generate the .sdf container using openeye.
    if toolkit_to_use is str("openeye"):
        
        for cell in SHEET['B']:
            # Shift to cell A in *.xlsx to read ligand name.
            cell_shift = cell.row - 1
            lig_name = SHEET['A'][cell_shift].value
            smiles = str(cell.value)
            
            ifs = oemolistream()
            ifs.SetFormat(oechem.OEFormat_SMI)
            ifs.openstring(smiles)
            ofs = oemolostream('{}.sdf'.format(lig_name))
            mol = OEMol()

            if (OEParseSmiles(mol, smiles) == 1):
                # Add an isomeric canonical smiles code to the bottom of the .sdf file.
                OESetSDData(mol, "$SMI", OECreateIsoSmiString(mol))
                OEWriteMolecule(ofs, mol)
            else:
                sys.stderr.write("SMILES string for {} was invalid\n".format(lig_name))
                sys.exit()
          



#the saves need to be moved to the correct folders
            
#def organize_sdfs(input_excel_file):

        













