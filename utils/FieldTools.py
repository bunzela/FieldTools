#!/usr/bin/env python
# coding: utf-8

import sys
import os
import subprocess
import numpy as np  
import pytraj as pt
import pickle
import datetime

def Field_of_ABC_at_XYZ_projected_onto_VEC(_ABC,_XYZ,_VEC,_Charge):
    
    # Field at _xyz from Atom _abc with _Charge 
    _Couloumb_Const = 1.00 #(A.U.)
    _Vector = _ABC-_XYZ        
    _Length = np.linalg.norm(_Vector)
    _Normalized_Vector = (-1.0)*_Vector/_Length
    _Length_Bohr = _Length/0.529177249           # Bohr
    _Field = (_Couloumb_Const*_Charge/(_Length_Bohr**2))*_Normalized_Vector
    _Field = np.dot(_Field,_VEC)
    ### 1 a.u. = 5.14x10^9 V/cm
    ### 1 a.u. = 5.14x10^3 MV/cm
    _Field=_Field*5.14*1000
    return _Field

def Field_of_ABC_at_XYZ(_ABC,_XYZ,_Charge):
    
    # Field at _xyz from Atom _abc with _Charge 
    _Couloumb_Const = 1.00 #(A.U.)
    _Vector = _ABC-_XYZ        
    _Length = np.linalg.norm(_Vector)
    _Normalized_Vector = (-1.0)*_Vector/_Length
    _Length_Bohr = _Length/0.529177249           # Bohr
    _Field = (_Couloumb_Const*_Charge/(_Length_Bohr**2)) ### For Field only*_Normalized_Vector
    ### 1 a.u. = 5.14x10^9 V/cm
    ### 1 a.u. = 5.14x10^3 MV/cm
    _Field=_Field*5.14*1000
    return _Field

def fkt_Update_Charges(QM_Charges,QM_Dict,Frame_i,Names,Charges):
    stop="X"
    for QM_Charge_i in range(0,len(QM_Charges),1):
        if len(QM_Charges[QM_Charge_i])<>1: continue
        if QM_Charges[QM_Charge_i][0]==str(Frame_i+1):
            start = QM_Charge_i+1
        if QM_Charges[QM_Charge_i][0]==str(Frame_i+2):
            stop = QM_Charge_i
            QM_Charges = QM_Charges[start:stop]
            break
    if stop == "X":
        QM_Charges = QM_Charges[start:]   
            
    for QM_Atom_i in range(0,len(QM_Dict),1):
        for Atom_i in range(0,len(Names),1):
            if Names[Atom_i][2] == QM_Dict[QM_Atom_i][2]:      #QM_Atom_i index in QM_Charges to be changed from   
                if Names[Atom_i][0] == QM_Dict[QM_Atom_i][0]:  #Atom_i index in Charges to be changed 
                    Charges[Atom_i]=QM_Charges[QM_Atom_i][2]
                    break
    return Charges

def fkt_Load_QM_Charges(arg_qm_charges,arg_qm_dict):
    with open(arg_qm_charges) as f:
        QM_Charges = [i.split() for i in f.readlines()]
    with open(arg_qm_dict) as f:
        QM_Dict = [i.split() for i in f.readlines()]
    return QM_Charges,QM_Dict

def fkt_calc_fields(FIELDS,FIELD,XYZ,VEC,Frame,Charges,Names,Field_Components,arg_solvent,Self_i,arg_exclude_atoms):

### VARIABLES STRUCTURE
    Absolut_Fields = {}
    for Field_Component in Field_Components:
        Absolut_Fields[Field_Component] = 0.0
        
    for Atom_Index in range(0,len(Frame),1):   

#### Exclude Self-residue!!! 
        if arg_exclude_atoms <> False:
            if Names[Atom_Index][2] == Self_i: continue        
        else:
            if Atom_Index in arg_exclude_atoms[FIELD]: continue
        
        if len(FIELD.split("_")) == 1:   #POINT CALCULATION
            Absolut_Field = Field_of_ABC_at_XYZ(Frame[Atom_Index],XYZ,Charges[Atom_Index])
        if len(FIELD.split("_")) == 2:   #VECTOR CALCULATION
            Absolut_Field = Field_of_ABC_at_XYZ_projected_onto_VEC(Frame[Atom_Index],XYZ,VEC,Charges[Atom_Index])

        #### Total Field
        Absolut_Fields["Total"]+=Absolut_Field
        #### Solvent Field
        if Names[Atom_Index][1] in arg_solvent:
            Absolut_Fields["Solvent"]+=Absolut_Field
            Absolut_Fields[Names[Atom_Index][1]]+=Absolut_Field
        else:
            Absolut_Fields["Protein"]+=Absolut_Field
            Absolut_Fields[Names[Atom_Index][2]]+=Absolut_Field
    
    for Field_Component in Field_Components:
        FIELDS[FIELD][Field_Component]+=[Absolut_Fields[Field_Component]]
    return FIELDS

def fkt_Field_Components(FIELDS,Names,arg_solvent):
    
    Field_Components = ["Total","Protein","Solvent"]+arg_solvent
    
    for Name in Names:
        if Name[1] not in arg_solvent:
            if Name[2] not in Field_Components: Field_Components += [Name[2]] 

    for Field in FIELDS:
        for Field_Component in Field_Components:
            FIELDS[Field][Field_Component]=[]
    
    return FIELDS,Field_Components

def fkt_get_Target_VEC(Target_XYZ):
    
    Target_VEC = {}

    for Target in Target_XYZ:
        if len(Target_XYZ[Target])==2:
            Target_VEC[Target]={}
            XYZ_START = Target_XYZ[Target][0]
            XYZ_STOP  = Target_XYZ[Target][1]
            Target_VEC[Target]["Center"] = (XYZ_START+XYZ_STOP)/2 
            Target_VEC[Target]["Vector"] = (XYZ_STOP-XYZ_START)/np.linalg.norm(XYZ_STOP-XYZ_START)
            
    return Target_VEC

def fkt_Target_Index(Names,arg_target,FIELDS):
    
    Target_Index={}
    
    with open(arg_target) as f:
        Targets = [i.split() for i in f.readlines()]
        
    for Target_i in range(0,len(Targets),1):
        Target_Name = "_".join(Targets[Target_i])
        FIELDS[Target_Name]={}
        for Atom_i in range(0,len(Targets[Target_i]),1): #If this is a vector, it contains two Atoms!
            Atom_ResidueNr = Targets[Target_i][Atom_i].split("@")[0][1:]
            Atom_Name      = Targets[Target_i][Atom_i].split("@")[1]
            for PDB_i in range(0,len(Names),1):
                if Atom_ResidueNr == Names[PDB_i][2]:      #Check if Residue Number matches
                    if Atom_Name == Names[PDB_i][0]:       #Check if Atom Name matches
                        Targets[Target_i][Atom_i] = PDB_i  #Update Number in Target
                        break
        Target_Index[Target_Name] = Targets[Target_i]
    return Target_Index,FIELDS #Passed to Target_Index in main!

def fkt_get_Target_XYZ(Frame,Target_Index):

#Make empty library
    Target_XYZ = {}

    for Target in Target_Index:
        Target_XYZ[Target]=[]
        for Atom_i in range(0,len(Target_Index[Target]),1): #If this is a vector, it contains two Atoms!
            Target_XYZ[Target]=Target_XYZ[Target]+[Frame[Target_Index[Target][Atom_i]]]
                
    return Target_XYZ #Passed to Target_Index in main!

def fkt_Load_Trajectory(arg_nc,arg_param,arg_pdb):

    traj = pt.load(arg_nc, arg_param) 
    Trajectory = traj.xyz
    Charges = traj.topology.charge

    with open(arg_pdb) as f:
        Names = [i for i in f.readlines() if len(i)==81] # PDB lines of Atoms and their xyz coords have 81 characters
    for i in range(0,len(Names),1):
        temp = Names[i].split()
####################AtmName , ResName , ResiNR
####################0       , 1       , 2        
        Names[i] = [temp[2], temp[3], temp[4]]   
    return (Trajectory, Charges, np.array(Names))
 
def make_pdb(arg_tmp_dir,arg_nc,arg_param,arg_pdb):
    cpptraj = """parm """+arg_param+"""
trajin """+arg_nc+"""
autoimage    
outtraj """+arg_pdb+""" onlyframes -1
strip !("""+arg_qm_mask+""")
outtraj """+".".join(arg_out.split(".")[:-1])+""".pdb onlyframes -1
"""
    with open(arg_tmp_dir+"cpptraj.in", "w") as f:
        f.write(cpptraj)  
    with open(arg_tmp_dir+"/cppraj.out", 'w') as f:
        process = subprocess.call(['cpptraj', '-i', arg_tmp_dir+'cpptraj.in'], stdout=f)
    
def make_tmpdir(arg_tmp_dir):
    if os.path.isdir(arg_tmp_dir):
        for f in os.listdir(arg_tmp_dir):
            os.remove(os.path.join(arg_tmp_dir, f))
    else:
        os.mkdir(arg_tmp_dir)
        
def load_qm_mask(arg_qm_mask):
    with open(arg_qm_mask) as f:
        arg_qm_mask = f.read()[:-1]
    return arg_qm_mask

def usage():
    print "                                                                                                               "
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-------------                                 FieldTools Usage:                                 -------------|"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-nc             : Specify trajectory file                                                                    |"
    print "|-parm           : Specify parameter file                                                                     |"
    print "|-out            : Specify output file                                                                        |"
    print "|-tmp_dir        : Optional: Specify temporary directory                                [default ./tmp_field/]|"
    print "|-target         : Specify target file                                                                        |"
    print "|-solvent        : Optional: Comma seperated list of non-protein residues               [default WAT,Na+,Cl- ]|"
    print "|-use_qm_charges : Optional: Specify if qmcharges should be read                        [default false       ]|"
    print "|-exclude_atoms  : Optional: Specify what atoms will be excluded from Field calculation                       |"
    print "|                : If not set, Field calculation will exclude the full residue of the first atom defined      |"
    print "|                : for each target.                                                                           |"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-------------         If -use_qm_charges = true, the following options have to be set:          -------------|"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-qm_mask        : Specify qm region (amber selection mask)                                                   |"
    print "|-qm_charges     : Specify qmcharges file (qm partial charges of each frame, made with QMChargesTools)        |"
    print "|-qm_dict        : Specify qmdict file (containing qmatom names in .pdb and gaussian)                         |"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-------------                               How to define -target:                              -------------|"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-target should point to a file, each line in that file defines target to calculate the Field                 |"
    print "|Targets are selected using the amber selection syntax                                                        |"
    print "|To calculate the magnitude of the field at a specific atom:     :47@OG                                       |"
    print "|To calculate the magnitude of the field projected along a bond: :47@OG :47@OH                                |"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|-------------                                       INFO                                        -------------|"
    print "|-------------------------------------------------------------------------------------------------------------|"
    print "|Script recognizes if total field at one point or projection onto a bond is to be calculated                  |"
    print "|based on the input in -target                                                                                |"
    print "|-------------                                                                                                |"
    print "|Script excludes all Field contributions of the Residue at which the field is calculated.                     |"
    print "|For Vector Projections, this is only the first residue number!                                               |"
    print "|-------------                                                                                                |"
    print "|Solvent = All atoms as defined in -solvent,                                                                  |"
    print "|          for these residues contributions will also be grouped based on the residue name                    |"
    print "|Protein = All atoms that are not part of -solvent,                                                           |"
    print "|          for the protein, per-reside fields will be calculated                                              |"
    print "|-------------                                                                                                |"
    print "|Linker Atom charge is ignored!                                                                               |"
    print "|-------------------------------------------------------------------------------------------------------------|" 
    
def inputParser(arg):
    
    arg_nc            = False                                                                         ### Add new flag here
    arg_param         = False
    arg_out           = False
    arg_tmp_dir       ="tmp_field/"   
    arg_target        = False
    arg_solvent       = "WAT,Na+,Cl-"
    arg_exclude_atoms = False
    arg_use_qmcharges = False
    arg_qm_mask       = False
    arg_qm_charges    = False
    arg_qm_dict       = False
    
    if(len(arg) == 1) or arg[1] in ["-help","--help","-h"]:
        usage()
        quit()
    
    for index in range(1,len(arg)-1,2):
        param = arg[index]
        if param=="-nc"             : arg_nc            = os.getcwd()+"/"+str(arg[index+1])           ### Add new flag here
        if param=="-parm"           : arg_param         = os.getcwd()+"/"+str(arg[index+1])
        if param=="-out"            : arg_out           = os.getcwd()+"/"+str(arg[index+1])
        if param=="-tmp_dir"        : arg_tmp_dir       = os.getcwd()+"/"+str(arg[index+1])+"/"
        if param=="-target"         : arg_target        = os.getcwd()+"/"+str(arg[index+1])
        if param=="-solvent"        : arg_solvent       = str(arg[index+1])
        if param=="-exclude_atoms"  : arg_exclude_atoms = os.getcwd()+"/"+str(arg[index+1])
        if param=="-use_qm_charges" : arg_use_qmcharges = str(arg[index+1])
        if param=="-qm_mask"        : arg_qm_mask       = os.getcwd()+"/"+str(arg[index+1])
        if param=="-qm_charges"     : arg_qm_charges    = os.getcwd()+"/"+str(arg[index+1])
        if param=="-qm_dict"        : arg_qm_dict       = os.getcwd()+"/"+str(arg[index+1])

                                                                                                      ### Add new flag here
    if arg_nc        == False: print "ERROR. Please specify the name of the trajectory file with -nc"  
    if arg_param     == False: print "ERROR. Please specify the name of the parameter file with -parm"
    if arg_out       == False: print "ERROR. Please specify the name of the output file with -out"
    if arg_target    == False: print "ERROR. Please specify the name of the target file with -target.                                       \n -> See examples with --help"
    if arg_use_qmcharges == "True":
        if arg_qm_mask    == False: print "ERROR. Please specify qm mask with -qm_mask"
        if arg_qm_charges == False: print "ERROR. Please specify charge of the qm region with -qm_charge"
        if arg_qm_dict    == False: print "ERROR. Please specify qm_dict with -qm_dict"
        
    if False in [arg_nc,arg_param,arg_out,arg_tmp_dir,arg_target]:                                   ### Add new flag here
        usage()
        quit()
    if arg_use_qmcharges == "True":
        if False in [arg_qm_mask,arg_qm_charges,arg_qm_dict]:                                        ### Add new flag here
            usage()
            quit()
    
    print "-nc              Trajectory file      : ",arg_nc                                          ### Add new flag here
    print "-parm            Parameter file       : ",arg_param
    print "-out             Output file          : ",arg_out
    print "-tmp_dir         Temporary dir        : ",arg_tmp_dir
    print "-target          Field target         : ",arg_target
    print "-arg_solvent     Non-protein residues : ",arg_solvent
    if arg_exclude_atoms <> False:
        print "-exclude_atoms   Atoms excluded       : ",arg_exclude_atoms
    else:
        print "-exclude_atoms   Atoms excluded       : All atoms of the residue defined in the first Atom for each Target."
    if arg_use_qmcharges == "True":
        print "-use_qm_charges  QM_theory            :  True"
        print "-qm_mask         QM region            : ",arg_qm_mask
        print "-qm_charges      QM charges           : ",arg_qm_charges
        print "-qm_dict         QM dict              : ",arg_qm_dict 
    
    return arg_nc,arg_param,arg_out,arg_tmp_dir,arg_target,arg_solvent,arg_exclude_atoms,           arg_use_qmcharges,arg_qm_mask,arg_qm_charges,arg_qm_dict                                  ### Add new flag here
    
########################
### THE FINAL SCRIPT ###
########################
### Read in arguments                                                                             
arg_nc,arg_param,arg_out,arg_tmp_dir,arg_target,arg_solvent,arg_exclude_atoms,           arg_use_qmcharges,arg_qm_mask,arg_qm_charges,arg_qm_dict=inputParser(sys.argv)            ### Add new flag here
arg_solvent = arg_solvent.split(",")

print "\nField calculation RUNNING : ",datetime.datetime.now() 

#Load qm_mask.in
arg_qm_mask = load_qm_mask(arg_qm_mask)

### Empty dictionary storing all Fields
FIELDS = {}

### Make tmpdir
make_tmpdir(arg_tmp_dir)

### Make pdb of system and qm_mm region
arg_pdb = ".".join(arg_out.split(".")[:-1])+"_full.pdb"
make_pdb(arg_tmp_dir,arg_nc,arg_param,arg_pdb)

### Load QM charges
if arg_use_qmcharges == "True":
    QM_Charges,QM_Dict = fkt_Load_QM_Charges(arg_qm_charges,arg_qm_dict) 
    #print "QM_Dict :", QM_Dict,"\n"      
    #print "QM_Charges :", QM_Charges[:50],"\n"
     
### Load Trajectory
Trajectory, Charges, Names = fkt_Load_Trajectory(arg_nc,arg_param,arg_pdb)
#print "Trajectory :", Trajectory,"\n"
#print "Charges :", Charges,"\n"
#print "Names :", Names,"\n"
 
### Load Index of Target Atoms for Field calculation
Target_Index,FIELDS = fkt_Target_Index(Names,arg_target,FIELDS)
#print "Target_Index:", Target_Index,"\n"
#print "Field_Dict:", FIELDS,"\n"

### Load Index of Excluded Atoms
if arg_exclude_atoms <> False: arg_exclude_atoms,_ = fkt_Target_Index(Names,arg_exclude_atoms,FIELDS)
#print "Excluded_Atoms:", arg_exclude_atoms,"\n"

### Define Field Components (Total,Protein,Solvent,Each residue,Each solvent type) to calculate
FIELDS,Field_Components = fkt_Field_Components(FIELDS,Names,arg_solvent)
#print "Field_Dict:", FIELDS,"\n"

### Loop through Trajectory to calculate Fields
for Frame_i in range(0,len(Trajectory),1):

### Get XYZ coordinates of the Target Atoms
    Target_XYZ=fkt_get_Target_XYZ(Trajectory[Frame_i],Target_Index)
    #print "Target_XYZ:",Target_XYZ,"\n"
    Target_VEC=fkt_get_Target_VEC(Target_XYZ)
    #print "Target_VEC:",Target_VEC,"\n"    
### Update QM_Charges
    if arg_use_qmcharges == "True":
        Charges = fkt_Update_Charges(QM_Charges,QM_Dict,Frame_i,Names,Charges)  
        #print "Charges :", Charges,"\n"    
        
### Loop through FIELDS to calculate each Target
    for FIELD in FIELDS:
        if len(FIELD.split("_")) == 1:   #POINT CALCULATION
            XYZ = Target_XYZ[FIELD][0]
            VEC = ""
        if len(FIELD.split("_")) == 2:   #VECTOR CALCULATION
            XYZ = Target_VEC[FIELD]['Center']
            VEC = Target_VEC[FIELD]['Vector']
        Self_i = FIELD.split("@")[0][1:]
        
### Calculate Field Vectors at position XZY        
        FIELDS = fkt_calc_fields(FIELDS,FIELD,XYZ,VEC,Trajectory[Frame_i],Charges,Names,                                 Field_Components,arg_solvent,Self_i,arg_exclude_atoms)
        #for Field_Component in FIELDS[FIELD]:
        #    print Field,Field_Component,FIELDS[FIELD][Field_Component]    
        
print "\nField calculation DONE : ",datetime.datetime.now()         
pickle.dump(FIELDS,open(arg_out,"wb"))

