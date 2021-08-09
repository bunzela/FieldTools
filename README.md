# FieldTools
FieldTools calculated electric Fields in MD trajectories. It requires an Amber trajectory and parameter file as input. FieldTools can optionally read in user-defined charges, e.g. from QM/MM calculations. QMChargesTools can be used to calculate ChelpG charges for QM regions over the trajectory. 

An example trajectory containing only two frames is provided in ./4bs0/
The code to calculted the magnitued of the field vector projected along another vector was taken from Stuyver et al., J Comp Chem 2019. 

## Content
- QMChargesTools
- FieldTools
- Example data plotting

## QMChargesTools
The helpfile can be read with -h -help --help or without providing any flags
```
!python utils/QMChargesTools.py -help
```
QMChargesTools requires a trajectory and a parameter file as input. Furthermore, you need to specifiy the output file at which the charges are saved.

Options for the QM method are provided with -qm_mask, -qm_charge and -qm_theory.

The script creates a submission file in the temporary directory. Calculations will be started by running bash submit.sh, which runs submit.in using the option defined in -submit. With -submit, you can control how (and if) the calculations are executed.

The next cell runs QMCHargesTools. Note that the recommended level of theory for QM/MM is M062X/6-31++g(d,p).
```
!python utils/QMChargesTools.py -nc 4bs0/4bs0.nc -parm 4bs0/4bs0.parm7 -out 4bs0/4bs0.charge \
                                -tmp_dir 4bs0/tmp_chrg  \
                                -qm_charge 0 -qm_mask  4bs0/qm_mask.in -qm_theory lsda/sto-3g -submit bash
```

## FieldTools
The helpfile can be read with -h -help --help or without providing any flags
```
!python utils/QMChargesTools.py -help
```

```
!python utils/FieldTools.py -nc 4bs0/4bs0.nc -parm 4bs0/4bs0.parm7 -out 4bs0/4bs0_field.pkl \
                            -tmp_dir 4bs0/tmp_field -exclude_atoms 4bs0/field_target.dat \
                            -target 4bs0/field_target.dat -solvent WAT,5NO -use_qm_charges False \
                            -qm_mask  4bs0/qm_mask.in -qm_charges 4bs0/4bs0.charge -qm_dict 4bs0/4bs0.dict
```

## Example data plotting
```
import pickle
import matplotlib.pyplot as plt 
import numpy as np

FIELDS = pickle.load(open("4bs0/4bs0_field.pkl","rb"))
plt.rcParams["figure.figsize"]=[16,3]
for Target in FIELDS:
    x_label = sorted(FIELDS[Target])
    x = range(1,301,1)
    y = []
    for x_i in x:
        y+=[np.mean(FIELDS[Target][str(x_i)])]
    plt.bar(x,y)
    plt.title(Target)
    plt.xlabel("Residue Number")
    plt.ylabel("Field (MV/cm)")
    plt.xlim(0,302)
    plt.ylim(-250,250)
    plt.plot([0,302],[0,0],c="black",linewidth=1.0)
    plt.xticks(range(1,302,20))
    plt.show()
```
