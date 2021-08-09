# FieldTools
FieldTools calculated electric Fields in MD trajectories. It requires an Amber trajectory and parameter file as input. FieldTools can optionally read in user-defined charges, e.g. from QM/MM calculations. QMChargesTools can be used to calculate ChelpG charges for QM regions over the trajectory. 

An example trajectory containing only two frames is provided in ./4bs0/
The code is based on Stuyver et al., J Comp Chem 2019.

## Content
QMChargesTools
FieldTools

## QMChargesTools
Helpfile can be read with -h -help --help or without providing any flags
```
!python utils/QMChargesTools.py -help
```
QMChargesTools requires a trajectory and a parameter file as input. Furthermore, you need to specifiy the output file at which the charges are saved.

Options for the QM method are provided with -qm_mask, -qm_charge and -qm_theory.

The script creates a submission file in the temporary directory. Calculations will be started by running bash submit.sh, which runs submit.in using the option defined in -submit. With -submit, you can control how (and if) the calculations are executed.

The next cell runs QMCHargesTools. Note that the recommended level of theory for QM/MM is M062X/6-31++g(d,p).
```
!python QMChargesTools.py -nc 4bs0/4bs0.nc -parm 4bs0/4bs0.parm7 -out 4bs0/4bs0.charge \
                          -tmp_dir 4bs0/tmp_chrg  \
                          -qm_charge 0 -qm_mask  4bs0/qm_mask.in -qm_theory lsda/sto-3g -submit bash
```
