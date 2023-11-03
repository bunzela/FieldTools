# FieldTools

FieldTools calculates electric fields in MD trajectories. The script requires a trajectory and parameter file as input. A -target file needs to be provided that specifies the positions at which the field should be calculated.

To test field tools, click on: <a target="_blank" href="https://colab.research.google.com/github/bunzela/FieldTools/blob/main/FieldTools.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>
---

### Requirements
- python3.X
- MDanalysis (install using pip install mdanalysis)


### Usage
python utils/FieldTools.py -nc <trajectory file> 
                           -parm <parameter file> 
                           -target <target file> 
                           [-solvent <non-protein residues>]
                           [-exclude_atoms <exclusion list>] 
                           [-TIP4P <True|False>] 
                           [-verbose <True|False>] 
                           -out <output file> 

#### FieldTools can calculate QM/MM point charges for a more refined field calcualtion (still experimental!). Please contact [adrian.bunzel@bsse.ethz.ch](mailto:adrian.bunzel@bsse.ethz.ch) fur early access.

### Citation
Please cite the following paper when using FieldTools:
H. Jabben et al., bioRxiv 2023. 

### Contact
For questions, help, or to report any bugs, please feel free to reach out to me at [adrian.bunzel@bsse.ethz.ch](mailto:adrian.bunzel@bsse.ethz.ch).