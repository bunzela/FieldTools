# FieldTools

FieldTools calculates electric fields in MD trajectories. The script requires standard MD **trajectory** and **parameter** files as input.
Furthermore a **target** file needs to be provided that specifies the positions at which the field will be calculated.

<font color="red">**Attention: Trajectories must be immaged!**</font>

<font color="red">FieldTools calculates the fields from the exact location of all atoms in the simulation without considiring periodicity.</font>

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

**Note:** FieldTools can also calculate QM/MM point charges for a more refined field analysis (still experimental!).

Please contact [adrian.bunzel@bsse.ethz.ch](mailto:adrian.bunzel@bsse.ethz.ch) fur early access.

### Citation
Please cite the following paper when using FieldTools:
H. Jabben et al., bioRxiv 2023. 

### Contact
For questions, help, or to report any bugs, please feel free to reach out to me at [adrian.bunzel@bsse.ethz.ch](mailto:adrian.bunzel@bsse.ethz.ch).