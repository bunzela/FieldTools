# FieldTools

FieldTools.py can be used to calculate electric fields from MD trajectories. The script requires standard MD **trajectory** and **parameter** files as input.
Furthermore a **target** file needs to be provided that specifies the positions at which the field will be calculated.

To test field tools, click on: <a target="_blank" href="https://colab.research.google.com/github/bunzela/FieldTools/blob/main/FieldTools.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>
---

> [!NOTE]  
> Fields are defined in the **target** file and are either calcualted at an atom or along a bond. 
> To define a target, use amber selection masks to select either one (atom) or two (bond) atoms. 
> Several targets can be calculated in parallel, by adding additional lines to the **target** file

> [!WARNING]  
> FieldTools calculates the fields from the exact location of all atoms in the system without considering periodicity. 
> Trajectories must thus be imaged!

### Requirements
- Python3.X
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

> [!NOTE]  
> FieldTools can also calculate QM/MM point charges for a more refined field analysis (still experimental!). <br />
> Contact [adrian.bunzel@bsse.ethz.ch](mailto:adrian.bunzel@bsse.ethz.ch) for early access.

### Citation
Please cite the following paper when using FieldTools:
**H. Jabben et al., bioRxiv 2023**. 

> [!NOTE]  
> Compare to the published version, the FieldTools.py script provided here uses MDtraj instead of pytraj, because pytraj cannot readily installed in Google Colab

### Contact
For questions, help, or to report any bugs, please feel free to reach out to me at [adrian.bunzel@bsse.ethz.ch](mailto:adrian.bunzel@bsse.ethz.ch).