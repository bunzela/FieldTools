import sys
import os
import subprocess
import datetime

def fkt_submit(arg_tmp_dir,arg_qm_submit):
    if arg_qm_submit == "False": return
    with open(arg_tmp_dir+"submit.sh", "w") as f:
        f.write(arg_qm_submit+" "+arg_tmp_dir+"submission.in")
    print ("QMMM calculation stated with "+arg_qm_submit+" "+arg_tmp_dir+"submission.in")
    with open(arg_tmp_dir+"submit.out", 'w') as f:
        process = subprocess.call(['bash', arg_tmp_dir+"submit.sh"], stdout=f)
    
def make_qm_dict(arg_out):
    with open('.'.join(arg_out.split('.')[:-1])+'_qm_region.pdb') as f:
        pdb = [i.split() for i in f.readlines()]
    pdb = [i[2:5] for i in pdb if i[0] in ["ATOM","HETATM"]]
    with open('.'.join(arg_out.split('.')[:-1])+'_qm.dict', "w") as f:
        f.write("\n".join(["_".join(i) for i in pdb]))
        
def convert_charges_script(arg_out,arg_tmp_dir):
    rstfiles = [i for i in os.listdir(arg_tmp_dir) if i.split(".")[-2] == "rst"]
    script = """#Read in gaussian log and extract only charges into file
QM_Charges=[]
for i in range(1,"""+str(len(rstfiles)+1)+""",1):
    with open('"""+arg_tmp_dir+"""gaussian_'+str(i)+'.log') as f:
        charges = f.readlines()

    for out_i in range(0,len(charges),1):
        if charges[out_i][:13] in ' ESP charges:':
            i_start=out_i+2                
        if charges[out_i][:19]==' Sum of ESP charges':
            i_stop=out_i
            QM_Charges=QM_Charges+['Frame '+str(i)+'\\n']
            QM_Charges=QM_Charges+charges[i_start:i_stop]
        
with open('"""+arg_out+"""', 'w') as f:
    f.writelines(QM_Charges) 
"""
    with open(arg_tmp_dir+"convert_charges.py", "w") as f:
        f.write(script)
        
def make_submission_script(arg_tmp_dir,arg_param,arg_out,arg_qm_theory,arg_rm_temp,arg_qm_submit):
    rstfiles = [i for i in os.listdir(arg_tmp_dir) if i.split(".")[-2] == "rst"]
    if arg_qm_submit == "False": 
        submission=""
    if arg_qm_submit == "sbatch":
        submission="""#!/bin/env bash
#SBATCH --job-name=Chrg
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:0
#SBATCH --mem=10000M

# 1. Load module(s)
module load  apps/amber18/tools19-packmol-mpi-gnu-7.2.0
module load apps/gaussian/16
"""        
    submission+="""
# 2. Set directories
cd """+arg_tmp_dir+"""
rm *.o* *.e* 2>1
cp tmp.rst tmp.rst.1 2>1
"""
    for i in range(1,len(rstfiles)+1,1):
        submission += """
echo Frame """+str(i)+""" gaussian STARTED >> submission.out
sander -O -i amber.in -p """+arg_param+""" -o sander.out -c tmp.rst."""+str(i)+""" &>2
cat  gaussian.in                       > gaussian_"""+str(i)+""".gjf
tail -n +5 gau_job.inp                >> gaussian_"""+str(i)+""".gjf
echo ""                               >> gaussian_"""+str(i)+""".gjf
g16  gaussian_"""+str(i)+""".gjf      >> gaussian_"""+str(i)+""".log
echo Frame """+str(i)+""" gaussian DONE >> submission.out
"""
    submission += """
python """+arg_tmp_dir+"""convert_charges.py """+arg_out+"""
"""
    if arg_rm_temp == "True": submission += """
rm *
"""
    with open(arg_tmp_dir+"submission.in", "w") as f:
        f.write(submission)
        
def make_input_files(arg_tmp_dir,arg_qm_mask,arg_qm_charge,arg_qm_theory,arg_out,arg_param,arg_nc):
    ### Input file for sander
    amber = """Create Gaussian input file. Designed to crash sander after files are created
&cntrl
    ntpr = 1, ntwx = 0,
    imin = 1, maxcyc = 0,           !Single-point energy calculation
    ntb = 0,                        !Non-periodic
    cut = 9999.,                    !Calculate all interactions
    ifqnt = 1,                      !Switch on QM/MM coupled potential
/
&qmmm 
    qmmask = '"""+arg_qm_mask+"""', 
    qmcharge = """+arg_qm_charge+""", 
    spin = 1, 
    qm_theory = 'EXTERN',
    qmcut = 999.0,
    itrmax = 10000000,
    printcharges = 1, 
&end     
&gau
    method = 'XXX',  
    basis = 'XXX',    
&end        
"""
    with open(arg_tmp_dir+"amber.in", "w") as f:
        f.write(amber) 
        
    ### Input header for gaussian
    gaussian = """%chk="""+arg_tmp_dir+"""gaussian.chk 
%NProcShared=1 
%mem=8192MB 
# """+arg_qm_theory+""" SCF=(Conver=8,verytight,Maxcyc=1000) Integral=Ultrafine NoSymm POP=CHelpG Charge   
"""
    with open(arg_tmp_dir+"gaussian.in", "w") as f:
        f.write(gaussian) 

    cpptraj = """parm """+arg_param+"""
trajin """+arg_nc+"""
autoimage    
outtraj """+arg_tmp_dir+"""tmp.rst multi nobox
outtraj """+arg_tmp_dir+"""full.pdb onlyframes -1
strip !("""+arg_qm_mask+""")
outtraj """+".".join(arg_out.split(".")[:-1])+"""_qm_region.pdb onlyframes -1
"""
    with open(arg_tmp_dir+"cpptraj.in", "w") as f:
        f.write(cpptraj)  
        
def make_rst(arg_tmp_dir,arg_nc,arg_param):
    with open(arg_tmp_dir+"/cppraj.out", 'w') as f:
        process = subprocess.run(['cpptraj', '-i', arg_tmp_dir+'cpptraj.in'], stdout=f)
    if os.path.isfile(arg_tmp_dir+"/tmp.rst"):
        with open(arg_tmp_dir+"/tmp.rst", "r") as f:
            temp_rst = f.read()
        with open(arg_tmp_dir+"/tmp.1.rst", "w") as f:
            f.write(temp_rst)          
    print ("Cpptraj DONE. Trajectory "+arg_nc+" split into "+arg_tmp_dir+"tmp.rst.X")
        
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
    print ("                                                                                                               ")
    print ("|-------------------------------------------------------------------------------------------------------------|")
    print ("|----------                                  QMChargesTools Usage:                                  ----------|")
    print ("|-------------------------------------------------------------------------------------------------------------|")
    print ("|-nc            : Specify trajectory file                                                                     |")
    print ("|-parm          : Specify parameter file                                                                      |")
    print ("|-out           : Specify output file                                                                         |")
    print ("|-tmp_dir       : Optional: Specify temporary directory                           [default ./tmp_chrg/       ]|")
    print ("|-qm_mask       : Specify qm region (amber selection mask)                                                    |")
    print ("|-qm_charge     : Specify charge of qm region                                                                 |")
    print ("|-qm_theory     : Optional: Specify qm theory                                     [default M062X/6-31++g(d,p)]|")
    print ("|-submit        : Optional: Define how calcuations should be run                  [default bash              ]|")
    print ("|               :           i.e submission command (qsub, sbatch, ...)                                        |")
    print ("|               :           i.e run in terminal (bash)                                                        |")
    print ("|               :           i.e don't run (False)                                                             |")
    print ("|-rm_temp       : Optional: Remove temporary data                                 [default True              ]|")
    print ("|-------------------------------------------------------------------------------------------------------------|")
    
def inputParser(arg):
            
    arg_nc        = False                                                                             ### Add new flag here
    arg_param     = False
    arg_out       = False
    arg_tmp_dir   = "tmp_chrg" 
    arg_qm_mask   = False
    arg_qm_charge = False
    arg_qm_theory ="M062X/6-31++g(d,p)"
    arg_qm_submit = "bash"
    arg_rm_temp = "True"
    
    if(len(arg) == 1) or arg[1] in ["-help","--help","-h"]:
        usage()
        quit()
    
    for index in range(1,len(arg)-1,2):
        param = arg[index]
        if param=="-nc"           : arg_nc        = os.getcwd()+"/"+str(arg[index+1])                 ### Add new flag here
        elif param=="-parm"       : arg_param     = os.getcwd()+"/"+str(arg[index+1])
        elif param=="-out"        : arg_out       = os.getcwd()+"/"+str(arg[index+1])
        elif param=="-tmp_dir"    : arg_tmp_dir   = os.getcwd()+"/"+str(arg[index+1])+"/"
        elif param=="-qm_mask"    : arg_qm_mask   = os.getcwd()+"/"+str(arg[index+1])
        elif param=="-qm_charge"  : arg_qm_charge = str(arg[index+1])  
        elif param=="-qm_theory"  : arg_qm_theory = str(arg[index+1])
        elif param=="-submit"     : arg_qm_submit = str(arg[index+1])
        elif param=="-rm_temp"    : arg_rm_temp   = str(arg[index+1])
        else:
            print ("Error! "+param+" not a valid option.")
            usage()
            quit()     
 
    if arg_nc        == False: print ("ERROR. Please specify the name of trajectory file with -nc")     ### Add new flag here
    if arg_param     == False: print ("ERROR. Please specify the name of parameter file with -parm")
    if arg_out       == False: print ("ERROR. Please specify the name of output file with -out")
    if arg_qm_mask   == False: print ("ERROR. Please specify qm mask with -qm_mask")
    if arg_qm_charge == False: print ("ERROR. Please specify charge of the qm region with -qm_charge")
        
    if False in [arg_nc,arg_param,arg_out,arg_qm_mask,arg_qm_charge]:                                 ### Add new flag here
        usage()
        quit()
    
    print ("-nc        Trajectory file : ",arg_nc)
    print ("-parm      Parameter file  : ",arg_param)
    print ("-out       Output file     : ",arg_out)
    print ("-tmp_dir   Temporary dir   : ",arg_tmp_dir)
    print ("-qm_mask   QM_mask         : ",arg_qm_mask)
    print ("-qm_charge QM_charge       : ",arg_qm_charge)
    print ("-qm_theory QM_theory       : ",arg_qm_theory)
    print ("-qm_submit Submit          : ",arg_qm_submit)
    print ("-rm_temp   rm tmp data     : ",arg_rm_temp)
    print ()
                                                                                                      ### Add new flag here
    return arg_nc,arg_param,arg_out,arg_tmp_dir,arg_qm_mask,arg_qm_charge,arg_qm_theory,arg_qm_submit,arg_rm_temp  
 
def main():
    ########################
    ### THE FINAL SCRIPT ###
    ########################
    #Read in arguments                                                                                ### Add new flag here
    arg_nc,arg_param,arg_out,arg_tmp_dir,arg_qm_mask,arg_qm_charge,arg_qm_theory,arg_qm_submit,arg_rm_temp    =inputParser(sys.argv) 
    print ("\nQM/MM calculation RUNNING : ",datetime.datetime.now() )
    #Make tmpdir
    make_tmpdir(arg_tmp_dir)
    #Load qm_mask.in
    arg_qm_mask = load_qm_mask(arg_qm_mask)
    #Prepare .in files
    make_input_files(arg_tmp_dir,arg_qm_mask,arg_qm_charge,arg_qm_theory,arg_out,arg_param,arg_nc)
    #Split .nc into .rst files
    make_rst(arg_tmp_dir,arg_nc,arg_param)
    #Make Scripts
    convert_charges_script(arg_out,arg_tmp_dir)
    make_submission_script(arg_tmp_dir,arg_param,arg_out,arg_qm_theory,arg_rm_temp,arg_qm_submit)
    make_qm_dict(arg_out)
    #Submit
    fkt_submit(arg_tmp_dir,arg_qm_submit)
    print ("\nQM/MM calculation DONE using : ",arg_qm_submit, datetime.datetime.now()  )
    
if __name__ == "__main__":
    main()
