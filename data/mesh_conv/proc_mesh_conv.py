import sys
import os
home_dir = os.path.expanduser('~')
sys.path.append(os.path.join(home_dir, 'data', 'tools'))
import sparc_tools as sparc
import numpy as np


# meshes = np.linspace(start=0.10, stop=0.90, num = 17)

# nstates = range(50,310,10)

# prefix of subfolders
prefix = "mesh_"
len_prfx = len(prefix)


# folders = []
# for folder in os.listdir('.'):
#     if folder.startswith(prefix):
#         folders.append(folder) # note that folders might not be sorted
#         # print(folder)
#         
# folders = sorted(folders,key = lambda x: float(x[len_prfx:])) # sort by the trailing numbers

folders = ["mesh_%.2f"%i for i in np.arange(0.20,0.81,0.05)]
# append ref folder to the end if necessary
# folders.append("mesh_0.10")
#print(folders)

folder_names = []
eatoms = []
forces = []
pressures = []
stresses = []
nscfs = []
wtimes = []
for folder_name in folders:
    # folder_name = "./ns_%.0f"%ns
    # static_fname = 'force.txt'
    static_fname = 'sprc-calc.static'
    out_fname = 'sprc-calc.out'
    # print(folder_name)

    os.chdir(folder_name)
    
    for filename in os.listdir('.'):
        if filename.endswith('.static'):
            static_fname = filename
        elif filename.endswith('.out'):
            out_fname = filename
        elif filename.endswith('.aimd'):
            aimd_fname = filename
   
    #print("filename: " + out_fname)

    # grep natom
    natom = sparc.grep_natom(out_fname)

    # grep energy
    eatom_i = sparc.grep_energy(out_fname)
    
    if not eatom_i:
        print("skipping " + folder_name)
        os.chdir("..")
        continue # skip current folder


    # grep forces
    #forces_i = sparc.read_forces_from_text_numpy(static_fname)
    forces_i = sparc.grep_forces(static_fname,natom)
    
    if not forces_i:
        print("skipping " + folder_name)
        os.chdir("..")
        continue # skip current folder
    
    # grep pressure
    press_i = sparc.grep_pressure(out_fname)
    #press_i = -1
    
    # grep stress
    stress_i = sparc.grep_stress(static_fname)


    folder_names.append(folder_name)
    eatoms.append(eatom_i)
    forces.append(np.array(forces_i))
    pressures.append(press_i)
    stresses.append(np.array(stress_i))

    # grep nscf
    nscf = sparc.grep_nscf(out_fname)
    nscfs.append(nscf)

    # grep walltime
    wtime = sparc.grep_walltime(out_fname)
    wtimes.append(wtime)

    os.chdir("..")


ref_ind = 0
eatoms_ref = eatoms[ref_ind]
forces_ref = forces[ref_ind]
press_ref = pressures[ref_ind]
stress_ref = stresses[ref_ind]

# abinit reference
# eatoms_ref = -5.5073192860E+01 / natom
# forces_ref = np.array([
# [  3.2286433553E-03, -1.8745295425E-02, -3.7631097838E-03],
# [ -3.1021405602E-02, -4.4746182584E-02,  1.1924336279E-02],
# [  3.6956131339E-02,  3.5929933012E-02,  8.7783496033E-03],
# [ -1.1845618846E-03,  3.6762017936E-02, -3.2821471852E-02],
# [  8.5161987562E-03,  3.2833635187E-02,  2.4946441398E-02],
# [  7.3458157768E-03, -1.6846471958E-02,  2.1522543270E-02],
# [ -8.7608157378E-03,  1.9617565978E-02, -1.7548310092E-02],
# [  1.3445288285E-02,  7.9631567026E-03,  7.9782938140E-03],
# [ -4.7906337474E-04, -9.8365483820E-03, -1.8216304127E-02],
# [ -1.2489898893E-02,  1.4003600746E-03,  6.7650295353E-03],
# [ -3.2121250775E-03, -4.7699996554E-03,  1.2249554607E-03],
# [  8.2594177533E-04,  3.7605724913E-03, -1.2785250535E-02],
# [  1.2807125506E-02, -4.9525051392E-03,  8.7809907080E-03],
# [ -4.3395500988E-02, -4.8162004828E-02,  8.4488268457E-03],
# [  5.6738538038E-02,  4.8711683227E-02, -1.3488372974E-02],
# [ -2.3726613449E-02, -5.7667073491E-03,  5.2064281008E-03],
# [ -8.9274007520E-03, -3.8623089744E-02,  4.6689790158E-02],
# [ -4.6066509143E-03, -4.9489701282E-03, -6.8495150402E-03],
# [ -2.8629480918E-04,  9.3550851448E-04,  9.7794602832E-04],
# [  1.6332432514E-02, -6.6789513878E-03, -3.9485054074E-05],
# [ -2.9146962462E-02,  2.2500699424E-02, -2.9291462571E-02],
# [  8.3120878099E-03, -2.6573274263E-03, -2.1074834834E-03],
# [  9.2351923161E-03, -2.7283518097E-02, -3.1393304328E-02],
# [ -6.5061015282E-03,  2.3602439558E-02,  1.5060138638E-02],
# ])
# 
# press_ref = 9.2200
# #- sigma(1 1)= -1.13115290E+01  sigma(3 2)=  6.32024746E-01
# #- sigma(2 2)= -1.84012877E+01  sigma(3 1)= -9.66045154E-01
# #- sigma(3 3)=  2.05269562E+00  sigma(2 1)= -5.80275184E-01
# stress_ref = np.array([[-1.13115290E+01, -5.80275184E-01, -9.66045154E-01],
#                        [-5.80275184E-01, -1.84012877E+01,  6.32024746E-01],
#                        [-9.66045154E-01,  6.32024746E-01,  2.05269562E+00]])


eatom_diffs = []
force_diffs = []
press_diffs = []
stress_diffs = []
for i in range(len(eatoms)):
    forces_i = forces[i]
    eatoms_i = eatoms[i]
    press_i = pressures[i]
    stress_i = stresses[i]
    nscfs_i = nscfs[i]
    eatom_diffs.append(abs(eatoms_i - eatoms_ref))
    force_diffs.append(np.max(np.abs(forces_i - forces_ref)))
    press_diffs.append(abs((press_i - press_ref)/press_ref))
    stress_diffs.append(np.max(np.abs(stress_i - stress_ref)))

print("%7s   eatom (Ha/atom)  deatom   dforce      press (GPa)   dp/p0     dStress (GPa)   nscf   walltime (s)"%prefix)
#for fdr,eatom,eatom_d,frc_d,pres,pres_d,nscf in zip(folders,eatoms,eatom_diffs,force_diffs,pressures,press_diffs,nscfs):
for i in range(len(eatoms)):
    fdr = folder_names[i]
    #ns = float(fdr[len_prfx:])
    ns = (fdr[len_prfx:])
    eatom = eatoms[i]
    deatom = eatom_diffs[i]
    df = force_diffs[i]
    pres = pressures[i]
    dp_p0 = press_diffs[i]
    dstress = stress_diffs[i]
    nscf = nscfs[i]
    wtime = wtimes[i]
    print("%7s %15.10f   %7.2e  %7.2e %15.10f  %7.2e     %7.2e   %4d     %.3f"%(ns,eatom,deatom,df,pres,dp_p0,dstress,nscf,wtime))


# for d in force_diffs:
#    print("%.2e"%d)


