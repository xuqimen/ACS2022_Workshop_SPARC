import sys
import re
import os
import numpy as np
import pandas as pd

def grep_out_info(fname):
    """
    read general info from .out, including:
    - cell sizes: cellsizes
    - lattice vectors: latvecs
    - meshes: hx, hy, hz
    - number of atoms: natom
    - energy: etot, eatom = etot/natom
    - pressure: pressure (not always exist)
    - wall time: walltime
    """
    info_dict = {}
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Free energy per atom .+: (-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        eatom = float(match)
        info_dict['eatom'] = eatom

    return info_dict



def grep_walltime(fname):
    """
    read walltime from SPARC .out file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Total walltime.+: (-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            return float(match.group(1))
        else:
            print("Energy not found!")



def grep_energy(fname):
    """
    read energy per atom from SPARC .out file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Free energy per atom .+: (-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            energy = float(match.group(1))
            return energy
        else:
            print("Energy not found!")


def grep_pressure(fname):
    """
    read pressure from SPARC .out file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Pressure .+: (-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            pressure = float(match.group(1))
            return pressure
        else:
            print("Pressure not found!")


def grep_natom(fname):
    """
    read number of atoms from SPARC .out file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Total number of atoms .+: (-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            natom = int(match.group(1))
            return natom
        else:
            print("Number of atoms not found!")
 
 
def grep_nscf(fname):
    """
    read number of SCFs from SPARC .out file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"number of SCF: (-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            nscf = int(match.group(1))
            return nscf
        else:
            print("Number of SCF not found!")



def grep_stress(fname):
    """
    read stress from SPARC .static file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    stress = []
    read_flag = False
    stress_count = 0
    with open(fname, 'r') as f:
        for line in f:
            if read_flag:
                stress_i = [float(x) for x in line.split()]
                stress.append(stress_i)
                stress_count += 1
                if stress_count == 3:
                    break
            elif re.search(r"Stress.*GPa.*", line):
                read_flag = True

    if not read_flag:
        pass
        #print("stress not found!")

    return stress



def grep_forces(fname, natom):
    """
    read pressure from SPARC .static file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    forces = []
    read_flag = False
    force_count = 0
    with open(fname, 'r') as f:
        for line in f:
            if read_flag:
                force_i = [float(x) for x in line.split()]
                forces.append(force_i)
                force_count += 1
                if force_count == natom:
                    break
            elif re.search(r"Atomic forces", line):
                read_flag = True

    if not read_flag:
        print("Forces not found!")

    return forces



def grep_forces_aimd(fname, natom):
    """
    read pressure from SPARC .aimd file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    forces = []
    read_flag = False
    force_count = 0
    with open(fname, 'r') as f:
        for line in f:
            if read_flag:
                force_i = [float(x) for x in line.split()]
                forces.append(force_i)
                force_count += 1
                if force_count == natom:
                    break
            elif re.search(r":F:", line):
                read_flag = True

    if not read_flag:
        print("Forces not found!")

    return forces




def read_number_from_text(fname):
    """
    read one number from .txt file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    with open(fname) as f:
        return float(f.read())



def read_forces_from_text(fname):
    """
    read force matrix from .txt file, without numpy
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    forces = []
    with open(fname, 'r') as f:
        for line in f:
            force_i = [float(x) for x in line.split()]
            forces.append(force_i)
            
    return forces


def read_forces_from_text_numpy(fname):
    """
    read force matrix from .txt file, using numpy
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    forces = np.loadtxt(fname)
    return forces



def first_match_line(string, filename):
    """
    This function finds the first line in filename that contains
    the given string, and returns the line number
    """
    with open(filename) as myFile:
        for num, line in enumerate(myFile, 1):
            if string in line:
                # print("found at line:", num)
                break

    return num


def read_mesh_conv(filename, delimiter='\s+'):
    """
    This function reads the mesh convergence result into a pandas
    dataframe.
    """
    headline = first_match_line('dforce', filename)
    df = pd.read_csv(filename, delimiter, skiprows=headline-1)

    return df


def select_mesh_from_dataframe(df_mesh_conv, criteria=None, tols=None):
    """
    Select mesh based on the given criteria and tols. The mesh convergence
    data is assumed to be stored in a Pandas dataframe
    E.g., criteria=["dforce"], tols=[1e-3] finds the mesh required
    for achieving < 1e-3 Ha/bohr accuracy in forces.
    """
    
    meshes = df_mesh_conv.iloc[:,0]

    # loop over meshes in reverse order
    for i,mesh in reversed(list(enumerate(meshes))):
        # check if current mesh meets the criteria
        isOk = True
        for ci,ct in enumerate(criteria):
            ct_data = df_mesh_conv[ct][i]
            tol = tols[ci]
            if ct_data > tol:
                isOk = False
                break

        if isOk:
            return mesh, i


def select_mesh(fname, criteria=None, tols=None):
    """
    Select mesh based on the given criteria and tols. This routine
    returns the final selected mesh, the corresponding index of the
    mesh, and the whole mesh convergence dataframe.

    E.g., criteria=["deatom","dforce"], tols=[1e-3, 1e-3] finds the
    mesh required for achieving < 1e-3 accuracy in both energy and
    forces.
    """

    # fist read mesh convergence result into a Pandas dataframe
    df = read_mesh_conv(fname)

    # select mesh based on the given criteria and tolerances
    mesh, ind = select_mesh_from_dataframe(df, criteria, tols)

    return mesh, ind, df


def grep_results(folders):
    """
    This function greps results from a set of folders, and reads
    Eatoms, forces, pressures, stresses, nscfs, and walltimes
    from the .out and .static files available in each folder.
    The result is returned as a bunch of lists
    """
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
        cwd = os.getcwd() # save current working director (pwd)
        
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
        natom = grep_natom(out_fname)

        # grep energy
        eatom_i = grep_energy(out_fname)
        
        if not eatom_i:
            print("skipping " + folder_name)
            os.chdir(cwd)
            continue # skip current folder


        # grep forces
        #forces_i = read_forces_from_text_numpy(static_fname)
        forces_i = grep_forces(static_fname,natom)
        
        if not forces_i:
            print("skipping " + folder_name)
            os.chdir(cwd)
            continue # skip current folder
        
        # grep pressure
        press_i = grep_pressure(out_fname)
        
        stress_i = grep_stress(static_fname)
        if not stress_i:
            stress_i = grep_stress(out_fname)
        
        folder_names.append(folder_name)
        eatoms.append(eatom_i)
        forces.append(np.array(forces_i))
        pressures.append(press_i)
        stresses.append(np.array(stress_i))

        # grep nscf
        nscf = grep_nscf(out_fname)
        nscfs.append(nscf)

        # grep walltime
        wtime = grep_walltime(out_fname)
        wtimes.append(wtime)

        os.chdir(cwd)
    
    return folder_names,eatoms,forces,pressures,stresses,nscfs,wtimes


def reject_outliers(data, m = 5.189):
    """
    data must be a numpy array
    """
    if len(data) == 1: return data

    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]



def grep_timing_component(search_regex,content,verbose=True,includeOutliers=False,excludeList=None):
    """
    grep timing of a specific component from the content read from log file
    """
    timings = re.findall(search_regex, content, re.IGNORECASE)
    timings = [float(x) for x in timings]
    if timings:
        # exclude outliers
        # print(search_regex)
        # print("with outliers")
        # print(timings)
        if excludeList:
            timings = [j for i, j in enumerate(timings) if i not in excludeList]

        if includeOutliers:
            return timings, np.mean(np.array(timings))
        else:
            timings = reject_outliers(np.array(timings))
            # print("without outliers")
            # print(timings)
            t_avg = np.mean(timings)
            timings = timings.tolist() # convert numpy array back to python list
            return timings,t_avg
    else:
        if verbose:
            print("Timings for the expression " + search_regex[:40] + "... not found!")
        # return None,None
        return np.array([-1000.0]),np.array([-1000.0])

    


def grep_timing_components(fname,verbose=True):
    """
    grep timing components (time split) from .log file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return

    numbers_regex = r"-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?"
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Total number of processors:\s*(-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            nproc = int(match.group(1))
        else:
            match = re.search(r"nproc = ({0})".format(numbers_regex), content)
            if match:
                nproc = int(match.group(1))
            else:
                print("Number of proc not found!")

        # FORMAT
        # search_regex = r"Total time for Chebyshev filtering.* ({0})\s?ms".format(numbers_regex)
        # all_timings,average = grep_timing_component(search_regex,content,verbose)
        search_regex = r"Total time for Chebyshev filtering.* ({0})\s?ms".format(numbers_regex)
        cheb_filts,t_cheb_filts = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Total time for projection:.* ({0})\s?ms".format(numbers_regex)
        projections,t_proj = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Total time for solving.*eigenvalue problem:.* ({0})\s?ms".format(numbers_regex)
        eigsolve,t_eigsolve = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Total time for subspace rotation:.* ({0})\s?ms".format(numbers_regex)
        rot,t_rot = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Solving Poisson took.* ({0})\s?ms".format(numbers_regex)
        poisson,t_poisson = grep_timing_component(search_regex,content,verbose)
        if not poisson: # for old version of SPARC
            search_regex = r"AAR took.* ({0})\s?ms".format(numbers_regex)
            poisson,t_poisson = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Mixing.* ({0})\s?ms".format(numbers_regex)
        mixing,t_mixing = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Transfering Veff from phi-domain to psi-domain took.* ({0})\s?ms".format(numbers_regex)
        transferVeff,t_transferVeff = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Time for calculating local force components:.* ({0})\s?ms".format(numbers_regex)
        floc,t_floc = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Time for calculating nonlocal force components:.* ({0})\s?ms".format(numbers_regex)
        fnloc,t_fnloc = grep_timing_component(search_regex,content,verbose)

        search_regex = r"Calculating b & b_ref took.* ({0})\s?ms".format(numbers_regex)
        b_bref,t_b_bref = grep_timing_component(search_regex,content,verbose)


    return nproc,(t_cheb_filts,t_proj,t_eigsolve,t_rot,t_poisson,t_mixing,t_transferVeff,t_floc,t_fnloc,t_b_bref)



def pyTimeSplit(fnames,verbose=True):
    # if verbose:
    #     #print("npro   Cheb_filt   projection  densKernel  subsp_rot   poisson     mixing      trans_veff  force_loc   force_nloc  b_bref")
    #     print("npro   Cheb_filt   projection  eigsolve    subsp_rot    poisson     mixing      trans_veff  force_loc   force_nloc  b_bref")
    print("npro   Cheb_filt   projection  eigsolve    subsp_rot    poisson     mixing      trans_veff  force_loc   force_nloc  b_bref")
    #(t_cheb_filts,t_proj,t_chebMats,t_densMats,t_distrDs,t_rotDs,t_poisson,t_transferVeff,t_floc,t_fnloc,t_b_bref)
    for fname in fnames:
        nproc,timings_millsec = grep_timing_components(fname,verbose)
        # convert timings to seconds
        timings = [x/1000.0 for x in timings_millsec]
        print("%-4d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f"%
            (nproc,timings[0],timings[1],timings[2],timings[3],timings[4],timings[5],
                timings[6],timings[7],timings[8],timings[9]))


def grep_timing_scf(fname):
    """
    grep average SCF timings (excluding the first SCF)
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return
    
    numbers_regex = r"-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?"
    
    with open(fname, 'r') as f:
        content = f.read()
        
        match = re.search(r"Total number of processors:\s*(-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            nproc = int(match.group(1))
        else:
            print("Number of proc not found!")
            nproc = None

        search_regex = r"This SCF took.* ({0})\s?ms".format(numbers_regex)
        scfs,t_scf_avg = grep_timing_component(search_regex,content)
        t_scf_avg = np.mean(scfs[1:]) # find mean, excluding 1st scf

        return (nproc, t_scf_avg, scfs) # timings in ms


def pySCFtiming(fnames, nscf_per_MD=0, verbose=True):
    if verbose:
        if nscf_per_MD <= 0:
            print("nproc   timing (s/SCF)")
        else:
            print("nproc     t (s/SCF)   t (s/MD)")
    
    for fname in fnames:
        nproc,tscf_ms,timings_millsec = grep_timing_scf(fname)
        # convert timings to seconds
        # timings = [x/1000.0 for x in timings_millsec]
        tscf = tscf_ms / 1000.0
        if nscf_per_MD <= 0:
            print("%5d  %10.3f"%(nproc,tscf))
        else:
            print("%5d  %10.3f   %10.3f"%(nproc,tscf,tscf*nscf_per_MD))
        #print("%4d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f"%
        #    (nproc,timings[0],timings[1],timings[2],timings[3],timings[4],timings[5],
        #        timings[6],timings[7],timings[8],timings[9],timings[10],timings[11]))






def grep_timing_components_sq3(fname):
    """
    grep timing components (time split) from .log file
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return

    numbers_regex = r"-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?"
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Total number of processors:\s*(-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            nproc = int(match.group(1))
        else:
            nproc = None
            print("Number of proc not found!")

        # FORMAT
        # search_regex = r"Total time for Chebyshev filtering.* ({0})\s?ms".format(numbers_regex)
        # all_timings,average = grep_timing_component(search_regex,content)
        search_regex = r"Total time for Chebyshev filtering.* ({0})\s?ms".format(numbers_regex)
        cheb_filts,t_cheb_filts = grep_timing_component(search_regex,content)

        search_regex = r"Orthogonalization of orbitals took:.* ({0})\s?ms".format(numbers_regex)
        orths1,t_orth1 = grep_timing_component(search_regex,content)
        search_regex = r"Updating orbitals took:.* ({0})\s?ms".format(numbers_regex)
        orths2,t_orth2 = grep_timing_component(search_regex,content)
        # orths = [x+y for x,y in zip(orths1, orths2)]
        t_orth = t_orth1 + t_orth2

        search_regex = r"Total time for projection:.* ({0})\s?ms".format(numbers_regex)
        projections,t_proj = grep_timing_component(search_regex,content)
        t_proj = t_proj - t_orth # exclude the timing for orth from projection

        search_regex = r"Total time for calculating Chebyshev matrix-vector component:.* ({0})\s?ms".format(numbers_regex)
        chebMats,t_chebMats = grep_timing_component(search_regex,content)

        search_regex = r"Distributing and broadcasting Hp took.* ({0})\s?ms".format(numbers_regex)
        bcastHps,t_bcastHp = grep_timing_component(search_regex,content)

        search_regex = r"Calculating Density matrix:.* ({0})\s?ms".format(numbers_regex)
        densMats,t_densMats = grep_timing_component(search_regex,content)

        search_regex = r"Redistribute Density matrix:.* ({0})\s?ms".format(numbers_regex)
        distrDs,t_distrDs = grep_timing_component(search_regex,content)

        search_regex = r"Total time for updating rho and psi:.* ({0})\s?ms".format(numbers_regex)
        rotDs,t_rotDs = grep_timing_component(search_regex,content)

        search_regex = r"Solving Poisson took.* ({0})\s?ms".format(numbers_regex)
        poisson,t_poisson = grep_timing_component(search_regex,content)
        if not poisson: # for old version of SPARC
            search_regex = r"AAR took.* ({0})\s?ms".format(numbers_regex)
            poisson,t_poisson = grep_timing_component(search_regex,content)

        search_regex = r"Mixing.* ({0})\s?ms".format(numbers_regex)
        mixing,t_mixing = grep_timing_component(search_regex,content)

        search_regex = r"Transfering Veff from phi-domain to psi-domain took.* ({0})\s?ms".format(numbers_regex)
        transferVeff,t_transferVeff = grep_timing_component(search_regex,content)

        search_regex = r"Time for calculating local force components:.* ({0})\s?ms".format(numbers_regex)
        floc,t_floc = grep_timing_component(search_regex,content)

        search_regex = r"Time for calculating nonlocal force components:.* ({0})\s?ms".format(numbers_regex)
        fnloc,t_fnloc = grep_timing_component(search_regex,content)

        search_regex = r"Calculating b & b_ref took.* ({0})\s?ms".format(numbers_regex)
        b_bref,t_b_bref = grep_timing_component(search_regex,content)


    return nproc,(t_cheb_filts,t_orth,t_proj,t_bcastHp, t_chebMats+t_densMats,t_distrDs,t_rotDs,t_poisson,t_mixing,t_transferVeff,t_floc,t_fnloc,t_b_bref)



def pyTimeSplit_sq3(fnames,verbose=True):
    if verbose:
        #print("npro   Cheb_filt   projection  densKernel  subsp_rot   poisson     mixing      trans_veff  force_loc   force_nloc  b_bref")
        print("npro   Cheb_filt     orthPsi   projection   bcastHp    densKernel   distrDs    subsp_rot    poisson     mixing      trans_veff  force_loc   force_nloc  b_bref")
    
    #(t_cheb_filts,t_proj,t_chebMats,t_densMats,t_distrDs,t_rotDs,t_poisson,t_transferVeff,t_floc,t_fnloc,t_b_bref)
    for fname in fnames:
        nproc,timings_millsec = grep_timing_components_sq3(fname)
        # convert timings to seconds
        timings = [x/1000.0 for x in timings_millsec]
        print("%-4d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f"%
            (nproc,timings[0],timings[1],timings[2],timings[3],timings[4],timings[5],
                timings[6],timings[7],timings[8],timings[9],timings[10],timings[11],timings[12]))



def grep_timing_components_ddbp(fname, excludeList=None):
    """
    grep timing components (time split) from .log file
    excludeList: indices (starting from 0) of DDBP basis reltated timings to exclude
    """
    if not os.path.isfile(fname):
        print("file " + fname + " doesn't exist!")
        return

    numbers_regex = r"-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?"
    with open(fname, 'r') as f:
        content = f.read()
        match = re.search(r"Total number of processors:\s*(-?\s*\d+\.?\d*(?:[Ee]\s*[+-]?\s*\d+)?)", content)
        if match:
            nproc = int(match.group(1))
        else:
            nproc = None
            print("Number of proc not found!")

        # FORMAT TEMPLATE
        # search_regex = r"Total time for Chebyshev filtering.* ({0})\s?ms".format(numbers_regex)
        # all_timings,average = grep_timing_component(search_regex,content)
        search_regex = r"Calculating DDBP basis took.* ({0})\s?ms".format(numbers_regex)
        calc_basis,t_calc_basis = grep_timing_component(search_regex,content,True,True,excludeList)
        # calc_basis = calc_basis[1:]
        t_calc_basis = np.mean(calc_basis)

        search_regex = r"Calculating DDBP basis overlap matrix took.* ({0})\s?ms".format(numbers_regex)
        basis_overlap,t_basis_overlap = grep_timing_component(search_regex,content,True,True,excludeList)
        
        search_regex = r"Aligning orbitals with current basis took.* ({0})\s?ms".format(numbers_regex)
        align_basis,t_align_basis = grep_timing_component(search_regex,content,True,True,excludeList)

        search_regex = r"Calculating DDBP Hamiltonian took.* ({0})\s?ms".format(numbers_regex)
        h_ddbp,t_h_ddbp = grep_timing_component(search_regex,content,True,True,excludeList)
        #print(h_ddbp)
        if not h_ddbp:
            search_regex = r"Construction of H_DDBP_loc:.* ({0})\s?ms".format(numbers_regex)
            h_ddbp_loc,t_h_ddbp_loc = grep_timing_component(search_regex,content,True,True)
            search_regex = r"Construction of Vnl_DDBP:.* ({0})\s?ms".format(numbers_regex)
            h_ddbp_nloc,t_h_ddbp_nloc = grep_timing_component(search_regex,content,True,True)
            t_h_ddbp = t_h_ddbp_loc + t_h_ddbp_nloc

        search_regex = r"Total time for Chebyshevfilter_constants_DDBP:.* ({0})\s?ms".format(numbers_regex)
        cheb_const,t_cheb_const = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Total time for DDBP Chebyshev filtering.* ({0})\s?ms".format(numbers_regex)
        Cheb_filt,t_cheb_filt = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Total time for projection.*: ({0})\s?ms".format(numbers_regex)
        proj,t_proj = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Total time for solving subspace eigenvalue problem for DDBP:.* ({0})\s?ms".format(numbers_regex)
        eig,t_eig = grep_timing_component(search_regex,content,True,False)
        # print(eig)
        # print(t_eig)
        search_regex = r"DDBP subpsace eigenproblem: bcast eigvals took.* ({0})\s?ms".format(numbers_regex)
        bcast_eig,t_bcast_eig = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Total time for subspace rotation \(DDBP\):.* ({0})\s?ms".format(numbers_regex)
        rot,t_rot = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Total time for calculating psi on the grid \(DDBP\).* ({0})\s?ms".format(numbers_regex)
        recover_psi,t_recover_psi = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Transfering density took.* ({0})\s?ms".format(numbers_regex)
        transferRho,t_transferRho = grep_timing_component(search_regex,content,True,False)

        search_regex = r"DDBP transfer orbitals E2D took.* ({0})\s?ms".format(numbers_regex)
        transferPsi,t_transferPsi = grep_timing_component(search_regex,content,True,False)        

        search_regex = r"Solving Poisson took.* ({0})\s?ms".format(numbers_regex)
        poisson,t_poisson = grep_timing_component(search_regex,content,True,False)
        if not poisson: # for old version of SPARC
            search_regex = r"AAR took.* ({0})\s?ms".format(numbers_regex)
            poisson,t_poisson = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Mixing.* ({0})\s?ms".format(numbers_regex)
        mixing,t_mixing = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Transfering Veff from phi-domain to psi-domain took.* ({0})\s?ms".format(numbers_regex)
        transferVeff,t_transferVeff = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Time for calculating local force components:.* ({0})\s?ms".format(numbers_regex)
        floc,t_floc = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Time for calculating nonlocal force components:.* ({0})\s?ms".format(numbers_regex)
        fnloc,t_fnloc = grep_timing_component(search_regex,content,True,False)

        search_regex = r"Calculating b & b_ref took.* ({0})\s?ms".format(numbers_regex)
        b_bref,t_b_bref = grep_timing_component(search_regex,content,True,False)

    return nproc,(t_calc_basis,t_basis_overlap,t_align_basis,t_h_ddbp,t_cheb_const,t_cheb_filt,t_proj,t_eig,t_bcast_eig,t_rot,t_recover_psi,t_transferRho,t_poisson,t_mixing,t_transferVeff,t_floc,t_fnloc,t_b_bref,t_transferPsi)



def pyTimeSplit_ddbp(fnames,verbose=True,BasisExcludeList=None):
    if verbose:
        #print("npro   Cheb_filt   projection  densKernel  subsp_rot   poisson     mixing      trans_veff  force_loc   force_nloc  b_bref")
        #print("npro   calc_basis   align_basis   H_DDBP  cheb_const   Cheb_filt   projection    eig     subsp_rot    recover_psi   transferRho  poisson     mixing    trans_veff    force    b_bref")
        print("npro   calc_basis  align_basis   H_DDBP    cheb_const   Cheb_filt   projection    eig      subsp_rot   recover_psi  transferRho  poisson    mixing    trans_veff    force      b_bref      E2D_Psi")
    
    #(t_cheb_filts,t_proj,t_chebMats,t_densMats,t_distrDs,t_rotDs,t_poisson,t_transferVeff,t_floc,t_fnloc,t_b_bref)
    for fname in fnames:
        if not BasisExcludeList:
            BasisExcludeList = [0]
        nproc,timings_millsec = grep_timing_components_ddbp(fname,BasisExcludeList)
        # convert timings to seconds
        timings = [x/1000.0 for x in timings_millsec]
        print("%-4d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f %10.6f  %10.6f  %10.6f  %10.6f  %10.6f"%
            (nproc,timings[0],timings[1]+timings[2],timings[3],timings[4],timings[5],timings[6],timings[7]+timings[8],timings[9],timings[10],timings[11],timings[12],timings[13],timings[14],timings[15]+timings[16],timings[17],timings[18]))
        #print("%-4d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f %10.6f  %10.6f  %10.6f  %10.6f  %10.6"%
        #    (nproc,timings[0],timings[1]+timings[2],timings[3],timings[4],timings[5],timings[6],timings[7]+timings[8],timings[9],timings[10],timings[11],timings[12],timings[13],timings[14],timings[15]+timings[16],timings[17],timings[18]))





def main():
    fname = 'Si8-ONCV-0.4.refout'
    eatom = grep_energy(fname)
    pressure = grep_pressure(fname)

    static_fname = 'Si8-ONCV-0.4.static'
    forces = grep_forces(static_fname,8)
    print(forces)


if __name__ == '__main__':
    # main()
    fname = "np_1000/sprc-calc-1000.log"
    pyTimeSplit(fname)


