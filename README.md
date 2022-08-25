## SPARC tutorial session

### 0. Accessing the PACE-ICE cluster
Everyone attending this session has been given temporary access to the [PACE-ICE cluster](https://docs.pace.gatech.edu/ice_cluster/ice-guide/) at Georgia Tech. In order to do that, we have created a guest GT account for everyone. You should have received your GT username and link to set up VPN through email.

To access the PACE-ICE cluster, you will need: 1) your GT username and password, and 2) VPN. Follow the following instructions to log in to the PACE-ICE cluster through ssh.
* Connect to Georgia Tech VPN using the GlobalProtect VPN client. (If you haven't already downloaded it, it can be downloaded here: `https://vpn.gatech.edu/`. For Linux, click on the `Download Linux VPN Client` icon. For other OS's, click on the `Getting Started With VPN FAQ` icon and find the corresponding link there).
* Open a terminal, run  
`ssh <username>@login-pace-ice.pace.gatech.edu`  
where `username` is your GT username. You will be asked to type in the password. Once you log in successfully, your bash prompt will look like this  
`[qxu78@login-pace-ice-1 ~]$ `  


### 1. Setup
Download the materials to the cluster:
```bash
git clone https://github.com/xuqimen/ACS2022_Workshop_SPARC.git
```
Go to `setups/` and run the following command
```bash
source setup.sh
```
to complete the setup. Note that this is not required for installing SPARC. It simply changes some settings of the shell enviroment for your convenience, and copies the `data` folder to your home directory.


### 2. Installing SPARC
#### 2.1 Download and install
SPARC is an open-source code released on GitHub (https://github.com/SPARC-X/SPARC). Installation guides can be found on the GitHub page (README.md or the documentation in doc/). There are several optional ways you can choose to install SPARC on your cluster, depending on the dependencies you want to use. In this tutorial, we will use the default option (option 2: compiling with MKL).

On the PACE-ICE cluster, you would need two modules, 1) intel/19.0.5 and 2) mvapich2/2.3.2. These two modules are already loaded on PACE-ICE by default. You can verify this by running
```bash
module list
```
You should see something like this:
```bash

Currently Loaded Modules:
  1) xalt/2.8.4     3) mvapich2/2.3.2   5) git/2.25.0                7) anaconda3/2021.05
  2) intel/19.0.5   4) pace/2020.01     6) gcc-compatibility/8.3.0
```

Now you are ready to install SPARC! Download the source code from GitHub by
```bash
git clone https://github.com/SPARC-X/SPARC.git
```
Once the code is downloaded, go to the `SPARC/src` directory. Compile the code by
```bash
make clean; make -j 2
```
The `-j 2` option is just to speedup the compilation. If you have more processors available, you can increase this number (e.g. -j 8) to make it faster. On PACE-ICE, the installation will take around 1 minute with `-j 2`. Once the code is compiled successfully, a binary named `sparc` will be created in the `SPARC/lib` directory.


#### 2.2 Verify your installation
To verify your installation, go to `SPARC/tests`,
```bash
cd ../tests
```
and run
```bash
python SPARC_testing_script.py quick_run
```
It will run a few quick tests to verify the installation. If all tests passed, congratulations! You have installed SPARC correctly.


### 3. Running tests


#### 3.1 Running jobs on PACE-ICE
Where you log in, you're on **login nodes**. Login nodes are shared by all. You should not use the login nodes for any resource-intensive activities, as it prevents others from using the cluster. Instead, you should request a **compute nodes** for doing computations.

There are two ways to request compute nodes.

- **Interactive jobs**.
  - An interactive job reserves resources on compute nodes to use interactively.
  - To request 1 node with 12 cores for an interactive job on the pace-ice queue for 2 hours:  
  `qsub -I -q pace-ice -l nodes=1:ppn=12,pmem=7gb,walltime=2:00:00`  

- **Job submission through [PBS scripts](https://docs.pace.gatech.edu/scheduler/job_submission/)**.
  - We have prepared a sample PBS script in the `/data/tools` directory in the [materials repo](https://github.com/xuqimen/ACS2022_Workshop_SPARC). 
  - Once you prepare the PBS script, the next step is to submit the script to the scheduler using `qsub`:  
  `qsub helloexmple.pbs`  
  Once this step is complete, the job has been successfully submitted to the scheduler.

Let's use the **interactive job** way as a start. We'll try using the PBS script later.


#### 3.2 Input files
The required input files to run a simulation with SPARC are (with shared names)
- ".inpt" file -- User options and parameters.
- ".ion" file -- Atomic information.

It is required that the ".inpt" and ".ion" files are located in the same directory and share the same name. A detailed description of the input options is provided in the documentation located in SPARC/doc/. Examples of input files can be found in the SPARC/tests directory .

#### 3.3 Pseudopotentials
SPARC requires pseudopotential files in the [psp8 format](https://docs.abinit.org/developers/psp8_info/), which can be generated by D. R. Hamann's open-source pseudopotential code [ONCVPSP](http://www.mat-simresearch.com/). There are accurate and efficient pseudopotential sets available online. For example, the [SG15 ONCV pseudopotentials](http://www.quantum-simulation.org/potentials/sg15_oncv/), which can be downloaded in psp8 format [here](https://github.com/xuqimen/SG15_pseudopotentials_psp8_upf); the [pseudoDOJO ONCV potentials](http://www.pseudo-dojo.org/). We have also developed our own pseudopotentials and will be released soon.

#### 3.4 Some example tests
Let's run some example tests. There are many tests in the `SPARC/tests/`, we'll run some of them and explain the input files along the way.


### 3.5 Extra tests
- Mesh convergence tests.
- Strong scaling of a big system.
