## SPARC tutorial session

#### 0. Access the PACE-ICE cluster
Everyone attending this session has been given temporary access to the PACE-ICE cluster at Georgia Tech. In order to do that, we have created a guest GT account for everyone. You should have received your GT username and link to set up VPN through email.

To access the PACE-ICE cluster, you will need: 1) your GT username and password, and 2) VPN. Follow the following instructions to log in to the PACE-ICE cluster through ssh.
* Connect to Georgia Tech VPN using the GlobalProtec VPN client. (If you haven't already downloaded it, it can be downloaded here: `https://vpn.gatech.edu/`. For Linux, click on the `Download Linux VPN Client` icon. For other OS's, click on the `Getting Started With VPN FAQ` icon and find the corresponding link there).
* Open a terminal, run
```bash
$ ssh <username>@login-pace-ice.pace.gatech.edu
```
where `username` is your GT username. You will be asked to type in the password. Once you log in successfully, your bash prompt will look like this
```bash
[qxu78@login-pace-ice-1 ~]$ 
```


#### 1. Setup
Download the materials to the cluster (`git clone https://github.com/xuqimen/ACS2022_Workshop_SPARC.git`). Go to `setups/` and run the following command
```bash
$ source setup.sh
```
to complete the setup. Note that this is not required for installing SPARC. It simply changes some settings of the shell enviroment for your convenience, and copies the `data` folder to your home directory.


#### 2. Installing SPARC
SPARC is an open-source code released on GitHub (https://github.com/SPARC-X/SPARC). Installation guides can be found on the GitHub page (README.md or the documentation in doc/). There are several optional ways you can choose to install SPARC on your cluster, depending on the dependencies you want to use. In this tutorial, we will use the default option (option 2: compiling with MKL).

On the PACE-ICE cluster, you would need two modules, 1) intel/19.0.5 and 2) mvapich2/2.3.2. These two modules are already loaded on PACE-ICE by default. You can verify this by running
```bash
$ module list
```
You should see something like this:
```bash

Currently Loaded Modules:
  1) xalt/2.8.4   2) intel/19.0.5   3) mvapich2/2.3.2   4) gcc-compatibility/8.3.0   5) pace/2020.01
```

Now you are ready to install SPARC! Download the source code from GitHub by
```bash
$ git clone https://github.com/SPARC-X/SPARC.git
```
Once the code is downloaded, go to the `SPARC/src` directory. Compile the code by
```bash
$ make clean; make -j 2
```
The `-j 2` option is just to speedup the compilation. If you have more processors available, you can increase this number (e.g. -j 8) to make it faster. On PACE-ICE, the installation will take around 1 minute with `-j 2`. Once the code is compiled successfully, a binary named `sparc` will be created in the `SPARC/lib` directory.



