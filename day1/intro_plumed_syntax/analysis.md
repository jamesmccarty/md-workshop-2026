# A Brief Introduction to PLUMED Syntax and Making Histograms

The aim of this tutorial is provide an introduction to the PLUMED syntax. You will learn how to calculate some simple collective variables for analyzing existing trajectories.

Once this tutorial is completed students will be able to:

- Write a simple PLUMED input file and use it to analyze a trajectory
- Print collective variables such as distances, torsional angles, and coordination numbers using the PRINT action
- Be able to use PLUMED to calculate ensemble averages and histograms

**Files**
Files to complete this tutorial can be accessed here:
[tutorial files](coming soon)

These files are already located on bigzam:
/opt/workshop/intro_2_plumed/

## Getting Started

Use PuTTY to connect to bigzam as you did in the previous [tutorial](../lj_fluid/lj_fluid_tutorial.md). Open PuTTY from the Window Start menu and enter `bigzam.local` for the Host Name. Login using the terminal using your username and password. 

**Important**: Once connected to the workshop computer, set your environment variables by typing:

{% highlight git %}
source setup.sh
{% endhighlight %}

Copy the tutorial files by typing in the terminal:

In the terminal type:
{% highlight git %}
cp -r /opt/workshop/intro_2_plumed/ ~/
{% endhighlight %}

This will copy the necessary tutorial files to your home directory on bigzam.

**Tip**: You can press the Tab key to automatically complete file and directory names. This can save time and help avoid typing errors.

Move into the intro_2_plumed directory:

{% highlight git %}
cd ~/intro_2_plumed
{% endhighlight %}

## What is PLUMED? 

So far we have been using GROMACS to build and run MD simulations. GROMACS comes with some analysis tools (such as RMSD, RMSF, and radius of gyration that you used in the previous tutorial). However, the features in GROMACS for analyzing trajectories are somewhat limited. Furthermore, some researchers may use other MD codes and analysis codes developed for one code may not be compatible with other MD codes. To overcome this issue, the MD community of researchers has developed a set of analysis tools called the The community-developed PLUgin for MolEcular Dynamics ([PLUMED](https://www.plumed.org/)) that can be used with several different MD codes. 

PLUMED is a plug-in library that can be incorporated into many MD codes including GROMACS. Once it is incorporated you can use PLUMED to perform a variety of different analyses on-the-fly and to bias the sampling in MD simulations as we will do in Day 2 of this workshop! Additionally, PLUMED can be used as a standalone code for analyzing trajectories. If you are using the code in this way you can run the PLUMED executable simply by issuing the command:

{% highlight git %}
plumed
{% endhighlight %}

To see a list of all option that PLUMED can do, type the following into the terminal: 

{% highlight git %}
plumed --help
{% endhighlight %}

The output of this command is a list of tasks that PLUMED can perform. In this tutorial we will use PLUMED as a post-processing tool to analyze MD simulation trajectories. The PLUMED tool for analyzing existing trajectories is the `driver` tool. Let's look at the options of PLUMED driver by issuing the following command:

{% highlight git %}
plumed driver --help
{% endhighlight %}

This shows you a list of all the things we can do with the PLUMED driver. For all of these options, however, we are going to need to create a PLUMED input file. 
This tutorial will introduce you to the syntax of the PLUMED input file. We will use this syntax later in the second day of the workshop to run enhanced sampling simulations on-the-fly during a MD simulation, so all the things that you will learn now will be useful later when you will run PLUMED coupled to GROMACS.

## Computing and printing simple collective variables

Chemical systems contain an enormous number atoms, which, in most cases makes it simply impossible for us to understand anything by monitoring the atom positions directly. Consequently, we introduce Collective variables (CVs) that describe the chemical processes we are interested in and monitor these simpler quantities instead.

A collective variable is any mathematical function of the atomic coordinates that can be useful for tracking structural changes in a complex molecular system containing hundreds or thousands of atoms. Instead of tracking and visualizing each individual atom's motion, it is often more informative to monitor a few judiciously chosen collective variables (CVs) that can tell us about what structural state the molecule is in. 

In this tutorial, you will prepare input files to calculate and print some common simple collective variables on a pre-calculated trajectory. The trajectory is of the folding of the GB1 protein (Immunoglobulin-binding protein G, domain B1). This is a small 56-amino acid protein. The tutorial contains the following files:
- GB1_native.pdb: A PDB file of the native GB1 protein
- traj-whole.xtc: A GROMACS trajectory in .xtc format. GB1 has already been make whole by fixing the periodic boundary conditions using the `trjconv` code in GROMACS.
- traj-broken.xtc: The same GROMACS trajectory as it was originally produced by GROMACS. Here the protein molecule is artificially broken by the periodic box convention. 
- traj-whole-skip10.pdb: A pdb file of every 10th frame of the trajectory that you can use to view the folding simulation in PyMOL.

It is always a good idea to view the movie of the trajectory to try and understand the simulation before computing any collective variables. Use the WinSCP app on your Windows machine to transfer the traj-whole-skip10.pdb file to your local machine and load this into PyMOL. You should be able to see the protein folding around frame 240. 

Folding of GB1 around frame 240:

![GB1 folded](../../images/GB1_folded.png)

Recall that we skipped every 10 frames, so frame 240 corresponds to frame 2400 in the traj-whole.xtc trajectory. 

Now that you have a sense of the trajectory, we can compute and print some collective variables. One of the simplest collective variable (CV) we can monitor is the distance between two atoms.

On bigzam, have a look at the first ten lines of the pdb file by typing:

{% highlight git %}
head GB1_native.pdb
{% endhighlight %}

You should see the following:
{% highlight git %}
ATOM      1  N   MET     1      12.970  18.510  30.950  1.00  1.00            
ATOM      2  CA  MET     1      13.940  18.530  29.840  1.00  1.00            
ATOM      3  C   MET     1      13.140  18.690  28.520  1.00  1.00            
ATOM      4  O   MET     1      12.010  18.220  28.400  1.00  1.00            
ATOM      5  CB  MET     1      14.730  17.220  29.880  1.00  1.00            
ATOM      6  CG  MET     1      15.740  16.980  28.740  1.00  1.00            
ATOM      7  SD  MET     1      17.380  17.020  29.360  1.00  1.00            
ATOM      8  CE  MET     1      17.170  16.060  30.820  1.00  1.00            
ATOM      9  N   THR     2      13.720  19.410  27.570  1.00  1.00            
ATOM     10  CA  THR     2      13.090  19.660  26.280  1.00  1.00            
{% endhighlight %}

In the pdb file atoms are indexed starting from 1 in the second column. Keep in mind that the numbering scheme in PLUMED also starts from 1. If you want to see that atom index for the alpha carbons for example, you can type:

{% highlight git %}
grep 'CA' GB1_native.pdb
{% endhighlight %}

Here we see a list of just the alpha carbons (atom type CA):
{% highlight git %}
ATOM      2  CA  MET     1      13.940  18.530  29.840  1.00  1.00            
ATOM     10  CA  THR     2      13.090  19.660  26.280  1.00  1.00            
ATOM     17  CA  TYR     3      12.730  17.030  23.610  1.00  1.00            
ATOM     29  CA  LYS     4      12.180  17.660  19.890  1.00  1.00            
ATOM     38  CA  LEU     5      10.250  15.790  17.220  1.00  1.00            
ATOM     46  CA  ILE     6      11.080  16.100  13.480  1.00  1.00            
ATOM     54  CA  LEU     7       8.010  15.160  11.390  1.00  1.00            
ATOM     62  CA  ASN     8       8.630  13.970   7.910  1.00  1.00            
ATOM     70  CA  GLY     9       5.210  12.640   6.970  1.00  1.00            
ATOM     74  CA  LYS    10       3.590  12.600   3.500  1.00  1.00            
{% endhighlight %}

Suppose we want to monitor the distance between CA on residue 1 and CA on residue 10 during the simulation. To do this, we can create a plumed file called plumed.dat with the following lines:

{% highlight git %}
d: DISTANCE ATOMS=2,74

PRINT ARG=d FILE=distance.dat STRIDE=1
{% endhighlight %}

In the above input file, we define the variable d as the distance between atom 2 (CA on residue 1) and atom 74 (CA on residue 74). Then we use the PRINT command to print the variable d to the file `distance.dat`. The STRIDE keyword sets how frequently to print to the output file. (STRIDE=1 means print every frame). 



