# Simulation of miniprotein in water

In this tutorial you will use GROMACS in the Unix command line to build and equilibrate a solvated protein system for the mini-protein [YYDPETGTWY](https://pubs.acs.org/doi/10.1021/ja8030533). The goal is to understand the MD simulation parameters that are chosen and the details about how the system was built.

Once this tutorial is completed students will be able to:

- Set up a protein in explicit solvent for MD simulation
- Understand GROMACS output files and be able to transfer files for analysis
- Calculate RMSD and monitor structural stability
- Analyze protein dynamics

**Files**
Files to complete this tutorial can be accessed here:
[tutorial files](coming soon)

These files are already located on bigzam:
/opt/workshop/miniprotein/

## Getting Started

Use PuTTY to connect to bigzam as you did in the previous [tutorial](../lj_fluid/lj_fluid_tutorial.md). Open PuTTY from the Window Start menu and enter `bigzam.local` for the Host Name. Log in using the terminal using your username and password:

![terminal window](../../images/LoginScreen_screenshot.png)

Once connected to the workshop computer, set your environment variables:

{% highlight git %}
source setup.sh
{% endhighlight %}

Copy the miniprotein tutorial files by typing in the terminal:

In the terminal type:
{% highlight git %}
cp -r /opt/workshop/miniprotein/ ~/
{% endhighlight %}

This will copy the necessary tutorial files to your home directory on bigzam.

**Tip**: You can press the Tab key to automatically complete file and directory names. This can save time and help avoid typing errors.

Move into the miniprotein directory:

{% highlight git %}
cd ~/miniprotein
{% endhighlight %}

From this directory, you will run the simulation. 

## Convert a pdb structure file (.pdb) to a GROMACS topology (.top) and structure (.gro) file

The Protein Data Bank (PDB) file [2RVD.pdb](https://www.rcsb.org/structure/2RVD) was obtained from NMR experiment and contains 20 conformers. We will only need to first conformer for starting our simulation. 

We first need to convert this pdb file to a format that GROMACS uses. The GROMACS command `pdb2gmx` will generate a GROMACS topology file (.top), a position restraint file (posre.itp), and a GROMACS structure file (.gro). Generating a topology file can sometimes be tricky if there are missing atoms, non-standard amino acids, or small molecules in the pdb file. A basic usage of `pdb2gmx` is 

{% highlight git %}
gmx pdb2gmx -f 2RVD.pdb -o conf.gro -p topol.top -ignh
{% endhighlight %}   

where -f signals the input pdb file, here called `2RVD.pdb` and -o is the output structure file and -p is the output topology file. The flag -ignh means to ignore H atoms in the PDB file. The -ignh flag is useful to deal with non-conventional naming of H atoms in pdb files; however, if you need to persevere the exact position of the H coordinates in the structure file, then you should not use the -ignh flag.

Typing the above command will prompt you to select a force field. The force field is the set of parameterized potential energy functions that are used to calculate all the forces acting on the atoms. A number of force field models are provided within GROMACS. Select the force field model by typing the appropriate number and hitting the Enter/Return key. You will also be prompted to select a water model. For this tutorial select the AMBER99SB-ILDN force fields from [Lindorff-Larsen et al., Proteins 78, 1950-58, 2010](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.22711), and when prompted, select the TIP3P water model. You should now have generated a structure file conf.gro, a topology file topol.top, and a position restraint file posre.itp.

The topology file `topol.top` contains all the information necessary to define the molecule within a simulation. This information includes nonbonded parameters (atom types and charges) as well as bonded parameters (bonds, angles, and dihedrals).

## Defining a simultion box and adding solvent molecules to the box

Now that we have converted the protein pdb structure file into a GROMACS topology and structure file, we need to define a simulation box. To create a simulation box of specified dimensions, use the `editconf` command by typing in the terminal:


{% highlight git %}
gmx editconf -f conf.gro -o conf_box.gro -c -d 1.0 -bt cubic
{% endhighlight %}

The -f flag signals the input GROMACS structure file (conf.gro) which we created above. The -o signals a new output structure file (conf_box.gro) which contains the box definition. The -c flag centers the protein in the box, and the -d specifies a minimum distance from the protein to the edge of the box (in nm). Specifying a solute-box distance of 1.0 nm means that there are at least 2.0 nm between any two periodic images of a protein. This distance will be sufficient for just about any cutoff scheme commonly used in simulations. The -bt cubic flag indicates to use a cubic box shape. 

Now that we have generated a simulation box, which we must fill with water. We add solvent by typing

{% highlight git %}
gmx solvate -cp conf_box.gro -cs spc216.gro -o conf_solv.gro -p topol.top
{% endhighlight %}

The input structure file is signaled by the -cp flag and is the structure file that was generated by the previous command (conf_box.gro). The -cs flag is the solvent configuration structure file that comes standard with a GROMACS installation. We are using spc216.gro, which is a generic equilibrated 3-point solvent model. You can use spc216.gro as the solvent configuration for SPC, SPC/E, or TIP3P water, since they are all three-point water models. The -o flag indicates a new output structure file that contains both the protein and the solvent (conf_solv.gro). The -p flag is our topology file generated above, which will be modified to include the solvent force fields as well.

**WARNING** When you run the `gmx solvate` command, the topology file is modified. This means that if you run this command more than once on the same topology file, you will get a mismatch in atoms between your structure file and topology file that can lead to errors later on. Whenever GROMACS modifies an existing file, a backup of the original file is created starting with a hashtag (#topol.top). 

## Adding ions

Usually, we will want to add small molecule ions (NaCl) in order to neutralize the system (if the protein has a nonzero net charge) and to perform the simulation at physiological or experimental salt conditions. The GROMACS command `genion` will randomly replace water molecules with ions. In order to run `genion` we need to generate a compiled GROMACS input file (.tpr) file. A .tpr file contains information from both the structure file (.gro) and topology file (.tpr) and is generated using the grompp (GROMACS pre-processor) command, which will also be used later when we run our first simulation. 

To produce a processed input (.tpr) file, we need an additional input file with extension .mdp (molecular dynamics parameter file) that contains specified parameters for running a simulation. For this tutorial, you can use the `ions.mdp` file provided. Generate the .tpr file by typing the following in the terminal:

{% highlight git %}
gmx grompp -f ions.mdp -c conf_solv.gro -p topol.top -o ions.tpr
{% endhighlight %}

where -o specifies the output .tpr file which we have called `ions.tpr`. Now that you have generated the `ions.tpr` file, you can add small ions with the command:

{% highlight git %}
gmx genion -s ions.tpr -o conf_solv_ions.gro -p topol.top -conc 0.15 -neutral
{% endhighlight %}

where the -conc signals to make the final ion concentration 0.15 M and -neutral signals to neutralize the system so that the net charge in the system is zero. GROMACS will randomly replace certain molecules with small ions to reach our specified concentration and charge balance. When prompted select group 13 for SOL, which is the solvent group. This will tell GROMACS to replace some of the solvent molecules with small ions. 

Congratulation, you have now built a simulation box containing a small protein, water, and ions. To run an MD simulation of this model, we will need the generated GROMACS structure file (conf_solv_ions.gro) that contains the initial positions of all the atoms in the simulation box and the GROMACS topology file (topol.top), that contains information on how to handle the calculation of the forces on the atoms. Before performing an MD simulation of this model, we will first perform an energy minimization of this structure to relax any steric restraints. 

## Energy minimization 

