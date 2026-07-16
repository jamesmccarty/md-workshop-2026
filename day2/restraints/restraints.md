# Implementing Restraints

In this tutorial you will learn how to implement a restraint using PLUMED. By applying restraints on one or more collective variables, you can restrict the sampling to a desired region. This will introduce you to the idea of using a bias to manipulate a simulation on-the-fly. 

Once this tutorial is completed students will be able to:

- Apply a restraint on a simulations over one or more collective variables
- Understand the effect of a restraint on the acquired statistics
- Perform an out-of-equilibrium simulation with a moving restraint

 
![restraint_figure](../../images/)

In the [previous introduction to PLUMED syntax tutorial](../../day1/intro_plumed_syntax/analysis.md), you used PLUMED as a post-processing tool to  calculate properties of the system on a previously generated MD simulation trajectory. In addition to being a post-processing tool, PLUMED can also interface with the MD code "on-the-fly" during a MD simulation. The atomic coordinates and atoms at a given instant can be passed to the PLUMED code to manipulate a simulation on-the-fly. This allows one to bias the simulation by adding extra constraints in addition to the standard force field terms. The additional energy terms are usually referred as **Bias**. 

## The molecule of the day

Throughout today's tutorial you will play with the terminally capped alanine molecule (specifically, N-acetyl-L-alanine-N'-methylamide). The molecule has the total chemical formula C<sub>6</sub>H<sub>12</sub>N<sub>2</sub>O<sub>2</sub> and consists of a central alanine amino group, with flanking acetyl and methylamide groups to have neural (uncharged) end-groups. In the computational literature this molecule is frequently abbreviated as Ace-Ala-NMe. 

![Alanine_dipeptide1](../../images/belfast-2-ala.png)

While simple, this system is a nice example because it presents two metastable states separated by a high free-energy barrier. It is conventional use to characterize the two states in terms of Ramachandran dihedral angles, which are denoted with $$\phi$$ and $$\psi$$.

![Alanine_dipeptide2](../../images/belfast-2-transition_ala.png)

**Files**
Files to complete this tutorial can be accessed here:
[tutorial files](coming soon)

These files are already located on bigzam:
/opt/workshop/using_restraints/

## Getting Started 

Use PuTTY to connect to bigzam as you did in previous [tutorials](../lj_fluid/lj_fluid_tutorial.md). Open PuTTY from the Window Start menu and enter `bigzam.local` for the Host Name. Login using the terminal using your username and password.

**Important**: Once connected to the workshop computer, set your environment variables by typing:

{% highlight git %}
source setup.sh
{% endhighlight %}

Copy the tutorial files by typing in the terminal:

In the terminal type:
{% highlight git %}
cp -r /opt/workshop/using_restraints/ ~/
{% endhighlight %}

This will copy the necessary tutorial files to your home directory on bigzam.

**Tip**: You can press the Tab key to automatically complete file and directory names. This can save time and help avoid typing errors.

Move into the using_restraints directory:

{% highlight git %}
cd ~/using_restraints
{% endhighlight %}

Within this directory you will find the following files:

- dialaA.pdb: A reference PDB structure file of the molecule
- alanine_dipeptide.gro: A GROMACS structure file (.gro)
- topol.top: A GROMACS topology file (.top)
- vacuum.mdp: A GROMACS parameter file (.mdp) 

## Adding a constant bias potential 
In the following we will see how to apply a constant bias potential. We can add a harmonic restraint on any variable according to the formula:

$$V_{bias}(x)=\frac{1}{2}\kappa (x-x_0)^2 $$

where $$x_0$$ is the value around which we want to put a restraint and $$\kappa$$ determines the strength (spring constant) of the restraint. Notice that if we implement this restraint, there will be a linear restoring force proportional to $$\kappa$$ to keep $$x$$ near the value of $$x_0$$. 

An example PLUMED input file called plumed_example1.dat is provided with the workshop materials. Have a loot at this file with:

{% highlight git %}
cat plumed_example1.dat 
{% endhighlight %}  

Here we see the plumed_example1.dat file contents:

{% highlight git %}
# set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

# Set up harmonic restraint 
restraint: RESTRAINT ARG=phi KAPPA=10 AT=-1.5

# Print output
PRINT FILE=dihedrals_weak_restraint.dat ARG=phi,psi,restraint.bias STRIDE=100 
{% endhighlight %}

The first two lines tell PLUMED to compute the $\phi$ and $\psi$ angles and store these values as variables `phi` and `psi`. The next line specifies that we are putting a restraint on the `phi` variable centered around a values of -1.5 radians with a harmonic spring constant of $$\kappa=10$$ kJ/mol specified by the `KAPPA` keyword. Finally, we are printing both the dihedral angles and the value of the bias potential to the file called `dihedrals_weak_restraint.dat`. The frequency of writing to the output file is specified by the STRIDE keyword. Here we are printing every 100 steps (0.2 ps). After looking at this file, you can set up a **biased** simulation in GROMACS by typing:

{% highlight git %}
gmx grompp -f vacuum.mdp -c alanine_dipeptide.gro -p topol.top -o run1.tpr
{% endhighlight %}  

and run the simulation by typing:

{% highlight git %}
gmx mdrun -v -deffnm run1 -plumed plumed_example1.dat -nt 1
{% endhighlight %}

Notice here the command for running the biased simulation is identical to running an unbiased simulation using GROMACS `mdrun` except that we must specify a plumed input file with the `-plumed` flag in the `mdrun` command. This tells GROMACS to run with PLUMED compiled on-the-fly. 

Note: If you run the simulation more than once, PLUMED will overwrite the previous output file and create a new `dihedrals_weak_restraint.dat`, but store the previous file as `bck.0.dihedrals_weak_restraint.dat`.

Looking at the output file by typing:

{% highlight git %}
head dihedrals_weak_restraint.dat
{% endhighlight %}

we see that we have written to the file the time (in ps), the $\phi$ and $\psi$ angles (in radians) and the value of the bias potential (in kJ/mol). 




## Adding a moving restraint 




## (Optional) Umbrella Sampling in PLUMED


