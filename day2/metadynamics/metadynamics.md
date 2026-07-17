
# Metadynamics

In this tutorial you will learn how to perform a metadynamics enhanced sampling simulation using PLUMED. We will use the alanine dipeptide molecule as a simple model system for exploring conformational transitions and free energy calculations.

Once this tutorial is completed students will be able to:

- Define and assess collective variables for enhanced sampling
- Run a metadynamics simulation using PLUMED
- Calculate free energies from a metadynamics simulation

![figure_goes_here](../../images/)

**Files**
Files to complete this tutorial can be accessed here:
[tutorial files](coming soon)

These files are already located on bigzam:
/opt/workshop/metadynamics/

## Getting Started

Use PuTTY to connect to bigzam as you did in previous [tutorials](../../day1/lj_fluid/lj_fluid_tutorial.md). Open PuTTY from the Window Start menu and enter `bigzam.local` for the Host Name. Login using the terminal using your username and password.

**Important**: Once connected to the workshop computer, set your environment variables by typing:

{% highlight git %}
source setup.sh
{% endhighlight %}

Copy the tutorial files by typing in the terminal:

In the terminal type:
{% highlight git %}
cp -r /opt/workshop/metadynamics/ ~/
{% endhighlight %}

This will copy the necessary tutorial files to your home directory on bigzam.

**Tip**: You can press the Tab key to automatically complete file and directory names. This can save time and help avoid typing errors.

Move into the using_restraints directory:

{% highlight git %}
cd ~/metadynamics
{% endhighlight %}

Within this directory you will find the following files:

- dialaA.pdb: A reference PDB structure file of the molecule
- alanine_dipeptide.gro: A GROMACS structure file (.gro)
- topol.top: A GROMACS topology file (.top)
- vacuum.mdp: A GROMACS parameter file (.mdp)

## Metadynamics

In the previous tutorial on [using restraints](../day2/restraints/restraints.md), you learned how to use PLUMED to bias an MD simulation on-the-fly. In this tutorial we will use metadynamics to accelerate the sampling and reconstruct the unbiased free energy surface. Here you can find a brief recap of the metadynamics theory. 

[Coming soon]

## Bias along a single coordinate

In this exercise we will setup and perform a well-tempered metadynamics run using the backbone dihedral $$\phi$$ as collective variable. During the calculation, we will also monitor the behavior of the other backbone dihedral $$\psi$$.

In the workshop directory  you can find a sample PLUMED input called `plumed_metad_phi.dat` that you can use as a template. Whenever you see an highlighted `__FILL__` string, this is a string that you must replace.

Have a look at the contents of `plumed_metad_phi.dat` with:

{% highlight git %}
nano plumed_metad1.dat
{% endhighlight %}

The first several lines define the collective variables $$\phi$$ and $$\psi$$. 
For help defining the $$\phi$$ and $$\psi$$ angles, you might need to refer back to the previous tutorial on [Introduction to PLUMED syntax](../../day1/intro_plumed_syntax/analysis.md).

{% highlight git %}
MOLINFO STRUCTURE=dialaA.pdb
# Compute the backbone dihedral angle phi, defined by atoms C-N-CA-C
# you might want to use MOLINFO shortcuts
phi: TORSION ATOMS=__FILL__
# Compute the backbone dihedral angle psi, defined by atoms N-CA-C-N
# here also you might want to use MOLINFO shortcuts
psi: TORSION ATOMS=__FILL__
{% endhighlight %}

The next lines initiate metadynamics with the [METAD](https://www.plumed.org/doc-v2.9/user-doc/html/_m_e_t_a_d.html) action. Read this section carefully; it is important to understand what each of these input options means.

{% highlight git %}
# Activate well-tempered metadynamics in phi
metad: METAD ...
    ARG=phi
   # Deposit a Gaussian every 500 time steps
    PACE=__FILL__ 
    # set initial Gaussian height equal to 1.2 kJ/mol  
    HEIGHT=__FILL__
     # set Gaussian width (sigma)
    SIGMA=0.1
    # The bias factor should be wisely chosen 
    BIASFACTOR=8
   # Gaussians will be written to HILLS file and stored on grid 
   FILE=HILLS GRID_MIN=-pi GRID_MAX=pi
...
{% endhighlight %}

The above lines tell PLUMED that we will run metadyanics on the $$\phi$$ angle. The `PACE` sets the Gaussian deposition rate (how frequently a new Gaussian is deposited). Set this to `PACE=500`, meaning we will add a new Gaussian every 500 MD steps. The `SIGMA` and `HEIGHT` section specify the Gaussian width and height, respectively. The units here for `SIGMA` will be in radians (since $$\phi$$ is in radians) and for `HEIGHT` is in energy units (kJ/mol). In this example, set `SIGMA=0.1` radians and `HEIGHT=1.2` kJ/mol. The `BIASFACTOR` dictates the level of sampling. In so-called **well-tempered metadynamics** the biasfactor determines the rescaling of the initial Gaussian height so that the bias potential smoothly converges to a fixed level in long time limit. 

A good rule of thumb is that the `BIASFACTOR` should be large enough to accelerate barrier crossing, but not so large that you overfill the free energy landscape or require prohibitively long convergence times. `BIASFACTOR=10` is usually a  good starting point. A useful way to think about the BIASFACTOR is in terms of the free energy barrier height you expect ($$\Delta G^{\dagger}$$). Suppose your free energy barrier is about $$\Delta G^{\dagger} = 10 RT$$. With a biasfactor = 10, the effective barrier becomes $$\sim 1 RT$$, making transitions frequent while still allowing the bias to converge smoothly. 

It is computationally more efficient to use a grid to store the accumulated bias. The following lines within the METAD action tell PLUMED to store the Gaussian bias kernels in a file called `HILLS` that will be constructed on a grid with range $$-\pi$$ to $$\pi$$. 

{% highlight git %}
FILE=HILLS GRID_MIN=-pi GRID_MAX=pi
{% endhighlight %}

Finally, we print everything to the file `COLVAR_phi` with the PRINT action. 

{% highlight git %}
PRINT ARG=phi,psi,metad.bias FILE=COLVAR_phi STRIDE=100
{% endhighlight %} 

**Important**: When running a metadynamics simulation, don't forget to PRINT the metadynamics bias (metad.bias) because you will need this for post-processing.

  
