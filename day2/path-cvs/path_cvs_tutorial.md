# Path CVs

This tutorial will introduce you to the concept of a path collective variable (CV). Biasing a path CV, you can enhance the sampling along a chosen path along the free energy surface. This is a powerfull approach for simulating complex reaction processes or structural transitions between two states where the endpoints are known.   

Once this tutorial is completed students will be able to:

- Generate a putative reaction path using `pathtools`.
- Prepare a PLUMED input file implementing a path collective variable (Path CV). 
- Perform a metadynamics simulation biasing the progress along the path and assess the sampling efficiency.
- Improve the sampling by performing two-dimensional metadynamics biasing both the progress along the path and the distance from the path.

![restraint_figure](../../images/Image_pathCV.png)

## Background 

In the previous [metadynamics tutorial](../metadynamics/metadynamics.md), we described the conformation of alanine dipeptide using the two backbone dihedral angles ($$\phi$$) and ($$\psi$$). Although these two variables provide a useful description of the conformational free-energy landscape for this simple molecule, many molecular transitions involve coordinated changes in a much larger number of atomic coordinates. When you deal with a complex conformational transition that you want to analyze (or bias), very often you cannot just describe it with a single order parameter (or collective variable).

As an example, consider a large conformational transition like the activatation of a kinase via an open-close transition of the activation loop. In the sticks representation, you see the part of the molecule involved in the large conformational change, whereas the structure of the rest of molecule (shown only in backbone cartoon) remains mostly unchanged between the X-ray structures (PDB: 2C5X and 2C5Y). This is a complex transition and it is hard to tell which is the collective variable (or order parameter) that best describes the transition.

![Kinase figure](../../images/kinase_overlay.png)

We would like to completely divide these two configurations with a colletive variable that discriminates between what intuitively one would think as the open and closed conformation. We would also like to capture the progress from open to closed state through a hypothetical transition state to see certain crucial interactions forming/breaking so to better explain what is really happening during this process.  

A **path collective variable** provides a convenient way to describe such a transition. Rather than selecting one or two individual structural variables, we define a series of reference configurations that connect an initial state to a final state. Any instantaneous molecular configuration can then be described by two quantities: its progress along the reference path and its distance away from that path.

In a nutshell, your reaction might be very complex and involving many simultaneous degree of freedom, but intuitively can be tracked along a "reaction coordinate" along which the reaction proceeds. So what we need is a coordinate that, given a conformation, just tells which point along the "reaction coordinate" is closest.

![Path CV figure](../../images/PathCV_reaction_coord1.png)

For example, in the above figure, you see a typical chemical reaction (hydrolysis of methyl phosphate) with the two end-points denoted by A and B. If you are given a third point, just by looking at it, you might find that this resembles the reactant more than the product. Hypothetically, if we create a parameter that describes progress along the reaction coordinate, $$\xi$$, and set $$\xi=1$$ for a configuration in the A state and $$\xi=2$$ for a configuration in the B state, then the third point along the path, $$\xi$$, would probably be something like $$\xi=1.3$$.

Path collective variables are the extension of this concept to the case where you have multiple conformations that describe a path, starting from the reactant state A and moving to the product state B. Therefore, instead of an index that goes from 1 to 2 you have an index that goes from 1 to N , where N is the number of conformations that you use to describe your path.

Mathematically, the progress along the path is calculated with the following equation:

$$s = \frac{\sum_{i=1}^N i \exp \left(-\lambda \|X-X_i\| \right)}{\sum_{i=1}^N \exp \left(-\lambda \|X-X_i\| \right)}$$

where $$\|X−X_i\|$$ represents a distance between any instantaneous configuration $$X$$ tp be analyzed and another reference configuration $$X_i$$ from a set that compose the path made of $$N$$ configurations. The parameter $$\lambda$$ is a positive value that is tuned in a way explained later.

The negative exponential function is something that is 1 whenever the value of the exponent is zero, and gets progressively smaller when the value is larger than zero. (Trivially, the value of the exponential function can never be larger than 1 since lambda is a positive quantity and the distance $$\|X-X_i\|$$ is by definition positive). Whenever you sit exactly on a specific images $$X_j$$ then all the other terms in the sum disappear (if $$\lambda$$ is large enough) and only the value $$j$$ survives, returning exactly $$s=j$$. The index $j$ will have integer values from 1 to N where N is the number of configurations in your path. Therefore, the quantity $$s$$ will have values between 1 and N and descibes the progress along the path.

At this point, we should be more specific about what we mean by the "distance" between configurations, since a distance between two conformations can be calculated in several ways. In this tutorial, we will use the RMSD distance for a subset of atoms after optimal alignment. This means that at each step in which the analysis is performed, a number N of optimal alignments must be performed.

As described above, $$s$$ measures the progress along the path from state A to state B. Another useful variable is the distance away from the closest point along the path, which is denoted with the $$z$$. This is defined as 

$$z = -\frac{1}{\lambda} \ln \left[\sum_{i=1}^N \exp \left(-\lambda \|X-X_i\|\right) \right]$$

In the above equation, we see that in case of perfect match of $$X=X_i$$ this equation gives a value of $$z=0$$. The two variables $$s$$ and $$z$$, put together can be visualized as:

![PathCV_version2](../../images/PathCV_reaction_coord2.png) 

In the above figure, the $$s$$ variable can be thought as the length of the red segment, while the $$z$$ variable is the length of the green segment. 

Monitoring the $$z$$ variable in addition to $$s$$ is important because whenever your simulation is running close to the path (low Z values), then you know that you are reproducing reliably the path you provided, but if by chance you find some other path that goes from $$s=1$$ to $$s=N$$ via large $$z$$ values, then it might well be that you have just discovered a good alternative pathway. 

## Getting started

In this tutorial you will again study the C7eq to C7ax transition of the so-called alanine dipeptide molecule (N-acetyl-L-alanine-N'-methylamide):

![Alanine_dipeptide2](../../images/belfast-2-transition_ala.png)

Unlike the previous metadynamics tutorial, here we will demonstrate the use of path CVs to accomplish this transition. The example is particularly useful because we know that the Ramachandran space defined by the $$\phi$$ and $$\psi$$ angles will be a good descriptor of the conformational landscape. Therefore, we can assess how well are path is able to go from state A to state B by monitoring the $$\phi$$ and $$\psi$$ angles even though we will not be directly biasing them. 

**Files**
Files to complete this tutorial can be accessed here:
[tutorial files](coming soon)

These files are already located on bigzam:
/opt/workshop/path-CVs/ 

Use PuTTY to connect to bigzam as you did in previous [tutorials](../../day1/lj_fluid/lj_fluid_tutorial.md). Open PuTTY from the Window Start menu and enter `bigzam.local` for the Host Name. Login using the terminal using your username and password.

**Important**: Once connected to the workshop computer, set your environment variables by typing:

{% highlight git %}
source setup.sh
{% endhighlight %}

Copy the tutorial files by typing in the terminal:

In the terminal type:
{% highlight git %}
cp -r /opt/workshop/path-CVs/ ~/
{% endhighlight %}

This will copy the necessary tutorial files to your home directory on bigzam.

**Tip**: You can press the Tab key to automatically complete file and directory names. This can save time and help avoid typing errors.

Move into the path-CVs directory:

{% highlight git %}
cd ~/path-CVs
{% endhighlight %}

Within this directory you will find the following files:

- dialaA.pdb: A reference PDB structure file of the molecule
- alanine_dipeptide.gro: A GROMACS structure file (.gro)
- topol.top: A GROMACS topology file (.top)
- vacuum.mdp: A GROMACS parameter file (.mdp)


## Generating a path collective variable

Suppose we have limited information about the transition state or reaction path, but we know the starting and ending configuration. In our example molecule, this corresponds to the C7eq state and the C7ax state that we can visualize on a 2D Ramachandran plot:

![Ramachandran_plot_pathCV](../../images/Ramachandran_4PathCV.png)

In this case we might consider creating a path that starts at one basin and goes to the other along a straight line. To do this we need to generate equally spaced configurations that extrapolate between the two basins. A reasonable question to ask is how many frames do we need? The answer depends on the limiting scale in your reaction. For example, if in your process you have a torsion rotation as the smallest event that you want to capture with a path collective variable, then it is important that you mimic that torsion in the path and that this does not contain simply the initial and final point but also some intermediate. Similarly, if you have a concerted bond breaking, it might be that this takes place in the range of an Angstrom or so. In this case you should have intermediate frames that cover the sub-Angstrom scale. If you have both in the same path, then the smallest scale motion dominates.

In this example, you will use the `pathtools` feature in PLUMED to generate an initial linear path connecting two structure (pdb) frames representing the initial reactant and final product states. The path itself is an ordered set of equally-spaced, frames that interpolate between the reactant and product states.

In this case, I have provided a reference structure (pdb) file for the reactant state, `c7eq.pdb`, and product state `c7ax.pdb`.

In the following example, we will generate a linear path by typing in the terminal:

{% highlight git %}
plumed pathtools --start c7eq.pdb --end c7ax.pdb --nframes 10 --metric OPTIMAL --out linear_path.pdb
{% endhighlight %} 

  
