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

The negative exponential function is something that is 1 whenever the value of the exponent is zero, and gets progressively smaller when the value is larger than zero. (Trivially, the value of the exponential function can never be larger than 1 since lambda is a positive quantity and the distance $$\|X-X_i\\$$ is by definition positive). Whenever you sit exactly on a specific images $$X_j$$ then all the other terms in the sum disappear (if $$\lambda$$ is large enough) and only the value $$j$$ survives, returning exactly $$s=j$$. $j$ will have integer values from 1 to N where N is the number of configurations in your path. Therefore, the quantity $$s$$ will have values between 1 and N and descibes the progress along the path.

At this point, we should be more specific about what we mean by the "distance" between configurations, since a distance between two conformations can be calculated in several ways. In this tutorial, we will use the RMSD distance for a subset of atoms after optimal alignment. This means that at each step in which the analysis is performed, a number N of optimal alignments must be performed.

As described above, $$s$$ measures the progress along the path from state A to state B. Another useful variable is the distance away from the closest point along the path, which is denoted with the $$z$$. This is defined as 

$$z = -\frac{1}{\lambda} \ln \left[\sum_{i=1}^N \exp \left(-\lambda \|X-X_i\| \right]$$


## Generating a path 

