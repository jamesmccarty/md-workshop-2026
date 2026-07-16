# Implementing Restraints

In this tutorial you will learn how to implement a restraint using PLUMED. By applying restraints on one or more collective variables, you can restrict the sampling to a desired region. This will introduce you to the idea of using a bias to manipulate a simulation on-the-fly. 

Once this tutorial is completed students will be able to:

- Apply a restraint on a simulations over one or more collective variables
- Understand the effect of a restraint on the acquired statistics
- Perform an out-of-equilibrium simulation with a moving restraint

 
![restraint_figure](../../images/)

In the [previous introduction to PLUMED syntax tutorial](../../day1/intro_plumed_syntax/analysis.md), you used PLUMED as a post-processing tool to  calculate properties of the system on a previously generated MD simulation trajectory. In addition to being a post-processing tool, PLUMED can also interface with the MD code "on-the-fly" during a MD simulation. The atomic coordinates and atoms at a given instant can be passed to the PLUMED code to manipulate a simulation on-the-fly. This allows one to bias the simulation by adding extra constraints in addition to the standard force field terms. The additional energy terms are usually referred as **Bias**. 

## The molecule of the day



## Adding a constant bias potential 
In the following we will see how to apply a constant bias potential



## Adding a moving restraint 




## (Optional) Umbrella Sampling in PLUMED


