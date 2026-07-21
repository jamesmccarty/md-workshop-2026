
# Day 2 - Enhanced Sampling Methods and Statistical Analysis

The goals of today's sessions are to introduce enhanced sampling methods and statistical analysis techniques commonly used to study conformational changes and free energy landscapes in molecular simulations.  

**Learning Objectives**

By the end of today, participants should be able to:

- Explain why conventional molecular dynamics simulations can be insufficient for studying rare events.
- Describe the basic principles of enhanced sampling methods.
- Run a metadynamics simulation using PLUMED
- Construct and interpret free energy surfaces.
- Estimate statistical uncertainties of computed histograms and free energy surfaces.


## Session 1: Introduction to Enhanced Sampling

Slides coming soon!

## Session 2: Using Restraints 

In this tutorial you will learn how to implement a restraint in PLUMED. By applying restraints on one or more collective variables, you can restrict the sampling to a specified region. This will introduce you to the idea of using a bias to manipulate a simulation on-the-fly.

Once this tutorial is completed students will be able to: 

- Apply a restraint on a simulations over one or more collective variables
- Understand the effect of a restraint on the acquired statistics
- Perform an out-of-equilibrium simulation with a moving restraint 

**Tutorial**

[Using Restraints Tutorial](day2/restraints/restraints.md)

## Session 3: Metadynamics Tutorial

In this tutorial you will learn how to perform a metadynamics enhanced sampling simulation using PLUMED. We will use the alanine dipeptide molecule as a simple model system for exploring conformational transitions and free energy calculations.

Once this tutorial is completed students will be able to:

- Define and assess collective variables for enhanced sampling
- Run a metadynamics simulation using PLUMED
- Calculate free energies from a metadynamics simulation

**Tutorial**

[Metadynamics Tutorial](day2/metadynamics/metadynamics.md)

## Session 4: Path Collective Variables Tutorial

The primary goals of this tutorial are to show you how you might use path collective variables for describing and simulating activated molecular processes, and to provide hands-on experience with setting up, running, and analyzing biased path-CV simulations, such as path-metadynamics. We will study alanine dipeptide with path collective variables to demonstrate how one might perform enhanced sampling along a putative reaction coordinate.

Once this tutorial is completed students will be able to:

- Generate a putative reaction path using `pathtools`.
- Prepare a PLUMED input file implementing a path collective variable (Path CV).
- Perform a metadynamics simulation biasing the progress along the path and assess the sampling efficiency.
- Improve the sampling by performing two-dimensional metadynamics biasing both the progress along the path and the distance from the path.

**Tutorial**

[Path CVs Tutorial](day2/path-cvs/path_cvs_tutorial.md)

## Session 5: Error Analysis and Quantifying Uncertainty from MD simulations 

The goals of this tutorial are to The aim of this tutorial are for you to understand how to analyze statistical uncertainties in measured quanties from MD simulations.

Once this tutorial is completed students will be able to:

- Use bootstrapping to estimate the uncertainty and confidence intervals from computed histograms from an MD simulation.
- Compare distributions measured from separate independent MD simulations to see if the two trajectories are sampling the same distribution.
- Perform a block-analysis to calculate error bars on the free energy surface claculated from a metadynamics trajectory.

**Tutorial**

[Error Analysis and Quantifying Uncertainty from MD simulations](day2/error_analysis/error_analysis_plumed.md) 


