
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

## Session 2: Implementing Restraints 

In this tutorial you will learn how to implement a restraint using PLUMED. By applying restraints on one or more collective variables, you can restrict the sampling to a desired region. This will introduce you to the idea of using a bias to manipulate a simulation on-the-fly.

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

[Metadynamics Tutorial](coming soon)

## Session 4: Path Collective Variables Tutorial

The primary goals of this tutorial are to show you how you might use path collective variables for describing and simulating activated molecular processes, and to provide hands-on experience with setting up, running, and analyzing biased path-CV simulations, such as path-metadynamics. We will study alanine dipeptide with path collective variables to demonstrate how one might perform enhanced sampling along a putative reaction coordinate.

Once this tutorial is completed students will be able to:

- Build a Path-CV from a series of MD simulation snapshots.
- Write a PLUMED input file to run a path-metadynamics simulation.

[Path CVs Tutorial](coming soon)

## Session 5: Error Analysis of Free Energy Surfaces Tutorial 

The goals of this tutorial are to The aim of this tutorial are for you to understand how to analyze free energy surfaces, compute error bars, and understand statistical uncertainties from MD simulation results. 

Once this tutorial is completed students will be able to:

- Assess sampling and convergence of a metadynamics simulation
- Be able to use PLUMED to perform block analysis of trajectory data
- Be able to use PLUMED to calculate unbiased ensemble averages and histograms from biased trajectories
- Be able to perform bootstrap resampling to obtain error bars

[Analysis of Free Energy Surfaces Tutorial](coming soon) 
