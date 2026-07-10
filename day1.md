
# Day 1 - Foundations of Molecular Simulation  

The goals of today's sessions are to introduce the basic concepts of molecular dynamics (MD) simulations and provide hands-on experience running and analyzing simulations using the GROMACS MD package.  

**Learning Objectives**

By the end of today, participants should be able to:

- Explain the basic principles of molecular dynamics simulations.
- Setup and run a simple molecular simulation using GROMACS.
- Analyze simulation trajectories and interpret common quantites such as energy, RMSD, and radial distribution function.
- Become familiar with the basic syntax of PLUMED
- Be able to calculate simple collective variable use them to analyze existing trajectories

## Session 1: Introduction to Molecular Dynamics

Slides coming soon!

## Session 2: Simulating a Simple Fluid Tutorial

In this tutorial you will learn how to perform a molecular dynamics (MD) of liquid argon using the GROMACS MD package. This will introduce you to the basics of an MD simulation and how to analyze structural properties of a fluid.  

Once this tutorial is completed students will be able to:

- Understand the basic format of GROMACS structure, topology, and parameter files
- Prepare a simple fluid simulation and run a molecular dynamics simulation
- Calculate the radial distribution function for a fluid and connect to scattering experiments
- Calculate the mean square displacement and use the plot of MSD vs. time to estimate the diffusion coefficient

**Tutorial**

[Fluid Tutorial](day1/lj_fluid/lj_fluid_tutorial.md)

## Session 3: MD Simulation of a Mini-Protein Tutorial

In this tutorial you will use GROMACS in the Unix command line to build and equilibrate a solvated protein system for the mini-protein [YYDPETGTWY](https://pubs.acs.org/doi/10.1021/ja8030533). The goal is to understand the MD simulation parameters that are chosen and the details about how the system was built.  

Once this tutorial is completed students will be able to:

- Simulate a protein in explicit solvent
- Understand GROMACS output files and be able to transfer files for analysis 
- Calculate RMSD and monitor structural stability
- Analyze protein dynamics

**Tutorial**

[Mini-Protein Tutorial](day1/mini-protein/mini-protein_tutorial.md)

## (Optional) Trp-Cage System

For additional practice, you can try setting up and running an MD simulation of the larger [trp-cage protein](https://www.rcsb.org/structure/1L2Y).
   

## Session 4: Introduction to PLUMED Syntax and Histograms

The aim of this tutorial is to introduce the users to the PLUMED syntax. We will go through the writing of simple collective variable and we will use them to analyze existing trajectories.

Once this tutorial is completed students will be able to:

- Write a simple PLUMED input file and use it to analyze a trajectory
- Print collective variables such as distances, torsional angles, and coordination numbers using the PRINT action
- Be able to use PLUMED to calculate ensemble averages and histograms

**Tutorial**

[Analyzing Trajectories Using PLUMED Tutorial](day1/intro_plumed_syntax/analysis.md)   
