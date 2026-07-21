# Error Analysis and Quantifying Uncertainty from MD simulations

By now you know how to measure observables from an MD trajectory using PLUMED and and how to plot a histogram of observable that represents the probability distribution of that given observable (see Tutorial on [A Brief Introduction to PLUMED Syntax and Making Histograms](../../day1/intro_plumed_syntax/analysis.md)). Also, we have seen how to perform a metadynamics simulation to enhance the sampling of configuration space and how to calculate the reweighted free energy surface along a chosen order parameter (see Tutorial on [Metadynamics](../metadynamics/metadynamics.md)). In this tutorial, we now attempt to quantify how we know what we observe from an MD simulation trajectory is **statistically meaningful**? 

Once this tutorial is completed students will be able to:

- Use bootstrapping to estimate the uncertainty and confidence intervals from computed histograms from an MD simulation.
- Compare distributions measured from separate independent MD simulations to see if the two trajectories are sampling the same distribution. 
- Perform a block-analysis to calculate error bars on the free energy surface claculated from a metadynamics trajectory. 

**Files**
Files to complete this tutorial can be accessed here:
[tutorial files](coming soon)

Because we are analyzing trajectories, you not need bigzam for this tutorial. This tutorial will be completed using Python notebooks. 

## Estimating uncertainty from a MD trajectory

In this first example, I have run a 1-$$\mu$s MD simulation trajectory of the RC9 peptide, GGKGMGFGL, shown here:

![Figure_RC9](Figure_RC9.png)

The file `RC9_distance.dat` contains the distance between atom HA1 on Gly1 and atom HN on Leu9. For transparency, I created this data using the PLUMED driver with the following input:

{% highlight git %}
d1: DISTANCE ATOMS=6,98 # distance between atom HA1 on Gly1 and HN on Leu9

PRINT ARG=d1 STRIDE=80 FILE=distance.dat
{% endhighlight %}

This data is provided for download [here](https://drive.google.com/file/d/1j48wvCrlDujbNO4fM-r1iEl4w4065Agx/view?usp=sharing). Use this link to download the data file to your Windows machine. 



## Comparing histograms from two independent MD simulations 

## Error bars on free energy surfaces from metadynamics  
