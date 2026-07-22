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

In this first example, I have run a 1-$$\mu$$s MD simulation trajectory of the RC9 peptide, GGKGMGFGL, shown here:

![Figure_RC9](Figure_RC9.png)

The file `RC9_distance.dat` contains the distance between atom HA1 on Gly1 and atom HN on Leu9. For transparency, I created this data using the PLUMED driver with the following input:

{% highlight git %}
d1: DISTANCE ATOMS=6,98 # distance between atom HA1 on Gly1 and HN on Leu9

PRINT ARG=d1 STRIDE=80 FILE=distance.dat
{% endhighlight %}

This data is provided for download [here](https://drive.google.com/file/d/1j48wvCrlDujbNO4fM-r1iEl4w4065Agx/view?usp=sharing). Use this link to download the data file to your Windows machine. 

Next, we will use Python to generate a histogram with error bars. Upload your data file `distance.dat` to the [Colab code here](https://colab.research.google.com/drive/10qtu-z5JoIigB5oEkyZfaieOSw37oupp?usp=sharing)

The beginning of this document is similar to what you have done before in earlier tutorials. First you upload your data file and store each column as an array:

![Figure_bootstrap_notebook1](../../images/Figure_bootstrap_notebook1.png)

Then we calculate the mean and standard deviation as follows:

{% highlight git %}
d_mean_value = np.mean(d)
d_std_value = np.std(d)
{% endhighlight %}

What is the value of the mean distance between HA1 and HN in our simulation?

As before we plot the distance as a function of time across the 1-$$\mu$$s MD simulation (left) and the histogram of that distance (right):

![Figure_distance_hist_bootstrap2](../../images/Figure_distance_hist_bootstrap2.png)

The histogram gives an estimate of the underlying **probability distribution** of the distance. However, becauese the MD simulation contains only a finite number of frames, the histogram itself is subject to statistical uncertainty. Suppose we wanted to quantify the uncertainty of the height of each bar in the above histogram. How could we do this? 

One way to answer this question is through **bootstrap analysis**. [Bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_(statistics) is a technique widely used in statistics for estimating the uncertainty of a distribution. Bootstrapping estimates the uncertainty by repeatedly generating new datasets from the original trajectory. The basic idea is simple:

- Randomly select frames from our original trajectory **with replacement** until a new dataset containing the same number of frames has been created.
- Compute the histogram for this new dataset
- Repeat this process many times (typicaly hundreds or thousands of times).
- Use the variation among these histograms to estimate a confidence interval for the probability distribution. 

The advantage of this approach is that we do not make any assumptions about the probabily distribution that generated our original data set. 

A simple four lines of Python code will perform the bootstrapping analysis:

{% highlight git %}
for i in range(n_bootstrap):
    sample = rng.choice(d, size=n, replace=True) # sampling with replacement
    boot_hist[i], _ = np.histogram(sample, bins=bins, density=False) # generating new histogram 
    boot_hist[i] = boot_hist[i] / np.sum(boot_hist[i]) # normalize 

# confidence intervals (95% interval lies between the 2.5th percentile and 97.5th percentile)
lower, upper = np.percentile(boot_hist, [2.5, 97.5], axis=0)
{% endhighlight %}

In the Colab notebook provided, we choose 1000 bootstrap samples. Note that because each bootstrap sample is drawn from the original trajectory with replacement, some frames will randomly be chosen multiple times, while others may not be chosen at all. Repeating the resampling proceedure many times allows us to estimate how much the histogram would vary if we repeated the MD simulation under similar conditions. 

Use the Colab notebook to plot the distance distribution with a 95% confidence interval for the distance distribution as a shaded region:

![Figure_histogram_bootstrap](../../images/Figure_histgram_bootstrap.png)

**Main idea**: If you were to repeat the MD simulation 100 times under the same conditions, generating 100 similar distance distributions. The shaded region shows the range where you would expect the histogram from 95 of those simulations to lie.  
 
# Comparing histograms from independent MD simulations 

Because the equations of motion that generate a molecular dynamics trajectory are nonlinear, two simulations started from the same initial state but with different random initial velocities will diverge and follow different trajectories. If each simulation is sufficiently long and **samples the same equilibrium ensemble**, the distributions of observables (such as interatomic distances) should be statistically indistinguishable, even though the individual trajectories are different.  

For this reason, it is good practice to perform multiple independent MD simulations using different randomized initial velocies. These independent replicas provide a way to assess the **reproducibility** of the simulation results and determine if the sampled distributions are consistent with one another. 

A second common application is to **compare simulations performed under different conditions**. For example, you may wish to compare a wild-type protein with a mutant, compare a reaction under different environmental conditions, or compare ligand binding with a functional group modification. In these cases, the goal is to determine whether the underlying probability distribution of some observable has changed in a statistically meaningful way.    

In this tutorial, we will use the **Kolmogorov-Smirnov (KS) test** to compare distance distributions from three independent simulations. The KS test provides an objective measure of whether the two datasets are consistent with having been drawn from the same underlying probability distribution. 

We will again use the RC9 peptide from the previous example as our model system. I have performed three independent 1 $$\mu$$s MD simulations of the peptide. Two of the simulations were run under identical conditions, and therefore should produce statistically indistinguishable distance distributions. However, in the the third MD simulation, I implemented a harmonic restraint applied to the peptide end-to-end distance, altering the conformational sampling. 

Our goal is to determine whether we can identify the restrained simulation based solely on the distance distributions. Rather than visualzing the histograms, we will use the KS test to determine whether the observed differences are statistically significant. 

## The KS test

Instead of comparing two histograms directly, the KS test compares the **cumulative distribution functions (CDFs)** of two datasets. The cumulative distribution is the probability that an observable is less than or equal to a given value. The KS statistic, $$D$$, is simply the **largest vertical distance** between two cumulative distribution functions as shown in the figure below:

![Figure_CDF_KS](../../images/figure_KS_CDF_example.png)

Conceptually, we can understand that if two distrivutions are identical, their CDFs will overlap and $$D$$ will be small. Very different distributions will have larger separation in the CDF, giving larger values of $$D$$. 

After computing the KS statistic, the test calculates a $$p-value$$. The p-value tells you how likely it would be to observe a distance at least this large by random chance if the two datasets were actually sampled for the same underlying probability distribution. A general rule is that a p-value $$>$$ 0.05 means that we cannot reject the hypothesis that the two simulations sample the same distribution. On the othe hand, a p-value $$ <$$ 0.05 means that the two simulations sample different distributions. 

Performing the KS test is quite easy in Python using the `scipy.stats` library. Here we should use the two-sample Kolmogorov-Smirnov test, [ks_2samp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html) since we are comparing two different data sets. If we have two sets of data, `trialA` and `trialB`, we can run the KS test in Python with the lines:

{% highlight git %}
from scipy.stats import ks_2samp

result = ks_2samp(trialA, trialB)

{% endhighlight %} 

The result will have two components that we can print, the KS statistic ($$D$$) and the p-value:


{% highlight git %}
print(result.statistic)
print(result.pvalue)
{% endhighlight %}

We can perform this test on any pair of trials to decide if they sample the same distribution or are statistically different. 

In this tutorial, you can find a link to download the three trial datasets here:

[distance_trials.tar.gz](https://drive.google.com/file/d/18J2AnzSspJBqxU-M5GfqNxeIYNOp2xyQ/view?usp=sharing)

Download the compressed archive containing the three trials to your local Windows machine. Then upload your files to the Colab notebook here for running the KS test:

[KS test notebook](https://colab.research.google.com/drive/1GUF6SWx2ibZVbENVv-pmxqOzApnDeAUz?usp=sharing) 

We unzip the archived data files with the line: 
{% highlight git %}
!tar -xzf distance_trials.tar.gz
{% endhighlight %}

Then store the three datasets as `trialA`, `trialB`, and `trialC`. 

Performing the KS test on `trialA` and `trialB` gives the following:

{% highlight git %}
KS statistic = 0.7463
p-value = 3.6019e-272
{% endhighlight %}

Here we see the p-value is much smaller that 0.05 indicating the distributions are statistically different. We can confirm this by plotting the CDFs for the two datasets. Notice that the CDFs hardly overlap and the vertical difference between them is large. 

![Figure_CDF_AB](../../images/Figure_CDF_AB.png)

For comparison, perform the KS test on `trialB` and `trialC` gives a small KS statistic $$D=0.0265$$ and a p-value = 0.8 which is well above 0.05. We conclude that data from `trialB` and data from `trialC` are drawn from the same underlying distribution. This means that averages obtained from these two independent MD simulations should be statistically equivalent. 

Plotting the CDFs for the two trials, we see a large degree of overlap as expected for such a large p-value:

![Figure_CDF_BC](../../images/Figure_CDF_BC.png)

Based on the KS test analysis we should conclude the following: trial B and trial C represent samples drawn from the same underlying distribution; therefore, these trials must be the two identical unrestrained MD simulation samples. On the other hand, trial A is drawn from a different underlying distribution and must correspond to the restrained MD simulation. We conclude that the restraint shifted the underlying distribution of sampled distances.      

In fact we can readily see this just by looking at the histograms of the three trials, where we see that the restraint has shifted the sampled distances to much lower values:

![Figure_histogram_overlay](../../images/Figure_histogram_overlay.png)

 
## Error bars on free energy surfaces from metadynamics  
