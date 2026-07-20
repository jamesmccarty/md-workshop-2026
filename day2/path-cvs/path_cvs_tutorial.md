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

In this case, I have provided a reference structure (pdb) file for the reactant state, `c7eq.pdb`, and product state `c7ax.pdb`. Pay careful attention to the format of these pdb files:

{% highlight git %}
ATOM      2  CH3 ACE     1      -3.220   0.160   1.920  1.00  1.00
ATOM      5  C   ACE     1      -1.800  -0.210   1.630  1.00  1.00
ATOM      6  O   ACE     1      -1.100  -0.890   2.410  1.00  1.00
ATOM      7  N   ALA     2      -1.400   0.350   0.520  1.00  1.00
ATOM      8  H   ALA     2      -1.940   0.980  -0.030  1.00  1.00
ATOM      9  CA  ALA     2      -0.060   0.000   0.060  1.00  1.00
ATOM     10  HA  ALA     2       0.060  -1.060   0.230  1.00  1.00
ATOM     11  CB  ALA     2      -0.020   0.050  -1.500  1.00  1.00
ATOM     15  C   ALA     2       1.160   0.770   0.730  1.00  1.00
ATOM     16  O   ALA     2       1.790   1.670   0.150  1.00  1.00
ATOM     17  N   NME     3       1.470   0.480   1.990  1.00  1.00
ATOM     18  H   NME     3       0.910  -0.190   2.480  1.00  1.00
ATOM     19  CH3 NME     3       2.650   1.120   2.570  1.00  1.00
END
{% endhighlight %}

You do not need to include every atom in the pdb file to define the path, but the atom number *must match the full structure** for your simulation. Notice here I am not including hydrogen atoms on the ACE group or NME group, so the atom numbers run from 2,5,6,7,8,9,10,11,15,16,17,18,19. Confirm that these match the atom numbers for the whole molecule we will use from the simulation:

{% highlight git %}
cat alanine_dipeptide.gro 
{% endhighlight %}

The atom numbers in the third column of the `alanine_dipeptide.gro` file.

**Important**: In your reference structure (pdb) file the last column is the occupancy number and this must be all 1.00 to be included in the optimal alignment. If these occupancy numbers are zero, then PLUMED cannot do the alignment. Therefore, always check to make sure the last column is all 1.00's. 
  
In the following example, we will generate a linear path by typing in the terminal:

{% highlight git %}
plumed pathtools --start c7eq.pdb --end c7ax.pdb --nframes 10 --metric OPTIMAL --out linear_path.pdb
{% endhighlight %} 

Here we are telling PLUMED to extrapolate a path of 10 equally-spaced frames between the c7eq.pdb structure and c7ax.pdb structure. The output file, `linear_path.pdb`, will be a pdb file with 12 total frames (one for each the start and the end state and 10 frames in between). 

If you want to visualize the path, I have included a python code, `convert_path2_pymol.py` that will write the path frames to a format that can be read by PyMOL. In the terminal type:

{% highlight git %}
python convert_path2_pymol.py linear_path.pdb c7ax.pdb linear_path_4pymol.pdb
{% endhighlight %}   

Then transfer the output `linear_path_4pymol.pdb` to your local Windows machine using the WinSCP app and load your structure in PyMOL.

After loading the `linear_path_4pymol.pdb` file into PyMOL, type the following into the consol to align the frames:

{% highlight git %}
intra_fit linear_path_4pymol, 1
{% endhighlight %}  

![figure_pymol_consol](../../images/pymol_consol.png)

You can then play through the frames by clicking on the play button on the bottom right. To see the movie slower you can set in the PyMOL toolbar: Movie –> Frame Rate –> 5 FPS.

![figure_pymol_linear_path](../../images/pymol_linear_path.png)

You should observe that our path smoothly rotates about the dihedral angles. Next, we need to determine the suitable choice for the $$\lambda$$ parameter in our path collective variable definition. A general rule of thumb is to use the following formula:

$$\lambda = \frac{2.3 N}{\sum_{i=1}^N \| X_i - X_{i+1}\|}$$ 

which implies that one should calculate the average distance between consecutive frames composing the path. Also, the distance between successive frames should be more or less similar between all frames. I have included a python code called `get-lambda-est.py` to assist in estimated the $$\lambda$$ paramter by evaluating the above equation. In the terminal, type the following:

{% highlight git %}
 python get-lambda-est.py linear_path.pdb 
{% endhighlight %} 

The output should be something like:

{% highlight git %}
Number of path frames: 12
Number of neighboring pairs: 11

Pair       RMSD (Å)       RMSD (nm)
01-02       0.125693       0.012569
02-03       0.125609       0.012561
03-04       0.125721       0.012572
04-05       0.125519       0.012552
05-06       0.125851       0.012585
06-07       0.125592       0.012559
07-08       0.125727       0.012573
08-09       0.125598       0.012560
09-10       0.125616       0.012562
10-11       0.125511       0.012551
11-12       0.125884       0.012588

Mean neighboring RMSD = 0.012567 nm
Suggested LAMBDA      = 183.03
Relative variation    = 0.1%
{% endhighlight %}

Here we see the frames are evenly spaced (similar RMSD between frames, and the suggested $$\lambda$$ based on these distances is $$\lambda=183$$. We are now ready to define our path CV in PLUMED. 

In PLUMED it is straightforward to implement a path CV with the following line:

{% highlight git %}
path: PATH REFERENCE=linear_path.pdb TYPE=OPTIMAL LAMBDA=183
{% endhighlight %}    

This defines a variable called `path` that takes as input a reference PDB file that containts the path frames. The path variable will have two components: `path.spath` - the distance along the path and `path.zpath` - the distance from the path.

## Metadynamics using a path CV 

I have provided a template PLUMED input file called `plumed-linear-path.dat` that will perform metadynamics on the $$s$$ path variable. Have a look at this file by typing:

{% highlight git %}
cat plumed-linear-path.dat
{% endhighlight %}  

{% highlight git %}
# set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

path: PATH REFERENCE=linear_path.pdb TYPE=OPTIMAL LAMBDA=183.0

metad: METAD ...
  ARG=path.spath
  SIGMA=0.05 
  HEIGHT=1.5 
  TEMP=300 
  BIASFACTOR=12 
  PACE=500
  FILE=HILLS GRID_MIN=0 GRID_MAX=12.5
...

PRINT ARG=phi,psi,path.spath,path.zpath,metad.bias STRIDE=100 FILE=COLVAR_linear_path.dat
{% endhighlight %}  

The first lines define the usual $$\phi$$ and $$\psi$$ dihedral angles that we can use to assess the transition from reactant to product state. 

The next line defines the path CV with the variable name `path`. We then setup a metadynamics bias similar to the previous [metadynamics tutorial](../metadynamics/metadynamics.md). Here, the `ARG=path.spath` specifies that we will perform metadynamics on the $$s$$ component of the path (distance along the path). Finally, we print the output to a file called `COLVAR_linear_path.dat`. 

To run this metadynamics simulation type:

{% highlight git %}
gmx grompp -f vacuum.mdp -c alanine_dipeptide.gro -p topol.top -o path-mdrun1.tpr

gmx mdrun -v -deffnm path-mdrun1 -plumed plumed-linear-path.dat
{% endhighlight %}

When this job finishes, the output file will be `COLVAR_linear_path.dat`. Transfer this file from bigzam to your local Windows machine using WinSCP. The columns of this file are specified by the `#! FIELDS` line:

{% highlight git %}
#! FIELDS time phi psi path.spath path.zpath metad.bias
{% endhighlight %} 

The first column is the simulation time (ps), the second and third column are the $$\phi$$ and $$\psi$$ angles (used for monitoring the progress of the reaction, but not directly biased), and the fourth and fifth column are the $$s$$ and $$z$$ path variables. Finally, the sixth column is the metadynamics bias. Upload your `COLVAR_linear_path.dat` file to the following [Google Colab link](https://colab.research.google.com/drive/1N4rXPjY5-O4Hhe7dbcL2sH7SRZbHU7eq?usp=sharing).

First we plot the path $$s$$ variable vs. time to examine how the system moves along the chosen reaction path.

![Figure_linear_spath](../../images/linear_path_spath_progress.png)

The trajectory clearly transitions between low values of $$s$$ (approximately 3–4) and high values (approximately 9–10), demonstrating that metadynamics successfully drives the system between the two endpoint states. Long periods where $$s$$ fluctuates within a narrow range correspond to the system residing in one metastable basin, while the sharp transitions indicate barrier crossings promoted by the increasing metadynamics bias.

Next, we can compare this with the $$\phi$$ vs. time plot:

![Figure_linear_spath_phi](../../images/linear_path_phi_angle.png)

Notice that transitions in the $$s$$ variable occur at the same time as jumps in the $$\phi$$ angle, which also indicates that the system is changing from the C7eq state to the C7ax state. This is confirmed by plotting the 2D Ramachandran plot where we see two sampled regions corresponding to the two metastable states. 

![Figure_linear_spath_2D](../../images/linear_path_Ramachandran.png)

Notice, that although we biased along a linear path from the C7eq to the C7ax state, the actual transition pathway does not appear to be a straight line in $$\phi$$ and $$\psi$$ space. Our linear extrapolated path captures the overall transition but is only an approximation to the underlying transition path. 

Finally, it is instructive to examine the 2D plot of the distance from the path ($$z$$) vs. the progress along the path ($$s$$). Here we see a broad loop in the $$z$$ coordinate near the middle of the path and the spread in $$z$$ is larger in the middle of the path than near the endpoints. This means the system is going through paths betwen A and B that are different from our putative path, $$s$$.  

![Figure_linear_2D_s_z](../../images/linear_path_2D.png)

This is an important point that because of the way metadynamics builds the bias to encourage exploration of previously unvisited configurations, the system may sample other paths different from the chosen path 

