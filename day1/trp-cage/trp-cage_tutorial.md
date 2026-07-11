### Principal Component Analysis

A very common analysis method is to extract the “principal” or “essential” motions that have the largest amplitudes and involve the largest parts of the structure. Principal component analysis (PCA) of the trajectory, which is sometimes also referred to as ‘essential dynamics’ (ED), aims at identifying large scale collective motions of atoms and thus reveal the structures underlying the atomic fluctuations. The fluctuations of particles are correlated due to coupled interactions between particles. The degree of correlation will vary and notably particles which are directly connected through bonds or lie in the vicinity of each other will move in a concerted manner. The correlations between the motions of the particles give rise to collective motions in the system that is often directly related to its function or (bio)physical properties. The study of the structure of the atomic fluctuations can give valuable insight in the behavior of a macromolecule.

The first step in PCA is the construction of the covariance matrix, which captures the degree of collinearity of atomic motions for each pair of atoms. This matrix is subsequently diagonalized, yielding a matrix of eigenvectors and a diagonal matrix of eigenvalues. Each of the eigenvectors describes a collective motion of particles, where the components of the vector indicate how much the corresponding atom participates in the motion. The associated eigenvalue is a measure of the total motility associated with an eigenvector. Usually most of the motion in the system (>90%) is described by less than 10 eigenvectors or principal components. Since the covariance analysis produces a lot of files, the analysis is best performed in a subdirectory below the directory of the MD run:

{% highlight git %}
mkdir COVAR

cd COVAR
{% endhighlight %}

The program `covar` will construct the covariance matrix and perform the diagonalization. Type the following command:


The program `covar` will construct the covariance matrix and perform the diagonalization. Type the following command:

{% highlight git %}
gmx covar -s ../protein-reference.pdb -f ../protein_noPBC.xtc -o eigenvalues.xvg -v eigenvectors.trr
{% endhighlight %}

