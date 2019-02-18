# msc-thesis-code
Code for [Molecular movies and geometry reconstruction using Coulomb explosion imaging](https://uwspace.uwaterloo.ca/handle/10012/12190) (and thesis plots).

Other relevant repositories:
* [mmmgrubs](https://github.com/ali-ramadhan/mmmgrubs): MMMGRUBS: Molecular Motion Movies and Geometry Reconstruction Using Bayesian Statistics. This was my attempt at reconstructing molecular geometries using Bayesian inference and inverse modeling (in Stan). I couldn't get the Markov chain to sample at all, so it doesn't work right now.
* [multi-start-4ion](https://github.com/ali-ramadhan/multi-start-4ion): Creates molecular movies of 4-atom molecules from Coulomb explosion imaging experiments. This was an extension of the triatomic molecular geometry reconstruction to 4-atom molecules. It also uses MATLAB's MultiStart global optimization solver to locate multiple local solutions in hopes of finding the global solution. I tested it on acetylene (C2H2) and it was able to reconstruct some geometries but convergence was difficult. In hindsight this may have been because I used a logarithmic objective function (which comes with a singularity). A logarithmic objective function worked for triatomic molecules (3-dimensional parameter space) but might work badly in higher dimensional spaces (i.e. for 4+ atom molecules).
* [msc-thesis](https://github.com/ali-ramadhan/msc-thesis): Repo with LaTeX used to typeset thesis.
