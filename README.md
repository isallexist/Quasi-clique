# Quasi-clique
This repository contains codes and instances used in the computational study reported in the article "An Ellipsoidal Bounding Scheme for the Quasi-Clique Number of a Graph" that has been accepted for publication in INFORMS Journal of Computing (citation will be updated when DOI information and URL are available).

@article{MiaoBala2019quasiclique,
Author = {Zhuqi Miao and Balabhaskar Balasundaram},
Journal = {{INFORMS} Journal on Computing},
Month = {July},
Note = {Accepted for publication.},
Title = {An Ellipsoidal Bounding Scheme for the Quasi-Clique Number of a Graph},
Year = {2019}}

User instructions: Packages required to run codes-- not provided here, other licenses required; compile and run instructions-- input parameters, etc.

LDB.cpp: The code for generating the proposed ellipsoidal bound.

MIP.cpp: The code for generating the MIP-based bounds by solving MIP formulations for the quasi-clique problem.

MakeCommand_ldb: Command to compile LDB.cpp in a linux system.

Makefile_mip: Make file to compile MIP.cpp in a linux system.

Input instance format(s): DIMACS II instance format. Examples can be found at http://iridia.ulb.ac.be/~fmascia/maximum_clique/DIMACS-benchmark

