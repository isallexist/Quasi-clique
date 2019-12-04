# Quasi-clique
This repository contains codes and instances used in the computational study reported in the article *"An Ellipsoidal Bounding Scheme for the Quasi-Clique Number of a Graph"* that has been accepted for publication in INFORMS Journal of Computing (citation will be updated when DOI information and URL are available). If you wish to use or cite this code, please cite the paper:
```
@article{MiaoBala2019quasiclique,
Author = {Zhuqi Miao and Balabhaskar Balasundaram},
Journal = {{INFORMS} Journal on Computing},
Month = {July},
Note = {Accepted for publication.},
Title = {An Ellipsoidal Bounding Scheme for the Quasi-Clique Number of a Graph},
Year = {2019}}
```
User instructions: Packages required to run codes-- not provided here, other licenses required; compile and run instructions-- input parameters, etc.

LDB.cpp: The code for generating the proposed ellipsoidal bound.

MIP.cpp: The code for generating the MIP-based bounds by solving MIP formulations for the quasi-clique problem.

MakeCommand_ldb: Command to compile LDB.cpp in a linux system.

Makefile_mip: Make file to compile MIP.cpp in a linux system.

Input instance format(s): DIMACS II instance format. 

# Compiling the code
1. Configure the Linux environment with Gurobi and Intel MKL as described in the paper.
2. Create a folder LDB for runing LDB.cpp.
3. Copy the input folder and LDB.cpp to the LDB folder. Then create an LDB/output folder for outputs.
4. In terminal, go to the LDB directory. Copy and paste the command in the MakeCommand_ldb to the terminal, then tap enter to execute the command to compile. LDB.cpp
5. Create a folder MIP for runing MIP.cpp.
6. Copy the input folder, MIP.cpp, Makefile_mip to the MIP folder. Then create an MIP/output folder for outputs, and rename Makefile_mip to Makefile.
7. In terminal, go to the MIP directory. Type "make" and hit enter to compile MIP.cpp

# Execution
1. The code create an execution program named "prog". Use "qsub" command or double click to run it depending on the Linux environment 
2. The code runs algorithms on the instances listed in input/dimacs.txt or input/sparse.txt. Modify the two files to include the instances you'd like to run. 
3. input/bounds.csv gives the predefined bounds for F3, F4 MIP formulations. Modify as needed.
