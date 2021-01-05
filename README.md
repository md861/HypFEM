# HypFEM
A Fortran code to solve an Initial-Boundary-value scalar wave problem in 2D. The solver uses semidiscrete formulation with lagrange polynomials based Finite Element Method (FEM) for space and diagonally implicit Runge-Kutta family of methods for time discretizations, respectively. The solver also provides an option to enrich the solution space for better convergence rates [[1]](#1). The code is written for the Linux API, however, the system calls could be modified to be run over other operating systems as well. 

<img src="https://github.com/md861/HypFEM/blob/main/images/mesh_p2.png" width="400" height="450"> <img src="https://github.com/md861/HypFEM/blob/main/images/wave.gif" width="400" height="450">
## Features
* Direct import 2D meshes generated in [Gmsh](https://gmsh.info/) - an open source mesh generator.
* Dynamic allocation of all variables and arrays, depending on the imported mesh and order of FEM polynomials.
* Plot the mesh and numerical solutions in [Paraview](https://www.paraview.org/) - an open source alternative to Tecplot.
* Apply Neumann, Dirichlet and Robin boundary conditions on edges marked respectively in Gmsh. Automated calculation of edge normals for using Neumann boundaries.
* LAPACK libraries for matrix solutions.

A complete documentation for the description and usage of the package is under way. However, for the time-being, this read-me file should provide a succinct manual.

## Version system
The version number is based on the following format: *y*:*[m]m*:*n[s]*, for the *n*-th iteration of the code, in the *m*-th month of the *y*-th year in the current decade. The *s* denotes a stable version. Download the latest version available from the repository. The old(er) versions may not be archived for long periods in the future, depending on the number of iterations and the size of updates.

## Setting up the API
### Compiler 
Make sure that the compiler supports Fortran. For Linux based systems, generally this is provided by the GNU compiler. The latest version of HypFEM was compiled using gcc 9.3.0
### LAPACK
Download the latest version from [netlib.org](http://www.netlib.org/lapack/#_release_history). The latest version of HypFEM was compiled using LAPACK, version 3.9.0. Make sure the path of the library is resolved properly. In Linux, this can be easily checked with via running `ld -lblas -llapack --verbose` in the terminal.
### [Gmsh](https://gmsh.info/)
This open-source meshing tool is required to produce the requisite format of data file, that is fed to the HypFEM for computations. Read the instructions below on how to prepare the mesh and define the computational domain.
### [Paraview](https://www.paraview.org/)
This tool is not required to run the HypFEM code, although, it does provide a way to view the output numerical results (if plotted) using this open-source software. 

## Directory tree and output files
The output files are stored in a folder named *case_default* created at runtime. The following files are created:
* *logfile.txt* - This file logs the information about the problem solved such as total degrees of freedom, mesh coordinates, node and edge mappings, integration points used, etc.
* *Plots* - Paraview plot files that contain the numerical (and analytical if available) values over the mesh for each time step.
* *error_data* - If analytical solutions are given, then the normed errors in numerical solution are stored in these files.

Towards the end of the run, the HypFEM code collates all the output files generated and stores them inside the *case_default* directory. By default, at compilation, the code cleans any previous existing folder named *case_default*, including its contents, thus caution must be taken (by renaiming the *case_default* folder, or backing it up, from any previous runs)

## Mesh preparation
The HypFEM code searches for a file named *dat*, inside the current working directory, during runtime to get the mesh information. If this file is not present in the folder, the code would compile correctly, but would generate an error during runtime (NB: the type of error generated is API dependent, and is generally a segmentation fault. However future development of the code would consist a proper handler for this, and other potential, errors).

To prepare the mesh, create a 2D geometry in [Gmsh](https://gmsh.info/) and export the mesh (only meshed using quadrilateral elements) as the *SU2* format and name it as *dat*. A sample *dat* file is given in the [Example](https://github.com/md861/HypFEM/tree/main/Example) folder. Then make sure that this *dat* file (or a copy of it) resides in the same folder as the compiled HypFEM code. 

Make sure to mark the boundary curves/nodes for each type of boundary condition. The boundaries can be described in any sequence, however do make note of the sequence of their definition as this would provide the index (for *pellib.f90* file) to implement the corresponding conditions for each boundadry type. 

NB: At the moment, only non-homogenous Neumann and homogenous Dirichlet conditions may be applied. Although, the future versions would contain the provision for non-homogenous mixed boundaries (even with enriched solution basis).

## Implementing problem conditions
### Initial condition
Modify the *pellib_DIC.f90* file to specify the initial conditions for the wave amplitude (`FZ_Phi`) and its time derivative (`FZ_Vlcty`).
### Boundary condition(s)
* Non-homogeneous Neumann boundary: In *pellib.f90* under the "!Integrate over edges" section set the `NBC_pos` index as the same as the corresponding index of the boundary type defined during meshing (see [Mesh preparation](#mesh-preparation) for details).  


## Usage:
All the files should be in the same folder. Open a terminal in the code folder, and type the following for:
* Compilation: `./CleanNCompile`
* Run: `./femSolver`

 The terminal then outputs the time step currently being processed, with the error in numerical solution if the analytical solution is available. 
 
 ## Description of files:
 * *dat* - This is the "su2" format mesh file supplied by the user, generated from Gmsh. The file should be named as "dat" to be read by the solver.
 * *femSolver.f90* - The main code that coordinates the subroutines and functions. This is where you could change some numerical parameters e.g. 
    * the wavenumber and angular frequency of the problem, 
    * number of integration points to be used, 
    * step size in time for finite differences, 
    * total number of timesteps, 
    * number of plots to be stored, etc.
 * *pellib.f90* - This file allows to specify the boundary sources as well as any sources inside the domain.
 * *pellib_DIC.f90* - Modify this file to specify initial conditions.
 * *ln_norm.f90* and *proslib.f90* - These two files are used to specify the analytical solution (if available) for the computation of normed errors and plotting of analytical values over mesh, respectively.
 
 ## Example files:
 An example *dat* file that has a 2D mesh with 4-th order elements is located in the "Example" folder. The Paraview plots of the mesh and an example numerical solution for a progressive plane wave with Neuman boundaries solved over this mesh, are also available. 
 
## References
<a id="1">[1]</a> 
M. Drolia, *et al*. Enriched finite elements for initial-value problem of transverse electromagnetic waves in time domain. *Computers & Structures*, 182, 354-367, 2017.
