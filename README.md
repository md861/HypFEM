# HypFEM
A Fortran code to solve an Initial-Boundary-value scalar wave problem in 2D. The solver uses semidiscrete formulation with lagrange polynomials based Finite Element Method (FEM) for space and diagonally implicit Runge-Kutta family of methods for time discretizations, respectively. The solver also provides an option to enrich the solution space for better convergence rates [[1]](#1). The code is written for the Linux API, however, the system calls could be modified to be run over other operating systems as well. 

<img src="https://github.com/md861/HypFEM/blob/main/images/mesh_p2.png" width="400" height="450"> <img src="https://github.com/md861/HypFEM/blob/main/images/wave.gif" width="400" height="450">
## Features
* Directly import 2D meshes generated in [Gmsh](https://gmsh.info/) - an open source mesh generator.
* Seamless integration of user defined enrichment functions into the Finite Element solution space.
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
Download the latest version from [netlib.org](http://www.netlib.org/lapack/#_release_history). The latest version of HypFEM was compiled using LAPACK, version 3.9.0. Make sure the path of the library is resolved properly. In Linux, this can be easily checked via running `ld -lblas -llapack --verbose` in the terminal.
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

NB: At the moment, only Neumann conditions may be applied. Although, the future versions would contain the provision for non-homogenous mixed boundaries (even with enriched solution basis).

## Implementing problem sources
### Initial condition
Modify the *pellib_DIC.f90* file to specify the initial conditions for the wave amplitude (`FZ_Phi`) and its time derivative (`FZ_Vlcty`).
### Boundary condition(s)
* Non-homogeneous Neumann boundary: 
    * In *pellib.f90* under the "!Integrate over edges" section set the `NBC_pos` index as the same as the corresponding index of the boundary type defined during meshing (see [Mesh preparation](#mesh-preparation) for details). In this section, update the `GZ` variable that is used to implement the relevant condition for the given boundary type (chosen using `NBC_pos` index). Notice that `n(1)` and `n(2)` store the x and y components of the normal vector to a given edge. This normal vector is computed for each edge automatically during run-time. 
    * For each different boundary condition (of non-homogeneous Neumann type), copy and paste the entire "!Integrate over edges" section (i.e. all the 4 subsections integrating over the four sides of a quadrilateral) and set the `NBC_pos` index appropriately. 
* Homogeneous Neumann boundary:
    * Since these boundaries do not require integration, simply do not include any "!Integrate over edges" section for this `NBC_pos` index in the *pellib.f90* file.
### Source term(s)
In the *pellib.f90* file, under the "!Integrate inside domain" section, update the `FZ` variable that evaluates the source terms for the given problem. If there are no source terms present, then this section can simply be commented out. 

## Implementing numerical parameters
Numerical parameters, e.g. the quadrature points for integration, total number of time iterations, etc could be set in the *femSolver.f90* file. Some (most used) parameters are detailed below:
* `NGAUSS` - this variable defines the number of quadrature points to be used for line integrals, and the square of this for integrations inside the element.
* `NPOINTPLOT` - similar to `NGAUSS`, defines the number of quadrature points for plotting the numerical results. 
* `n_record_IStep` - sets the intervals at which the numerical results are plotted. To disable numerical plots, simply comment the section in the time loop.
* `dt` - step size in time.
* `n_istep` - number of time iterations to be run. Thus, total time = `dt*n_istep`.
* `RKM_Order` - determines the order of Runge-Kutta method for time integration. Make sure to initialize the corresponding Butcher tableau (`a_RKM`,`b_RKM`,`c_RKM`) in the "!Initialize" section (in the same file *femSolver.f90*).
* `Cp` - phase velocity of the wave.
* `NANGL` - number of enrichment functions per node. Set `NANGL = 1` for purely p-FEM solutions.
* `K_W` - wavenumber of the enrichment functions. NB: these are not used if `NANGL = 1`.

 ## Description of files:
 A short summary of some files (which are usually modified the most) is presented here. 
 * *dat* - This is the "su2" format mesh file supplied by the user, generated from Gmsh. The file should be named as "dat" to be read by the solver.
 * *femSolver.f90* - The main code that coordinates the subroutines and functions. This is where you could change some numerical parameters e.g. 
    * the wavenumber and angular frequency of the problem, 
    * number of integration points to be used, 
    * step size in time for finite differences, 
    * total number of timesteps, 
    * number of plots to be stored, etc.
 * *pellib.f90* - This file allows to specify the boundary sources as well as any sources (`FZ`) inside the domain.
 * *pellib_DIC.f90* - Modify this file to specify initial conditions.
 * *ln_norm.f90* and *proslib.f90* - These two files are used to specify the analytical solution (if available) for the computation of normed errors and plotting of analytical values over mesh, respectively.
 * *psflib.f90* - Modify the finite element shape functions used for the numerical solution. By default, Lagrange polynomials are used for `NANGL = 1` (pertaining to the order of elements defined by the mesh imported from *dat* file). For `NANGL > 1`, by default the shape functions are defined using plane wave enrichment functions [[1]](#1).
 
## Usage:
All the files should be in the same folder. Open a terminal in the code folder, and type the following for:
* Compilation: `./CleanNCompile`
* Run: `./femSolver`

The terminal then outputs the time step currently being processed, with the error in numerical solution if the analytical solution is available. 
 
## Example files:
An example *dat* file that has a 2D mesh with 2-nd order quadrilateral elements (generated with [Gmsh](https://gmsh.info/)) is located in the [Example/10pi_p2](https://github.com/md861/HypFEM/tree/main/Example/10pi_p2) folder. The Paraview plots of the mesh and an example numerical solution for a progressive plane wave with homogeneous Neuman boundaries and a non-zero (sinusoidal) source solved over this mesh, are shown as animated gifs at the beginning of this read-me file. 

## Citing this package
[![DOI](https://zenodo.org/badge/327013671.svg)](https://zenodo.org/badge/latestdoi/327013671)

## References
<a id="1">[1]</a> 
M. Drolia, *et al*. Enriched finite elements for initial-value problem of transverse electromagnetic waves in time domain. *Computers & Structures*, 182, 354-367, 2017.
