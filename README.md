# MemPhaseFlow: Membrane Phase-field Flow in 2D

 > Please cite as: A. F. Gallen, M. Castro and A. Hernandez-Machado, Soft Matter, 2021, DOI: 10.1039/D1SM00559F.

 > https://doi.org/10.1039/D1SM00559F

This repository is a Work in Progress and will be updated over time.

Code to simulate the temporal evolution of a biological membrane inside a fluid flow.

This is the code used in the peer-review article "Red Blood Cells in low Reynolds number flow: a vorticity-basedcharacterization of shapes in two dimensions" accepted for publication in the Soft Matter Journal. It is uploaded here for transparency and to make the use or implementation of our model easier for any person interested.

Straightfoward and simple application using a stream function formulation to solve the fluid flow. This way the incompressibility of the fluid is mantained for

Modifying the code as explained in the article you can simulate multiple flows: Poiseuille, Couette, temporal-dependent etc. Here we will provide a few versions of the code so you can simulate different systems without issue, but playing with the code to explore new possibilities is encouraged.

# Input

The program consists of the python program that implements the model and an input.txt file where you introduce the parametes that will define the system.
There is no need to understand all of them, we will start by explaining the most important ones.


* "simulation" will be the name of the folder where the data will be saved 

Depending on wether the simulation is a Poiseuille or Couette flow: 
* "Pspeed"/"wallspd" are Poiseuille average flow speed/the speed of the top wall for a Couette flow respectively in dx/dt units

The base codes start with an ellipse-shaped cell, as a Red Blood Cell shape can be obtained from an ellipse with a given relation between its enclosed volume vs surface area:
* "R1" semi-major axis of the ellipse for the initial conditions for the vesicle
* "R2" semi-minor axis of the ellipse for the initial conditions for the vesicle. For Red blood cells taken as 0.23*R1

Duration of the simulation
* "Tend" total number of steps for the simulation to end
* "t0" number of time steps to evolve the interface before starting running the fluid flow. If we start the simulation with flow right away we do not give time to the membrane to form and to measure the starting area and volume values.
* "tdump" data dump period

Size of the simulated system, it is a box of Nx*Ny points
* "Nx" lattice size X
* "Ny" lattice size Y

Variables that define the membrane and the phase field
* "Nx0" starting position of the cell in X axis
* "Ny0" starting position of the cell in Y axis
* "eps" interfacial width, to increase eps its necessary to decrease "dt" so that the simulation converges, reducing eps gives lower resolution for stream function and voritcity fields
* "C0" Spontaneous curvature of the membrane. Can take values in a (-0.4 , 0.4) range easily. Higher values might give problems. You can obtain the radius of curvature of the spontaneous curvature in dx units by simply R0=1/C0. 
* "M" mobility of the phase field. This value is important as the relaxation time of the phase field has to be much shorter than that of the fluid
* "kappa" bending modulus of the membrane
* "viscosityliq" viscosity value of the liquid of the channel (phi=-1)
* "viscositycell" viscosity value of the liquid inside the cell (phi=+1)
 
Resolution of the simulation
* "dt" temporal resolution. The bigger this is the faster the simulation goes but there is a limit at which the computation diverges and you get NaN results. Increasing some parameters require to reduce this value.
* "dx" spatial resolution in the x axis. This is the distance between 2 consecutive points in the x axis.
* "dy" spatial resolution in the y axis. This is the distance between 2 consecutive points in the y axis.
* "nit" number of iterations for the iterative Poisson solver. As the state of the simulation is always close to the next solution we dont need many iterations

Lagrange Multipliers
* "alpha" Strength of the lagrange mutiplier effect for AREA
* "beta" Strength of the lagrange mutiplier effect for VOLUME


# Output

The simulation will create a folder with the name that it is given in the input file. Inside this folder the code stores images of the state of the simulation every a given interval of iterations.
It also creates three folders, "phi", "stream", and "vorticity". In each of this folders the state of each of this fields is stored over time. 
In case there is interest on analysing more in depth the evolution of the system you can write a code that reads these output files and computes the desired values from the data without need of running the simulation all over again. This uses storage space but saves a lot of time if you want to compute a new value (for example the evolution of the center of masss of the cell (phi) over time) but you do not want to add the calculation to the code and run the simulation again.


# How does it work?

The details of the mathematical model can be read in the article. However here we will give a brief explanation of the numerical implementation in the near future.
