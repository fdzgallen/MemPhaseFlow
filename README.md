# MemPhaseFlow: Membrane Phase-field Flow in 2D

 > Please cite as: Gallen, A. F., Castro, M., & Hernandez-Machado, A. (2021). Red Blood Cells in low Reynolds number flow: a vorticity-basedcharacterization of shapes in two dimensions. Soft Matter, pre-print. 
 
 > https://doi.org/

This repository is a Work in Progress and will be updated over time.

Code to simulate the temporal evolution of a biological membrane inside a fluid flow.

This is the code used in the peer-review article "Red Blood Cells in low Reynolds number flow: a vorticity-basedcharacterization of shapes in two dimensions" accepted for publication in the Soft Matter Journal. It is uploaded here for transparency and to make the use or implementation of our model easier for any person interested.

Straightfoward and simple application using a stream function formulation to solve the fluid flow. This way the incompressibility of the fluid is mantained for

Modifying the code as explained in the article you can simulate multiple flows: Poiseuille, Couette, temporal-dependent etc. Here we will provide a few versions of the code so you can simulate different systems without issue, but playing with the code to explore new possibilities is encouraged.

# Input

The program consists of the python program that implements the model and an input.txt file where you introduce the parametes that will define the system.
There is no need to understand all of them, we will start by explaining the most important ones.

# Output

The simulation will create a folder with the name that it is given in the input file. Inside this folder the code stores images of the state of the simulation every a given interval of iterations.
It also creates three folders, "phi", "stream", and "vorticity". In each of this folders the state of each of this fields is stored over time. 
In case there is interest on analysing more in depth the evolution of the system you can write a code that reads these output files and computes the desired values from the data without need of running the simulation all over again. This uses storage space but saves a lot of time if you want to compute a new value (for example the evolution of the center of masss of the cell (phi) over time) but you do not want to add the calculation to the code and run the simulation again.


# How does it work?

The details of the mathematical model can be read in the article. However here we will give a brief explanation of the numerical implementation in the near future.
