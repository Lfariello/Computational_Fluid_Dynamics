Driven Cavity Flow in Slender Cavities with Unsteady Lid Velocity  
Code for simulating the 2D incompressible Navier-Stokes equations  
using the vorticity-streamfunction formulation ('psi-zita' model).  
The problem addressed in this code is the flow inside a customizable  
Rectangular (Variable Aspect Ratio) cavity with the top side moving parallel to itself  
with a constant, linear and sinusuidal wall velocity (u_wall).  

The model used is the 'psi-zita' model.  
The boundary conditions are Dirichlet for the elliptic equation for  
Psi (ψ), and 'Thom' conditions for Zita (𝜁), which estimate the vorticity  
production at the wall.  

The variable placement on the mesh is 'collocated,' meaning all variables  
are located at the grid nodes, which also lie directly on the boundaries.
