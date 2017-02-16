# Neonate Head Tank


 <img src="https://raw.githubusercontent.com/EIT-team/Tanks/master/doc/figures/nn_skull_final.png" height="400" alt="Neonatal tank"> <img src="https://raw.githubusercontent.com/EIT-team/Tanks/master/doc/figures/nn_tank_final.png" height="400" alt="Neonatal tank">

## I just want to print this!

The .stl files are in the `Printing` folder [here](./Construction/Printing)

## Whats in these folders?

- [Construction](./Construction/) has all the design files necessary to make the tank, change the electrode positions/diamaters, as well as changing the conductivity of the skull. This also has the printable `.stl` files
- [Forward_Solver](./Forward_Solver/) has the mesh files and code to generate the things necessary to use this tank in an EIT forward solver, most likely [PEITS](https://github.com/EIT-team/PEITS)
- [Meshing](./Meshing/) has the meshes and Matlab code to convert the solid models to FEMs. Also has the code to generate the Macro to create the holes in the skull.
