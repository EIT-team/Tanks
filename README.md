# Head tanks for Electrical Impedance Tomography
---
These head tanks are designed for Electrical Impedance Tomography (EIT) experiments. They are formed of two parts: a head shaped tank, and a skull. This repository contains both an Adult and a neonatal head tank. The geometry is derived from CT/MRI segmentations. The conductivity of the skull is controlled through surface perforations which allow saline to pass through.

<img src="https://raw.githubusercontent.com/EIT-team/Tanks/master/doc/figures/final_ad.png" height="300" alt="Adult tank with KHU connected"> <img src="https://raw.githubusercontent.com/EIT-team/Tanks/master/doc/figures/final_nn.png" height="400" alt="Neonatal tank">

### DESIGN FILES REQUIRE GIT LFS - CONTAINS LOADS OF SOLID MODELS ~2GB DOWNLOAD
Git lfs is *not* required to get printable `.stl` files, only if you want the design files or to edit the electrode positions or the skull conductivities.

### Making tanks
If you are interested in making one of these tanks, please create an __issue__ in the repo or email us, we would love to help you!

To make one of these tanks you will need to:

- Print the tank
- Print the skull
- Remove support material from both models
- Get electrodes from somewhere - either [stainless steel ones](./doc/electrodes) which we had made (We have spares too!), or Ag/AgCl electrodes for the neonatal tank from [biomed electrodes](http://www.biomedelectrodes.com/product/bmd-8/)
- Insert electrodes in tank (somewhat forcefully) and seal with silicon
- Wait to dry
- Insert skull into the tank
- Fill with 0.2% saline
- Collect EIT data
- Win

The `.stl` files are found [here for the neonate tank ](./Neonate/Construction/Printing) and [here for the adult tank ](./Adult/Construction/Printing)

### Using these tanks
If you are using the `.stl` files provided, without changing the conductivities, geometry or electrode positions, then you can use the Meshes and code to generate the `.dgf` files for use with [here](https://github.com/EIT-team/PEITS). The files for the Adult tank are [here](./Adult/Forward_Solver), and [here](./Neonate/Forward_Solver) for the neonatal tank.

### Citing this work
The accompanying journal article for these tanks is found [here](http://iopscience.iop.org/article/10.1088/1361-6579/aa6586)
```
Avery, J. et al., Reproducible 3D printed head tanks for electrical impedance tomography with realistic shape and conductivity distribution. Physiological measurement, 1116.

```

