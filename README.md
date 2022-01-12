# Charge Displacement Function
Small Fortran programs for the generation of a charge displacement function from gaussian cube files

## Workflow to generate the charge displacement function
1. Perform a geometry optimization of your molecule with a quantum chemistry software of your chose. 
2. Align the bond to be analyzed in the z-axis using for example with Avogadro. For terminal end-on dinitrogen complexes you can also use `align-xyz-to-z`.
3. Separate the molecule at the bond to be analyzed into two fragments.
4. Generate gaussian cube files of all occupied molecular orbitals for the complete molecule and the separated fragments in their unrelaxed geometries with the same 3D grid. Therefore, use a high resolution, e.g. 0.1 Å in x and y direction and 0.05 Å in the z direction and make the grid large enough, so not electron density will be cut off. 
5. Optional: Separate the molecular orbitals according to their irreducible representations into different groups. If the separation is not performed, the total electron densities of the molecule and the separated fragments can be directly generated with the quantum chemistry software. In this case, continue with step 8.
6. Transfer all molecular orbitals into electron densities with `mo-to-dens <occupation>`. The occupation is 1 for open-shell calculations and 2 for closed-shell calculations. 
7. Add all electron density cube files (of one group) together using `calc` to obtain the total electron density. 
8. Subtract the total electron densities (of each group) of the separated fragments from the those of the complete molecule using `calc`.
9. Calculate the charge displacement function from the difference electron density cube file (of each group) with `integrate-cube`. 
8. Determine the isodensity boundary on the z-axis between the separated fragments with 
`find-fragment-intersection <total-density-fragment1.cub> <total-density-fragment2.cub>`
9. Determine the value of the charge displacement function at the isodensity boundary on the z-axis with 
`read-function-at <cdf.dat> <z-value>`

## align-xyz-z
Aligns the metal dinitrogen bond of a dinitrogen complex in the z axis. The dinitrogen ligand has to be at position 1 and 2 of the xyz file and the metal atom on position 3.

Usage: `align-xyz-z <input.xyz> <output.xyz>`

## mo-to-dens
Transfers a cube file of a molecular orbital into the corresponding electron density by squaring and normalization. The occupation of the molecular orbital must be given. 

Usage: `mo-to-dens <mo.cub> <occupation> <density-output.cub>`

## calc
Allows simple arithmetic between the data points of two cube files at the same position, if they have the same grid. A calculation of each data point of one cube file with a number is also possible.

Usage: `calc <input1.cub> <+,-,*,/> <input2.cub> <output.cub>`

Examples: `calc density1.cub + density2.cub result.cub`

	`calc density.cub * 2 result.cub`

## integrate-cube
Calculates the charge displacement function from the difference electron density cube file by numerical integration according to: 
![cdf]( https://github.com/Manuel-Schmitt/charge-displacement-function/tree/main/pictures/cdf.png?raw=true)
An output of the following integral is also generated:
![integral]( https://github.com/Manuel-Schmitt/charge-displacement-function/tree/main/pictures/integral.png?raw=true)

Usage: `integrate-cube <input.cub>` 

## find-fragment-intersection
Determines the isodensity boundary on the z-axis between two separated fragments.

Usage: find-fragment-intersection <total-density-fragment1.cub> <total-density-fragment2.cub>

## read-function-at
Determines the value of the charge displacement function, generated with ` integrate-cube` at the given point on the z-axis. 

Usage: `read-function-at <cdf.dat> <z-value>`
