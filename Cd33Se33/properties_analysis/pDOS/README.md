# Plot average partial density of states (pDOS) for the molecular dynamics trajectory

Here we plot the pDOS over the whole trajectory. First we have defined a function called `convolve_pdos`. This function will use the `*.pdos` and `*.log` files produced by CP2K and weights them with [Gaussian fucntions](https://en.wikipedia.org/wiki/Gaussian_function).

This files combines this function with multiprocessing library of Python and weights the pDOS by Gaussian functions and plots them based on _orbital resolved_ or _atom resolved_ type. The main part starts by defining the `angular_momentum_cols` variable. This variable contains sets of columns in the `pdos` files. For example `list(range(3,12))` adopts all the columns in the `pdos` files and `list(range(4,7))` adopts all the columns related to the *p* orbitals i.e. *px, py, pz*. Here is an example on how to choose the columns in `pdos` files:

```
#===================== Orbital resolved columns
#                           Cd, total          Cd, s             Cd, p             Cd, d
angular_momentum_cols = [ [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ], \
#                           Se, total          Se, s             Se, p             Se, d
                          [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ] ]
```

If you want to weight the pDOS based on atom resolved you need to define the `angular_momentum_cols` as follows:

```
#===================== Atom resolved columns
# #                           Cd, total          Cd, s, p, d
# angular_momentum_cols = [ [ list(range(3,12)), list(range(3,12))],\
# #                           Se, total          Se, s, p, d
#                           [ list(range(3,12)), list(range(3,12))] ]
```

The `Cd, total` and `Se, total` are for the *Total* desity of states. The next ones are for each type.

We define the order of `labels` as the order defined in the `angular_momentum_cols`. The same is with the `colors` when we want to plot. For example:

`labels = ['Cd, s','Cd, p','Cd, d','Se, s','Se, p','Se, d']`

`colors = ['teal','red','orange','purple', 'blue', 'gray', 'green', 'olive', 'brown', 'pink']`

`outname`: This variable is the variable used to save the images. Since we are plotting the orbital resolved pDOS we use `outname = orbital`.


Now we want to start convolving them. So, we need to define the function variables as follows:

`time_step`: The time step is for MD in which the CP2K will append all the partial densities in one file for each type. Therefore, for static calculations we only use _0_.

`sigma`: The standard deviaton value for the Gaussian function. Here we use a vaule of _0.1_.

`coef`: The Gaussian coefficient. This value will weight the functions and you can choose any float number. Here we use _1_.

`npoints`: The number of points for the grid mesh vector. We use _2000_ here. This value should be larger than the number of orbitals in the pdos files.

`energy_conversion`: The energy conversion unit from Hartree to any other unit. Here we use *eV*. So the conversion unit is set to _27.211386_.

`nprocs`: The number of processors to read and convolve the `pdos` files. 


For the `pdos` files the k1 indicates Cd and k2 indicates Se. So, if we use `range(1,3)` or `[1,2]` it will start with the order of Cd and then Se. If we choose `[2,1]` it will first start with Se and then Cd angular momentum columns.

You need to specify the path to find all the `pdos` files using the `glob`. First you need to make a directory and copy all the `pdos` files from the `wd` folder. You can use the following command:

```
mkdir all_pdosfiles
for file in $(find path/to/wd -name '*.pdos'); do cp $file all_pdosfiles/.; echo $file; done
```

In the code you have to specify the path to the `all_pdosfiles`. After that, the code will find and append the names of the files in a variable called `DOS_files1` using the glob for each atom type (`k1` or `k2`). For example `DOS_files1 = glob.glob('Cd33Se33_all_pdosfiles/*k%d-1.pdos'%k)` as in the code.

The rest of the code will just plot the patial density of states first by _Total density of states_ and then the same as in `labels`.



**_NOTE_:** We highly recommend the user to use the latest inputs used in a recent project about the role of crystal symmetry in nonadiabatic dynamics of CsPbI3 perovskite. The new inputs have better functionality and can be found in [this link](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP).






