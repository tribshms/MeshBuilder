# MeshBuilder

MeshBuilder is a utility program required to generate teh `*.reach` graphfile for parallel [tRIBS](https://tribshms.readthedocs.io/en/latest/index.html) simulations.

## Obtaining MeshBuilder
Two options are available for using MeshBuilder: (1) use the [docker image](https://tribshms.readthedocs.io/en/latest/man/Docker.html) maintained by the tRIBS developers, or (2) Build your own from the source code in this repository. __We highly recommend the former for ease of use.__ Instructions for using the docker image are found on the tRIBS documentation [page](https://tribshms.readthedocs.io/en/latest/man/Docker.html).

If you decide to build MeshBuilder, note that is only supported by the GNU compiler collection. This version of MeshBuilder can be compiled and built using CMake as follows:

```bash
cmake -B path/to/build/directory -S path/to/source directory
cmake --build path/to/build/directory --target all
```
Additional information on working with CMake is available on the tRIBS documentation [page](https://tribshms.readthedocs.io/en/latest/man/Model_Execution.html#cmake).

## Working with MeshBuilder
MeshBuilder requires a .in file, similar to tRIBS, with at least the keywords:__POINTFILENAME__ and __OUTFILENAME__ provided. __OPTMESHINPUT__ should also be set to option 9. MeshBuilder can be run via the command line as follows.:

```meshBuilder example.in```

MeshBuilder will produce a number of output files in identifying how the model nodes are connected through runoff and sub-surface fluxes. Of these the most important one is the _connectivity.meshb_ file which provides a list of reaches to be partitioned. Once _connectivity.meshb_ is generated it must be converted into a readable format by [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview). METIS then takes the different reaches and partitions them based on three options:

1) SurfaceFlow (SF)
2) Surface-Subsurface Flow (SSF)
3) Surface-Subsurface Flow with Headwaters (SSF-H).

Once METIS is complete, the outputs are then converted back into a format readable by tRIBS denoted as _*.reach_. This is the graphfile required for the keyword __GRAPHFILE__ if __PARALLEMODE__ is turned on. To aid in these conversions we provide 2 pearl scripts _connectivity2metis.pl_ and  _metis_2tribs.pl_ utilized by the shell script. _run_metis.zsh_. These files are stored under _src/workflow_. Before you use this script you must build your own version of METIS as outlined [here](https://github.com/KarypisLab/METIS) and also ensure the _connectivity.meshb_ was generated from MeshBuilder. To run this basic utility script, your must provide the additional argument specifying the number of cores for partitioning, how the mesh should be partitioned based on the above options, and the base name of your simulation. For example to partition a mesh for the [Big Spring Benchmark](https://zenodo.org/records/10951574) onto three cores by Surface Flow (option 1) only, you would run:

```zsh
# ./run_metis.zsh NumberOfNodes OPT_Part basename
./run_metis.zsh 3 1 bigsp
```

Finally, we recommend that:

1) All executables, including MeshBuilder, METIS, and the perl/shell script are located in the same folder as your .in file.
2) That you review the output files from METIS denoted .out. These provide valuable information about whether the partitioning was successful on the requested number of cores. In some instances you can request that the mesh be partitioned finer then is possible. In this case you will over-allocate the number of cores required for the actual task of running your tRIBS simulation.
3) Finally, note that If the points file changes MeshBuilder must be run again.

To wrap up this section, some of the above steps can be bypassed by using the [docker image](https://tribshms.readthedocs.io/en/latest/man/Docker.html) of MeshBuilder. The image includes the METIS executable as well as the required pearl and shell scrips. All that is required in this case is to provide the input file as well as the points file to the appropriate volume (directory) that will be mounted to the image and then you can run _meshbuild_workflow.sh_ which will take you through the above steps. Further documentation is provided [here](https://tribshms.readthedocs.io/en/latest/man/Docker.html#meshbuilder).
