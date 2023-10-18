# MeshBuilder

Currently MeshBuilder is only supported by the GNU compiler collection. This version of meshBuilder can be compiled and installed using cmake as follows:

```bash
cmake -B path/to/build/directory -S path/to/source directory
cmake --build path/to/build/directory --target all
```

The meshBuilder runs off the .in file which used by tRIBS and in order to build a partitioned mesh the value 9 must be supplied for mesh buildin option.
This is demonstrated as follows:

OPTMESHINPUT:   Mesh input data option
9

Then run meshBuilder demonstrated by the below example:

```meshBuilder example.in```

You can always go back to running Option 8 by changing the OPTMESHINPUT to
verify that meshBuilder is getting the same answer.  Option 8 and Option 9,
serial and parallel should produce the same qout information.

You will know that you correctly got Option 9 because it will inform you that
it is loading files.

Read reach nodes: Pass 1
Read reach edges: Pass 1
Read reach nodes: Pass 2
Read reach edges: Pass 2

The meshBuilder operates only from a .points file at this time.  If the
points file changes the meshBuilder must be run again, but otherwise the
meshBuilder files are used to do the setup for a tRIBS run.

## Output
meshBuilder produces the following files:

### reach.meshb     
This is an ASCII directory file which gives each reach and the numbers of regular nodes and edges in that reach and the
the number of flux nodes and edges in the reach. The flux nodes are those which must communicate with the reach nodes but which might live on another processor.
The counts of nodes and edges are used to find the offset within the file for a particular reach.  The first thing that happens when option 9 is run, is that the partitionin of reaches among the processors is decided.  Then each processor finds the offset in the file for its information and reads in only things needed by its own reaches.

### nodes.meshb & edges.meshb
These are binary files with the data needed to build the mesh and flow net for any reach.  This will not contain all of the tCNode information but only the part of it needed for the flow net.  When option 9 has initialized the mesh is built, the flownet is built, and those variables needed by both are filled in.  The rest of the initialization is then done by the current simulation.  So tFlowNet will be initialized but tKinemat will not for instance.

### fluxnodes.meshb & fluxedges.meshb
These are binary files with the data needed for the fluxnodes and edges.  Because of runon, two nodes beyond thenode owned by a reach must be kept, and the edges which allow circling those nodes must be supplied.  These nodes are added as "inActive" on a processor.  Since it is possiblefor different reaches to have the same flux nodes, meshBuilder makes sure they are unique on a processor.

### flow.meshb
This binary file contains the member variables needed by tFlowNet which was already run on meshBuilder, but whichtRIBS will use.


