/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  bMeshBuilder.h: Header for bMeshBuilder.cpp 
**  	
**  bMeshBuilder Class builds the necessary mesh, flow network and
**  connectivity for very large datasets.  It runs one one processor
**  taking a *.points file as input.  It constructs the nodes, edges and
**  triangles by writing temporary files along the way to free memory.
**  It calls bMesh to build the mesh, bFlowNet to build the flow net
**  structure and bGraph to build the connectivity needed for partitioning.
**  bMeshBuilder writes information files about the structure which can
**  be read by any number of simulation processors before running the
**  simulation.
**
***************************************************************************/

#ifndef BMESHBUILDER_H
#define BMESHBUILDER_H

//=========================================================================
//
//
//                  Section 1: bMeshBuilder Include and Define Statements
//
//
//=========================================================================

#include "src/Headers/Inclusions.h"
#include "src/bInOut/bInputFile.h"
#include "src/bMesh/bMesh.h"
#include "src/bFlowNet/bFlowNet.h"
#include "src/bGraph/bGraph.h"

//=========================================================================
//
//
//                  Section 2: Mesh Builder Class Definition
//
//
//=========================================================================

class bMeshBuilder 
{
 public:
  bMeshBuilder(bInputFile& infile);
  ~bMeshBuilder();

  bMesh<bNode>        *baseMesh;    // Pointer to mesh for triangularization
  bMesh<bFlowNode>    *flowMesh;    // Pointer to mesh with flow information
  bFlowNet            *flowNet;     // Pointer to restart object
  bGraph              *graph;       // Pointer to partition structure
};

#endif 

//=========================================================================
//
//
//                             End of bMeshBuilder.h 
//
//
//=========================================================================
