/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**           and Los Alamos National Laboratory
**  
**
**  bGraph.h: Header for tConnectivity class and objects
**
**  bGraph Class used in tRIBS for the parallel version, providing 
**  information about upstream and downstream basins
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: bGraph Include and Define Statements
//
//
//=========================================================================

#ifndef BGRAPH_H
#define BGRAPH_H

#include <iostream>
#include <vector>
#include <set>

#include "src/bGraph/bGraphNode.h"
#include "src/bMesh/bMesh.h"
#include "src/bFlowNet/bFlowNet.h"

//=========================================================================
//
//
//                  Section 2: bGraph Class Definitions
//
//
//=========================================================================

class bGraph
{
public:
  bGraph();
  bGraph(bMesh<bFlowNode>*, bFlowNet*, bInputFile&);
  ~bGraph();

  void SetConnectivity();
  void OutputConnectivity();

  // Write the node and edge information by reach
  void SortNodesByReach();
  void WriteFlowMesh(); 
  void WriteFlowNode(fstream&, bFlowNode*) const;
  void WriteFlowEdge(fstream&, bEdge*) const;

private:
  bMesh<bFlowNode>*  mesh;            //!< Mesh
  bFlowNet*          flow;            //!< Kinematic flow

  int                numberOfReaches; //!< # of stream reaches
  int*               headID;          //!< Reach head node IDs
  int*               outletID;        //!< Reach outlet node IDs
  int*               nodeAboveID;     //!< Node above the outlet

  std::list<bFlowNode*>** ReachNodes; //!< Node ids in reach
  bGraphNode* connection;             //!< Reach connectivity
};

#endif

//=========================================================================
//
//
//                          End of bGraph.h 
//
//
//=========================================================================
