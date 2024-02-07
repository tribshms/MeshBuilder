// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bGraph.cpp: Functions for class bGraph (see bGraph.h)
**
***************************************************************************/

#include <cassert>
#include <map>

#include "src/bGraph/bGraph.h"

bGraph::bGraph()
{
   mesh = 0;
   flow = 0;
   numberOfReaches = 0;
}

/***************************************************************************
** 
** Constructor takes general information about all nodes and edges from
** bMesh and the flow reach information from bFlowNet and constructs
** information necessary for arranging nodes and edges by reach
**
***************************************************************************/

bGraph::bGraph(bMesh<bFlowNode>* m, bFlowNet* f, bInputFile& InputFile)
{
   mesh = m;
   flow = f;

   // Number of reaches are those identified by bFlowNet
   numberOfReaches = flow->getReachHeadList().size();

   // Allocate array of lists to hold nodes in each reach and connections
   connection = new bGraphNode[numberOfReaches];
   ReachNodes = new std::list<bFlowNode*>*[numberOfReaches];
   for (int reach = 0; reach < numberOfReaches; reach++)
      ReachNodes[reach] = new std::list<bFlowNode*>;

   // Create stream connectivity table
   cout << "\nCreating stream reach connectivity table..." << endl;
   SetConnectivity();
}

/***************************************************************************
** 
** Destructor
**
***************************************************************************/

bGraph::~bGraph()
{
   mesh = 0;
   flow = 0;
   numberOfReaches = 0;
   for (int reach = 0; reach < numberOfReaches; reach++) {
      std::list<bFlowNode*>* reachNodes = ReachNodes[reach];
      reachNodes->erase(reachNodes->begin(), reachNodes->end());  
   }
   delete [] ReachNodes;
   delete [] connection;
} 

/*************************************************************************
**
** Create stream reach connectivity table from bFlowNet object.
**
*************************************************************************/

void bGraph::SetConnectivity()
{
  // Get reach heads and outlets
  std::list<bFlowNode*>& hlist = flow->getReachHeadList();
  std::list<bFlowNode*>& olist = flow->getReachOutletList();
  
  // Figure out connectivity
  bFlowNode *chead, *coutlet, *cnext, *clast;
  std::_List_iterator<bFlowNode*> HeadIter;
  std::_List_iterator<bFlowNode*> OutletIter;

  int i;
  headID = new int[numberOfReaches];
  outletID = new int[numberOfReaches];
  nodeAboveID = new int[numberOfReaches];

  // First collect head and outlet IDs
  for (HeadIter = hlist.begin(), OutletIter = olist.begin(), i = 0;
       HeadIter != hlist.end();
       HeadIter++, OutletIter++, i++) {
    chead = (*HeadIter);
    coutlet = (*OutletIter);
    cout << "bGraph:connectivity Reach: " << i 
         << "\thead = " << chead->getID() 
         << "\toutlet = " << coutlet->getID() << endl;
    bGraphNode rnode(i);
    connection[i] = rnode;
    headID[i] = chead->getID();
    outletID[i] = coutlet->getID();

    // Get the last node before the outlet which is needed by simulator
    cnext = chead;
    clast = cnext;
    while (cnext != coutlet) {
      clast = cnext;
      cnext = cnext->getDownstrmNbr();
    }
    nodeAboveID[i] = clast->getID();
  }                  	                                                              
  // Figure out which are connected
  // For each reach, find all outlets == head, upstream
  //                 find all heads == outlet, downstream
  for (i = 0; i < numberOfReaches; i++) { 
    for (int j = 0; j < numberOfReaches; j++) {
      if (j != i) {
        if (outletID[j] == headID[i])
          connection[i].addUpstream(j);
        if (headID[j] == outletID[i])
          connection[i].addDownstream(j); 
      }
    }
  }

  // Set the ground flux nodes needed for connectivity.meshb file
  std::map<int, std::map<int, int> > reachFlux;
  bMeshList<bEdge*>* elist = mesh->getEdgeList();
  std::list<bEdge*>::iterator eiter;

  for (eiter = elist->begin(); eiter != elist->getLastActive(); eiter++) {

    bFlowNode* co = (bFlowNode*) ((*eiter)->getOriginPtrNC());
    bFlowNode* cn = (bFlowNode*) ((*eiter)->getDestinationPtrNC());

    int coReach = co->getReach();
    int cnReach = cn->getReach();

    if (coReach != cnReach)
      reachFlux[coReach][cnReach] = 1;
  }

  for (int i = 0; i < numberOfReaches; i++) {
    for (int j = 0; j < numberOfReaches; j++) {
      if ((i != j) && (reachFlux[i][j] == 1))
        connection[i].addFlux(j);
    }
  }
  reachFlux.clear();

  // Upstream and downstream connections per reach
  cout << "\nUpstream and Downstream Reach connection:" << endl;
  for (i = 0; i < numberOfReaches; i++) 
      cout << connection[i] << endl;

   // Graphviz representation
   ofstream graphvizStr("connectivity.dot");
   graphvizStr << "digraph abstract {" << endl;
   for (int i = 0; i < numberOfReaches; i++) {
      const std::vector<int>& downstream = connection[i].getDownstream();
      for (int j = 0; j < downstream.size(); j++)
         graphvizStr << "\t" << i << "->" << downstream[j] << ";" << endl;
   }
   graphvizStr << "}" << endl;
   graphvizStr.close();
}

/*************************************************************************
**
** Partition graph based on number of partitions previously set.
**
*************************************************************************/

void bGraph::OutputConnectivity()
{
  // Write out reach connectivity to a file
  ofstream connFile("connectivity.meshb");
  if (!connFile.good()) {
      cout << "connFile problem" << endl;
      exit(5);
  }

  connFile << "#ReachID PointCount HeadNodeID OutletNodeID ";
  connFile << "NumberDownstream [REACHID1 REACHID2 ...] ";
  connFile << "NumberFlux [FLUXID1 FlUXID2 ...]" << endl;

  connFile << numberOfReaches << endl;
  for (int i = 0; i < numberOfReaches; i++) {
    connFile << connection[i].getID() << " "
             << ReachNodes[i]->size() << " "
             << headID[i] << " "
             << outletID[i] << " ";

    std::vector<int> downReach = connection[i].getDownstream();
    connFile << downReach.size();
    for (int j = 0; j < downReach.size(); j++) {
      connFile << " " << downReach[j];
    }

    std::vector<int> flux = connection[i].getFlux();
    connFile << " " << flux.size();
    for (int j = 0; j < flux.size(); j++) {
        connFile << " " << flux[j];
    }
    connFile << "\n";
  }
  connFile.close();
}

/*************************************************************************
**
** Organize the nodes in the nodelist by reach number
** As a pointer to a node is put in the reach lists, remove it from the
** node list to save space
**
*************************************************************************/

void bGraph::SortNodesByReach()
{
  // Get complete node and edge list which is in flow order but not reach order
  bMeshList<bFlowNode*> *nlist = mesh->getNodeList();
  std::_List_iterator<bFlowNode*> niter;
  std::_List_iterator<bFlowNode*> nodeToMove;
  bFlowNode* cn;
  int reach;

  // Iterate on all active nodes moving them to the correct reach list
  for (niter = nlist->begin(); niter != nlist->getLastActive(); niter++) {
    cn = (*niter);
    reach = cn->getReach();
    ReachNodes[reach]->push_back(cn);
  }

  cout<<"\nDerived " << numberOfReaches << " stream reaches..."<<endl;
  std::list<bFlowNode*>::iterator iter;
  for (int i = 0; i < numberOfReaches; i++)
     cout << "Nodes in Reach " << i << ": " << ReachNodes[i]->size() << endl;
}

/*************************************************************************
**
** Write the node and edge information in reach order
** Node information was already calculated so each node knows its reach
** Each reach also has a list of node ids for nodes it contains
**
** For each reach node look at edges and decide if it is completely internal
** to the reach, or external going to a boundary or another reach.  If the
** edge is external then the destination might be a flux node so record
** that for the reach so that when reaches are partitioned between processors
** the flux nodes will be able to exist on both
**
** For each boundary (inactive) nodes record internal or external edges,
** plus the external nodes which are destination nodes in actual reaches,
** and record the edges of those destination reaches because they are needed
** to find the complement edge coming back to the boundary node which is
** needed for MakeBoundaryPolygon in the simulator
**
*************************************************************************/

void bGraph::WriteFlowMesh()
{
  // Collect the counts at the top of this method because after sorting
  // nodes by reach the nodeList no longer exists
  int totalNodeCount = mesh->getNodeList()->getTotalSize();
  int activeNodeCount = mesh->getNodeList()->getActiveSize();
  int totalEdgeCount = mesh->getEdgeList()->size();

  // Reach number is set for each node but we want to organize nodes by
  // reach for writing
  SortNodesByReach();

  // Write connectivity to a file, needed for figuring out partitions
  // Has to be done here because the above counts the nodes per reach
  OutputConnectivity();

  // Files for nodes and edges which contain one copy of each node and edge
  // Can be used in serial version, or parallel on one processor version
  fstream nodeStr("nodes.meshb", ios::out | ios::binary);
  fstream edgeStr("edges.meshb", ios::out | ios::binary);

  // Files with duplicate flux nodes and edges per reach for parallel version
  fstream fluxNodeStr("fluxnodes.meshb", ios::out | ios::binary);
  fstream fluxEdgeStr("fluxedges.meshb", ios::out | ios::binary);

  // Directory with reach, node and edge counts, and offsets within files
  fstream reachStr("reach.meshb", ios::out);

  // Save a place in the first word of edges and nodes for the final sizes
  BinaryWrite(nodeStr, totalNodeCount);
  BinaryWrite(edgeStr, totalEdgeCount);

  // Space used by each node and edge written to file
  int nodeBytes = 192;
  int edgeBytes = 56;

  // Write number of reaches and total number of nodes and edges
  reachStr << numberOfReaches << "\t" 
           << totalNodeCount << "\t" << nodeBytes << "\t"
           << totalEdgeCount << "\t" << edgeBytes << endl;

  bEdge *firstEdge, *curEdge;
  bFlowNode *origNode, *destNode;
  int origReach, destReach;

  // Nodes were previously stored in ReachNode vector
  std::list<bFlowNode*>::iterator niter;

  // Two edges between every pair of nodes
  std::list<bEdge*> internalEdge;
  std::list<bEdge*> externalEdge;
  std::list<bEdge*>::iterator eiter;

  // Flux edges and nodes are duplicates belonging to another reach but
  // necessary to the operation of the current reach when in parallel
  std::set<bEdge*> fluxEdge;
  std::set<bFlowNode*> fluxNode;
  std::set<bEdge*>::iterator feiter;
  std::set<bFlowNode*>::iterator fniter;

  // Keep offsets within files by reach for access in the parallel simulator
  int nodeOffset = 4;
  int edgeOffset = 4;
  int fluxNodeOffset = 0;
  int fluxEdgeOffset = 0;

  ////////////////////////////////////////////////////////////////////////
  //
  // Iterate on each reach writing the local node and edge information
  //
  for (int reach = 0; reach < numberOfReaches; reach++) {

    // Iterate on all nodes within a reach
    for (niter = ReachNodes[reach]->begin(); 
         niter != ReachNodes[reach]->end(); niter++) {
      origNode = (*niter);
      origReach = origNode->getReach();
      
      // Investigate all edges for this origin node
      firstEdge = origNode->getFlowEdg();
      curEdge = firstEdge;
      do {
        // Find destination node and reach for this edge
        destNode = (bFlowNode*) curEdge->getDestinationPtrNC();
        destReach = destNode->getReach();

        if (origReach == destReach) {
          // Internal edge has origin and destination in this reach
          internalEdge.push_back(curEdge);

        } else {
          // External edge has destination on another reach or in boundary area
          // Since destination will appear in file for another reach we save
          // it in the flux node file for this reach as a duplicate
          externalEdge.push_back(curEdge);
          fluxNode.insert(destNode);

          // For groundwater we need the slope of the flow edge of flux node
          // which means all the edges and nodes that come out of it also
          bEdge* destEdge = destNode->getEdg();
          bEdge* curDestEdge = destEdge;
          do {
            bFlowNode* curDestNode = 
              (bFlowNode*) curDestEdge->getDestinationPtrNC();
            fluxEdge.insert(curDestEdge);
            fluxNode.insert(curDestNode);
            curDestEdge = curDestEdge->getCCWEdg();
          } while (curDestEdge != destEdge);
        }
        curEdge = curEdge->getCCWEdg();
      } while (curEdge != firstEdge);

      WriteFlowNode(nodeStr, origNode);
    }

    // Write internal edges for this reach
    for (eiter = internalEdge.begin(); eiter != internalEdge.end(); eiter++)
       WriteFlowEdge(edgeStr, (*eiter));

    // Write external edges for this reach
    for (eiter = externalEdge.begin(); eiter != externalEdge.end(); eiter++)
       WriteFlowEdge(edgeStr, (*eiter));

    // Write flux nodes for this reach
    for (fniter = fluxNode.begin(); fniter != fluxNode.end(); fniter++)
       WriteFlowNode(fluxNodeStr, (*fniter));

    // Write flux edges for this reach
    for (feiter = fluxEdge.begin(); feiter != fluxEdge.end(); feiter++)
       WriteFlowEdge(fluxEdgeStr, (*feiter));

    // Write directory information
    reachStr << reach << "\t"

             << ReachNodes[reach]->size() << "\t"
             << nodeOffset << "\t"

             << fluxNode.size() << "\t"
             << fluxNodeOffset << "\t"

             << internalEdge.size() << "\t"
             << externalEdge.size() << "\t"
             << edgeOffset << "\t"

             << fluxEdge.size() << "\t"
             << fluxEdgeOffset << endl;

    // Accumulate number of nodes and edges and add to offset
    nodeOffset += ReachNodes[reach]->size() * nodeBytes;
    fluxNodeOffset += fluxNode.size() * nodeBytes;
    edgeOffset += (internalEdge.size() + externalEdge.size()) * edgeBytes;
    fluxEdgeOffset += fluxEdge.size() * edgeBytes;

    // Clear the temporary lists
    internalEdge.erase(internalEdge.begin(), internalEdge.end());
    externalEdge.erase(externalEdge.begin(), externalEdge.end());
    fluxNode.erase(fluxNode.begin(), fluxNode.end());
    fluxEdge.erase(fluxEdge.begin(), fluxEdge.end());
  }

  ////////////////////////////////////////////////////////////////////////
  //
  // Write all the boundary nodes which are different from a regular reach
  // in that the flux node must be recorded as well as all the CCW edges
  // for that flux node because returning edge must also be on processor
  // Boundary nodes are in the inactive part of the mesh node list
  //
  bMeshList<bFlowNode*> *nlist = mesh->getNodeList();
  std::_List_iterator<bFlowNode*> biter;
  for (biter = nlist->getLastActive(); biter != nlist->end(); biter++) {
    origNode = (*biter);
    origReach = origNode->getReach();

    // Investigate all edges for this node
    firstEdge = origNode->getEdg();
    curEdge = firstEdge;
    do {
        // Find destination node and reach for this edge
        destNode = (bFlowNode*) curEdge->getDestinationPtrNC();
        destReach = destNode->getReach();
        internalEdge.push_back(curEdge);

        if (origReach != destReach) {
          // External edge has destination in an actual reach
          // Must save the external node as a flux node for the boundary
          fluxNode.insert(destNode);

          // Must save all CCW edges for the destination as flux edges
          // For each flux edge save the destination node because we have to
          // know whether edge has flow allowed, but do not save the edges
          // from these ending nodes
          bEdge* startEdge = destNode->getEdg();
          bEdge* ccwEdge = startEdge;
          do {
             fluxEdge.insert(ccwEdge);
             fluxNode.insert((bFlowNode*) ccwEdge->getDestinationPtrNC());
             ccwEdge = ccwEdge->getCCWEdg();
          } while (ccwEdge != startEdge);
        }
      curEdge = curEdge->getCCWEdg();
    } while (curEdge != firstEdge);

    WriteFlowNode(nodeStr, origNode);
  }

  // Insert the head and outlet nodes of every reach as flux nodes
  // so that tKinemat can find them in the inactive list
  std::list<bFlowNode*>& hlist = flow->getReachHeadList();
  std::list<bFlowNode*>& olist = flow->getReachOutletList();
  std::_List_iterator<bFlowNode*> HeadIter;
  std::_List_iterator<bFlowNode*> OutletIter;

  for (HeadIter = hlist.begin(), OutletIter = olist.begin();
       HeadIter != hlist.end(); HeadIter++, OutletIter++) {
      fluxNode.insert((*HeadIter));
      fluxNode.insert((*OutletIter));
  }

  // Write internal edges for this reach
  for (eiter = internalEdge.begin(); eiter != internalEdge.end(); eiter++)
     WriteFlowEdge(edgeStr, (*eiter));

  // Write external edges for this reach
  for (eiter = externalEdge.begin(); eiter != externalEdge.end(); eiter++)
     WriteFlowEdge(edgeStr, (*eiter));

  // Write flux nodes for this reach
  for (fniter = fluxNode.begin(); fniter != fluxNode.end(); fniter++)
     WriteFlowNode(fluxNodeStr, (*fniter));

  // Write flux edges for this reach
  for (feiter = fluxEdge.begin(); feiter != fluxEdge.end(); feiter++)
     WriteFlowEdge(fluxEdgeStr, (*feiter));

  // Write directory information
  reachStr << "-1\t" 

           << (totalNodeCount - activeNodeCount) << "\t" 
           << nodeOffset << "\t"

           << fluxNode.size() << "\t"
           << fluxNodeOffset << "\t"

           << internalEdge.size() << "\t"
           << externalEdge.size() << "\t"
           << edgeOffset << "\t"

           << fluxEdge.size() << "\t"
           << fluxEdgeOffset << endl;

  // Clear the temporary lists
  internalEdge.erase(internalEdge.begin(), internalEdge.end());
  externalEdge.erase(externalEdge.begin(), externalEdge.end());
  fluxNode.erase(fluxNode.begin(), fluxNode.end());
  fluxEdge.erase(fluxEdge.begin(), fluxEdge.end());

  // Write connectivity information to the directory
  // Write reach head node IDs and reach outlet node IDs
  for (int i = 0; i < numberOfReaches; i++)
    reachStr << "Reach/Head/Outlet/Above\t\t"<< i << "\t" 
             << headID[i] << "\t" 
             << outletID[i] << "\t" 
             << nodeAboveID[i] << endl;

  // Write connectivity of stream reaches
  std::vector<int>::const_iterator viter;

  for (int reach = 0; reach < numberOfReaches; reach++) {
    reachStr << "Reach\t" << connection[reach].getID() << endl;

    const std::vector<int>& upstream = connection[reach].getUpstream();
    reachStr << "Upstream\t" << upstream.size() << "\t";
    for (viter = upstream.begin(); viter != upstream.end(); viter++)
      reachStr << (*viter) << "\t";
    reachStr << endl;

    const std::vector<int>& downstream = connection[reach].getDownstream();
    reachStr << "Downstream\t" << downstream.size() << "\t";
    for (viter = downstream.begin(); viter != downstream.end(); viter++)
      reachStr << (*viter) << "\t";
    reachStr << endl;
  }

  nodeStr.close();
  edgeStr.close();
  fluxNodeStr.close();
  fluxEdgeStr.close();
  reachStr.close();
}

/***************************************************************************
**
** Write the flow node information
** 
***************************************************************************/

void bGraph::WriteFlowNode(fstream& nodeStr, bFlowNode* curnode) const
{
  // Size of flow node written to file
  // (int) 8 * 4 = 32 bytes + (double) 20 * 8 = 160 bytes = (total) 192 bytes
  BinaryWrite(nodeStr, curnode->getID());
  BinaryWrite(nodeStr, curnode->getBoundaryFlag());
  BinaryWrite(nodeStr, curnode->getX());
  BinaryWrite(nodeStr, curnode->getY());
  BinaryWrite(nodeStr, curnode->getZ());
  BinaryWrite(nodeStr, curnode->getVArea());
  BinaryWrite(nodeStr, curnode->getEdg()->getID());
  BinaryWrite(nodeStr, curnode->getFloodStatus());
  BinaryWrite(nodeStr, curnode->getTracer());
  BinaryWrite(nodeStr, curnode->getHillPath());
  BinaryWrite(nodeStr, curnode->getTTime());
  BinaryWrite(nodeStr, curnode->getSrf());
  BinaryWrite(nodeStr, curnode->getHsrf());
  BinaryWrite(nodeStr, curnode->getPsrf());
  BinaryWrite(nodeStr, curnode->getSatsrf());
  BinaryWrite(nodeStr, curnode->getSbsrf());
  BinaryWrite(nodeStr, curnode->getEvapoTrans());
  BinaryWrite(nodeStr, curnode->getSoilMoistureSC());
  BinaryWrite(nodeStr, curnode->getSoilMoistureUNSC());
  BinaryWrite(nodeStr, curnode->getRootMoistureSC());
  BinaryWrite(nodeStr, curnode->getContrArea());
  BinaryWrite(nodeStr, curnode->getNwtNew());
  BinaryWrite(nodeStr, curnode->getRain());
  BinaryWrite(nodeStr, curnode->getCurvature());
  BinaryWrite(nodeStr, curnode->getStreamPath());
  BinaryWrite(nodeStr, curnode->getReach());
      
  bEdge* flowEdge = curnode->getFlowEdg();
  int novalue = -999;
  if (flowEdge == 0)
    BinaryWrite(nodeStr, novalue);
  else
    BinaryWrite(nodeStr, curnode->getFlowEdg()->getID());

  bFlowNode* streamNode = curnode->getStreamNode();
  if (streamNode == 0)
    BinaryWrite(nodeStr, novalue);
  else
    BinaryWrite(nodeStr, curnode->getStreamNode()->getID());
}

/***************************************************************************
**
** Write the flow edge information
** 
***************************************************************************/

void bGraph::WriteFlowEdge(fstream& edgeStr, bEdge* curedge) const
{
  // Size of flow edge written to file
  // (int) 4 * 4 = 16 bytes + (double) 5 * 8 = 40 bytes = (total) 56 bytes
  std::vector<double> rvtx(2);
  rvtx = curedge->getRVtx();
  BinaryWrite(edgeStr, curedge->getID());
  BinaryWrite(edgeStr, rvtx[0]);
  BinaryWrite(edgeStr, rvtx[1]);
  BinaryWrite(edgeStr, curedge->getLength());
  BinaryWrite(edgeStr, curedge->getSlope());
  BinaryWrite(edgeStr, curedge->getVEdgLen());
  BinaryWrite(edgeStr, curedge->getOriginPtrNC()->getID());
  BinaryWrite(edgeStr, curedge->getDestinationPtrNC()->getID());
  BinaryWrite(edgeStr, curedge->getCCWEdg()->getID());
}

//=========================================================================
//
//
//                        End of bGraph.cpp
//
//
//=========================================================================
