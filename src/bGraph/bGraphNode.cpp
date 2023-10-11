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
**  bGraphNode.cpp: Functions for class bGraphNode (see bGraphNode.h)
**
***************************************************************************/

#include "src/bGraph/bGraphNode.h"


/**************************************************************************
**
** Constructor that only creates object.
**
**************************************************************************/

bGraphNode::bGraphNode() : id(-1) {
}

/*************************************************************************
**
** Constructor that takes an id.
**
*************************************************************************/

bGraphNode::bGraphNode(int nid) : id(nid) {
}

/*************************************************************************
**
** Constructor that takes an id, upstream nodes, downstream nodes.
**
*************************************************************************/

bGraphNode::bGraphNode(int nid, std::vector<int>& nup, 
  std::vector<int>& ndown) : id(nid) {

  for (int i = 0; i < nup.size(); i++) {
    upstream.push_back(nup[i]);
  }
  for (int i = 0; i < ndown.size(); i++) {
    upstream.push_back(ndown[i]);
  }
}

/*************************************************************************
**
** Destructor
**
*************************************************************************/

bGraphNode::~bGraphNode(){
  id = -1;
  upstream.erase(upstream.begin(), upstream.end());
  downstream.erase(downstream.begin(), downstream.end());
} 

/*************************************************************************
**
** Write out bGraphNode
**
*************************************************************************/

std::ostream& operator<<(std::ostream& out, const bGraphNode& n)
{
  out << "bGraphNode:Reach: " << n.getID();
  std::vector<int> down = n.getDownstream();
  std::vector<int> up = n.getUpstream();
  out << "\tupstream = ";
  for (int i = 0; i < up.size(); i++)
    out << up[i] << ",";
  out << "\tdownstream = ";
  for (int i = 0; i < down.size(); i++) 
    out << down[i] << ",";
  return out;
}

//=========================================================================
//
//
//                        End of bGraphNode.cpp
//
//
//=========================================================================
