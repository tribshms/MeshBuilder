// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bGraphNode.h: Header for tBasin class and objects
**
**  bGraphNode Class used in tRIBS for the parallel version, providing 
**  information about upstream and downstream basins as graph nodes
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: bGraphNode Include and Define Statements
//
//
//=========================================================================

#ifndef BGRAPHNODE_H
#define BGRAPHNODE_H

#include <fstream>
#include <iostream>
#include <vector>


//=========================================================================
//
//
//                  Section 2: bGraphNode Class Definitions
//
//
//=========================================================================

class bGraphNode {

public:
  /// Default Constructor
  bGraphNode(); 
  /// Constructor that takes an ID
  bGraphNode(int nid);
  /// Constructor that takes an ID, upstream nodes, downstream nodes
  bGraphNode(int nid, std::vector<int>& nup, std::vector<int>& ndown);
  /// Destructor
  ~bGraphNode();

  /// Does another node flow into this node?
  bool hasUpstream() { if (upstream.size() > 0) return true;
                       else return false; }
  /// Does this node flow into another?
  bool hasDownstream() { if (downstream.size() > 0) return true;
                         else return false; }

  /// Return ID
  int getID() const { return id; }

  /// Return list of upstream nodes
  const std::vector<int>& getUpstream() const { return upstream; }
  /// Return list of downstream nodes
  const std::vector<int>& getDownstream() const { return downstream; }

  /// Return list of flux nodes
  std::vector<int> getFlux() const { return flux; }

  /// Return number of upstream nodes
  int getUpstreamCount() const { return upstream.size(); }
  // Return number of downstream nodes
  int getDownstreamCount() const { return downstream.size(); }
  /// Return number of flux nodes
  int getFluxCount() const { return flux.size(); }

  /// Add upstream node
  void addUpstream(int n) { upstream.push_back(n); }
  /// Add downstream node
  void addDownstream(int n) { downstream.push_back(n); }

  /// Add flux node
  void addFlux(int n) { flux.push_back(n); }

  /// Set ID
  void setId(int nid) { id = nid; }

private:

  int id;                        //!< Node ID
  std::vector<int> upstream;     //!< List of upstream nodes 
  std::vector<int> downstream;   //!< List of downstream nodes
  std::vector<int> flux;        //!< List of flux nodes
};

//! Write the bGraphNode to a given stream
std::ostream& operator<<(std::ostream& out, const bGraphNode& n);

#endif

//=========================================================================
//
//
//                          End of bGraphNode.h 
//
//
//=========================================================================
