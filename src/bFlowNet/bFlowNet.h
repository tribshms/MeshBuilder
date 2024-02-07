// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bFlowNet.h: Header for class bFlowNet (see bFlowNet.cpp) based on CHILD
**             and RIBS routines for Stream Net Routing
**
***************************************************************************/

#ifndef BFLOWNET_H
#define BFLOWNET_H

//=========================================================================
//
//
//                  Section 1: bFlowNet Include Statements
//
//
//=========================================================================

#include "src/bMesh/bMesh.h"

//=========================================================================
//
//
//                  Section 2: bFlowNet Class Definitions
//
//
//=========================================================================

class bFlowNet
{
public:
  bFlowNet();
  bFlowNet(bMesh< bFlowNode> *, bInputFile &);
  ~bFlowNet();

  int FindLakeNodeOutlet(bFlowNode*);
  bFlowNode* FindLowNeighbor(bFlowNode*);
  int IsConnected(bFlowNode*, bFlowNode*);  
  int IsConfluence(bFlowNode*, bFlowNode*);           
  int IsStreamHead(bFlowNode*);                
  int IsStreamDisjoint(bEdge*, bFlowNode*, 
                       std::list<bEdge*> &,
                       std::list<bEdge*> &);
  int IsToEliminate(bFlowNode*);

  int CheckNeighbor(bFlowNode*, bEdge*,
                    std::list<bFlowNode*> &,
                    std::list<bEdge*> &);
  int FindStreamDisjoints(bFlowNode*, int, std::list<bFlowNode*> &);
  int FindConfluenceDisjoints(std::list<bFlowNode*> &);

  void SetFlowVariables(bInputFile &);
  void SetBasinOutlet();
  void InitFlowDirs();
  void FlowDirs();
  void SortNodesByNetOrder();
  void FillLakes();
  void initializeTravelTime();
  void initializeTravelTimeOnly();     
  void setTravelTime();
  void setTravelVelocity(double);
  void SurfaceFlow();
  void DrainAreaVoronoi();
  void RouteFlowArea(bFlowNode *, double);
  void DeriveStreamReaches(bInputFile &); 
  void SortStreamNodes();                 
  void TellAboutNode(bFlowNode *cn);       
  void PrintArcInfoLinks(bInputFile &);  
  void WeightedShortestPath(); 
  void AddUnsettledNeighbors(bFlowNode*,std::list<bFlowNode*> &,
                                        std::list<bEdge*> &);
  void UpdatePathVariable(bFlowNode*);
  void CheckVDrainageWidths();
  void FixVoronoiEdgeWidth(bFlowNode *);
  void DeriveCurvature();

  int IsBetweenEndPnts(std::vector<double>&, std::vector<double>&, 
                       std::vector<double>&, std::vector<double>&, 
                       double, double);
  int AreSegmentsParallel(std::vector<double>&, std::vector<double>&, 
                          std::vector<double>&, std::vector<double>&);
  int IsInTriangle(std::vector<double>&, 
                   std::vector<double>&, 
                   std::vector<double>&, 
                   double, double);
  double FindAngle(bFlowNode*, bFlowNode*, bFlowNode*);     
  double FindDistance(double, double, double, double); 
  double CompElevationTendency(bFlowNode*, bFlowNode*, int*, int);   
  double polygonArea(double **, int);
  double getCurrDischarge(int);

  void SetReachInformation();
  void WriteFlowNet() const;
  void WriteGeometry() const;

  bFlowNode* getOutletPtr()			{ return OutletNode;}
  std::list<bFlowNode*>& getReachHeadList()	{ return NodesLstH; }
  std::list<bFlowNode*>& getReachOutletList()	{ return NodesLstO; }
  std::list< int >& getReachSizeList()		{ return NNodes; }

protected:
  bMesh<bFlowNode>*     gridPtr;
  bFlowNode*            OutletNode;    // Ptr to the outlet node, 1 so far...
  std::list<bFlowNode*> HeadsLst;      // List of basin stream heads
  std::list<bFlowNode*> NodesLstH;     // Heads of stream reaches 
  std::list<bFlowNode*> NodesLstO;     // Outlets of stream reaches 
  std::list<int>        NNodes;        // # of nodes in each reach 

  long flowboxes; 	        // Size of current discharge array
  double hillvel;   	        // Hillslope velocity, [m/sec]
  double streamvel; 	        // Stream velocity, [m/sec]
  double velratio;  	        // Ratio Stream/Hillslope
  double velcoef;   	        // Coeffcient velocity-discharge
  double flowexp;   	        // Power in velocity-discharge relationship
  double baseflow;  	        // Baseflow: leftover from event-based
  double flowout;   	        // Discharge at the outlet, [m^3/sec]
  double maxttime;  	        // MAX travel time defined for watershed, [hour]
  double dist_hill_max;  	// MAX distance on hillslope, [m] 
  double dist_stream_max;	// MAX distance in stream, [m]
  double BasArea;               // Total Basin Area, [m^2]
};

#endif
	
//=========================================================================
//
//
//                          End bFlowNet.h
//
//
//=========================================================================
