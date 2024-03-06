// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bFlowNet.cpp: Functions for class bFlowNet (see bFlowNet.h) based on
**             CHILD and RIBS routines for Stream Net Routing
**
***************************************************************************/

#include "src/bFlowNet/bFlowNet.h"
#include <ctime>
#include <cstring>

//=========================================================================
//
//
//                  Section 1: bFlowNet Constructors/Destructors
//
//
//=========================================================================

bFlowNet::bFlowNet() 
{
  gridPtr = 0;
}

bFlowNet::bFlowNet(bMesh<bFlowNode> *gridRef, bInputFile &infile) 
{
  gridPtr = gridRef;
  BasArea = 0.0;

  time_t clock;
  clock = time(NULL);
  cout <<"\nSetFlowVariables..." << ctime(&clock) << endl;
  SetFlowVariables( infile );
    
  clock = time(NULL);
  cout <<"\nInitFlowDirs..." << ctime(&clock) <<endl;
  InitFlowDirs();

  clock = time(NULL);
  cout <<"\nFlowDirs..." << ctime(&clock) <<endl;
  FlowDirs();

  clock = time(NULL);
  cout <<"\nFillLakes..." << ctime(&clock) <<endl;
  FillLakes();

  clock = time(NULL);
  cout <<"\nSortNodesByNetOrder..." << ctime(&clock) <<endl;
  SortNodesByNetOrder();

  clock = time(NULL);
  cout <<"\nSetBasinOutlet..." << ctime(&clock) <<endl;
  SetBasinOutlet();

  clock = time(NULL);
  cout <<"\nWeightedShortestPath..." << ctime(&clock) <<endl;
  WeightedShortestPath();  

  clock = time(NULL);
  cout <<"\nSortStreamNodes..." << ctime(&clock) <<endl;
  SortStreamNodes();

  clock = time(NULL);
  cout <<"\nDrainAreaVoronoi..." << ctime(&clock) <<endl;
  DrainAreaVoronoi();

  clock = time(NULL);
  cout <<"\nDeriveCurvature..." << ctime(&clock) <<endl;
  DeriveCurvature();

  clock = time(NULL);
  cout <<"\nDeriveStreamReaches..." << ctime(&clock) <<endl;
  DeriveStreamReaches( infile );

  clock = time(NULL);
  cout <<"\nInitializeTravelTime..." << ctime(&clock) <<endl;
  setTravelVelocity(0.0);  
  initializeTravelTime();

  clock = time(NULL);
  cout <<"\nCheckVDrainageWidths..." << ctime(&clock) <<endl;
  CheckVDrainageWidths();

  cout <<"\nSet reach numbers..."<<endl;
  SetReachInformation();

  clock = time(NULL);
  cout <<"\nDone..." << ctime(&clock) <<endl;
}

bFlowNet::~bFlowNet() 
{
   gridPtr = NULL;
   cout<<"bFlowNet Object has been destroyed..."<<endl<<flush;
}

//=========================================================================
//
//
//                  Section 2: bFlowNet Slopes/Direction Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  SetFlowVariables()
**  
**  Reads from InputFile and assigns basic varibales used in bFlowNet
**  As opposed to work conducted on mesh by other functions in constructor
**
*****************************************************************************/
void bFlowNet::SetFlowVariables(bInputFile &infile)
{
  velratio = infile.ReadItem(velratio, "VELOCITYRATIO");
  baseflow = infile.ReadItem(baseflow, "BASEFLOW");
  velcoef  = infile.ReadItem(velcoef,  "VELOCITYCOEF");
  flowexp  = infile.ReadItem(flowexp,  "FLOWEXP");

  flowout   = 0.0;
  streamvel = 0.0;
  return;
}

/*****************************************************************************
**  
**  bFlowNet::InitFlowDirs()
**  
**  Defines flow direction of flow for each Node  
**
*****************************************************************************/
void bFlowNet::InitFlowDirs()
{
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter = nodeList->begin();
  bFlowNode * curnode;
  bEdge  * flowedg;
  std::_List_iterator<bFlowNode*> nodeToMove;
  int ctr;
  int kMaxSpokes = 100; 

  // For every active (non-boundary) node, initialize it to flow to a
  // _non-boundary_ node (i.e., to a hillslope or stream node )

  curnode = (*nodIter);
  while (nodIter != nodeList->getLastActive()) {
    // Start with the node's default edge
    flowedg = curnode->getEdg();

    // As long as the current edge is a no-flow edge, 
    // advance to the next one, counter-clockwise
    ctr = 0;
    while ( (!flowedg->FlowAllowed() || 
	     flowedg->getVEdgLen() < THRESH)) {

      flowedg = flowedg->getCCWEdg();
      ctr++;
      if (ctr > kMaxSpokes) { //Make sure to prevent endless loops
	cout << "\nError in InitFlowDirs(): Node " << curnode->getID()
	     << " surrounded by closed boundary nodes ..." << endl;
	cout<<"\nNode boundary will be changed to kClosedBoundary..."<<endl;
	cout<<curnode->getID()
	    <<"\t"<<curnode->getX()<<"\t"<<curnode->getY()
	    <<"\t"<<curnode->getZ()<<"\t"<<curnode->getBoundaryFlag()
	    <<endl<<endl<<flush; 

	// Changing the boundary code and going to the next node
	curnode->setBoundaryFlag( kClosedBoundary );
	nodeToMove = nodIter;
        nodIter++;
	curnode = (*nodIter);
	flowedg = curnode->getEdg();
	ctr = 0;

	// Now, move that node to the end of the list
	nodeList->moveToBack( nodeToMove );
      }
    }
    curnode->setFlowEdg( flowedg ); //Flowedge for current node is set
    nodIter++;
    curnode = (*nodIter);
  }
  return;
}

/*****************************************************************************
**  
**  bFlowNet::FlowDirs()
**  
**  Defines flow direction of flow for each Node  
**
*****************************************************************************/
void bFlowNet::FlowDirs() 
{
  int ctr;
  double slp;        // steepest slope found so far
  int kLargeNegative  = -1000;
  int kMaxSpokes = 100;

  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter = nodeList->begin();
  bFlowNode *curnode;  // ptr to the current node
  bEdge * firstedg; // ptr to first edge
  bEdge * curedg;   // pointer to current edge
  bEdge * nbredg;   // steepest neighbouring edge so far

  cout.setf( ios::fixed, ios::floatfield);

  // Find the connected edge with the steepest slope 
  curnode = (*nodIter);
  while (nodIter != nodeList->getLastActive()) {

    firstedg = curnode->getFlowEdg();
    slp    = firstedg->getSlope();
    nbredg = firstedg;
    curedg = firstedg->getCCWEdg();
    ctr = 0;         
         
    // Check each of the various spokes, stopping when
    // we've gotten back to the beginning
    while ( curedg != firstedg )  {
      if (curedg->getSlope() > slp && 
	  curedg->getDestinationPtrNC()->getBoundaryFlag() != kClosedBoundary &&
	  curedg->getVEdgLen() > THRESH)
      {
	 slp    = curedg->getSlope();
	 nbredg = curedg;
      }
      curedg = curedg->getCCWEdg();
      ctr++;

      if (ctr > kMaxSpokes) { // Make sure to prevent endless loops
        cout<<"\nError in FlowDirs(): Node "<<curnode->getID()
	    << " surrounded by closed boundary nodes" << endl;
      }
    }
    curnode->setFlowEdg( nbredg ); // The steepest flowedge

    if ((slp <= 0.0) && (curnode->getBoundaryFlag() == kStream )) {
      //TellAboutNode(curnode);
    }

    // This is needed to check for pits
    if ( (slp > 0.0) && (curnode->getBoundaryFlag() != kClosedBoundary))
      curnode->setFloodStatus( kNotFlooded );
    else {
      curnode->setFloodStatus( kSink );
    }
    nodIter++;
    curnode = (*nodIter);
  }
  return;
}

//=========================================================================
//
//
//                  Section 3: bFlowNet Sorting Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  bFlowNet::SortNodesByNetOrder()
**  
**  This function sorts the list of nodes according to their order in the
**  network (upstream to downstream). The sorting algorithm is based on the 
**  "cascade" algorithm of Braun and Sambridge. 
**
**     The single-direction sorting algorithm works by initially assigning a 
**  tracer (like a packet of water) to each node. At each iteration, a tracer 
**  from each node is sent downstream. Any nodes that have zero tracers left 
**  are moved to the bottom of the list (a FIFO stack), so that for example 
**  the very first node moved will be the first node on the list when the 
**  sorting is completed. The process continues until no unsorted nodes remain.
**
**    The multi-flow option was added to allow for multiple flow directions
**  and kinematic-wave routing. The algorithm is slightly different. At
**  each pass, the unsorted nodes are "de-flagged" by setting their
**  tracer variables to zero. Then each unsorted node flags _all_ of the
**  adjacent nodes that are downhill by setting their tracer to 1 (this
**  is accomplished through a call to ActivateSortTracer). Finally, any
**  unflagged nodes are moved to the back of the list, and the process
**  is repeated until all nodes have been sorted.
**
**  It is possible that the "flagging" method, as opposed to the
**  "tracer movement" method, is more efficient even for single-flow
**  routing. The added cost lies in unflagging the unsorted nodes at
**  each pass. However, the gain comes from the fact that you never have
**  multiple "tracers" at a node that need to be removed one by one.
**  The two methods should be tested and compared. Yes - do the tracer reset.
**
*****************************************************************************/
void bFlowNet::SortNodesByNetOrder()
{
  int nThisPass;                   // Number moved in current iteration
  int i;
  int done=0;

  bFlowNode * cn;
  bMeshList<bFlowNode*> *nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> listIter = nodeList->begin();
  int nUnsortedNodes = nodeList->getActiveSize(); // Number not yet sorted

  // Create a temporary node list to hold nodes in flow order
  std::list<bFlowNode*> flowList;

  // Assign initial tracers, one per node
  for (listIter = nodeList->begin(); 
       listIter != nodeList->getLastActive(); listIter++)
       (*listIter)->ActivateSortTracer();

  // Move tracers downstream and remove until no nodes with tracers remain
  do {
    // Send tracers downstream
    listIter = nodeList->begin();
    cn = (*listIter);
    for ( i=0; i < nUnsortedNodes; i++ ) {
      assert( cn != 0 );
      cn->MoveSortTracerDownstream();
      listIter++;
      cn = (*listIter);
    }

    // Scan for nodes with no tracers, and move to temporary flow list in order
    std::_List_iterator<bFlowNode*> nodeToMove;
    nThisPass = 0;
    done = true;

    listIter = nodeList->begin();
    cn = (*listIter);
    for ( i=0; i < nUnsortedNodes; i++ ) {
      if ( cn->NoMoreTracers() ) {
        nodeToMove = listIter;
        listIter++;
        cn = (*listIter);

        // Add the node without a tracer to temporary flow list
        flowList.push_back((*nodeToMove));

        // Remove it from the current node list
        nodeList->removeActive(nodeToMove);
        nThisPass++;
      }
      else {
        listIter++;
        cn = (*listIter);
        done = false;
      }
    }

    nUnsortedNodes -= nThisPass;

      if((nThisPass ==  0) & (nUnsortedNodes>0)) {
          throw runtime_error("No more tracers identified, but unsortedNodes remain.\n Check Points File");
      }


      // Reset all active node tracers to value 1
    for (listIter = nodeList->begin(); 
         listIter != nodeList->getLastActive(); listIter++) {
         (*listIter)->ActivateSortTracer();
    }

  } while (done == false);

  // Move temporary flow list in order to the active part of the node list
  // which should have only inactive nodes at this point
  std::list<bFlowNode*>::iterator flowIter;
  for (flowIter = flowList.begin(); flowIter != flowList.end(); flowIter++)
    nodeList->insertAtActiveBack((*flowIter));

  // Changed to make edge IDs consistent
  bEdge     * ce;
  bTriangle * ct;  
  nodeList = gridPtr->getNodeList();
  bMeshList<bEdge*>* edgeList = gridPtr->getEdgeList();
  std::list<bTriangle*>* triList = gridPtr->getTriList();

  std::_List_iterator<bFlowNode*> niter;
  std::_List_iterator<bEdge*> eiter;
  std::_List_iterator<bTriangle*> titer;

  cout <<"\nRenumbering nodes, edges, triangles..."<<endl;
  int id;
  for (niter = nodeList->begin(), id = 0; 
       niter != nodeList->end(); niter++, id++)
    (*niter)->setID(id);
  for (eiter = edgeList->begin(), id = 0; 
       eiter != edgeList->end(); eiter++, id++)
    (*eiter)->setID(id);
  for (titer = triList->begin(), id = 0; 
       titer != triList->end(); titer++, id++)
    (*titer)->setID(id);

  return;
}

//=========================================================================
//
//
//                  Section 4: bFlowNet Velocity Functions
//
//
//=========================================================================

/*****************************************************************************
**  
**  bFlowNet::setTravelVelocity()
**  
**  Sets travel velocities depending on the discharge at the oulet
**
*****************************************************************************/
void bFlowNet::setTravelVelocity(double curr_discharge) 
{
  // Get the node velocities first 
  if ( !curr_discharge )  { 
     // If it's zero, a baseflow value from the input
     // file is used to define the stream velocity     
     if ( !baseflow && flowexp ) {
       cout<<"\nWarning: Baseflow is zero and thus the lower "
	   <<"limit of stream velocity is undefined --> Set to 0.001"<<endl;
       baseflow = 0.001; 
     }
     flowout = baseflow; //velocity will correspond to baseflow
  }
  else 
     flowout = curr_discharge;

  streamvel = velcoef;
  hillvel   = streamvel/velratio;
  
  return;
}

/*****************************************************************************
** 
**  bFlowNet::initializeTravelTime()
**  
**  Needed to initialize some of the class members  
**  Sets the stream node for every node which is used for deciding reach
**
*****************************************************************************/
void bFlowNet::initializeTravelTime() 
{
  bFlowNode *cn;
  bFlowNode *ctimer;
  bEdge  *ce;
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter;
  double hill, stream, tt;

  // Loop through the nodes and set velocity
  // Define maximum travel time        
  BasArea = 0.0;
  maxttime = 0.0;     		// SECONDS
  dist_hill_max = 0.0;   	// METERS
  dist_stream_max = 0.0; 	// METERS

  for (nodIter = nodeList->begin();
       nodIter != nodeList->getLastActive(); nodIter++) {
    cn = (*nodIter);
    BasArea += cn->getVArea();

    // Initialize with distance from centroid to node
    cn->getCentroidX();
    cn->getCentroidY();
    hill = sqrt(pow(cn->getX() - cn->getCentroidX(), 2.0) + 
                pow(cn->getY() - cn->getCentroidY(), 2.0));
    stream = tt = 0.0;

    ctimer = cn;
    assert(ctimer != 0);

    // If it is an hillslope, it flows at a hillslope velocity
    if (ctimer->getBoundaryFlag() == kNonBoundary) {
      while (ctimer->getBoundaryFlag() == kNonBoundary) {
	ce = ctimer->getFlowEdg(); //Get the steepest edge
	hill += ce->getLength();
	ctimer = ctimer->getDownstrmNbr();
      }
      // Set a stream node to which it contributes flow
      assert(ctimer != 0);
      cn->setStreamNode(ctimer);
    }
    
    // If it is a streamnode, it flows at a stream velocity
    // Go downstream to the outlet
    while ( ctimer->getDownstrmNbr() ) {
      ce = ctimer->getFlowEdg();  // Get the steepest edge
      if (ctimer->getBoundaryFlag() == kStream || 
	  ctimer->getBoundaryFlag() == kOpenBoundary)
	stream += ce->getLength();            // METERS
      ctimer = ctimer->getDownstrmNbr();
    }

    cn->setHillPath(hill);     		// Set the HILLSLOPE path for node
    cn->setStreamPath(stream); 		// Set the STREAM path for node

    tt = hill/hillvel + stream/streamvel; // SECONDS 
    tt /= 3600.;                          // HOURS  

    // Set travel time: travel time is only suitable for HYDROLOGIC routing 
    // The KINEMATIC scheme uses independent calculation of travel time

    cn->setTTime( tt );        		// Set Travel Time in HOURS

    if (tt > maxttime) {  //MAX values of paths with maxttime
      maxttime = tt;
      dist_hill_max = hill;       	// METERS
      dist_stream_max = stream;   	// METERS
    }
  }
  cout.setf( ios::fixed, ios::floatfield);
  cout<<endl;
  cout<<"Flow Characteristics: "<<endl;
  cout<<"The  Total  Basin  Area is   \t"<<BasArea<<" M^2"<<endl;
  cout<<"Maximum travel time: \t\t"<<maxttime<<" hours"<<endl;
  cout<<"Maximum hillslope path: \t"<<dist_hill_max<<" meters"<<endl;
  cout<<"Maximum stream path: \t\t"<<dist_stream_max<<" meters"<<endl;
  cout<<endl;

  return;
}

/*****************************************************************************
**  
**  initializeTravelTimeOnly()
**  
**  Difference with the preceeding function is that this one does  
**  not compute distances and deals only with the travel times
** 
*****************************************************************************/
void bFlowNet::initializeTravelTimeOnly() 
{
  bFlowNode *cn;
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter;
  double hill, stream, tt;

  // Loop through the nodes and set velocity
  // Define maximum travel time 
  maxttime = 0.0;        // HOURS
  dist_hill_max = 0.0;   // METERS
  dist_stream_max = 0.0; // METERS

  for (nodIter = nodeList->begin();
       nodIter != nodeList->getLastActive(); nodIter++) {
    cn = (*nodIter);
    hill   = cn->getHillPath();     // Get the HILLSLOPE path for node
    stream = cn->getStreamPath();   // Get the STREAM path for node
    tt = hill/hillvel + stream/streamvel;   //SECONDS
    tt /= 3600.0;                           //HOURS
    cn->setTTime( tt );             // Set Travel Time in HOURS

    if (tt > maxttime) {            // MAX values of paths with maxttime
        maxttime = tt;
        dist_hill_max = hill;       // METERS
        dist_stream_max = stream;   // METERS
    }
  }

  cout<<"\nMaximum Travel Time: \t\t"<<maxttime<<" hours"<<endl;
  cout<<"Maximum Hillslope path: \t"<<dist_hill_max<<" meters"<<endl;
  cout<<"Maximum Stream path: \t\t"<<dist_stream_max<<" meters"<<endl;

  return;
}

//=========================================================================
//
//
//                  Section 5: bFlowNet FillLakes Function
//
//
//=========================================================================

/*****************************************************************************
**
**  bFlowNet::FillLakes 
**
**  Finds drainage for closed depressions. The algorithm assumes
**  that sinks (nodes that are lower than any of their neighbors)
**  have already been identified during the flow directions
**  procedure. For each sink, the algorithm creates a list of
**  nodes in the current lake, which initially is just the sink
**  itself. It then iteratively looks for the lowest node on the
**  perimeter of the current lake. That node is checked to see 
**  whether it can be an outlet, meaning that one of its
**  neighbors is both lower than itself and is not already
**  flooded (or is an open boundary). If the low node is not an 
**  outlet, it is added to the current lake and the process 
**  is repeated. If it is an outlet, then all of the nodes on the 
**  current-lake list are identified draining it. The list is then
**  cleared, and the next sink is processed. If during the search
**  for the lowest node on the perimeter a flooded node is found
**  that isn't already part of the current lake (i.e., it was
**  flagged as a lake node when a previous sink was processed),
**  then it is simply added to the current-lake list --- in other
**  words, the "new" lake absorbs any "old" ones that are encountered.
**
**  Once an outlet has been found, flow directions for nodes in the
**  lake are resolved in order to create a contiguous path through
**  the lake.
**
**    Calls: FindLakeNodeOutlet
**    Called by: MakeFlow
**    Modifies:  flow direction and flood status flag of affected nodes
**
*****************************************************************************/
void bFlowNet::FillLakes()
{
   bFlowNode *cn;		// Node on list: if a sink, then process
   bFlowNode *thenode;		// Node on lake perimeter
   bFlowNode *lowestNode;	// Lowest node on perimeter found so far
   bFlowNode *cln;		// Current lake node
   bFlowNode *node;              // Placeholder

   bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();

   std::vector<bFlowNode*> lakeArray;		// Nodes in current lake
   std::vector<bFlowNode*> lowNodeArray;	// Lowest of neighbors of node

   bEdge *ce;              
   double lowestElev;      // Lowest elevation found so far on lake perimeter
   int done;               // Flag indicating whether outlet has been found

   double lowNeighborElev;	// Keep lowest elevation of every lake node
   bFlowNode *lowNeighbor;	// Keep lowest neighbor of every lake node

   int lakeIndex = 0;		// Index into lakeArray

   // Check each active node to see whether it is a sink
   std::_List_iterator<bFlowNode*> niter = nodeList->begin();
   for (; niter != nodeList->getLastActive(); niter++) {
      cn = (*niter);

      // Find lakes around sink nodes (lower elevation than any neighbor)
      if (cn->getFloodStatus() == kSink)
      {
         // Create a new lake-list, initially containing just the sink node.
         lakeArray.push_back( cn );
         cn->setFloodStatus( kCurrentLake );
         
         // lowestNode is either the sink node or the lowest neighbor
         // of all the nodes on the lakeArray
         // As new nodes are added to the lakeArray, only those nodes are
         // searched for low neighbors.  Results from previously processed
         // nodes are stored in the lowNodeArray
         done = false;
         lowestNode = cn;
         lowestElev = kMaxElevation;
         while (done == false) {
            // Check the neighbors of current lake node for lowest and store
            // This is so subsequent searched for low nodes can be done fast
            // Also look for the lowest overall
            while (lakeIndex < lakeArray.size())
            {
               cln = lakeArray[lakeIndex];
               ce = cln->getEdg();
               bFlowNode* lowNeighbor = 0;
               double lowNeighborElev = kMaxElevation;
               do {
                  thenode = (bFlowNode *) ce->getDestinationPtrNC();
                  // Is it a potential outlet (not flooded, not basin boundary)
                  if (thenode->getBoundaryFlag() == kOpenBoundary ||
		     (thenode->getFloodStatus()  == kNotFlooded &&
		      thenode->getBoundaryFlag() != kClosedBoundary))
                  { 
		     // If neighbor is basin outlet OR lower than lowest so far
		     if (thenode->getBoundaryFlag() == kOpenBoundary || 
			 thenode->getZ() < lowestElev)
                     {
		         lowestNode = thenode;
                         lowestElev = thenode->getZ();
		     }
                     // Record the lowest neighbor even if not lowest overall
		     if (thenode->getBoundaryFlag() == kOpenBoundary || 
			 thenode->getZ() < lowNeighborElev)
                     {
		         lowNeighbor = thenode;
                         lowNeighborElev = thenode->getZ();
                     }
                  }
                  // If it's a previous lake node or a sink, add it to the list
                  else if ((thenode->getFloodStatus()  == kFlooded ||
                            thenode->getFloodStatus()  == kSink) && 
			   (thenode->getBoundaryFlag() == kStream ||
			    thenode->getBoundaryFlag() == kNonBoundary))  {

                     // Problem if no low node was found and we have to use sink
                    bool found = false;
                    for (int j = 0; j < lakeArray.size(); j++)
                       if (lakeArray[j] == thenode)
                          found = true;
                    if (found == true) {
                        cout << "FillLakes: add node already in lakelist " 
                             << thenode->getID()
                            << " x " << thenode->getX()
                            << " y " << thenode->getY()
                            << " z " << thenode->getZ()
                            << endl;
                     } else {
                        lakeArray.push_back(thenode);
                        thenode->setFloodStatus( kCurrentLake );
                     }
                  }
               }  while ( ( ce=ce->getCCWEdg() ) != cln->getEdg() );

               // After looking at all neighbors record the lowest
               lowNodeArray.push_back(lowNeighbor);
               lakeIndex++;
            } // Process next node in the lake array

            // Now we've found the single lowest point on the perimeter
            // Test to see whether it's an outlet and get out of the loop
	    // 1.) If it's an open boundary, it's an outlet
            if ( lowestNode->getBoundaryFlag() == kOpenBoundary ) 
	       done = true;

	    // 2.) it's also an outlet if it can drain to a "dry" location.
            else  { 
               if (FindLakeNodeOutlet(lowestNode)) {
		  done = true;

               // lowest node is not an outlet so add it to lake Array
               // Since it no longer can be the lowest node for the
               // lake nodes that contributed it, all lake nodes which
               // report it to be their lowest must be reprocessed
               } else {
		  if (lowestNode->getBoundaryFlag() == kStream ||
		      lowestNode->getBoundaryFlag() == kNonBoundary) {
                    // Verify that we aren't creating a loop
                    bool found = false;
                    for (int j = 0; j < lakeArray.size(); j++)
                       if (lakeArray[j] == lowestNode)
                          found = true;
                    if (found == true) {
                       cout << "FillLakes: add low node already in lakelist " 
                            << lowestNode->getID() 
                            << " x " << lowestNode->getX()
                            << " y " << lowestNode->getY()
                            << " z " << lowestNode->getZ()
                            << endl;
                       done = true;
                    } else {
		       lakeArray.push_back( lowestNode );
		       lowestNode->setFloodStatus( kCurrentLake );

                       // Iterate on all current low nodes to find the
                       // next lowest to start the next pass with
                       bFlowNode* newLakeNode = lowestNode;
                       lowestNode = 0;
                       lowestElev = kMaxElevation;
                       for (int j = 0; j < lowNodeArray.size(); j++)
                       {
                          // Every lakelist node that reports current lowestNode
                          // as its lowest must find a new lowest
                          if (lowNodeArray[j] == newLakeNode) {
                             lowNodeArray[j] = 
                                FindLowNeighbor(lakeArray[j]);
                          }

                          // Locate the new lowest node overall
                          if (lowNodeArray[j] != 0)
                            if (lowNodeArray[j]->getZ() < lowestElev) {
                             lowestNode = lowNodeArray[j];
                             lowestElev = lowNodeArray[j]->getZ();
                          }
                       }
                       // If we find no lowest node use the sink node
                       if (lowestNode == 0) {
                          lowestNode = lakeArray[0];
                          lowestElev = kMaxElevation;
                       }
                    }
		  }
               }
            }
	    assert(lakeArray.size() <= gridPtr->getNodeList()->getActiveSize());
         } while ( !done );

         // Part 2:
         // Now we've found an outlet for the current lake.
         // This next bit of code assigns a flowsTo for each node so there's
         // a complete flow path toward the lake's outlet. This isn't strictly
         // necessary --- the nodes could all point directly to the outlet,
         // skipping anything in between --- but it prevents potential problems
         // in ordering the list by network order. This also works by pointing
         // each node toward the first neighboring node they happen to find
         // that has been flagged as having its flow direction resolved. 
         // Initially, the low node is thus flagged, and the algorithm repeats
         // until all the lake nodes are flagged as having a flow direction.
         // The algorithm isn't unique---there are many paths that could be
         // taken; this simply finds the most convenient one.

         lowestNode->setFloodStatus( kOutletFlag );
         do
         {
            done = true;  // assume done until proven otherwise
            // Iterate over entire lake array starting with the sink looking
            // looking for the nodes that touch the just found outlet node
            // Could we just start with the outlet node and back up instead??
            for (int j = 0; j < lakeArray.size(); j++) {
               cln = lakeArray[j];
               if ( cln->getFloodStatus() != kOutletFlag ) {
                  done = false;
                  ce = cln->getEdg();		// Check each neighbor
                  do {
                     node = (bFlowNode *) ce->getDestinationPtrNC();
                     if ( node->getFloodStatus() == kOutletFlag ) {
                        cln->setFloodStatus( kOutletPreFlag );
                        cln->setFlowEdg( ce );
                        
                     }
                  } while ( cln->getFloodStatus() != kOutletFlag
                           && ( ce=ce->getCCWEdg() ) != cln->getEdg() );
               } 
	    } // END for each lake node

            // Now flag all the "preflagged" lake nodes as outlets
            for (int j = 0; j < lakeArray.size(); j++) {
               cln = lakeArray[j];
               if ( cln->getFloodStatus() == kOutletPreFlag )
                  cln->setFloodStatus( kOutletFlag );
            }
         } while ( !done );

         lowestNode->setFloodStatus( kNotFlooded );
         
         // Finally, flag all of the nodes in it as "kFlooded" 
	 // and clear the list so we can move on to the next sink
         for (int j = 0; j < lakeArray.size(); j++) {
            lakeArray[j]->setFloodStatus(kFlooded);
         }
         lakeArray.erase(lakeArray.begin(), lakeArray.end());
         lowNodeArray.erase(lowNodeArray.begin(), lowNodeArray.end());
         lakeIndex = 0;
      } // END if Sink
   } // END Active Nodes

   return;  
}

/*****************************************************************************
**
**  FindLowestNeighbor
**
**  To make FillLakes() more efficient, when processing a lakeList node
**  for the first time, save the lowest neighbor.  Also remember which
**  lakeList node produced the lowest of all nodes, because on the next
**  pass through the lake list (looking for an outlet) the node which
**  contributed the lowest in the previous pass must now calculate a
**  new lowest.  Also and new node added to the lake list must calculate
**  its lowest node.
**
*****************************************************************************/
bFlowNode* bFlowNet::FindLowNeighbor( bFlowNode *cln )
{
   // Check all the neighbor nodes of perimeter nodes
   bEdge* ce = cln->getEdg();
   bFlowNode* lowNeighbor = 0;
   bFlowNode* thenode;
   double lowNeighborElev = kMaxElevation;

   do
   {
      thenode = (bFlowNode *) ce->getDestinationPtrNC();
      // Is it a potential outlet (not flooded, not basin boundary)
      if (thenode->getBoundaryFlag() == kOpenBoundary ||
		(thenode->getFloodStatus()  == kNotFlooded &&
		 thenode->getBoundaryFlag() != kClosedBoundary))
      { 
         // If it's the basin outlet OR lower than the lowest found so far?
	 if (thenode->getBoundaryFlag() == kOpenBoundary || 
	     thenode->getZ() < lowNeighborElev)
         {
	    lowNeighbor = thenode;
            lowNeighborElev = thenode->getZ();
	 }
      }
   }  while ( ( ce=ce->getCCWEdg() ) != cln->getEdg() );
   return lowNeighbor;
}

/*****************************************************************************
**
**  FindLakeNodeOutlet
**
**  This function is part of the lake-filling algorithm. It checks to see 
**  whether there is a valid outlet for the current node, and if so it 
**  assigns that outlet. An "outlet" essentially means a downhill neighbor 
**  that isn't already flooded to the level of the current node. The function 
**  performs basically the same operation as FlowDirs, but with stricter 
**  criteria. The criteria for a valid outlet are:
**
**  (1) It must be lower than the current node (slope > 0)
**  (2) It must not be part of the current lake (a lake can't outlet to itself)
**  (3) It must not be a closed boundary (_flowAllowed_ must be true)
**  (4) If the outlet is itself part of a different lake, the water surface
**      elevation of that lake must be lower than the current node.
**
**  Returns: true if a valid outlet is found, FALSE otherwise
**  Calls: (none)
**  Called by: FillLakes
**  Created: 6/97 GT
**  Updated: 12/19/97 SL; 1/15/98 gt bug fix (open boundary condition)
**
*****************************************************************************/
int bFlowNet::FindLakeNodeOutlet( bFlowNode *node )
{
   double maxslp = 0; // Maximum slope found so far
   bEdge * ce;        // Current edge
   bFlowNode *dn, *an;   // 'dn' - Potential outlet
                      // 'an' - Node ptr used to find outlet of a previously
                      // identified lake
   
   // Check all node's neighbors
   ce = node->getEdg();
   do
   {
      // If it passes this test, it's a valid outlet
      dn = (bFlowNode *) ce->getDestinationPtrNC();
      //assert( dn>0 );

      // If the node checked IS Outlet - the outlet has been found
      if (dn->getBoundaryFlag() == kOpenBoundary) {
	  maxslp = 99999.0;
          node->setFlowEdg( ce );
      }
      else if ( ce->getSlope() > maxslp &&
		dn->getBoundaryFlag() != kClosedBoundary && 
	        dn->getFloodStatus()  != kCurrentLake )
	{
         // Handle a very special and rare case: if the "target" node dn is
         // part of a previous lake, it's still a valid exit as long as its
         // water surface elevation is lower than the current lake (whose
         // wse, assuming an outlet is found, would be equal to _node_'s
         // elevation). It can sometimes happen that the target lake's wse is
         // exactly equal in elevation to _node_, in which case
         // the point is not considered an outlet---if it were, infinite loops
         // could result
         if ( dn->getFloodStatus() == kFlooded )
         {
            // Iterate "downstream" through the "old" lake until reaching the
            // outlet, then test its elevation. If the elevation is exactly
            // equal to _node_, skip the rest and go on to the next iteration.
	    an = dn;
            while ( an->getFloodStatus()  != kNotFlooded &&
		    an->getBoundaryFlag() != kOpenBoundary )
                    an = an->getDownstrmNbr();
            if ( an->getZ() == node->getZ() && 
		 an->getBoundaryFlag() != kOpenBoundary) continue;
         }
         // Assign the new max slope and set the flow edge accordingly
         maxslp = ce->getSlope();
         node->setFlowEdg( ce );
      }
   }  while ( ( ce=ce->getCCWEdg() ) != node->getEdg() );
   
   return( maxslp > 0 );
}

/*****************************************************************************
**
**  SetBasinOutlet
**
**  Reset drainage area for all active nodes to zero
**  and find basin outlet node
**
*****************************************************************************/
void bFlowNet::SetBasinOutlet()
{
  bFlowNode *cn;
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter;

  // Finding basin outlet node and setting zero contributing area 

  for (nodIter = nodeList->begin(); nodIter != nodeList->end(); nodIter++) {
    cn = (*nodIter);
    if (cn->getBoundaryFlag() != kClosedBoundary)
      cn->setContrArea( 0.0 );

    if (cn->getBoundaryFlag() == kOpenBoundary) {
      OutletNode = cn; //Outlet does not have Voronoi area
      cout<<"\nBasin outlet node determined..." 
          << OutletNode->getID() <<endl<<flush;
    }
  }
  return;
}

/*****************************************************************************
**
**  DrainAreaVoronoi
**
**  Computes drainage area for each node by summing the Voronoi areas of all
**  nodes that drain to it, using the following algorithm:
**
**    FOR each active node
**      Cascade downstream, adding starting node's Voronoi area to each
**      ownstream node's drainage area, until an outlet or sink is reached
**
**    Note that each node's drainage area includes its own Voronoi area.
**
**    Modifies:  node contrbuting area
**
*****************************************************************************/
void bFlowNet::DrainAreaVoronoi()
{
  bFlowNode * cn;
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter;

  // Send voronoi area for each node to the node at 
  // the other end of the flowedge and downstream  
  cout.setf( ios::fixed, ios::floatfield);

  for (nodIter = nodeList->begin();
       nodIter != nodeList->getLastActive(); nodIter++) {
    cn = (*nodIter);
    RouteFlowArea( cn, cn->getVArea() );
  }
  
  if (OutletNode->getContrArea() < THRESH)
    cout<<"\nOutlet contributing area:\t"
	<<(int)(OutletNode->getContrArea())<<" m^2"<<endl;
  else 
    cout<<"\nOutlet contributing area:\t"
	<<(OutletNode->getContrArea())*THRESH<<" km^2"<<endl;
  return;
}

/*****************************************************************************
**
**  RouteFlowArea
**
**  Starting from the current node 'cn', this routine increments 
**  the drainage area of the node and each node downstream by _addedArea_
**
*****************************************************************************/
void bFlowNet::RouteFlowArea( bFlowNode *cn, double addedArea )
{
   int niterations=0;  // Safety feature: prevents endless loops

   // As long as the current node is neither a boundary nor a sink, add
   // _addedArea_ to its total drain. area and advance to the next downstream

   do {
     cn->addContrArea( addedArea );
     cn = cn->getDownstrmNbr();
     niterations++;
     assert( niterations <= gridPtr->getNodeList()->getActiveSize() );
   } while ( cn != OutletNode);

   cn->addContrArea( addedArea ); //Add it to the Outlet node as well
   
   return;
}

/****************************************************************************
**
**  DeriveCurvature
**
**  Derives curvature for each element in the basin
**
*****************************************************************************/
void bFlowNet::DeriveCurvature() 
{
   double curv;
   bFlowNode *cnn, *cn;
   bEdge  *firstedg; 
   bEdge  *curedg;
   bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
   std::_List_iterator<bFlowNode*> nodIter;

   for (nodIter = nodeList->begin();
        nodIter != nodeList->getLastActive(); nodIter++) {
     cn = (*nodIter);
     curv = 0.0;
     firstedg = cn->getFlowEdg();
     curedg = firstedg->getCCWEdg();
     while (curedg != firstedg) {
       cnn = (bFlowNode*)curedg->getDestinationPtrNC();
       if (cnn->getBoundaryFlag() != kClosedBoundary &&
	   cnn->getBoundaryFlag() != kOpenBoundary) {

	 if (cnn->getFlowEdg()->getDestinationPtrNC() == (bNode*)cn )
	   curv += cnn->getFlowEdg()->getSlope();
       }
       curedg = curedg->getCCWEdg();
     }
     curv -= firstedg->getSlope();
     cn->setCurvature(curv);
   }
   return;
}

/*****************************************************************************
**  
**  DeriveStreamReaches
**
**  Define stream reaches in the proper order according to their position
**  in the drainage network
**
*****************************************************************************/
void bFlowNet::DeriveStreamReaches(bInputFile &infile)
{
  int cnt, flag, ll; // niterations;
  bFlowNode *cn;
  bFlowNode *cmove, *cprev;

  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> niter;
  std::_List_iterator<int> NNodesIter = NNodes.begin();

  cerr.setf( ios::fixed, ios::floatfield);
  cout.setf( ios::fixed, ios::floatfield);

  // Activating tracer and composing a list of river heads 

  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    (*niter)->ActivateSortTracer(); // tracer is assigned to '1' for each

  // The stream nodes must be ORDERED - that's why  
  // routine 'SortStreamNodes()' is called before.  
  // If the routine work correctly, then we should be
  // able to extract stream reaches quite correctly  
  // in terms of their hierarchical position in the drainage network.  

   ll = 0;
   for (niter = nodeList->begin();
        niter != nodeList->getLastActive(); niter++) {
     cn = (*niter);

     // If the node is a stream node and it has not been
     // included in any stream reach yet (it also can be
     // an outlet for some reach), take it as a new origin 

     if (cn->getBoundaryFlag() == kStream && !(cn->NoMoreTracers())) {
       NodesLstH.push_back( cn ); // Channel Head Node
       cprev = cmove = cn;
       flag = 0;
       cnt = 1;

       // Go downstream until you reach a confluence/outlet 
       // node, which we define in the code
       while ( !flag ) {
	 cprev = cmove; // Always keeps track of the previous node
	 cmove = cmove->getDownstrmNbr(); // Downstream node...
	 cmove->DeactivateTracer();
	 flag = IsConfluence(cmove, cprev);
	 if (flag)
	   cmove->ActivateSortTracer(); // Re-set tracer for confluence
	 cnt++;
       }

       //cout<<"SID = "<<ll<<"\t# of NODES = "<<cnt<<endl;
       NodesLstO.push_back( cmove ); // Channel Outlet (basin O included)
       NNodes.push_back( cnt );      //# of nodes in a stream reach
       ll++;
     }
   }
   // Now we have stream reaches between the points of confluence 
   // and origin. The first list (NodesLstH) contains the origin node. 
   // It is either starts at the river head OR at the node of confluence
   // Therefore, this node is always unique. The second list (NodesLstO)
   // contains outlet points, which are either the points of confluence 
   // for different stream reaches or the Outlet. There are various 
   // reaches having the same outlet node. This feature will be used to 
   // extract the hierarchical order of river network.

   PrintArcInfoLinks( infile ); 

   cout<<"\nDerived stream reaches..."<<endl;
   return;
}

/*****************************************************************************
**  
**  PrintArcInfoLinks
**
**  Prints to a file coordinates of stream links derived in the 
**  'DeriveStreamReaches' routine in ArcInfo input format
**  Use 'generate' command to obtain the coverage
**  
*****************************************************************************/
void bFlowNet::PrintArcInfoLinks(bInputFile &infile) 
{
  int cnt; // flag;
  char fullName[kMaxNameSize+20];
  ofstream ControlOut;

  bFlowNode *cn;
  bFlowNode *cmove, *cprev;
  std::_List_iterator<bFlowNode*> NodesIterO;
  std::_List_iterator<bFlowNode*> NodesIterH;
  std::_List_iterator<int> NNodesIter = NNodes.begin();

  infile.ReadItem(fullName, "OUTFILENAME" ); // basename 
  strcat( fullName, "_reach");

  ControlOut.open(fullName);
  if ( !ControlOut.good() ) {
    cerr<<"Unable to open: "<< fullName
	<<" exiting..."<<endl<<flush;
    exit(2);
  }
  ControlOut.setf( ios::fixed, ios::floatfield);

  NodesIterO = NodesLstO.begin();
  for (NodesIterH = NodesLstH.begin(), 
       NodesIterO = NodesLstO.begin(), 
       NNodesIter = NNodes.begin(), cnt=1;
       NodesIterH != NodesLstH.end();
       NodesIterH++, NodesIterO++, NNodesIter++,  cnt++) {
    cn = (*NodesIterH);
    cprev = (*NodesIterO); // Point to the current Outlet
    cmove = cn;                  // Assign to the current one

    ControlOut<<cnt<<endl;  
    ControlOut<<cmove->getX()<<","<<cmove->getY()<<endl;

    while (cmove != cprev ) {
      cmove = cmove->getDownstrmNbr();
      ControlOut<<cmove->getX()<<","<<cmove->getY()<<endl;
    }
    ControlOut<<"END"<<endl;
  }
  ControlOut<<"END"<<endl;
  ControlOut.close();

  return;
}

/*****************************************************************************
**  
**  WeightedShortestPath()
**
*****************************************************************************/

#define STEDGWEIGHT 2.5
#define SETTLED    3
#define INSTACK    2
#define RECURS     25

void bFlowNet::WeightedShortestPath()
{
  int flag = 1;
  int niterations = 0;
  bFlowNode *cn;

  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> niter;
  std::_List_iterator<bFlowNode*> HeadsIter;

  std::list<bFlowNode*> NodesLst; 
  std::_List_iterator<bFlowNode*> NodesIter;
  std::list<bEdge*> EdgeLst; 
  std::_List_iterator<bEdge*> EdgeIter;
  
  // Set tracer to '1' and use the contributing area variable 
  // to find the shortest path (summing weights for edge length) 
 
  for (niter = nodeList->begin();
       niter != nodeList->getLastActive(); niter++) {
    cn = (*niter);
    cn->DeactivateTracer();
    cn->setContrArea(1.0E+9); // Inf at the beginning
  }

  // "Settle" the Outlet node
  OutletNode->setTracer(SETTLED);
  OutletNode->setContrArea(0.0);

  // Start adding stream nodes to the least
  AddUnsettledNeighbors(OutletNode, NodesLst, EdgeLst);

  // The first 'while' loop is intended not to miss  
  // any possible disjoints in stream network of the 
  // basin. There might be many and so the routine  
  // searches for them. The second 'while' is for  
  // correction of flow directions between disjoints

  while ( flag ) { 
     // Loop through the nodes of the list 
     // Update the list as we move up in the basin
     NodesIter = NodesLst.begin();
     while ( !NodesLst.empty() ) {

       // Take the current node
       cn = (*NodesIter);
       AddUnsettledNeighbors(cn, NodesLst, EdgeLst);
   
       // Now, remove the current node from the list  
       if (!NodesLst.empty())
          NodesLst.erase(NodesLst.begin());
       if (!EdgeLst.empty())
          EdgeLst.erase(EdgeLst.begin());

       // If list is not empty, take the next node
       if (!NodesLst.empty())
	 NodesIter = NodesLst.begin();

       assert( niterations <= gridPtr->getNodeList()->getActiveSize());
       niterations++;
     }
     NodesLst.erase(NodesLst.begin(), NodesLst.end());
     EdgeLst.erase(EdgeLst.begin(), EdgeLst.end());

     // Now, define the basin stream heads among "settled" nodes

     for (niter = nodeList->begin();
          niter != nodeList->getLastActive(); niter++) {
       cn = (*niter);
       if (cn->getBoundaryFlag() == kStream && cn->getTracer() == SETTLED) {
	 if ( IsStreamHead(cn) ) {
	   HeadsLst.push_back( cn );
	 }
       }
     }

     // Once we have some stream heads, let's see if we have missed 
     // some stream nodes due to the disjoints in the TIN 
     flag = 0;
     for (HeadsIter = HeadsLst.begin();
          HeadsIter != HeadsLst.end(); HeadsIter++) {
       cn = (*HeadsIter);
       //cout<<"\n\n # Checking stream heads:"<<endl;
       //TellAboutNode(cn);
       flag += FindStreamDisjoints(cn, 1, NodesLst);
     }

     // If disjoint nodes have been found Flush the stream head list, it would
     // be in wrong order otherwise 
     if (flag)
        HeadsLst.erase(HeadsLst.begin(), HeadsLst.end());

     // If disjoint nodes have NOT been found by the above procedure 
     // (i.e. searching at a defined stream head), they still may exist, 
     // so use another procedure to figure out if they indeed exist 
     else {
        flag = FindConfluenceDisjoints(NodesLst);
     }
  }
  
  // Find STREAM nodes that have not been passed by during the  
  // preceding operations and make them 'HILLSLOPE' nodes
  int num_streams_nodes = 0;
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++) {
        cn = (*niter);
        if (cn->getBoundaryFlag() == kStream) {
            num_streams_nodes++;
        }
  }


  int count = 0;
  for (niter = nodeList->begin();
       niter != nodeList->getLastActive(); niter++) {
    cn = (*niter);
    if (cn->getBoundaryFlag() == kStream && cn->getTracer() < INSTACK ) {
      count++;
      cn->setBoundaryFlag( kNonBoundary );
    }
  }
  cout << "CHANGED " << count << " stream to nonboundary" << endl;

  if(num_streams_nodes == count){
      throw runtime_error("All stream nodes were converted to interior nodes, no reaches were identified. Increase stream node density.");
  }

  // Due to some peculiarities, we need to re-check stream heads again
  for (HeadsIter = HeadsLst.begin();
       HeadsIter != HeadsLst.end(); HeadsIter++) {
    cn = (*HeadsIter);
    AddUnsettledNeighbors(cn, NodesLst, EdgeLst);
  }
  HeadsLst.erase(HeadsLst.begin(), HeadsLst.end());

  // Re-define the basin stream heads one more time 
  for (niter = nodeList->begin();
       niter != nodeList->getLastActive(); niter++) {
    cn = (*niter);
    if (cn->getBoundaryFlag() == kStream && cn->getTracer() == SETTLED) {
      if ( IsStreamHead(cn) ) {
	HeadsLst.push_back( cn );
      }
    }
  }

  // Eliminate tributaries that have only two nodes: a stream head
  // and an outlet. Make stream head as a hillslope node
  std::_List_iterator<bFlowNode*> nodeToMove;
  niterations = HeadsLst.size();
  HeadsIter = HeadsLst.begin();
  cn = (*HeadsIter);

  for ( flag=1; flag <= niterations; flag++ ) {
    if ( IsToEliminate( cn ) ) {
      cn->setBoundaryFlag( kNonBoundary );
      nodeToMove = HeadsIter;
      HeadsIter++;
      HeadsLst.erase(nodeToMove);
      cn = (*HeadsIter);
      HeadsLst.push_front(cn);
      HeadsLst.erase(HeadsLst.begin());
      cn = (*HeadsIter);
    }
    else {
      HeadsIter++;
      cn = (*HeadsIter);
    }
  }

  // Set contributing area variable to back to '0.0' 
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    (*niter)->setContrArea(0.0);

  return;
}

/*****************************************************************************
**  
**  AddUnsettledNeighbors()
**  
**  The function considers neighboring to 'cn' nodes and sets the shortest 
**  path. It first take an edge to a node that has been "settled" already 
**  and searches for the other "settled" nodes to set the flowedge. All
**  "unsettled" nodes are added to the general list of nodes that are to 
**  be checked by the calling function
**
*****************************************************************************/
void bFlowNet::AddUnsettledNeighbors(bFlowNode *cn,
			             std::list<bFlowNode*> &NodesLst,
			             std::list<bEdge*> &EdgeLst)
{
   int cnt = 0;
   double ttt = 0;
   bFlowNode *cnn;
   bEdge  *firstedg;   // pointer to first edge
   bEdge  *curedg;     // pointer to current edge

   ttt = cn->getContrArea();

   // 'FlowEdge' should point to a downstream node
   //  in some optimized manner --> after FillLake()
   if (cn->getBoundaryFlag() == kNonBoundary ||
       cn->getBoundaryFlag() == kStream)
     firstedg = cn->getFlowEdg();
   else
     firstedg = cn->getEdg();

   // First, find any downslope "settled" node
   if (cn != OutletNode) {
     cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
     while (cnn->getTracer() != SETTLED && cnt < 100) {
       firstedg = firstedg->getCCWEdg();
       cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
       cnt++;
     }
     if (cnt >= 100) { 
       cout<<">>> Error in bFlowNet:: NO Settled neighbors found! "
	   <<"Exiting..."<<endl;
       exit(2);
     }
   }
   // Start checking the nodes... 
   cnt = CheckNeighbor( cn, firstedg, NodesLst, EdgeLst );

   curedg = firstedg->getCCWEdg();
   while (curedg != firstedg) {
       cnt += CheckNeighbor( cn, curedg, NodesLst, EdgeLst );
       curedg = curedg->getCCWEdg();
   }

   // Make all the necessary checks
   if (cn->getBoundaryFlag() == kNonBoundary ||
       cn->getBoundaryFlag() == kStream) {
     if (cn->getContrArea() >= ttt && cn->getTracer() != SETTLED ) {
       cout<<">>> bFlowNet: AddUnsettledNeighbors: "
	   <<"The path variable has not changed! Exiting..."<<endl;
       exit(2);
     }
     // Once the flow direction is chosen, mark this node as "settled"
     else
       cn->setTracer(SETTLED);
   }
   return;
}

/*****************************************************************************
**  
**  CheckNeighbor()
**  
**  The function takes as arguments a ptr to a current node 'cn', an edge
**  'curedge' that originates at 'cn' and a node list NodesLst. The routine
**  analyzes the node located on another end of 'curedg'. Depending to what 
**  its tracer equals to, it:
**      - Does nothing with it (tracer == INSTACK);
**      - Checks if 'cnn' is suitable for flowing into (tracer == SETTLED);  
**      - Adds the node to the list of nodes to be analyzed ((tracer < INSTACK);
**
*****************************************************************************/
int bFlowNet::CheckNeighbor(bFlowNode *cn, bEdge *curedg,
			    std::list<bFlowNode*> &NodesLst,
			    std::list<bEdge*> &EdgeLst)
{
  int cnt = 0;
  double tempo;
  double edgeWeight;
  bFlowNode *cnn;

  cnn = (bFlowNode*)curedg->getDestinationPtrNC();

  // Check STREAM nodes
  if (cnn->getBoundaryFlag() == kStream || 
      cnn->getBoundaryFlag() == kOpenBoundary) {

    // If this is a "Settled" node --> compare its weight
    // plus the weight of the connecting edge to the 
    // current weight of the node
    if (cnn->getTracer() == SETTLED) { 
      tempo = cnn->getContrArea();
      edgeWeight = pow((curedg->getLength() / 100.0), STEDGWEIGHT);
      tempo += edgeWeight;

      // If we change the flow direction for 'cn'
      // we need to check and update all its neighbors
      // if 'cn' is a better choice for them to flow to 
      // --> call 'UpdatePathVariable'
      if ( tempo < cn->getContrArea() ) {
	cn->setContrArea( tempo );
	cn->setFlowEdg( curedg );
	UpdatePathVariable( cn );
      }

      // If current weight assigned to 'cn' is 
      // less than the one assigned to 'cnn'
      // ---> change the flow direction for 'cnn'
      else {
	tempo = cn->getContrArea();
        edgeWeight = pow((curedg->FindComplement()->getLength() / 100.0), 
                          STEDGWEIGHT);
	tempo += edgeWeight;

	if ( tempo < cnn->getContrArea() ) {
	  cnn->setContrArea( tempo );
	  cnn->setFlowEdg( curedg->FindComplement() );
	  UpdatePathVariable( cnn );
	}
      }
    }
    // If the stream node is "Unsettled" --> put it in the stack 
    else if (cnn->getTracer() < INSTACK) {
      NodesLst.push_back( cnn );
      EdgeLst.push_back( curedg );
      cnn->setTracer(INSTACK);
      cnt++;
    }
  }
  // For the break points we also need to check HILLSLOPE nodes 
  else if (cnn->getBoundaryFlag() == kNonBoundary ) {
    ;
  }
  return cnt;
}

/*****************************************************************************
**  
**  IsStreamHead()
**  
**  Checks if node 'cn' is a stream head. For that it checks all neighbors of
**  'cn', if there is any that has a flow edge pointing to 'cn' --> 'cn' is 
**  NOT a stream head (return '0'). Otherwise, it is a stream head (return '1')
**  
*****************************************************************************/
int bFlowNet::IsStreamHead(bFlowNode *cn)
{
  int cnt = 0;
  bFlowNode *cnn;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge

  if (cn->getBoundaryFlag() != kStream) 
    return cnt;

  // The stream node which is the destination node of current edge
  // may have FLOWEDGE pointed to the 'cn' node (the node being tested)

  firstedg = cn->getFlowEdg();
  cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
  if (cnn->getBoundaryFlag() == kStream ) {
    if (cnn->getFlowEdg()->getDestinationPtrNC() == (bNode*)cn )
      cnt++;
  }

  curedg = firstedg->getCCWEdg();
  while (curedg != firstedg) {
    cnn = (bFlowNode*)curedg->getDestinationPtrNC();
    if (cnn->getBoundaryFlag() == kStream ) {
      if (cnn->getFlowEdg()->getDestinationPtrNC() == (bNode*)cn )
	cnt++;
    }
    curedg = curedg->getCCWEdg();
  }
  if (cnt > 0)
    cnt = 0;
  else if (cnt == 0)
    cnt = 1;

  return cnt;
}

/*****************************************************************************
**  
**  IsStreamHead()
**  
**  The function checks if the tributary starting at the stream head node 
**  'chead' has to be eliminated (the length of the tributary is 2 nodes:
**  the head 'chead' and the outlet which is the confluence node).  
**  
*****************************************************************************/
int bFlowNet::IsToEliminate(bFlowNode *chead)
{
  int cnt = 0;
  bFlowNode *cn, *cnn;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge

  if (chead->getBoundaryFlag() != kStream) 
    return 0;

  // Get the downstream node of a stream head node
  cn = chead->getDownstrmNbr();
   
  // The stream node which is the destination node of an edge
  // may have FLOWEDGE pointed to the 'cn' node (the node being tested)

  if ( cn != OutletNode ) {
    firstedg = cn->getFlowEdg();
    cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
    if (cnn->getBoundaryFlag() == kStream && cnn != chead) {
      if (cnn->getFlowEdg()->getDestinationPtrNC() == (bNode*)cn )
	cnt++;
    }

    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg) {
      cnn = (bFlowNode*)curedg->getDestinationPtrNC();
      if (cnn->getBoundaryFlag() == kStream && cnn != chead) {
	if (cnn->getFlowEdg()->getDestinationPtrNC() == (bNode*)cn )
	  cnt++;
      }
      curedg = curedg->getCCWEdg();
    }
  }
  return cnt;
}

/*****************************************************************************
**  
**  FindConfluenceDisjoints()
**  
**  Disjoints may exist at the confluence points so here we use another
**  procedure to figure out if they indeed exist. Return '0' if there are 
**  no disjoint tributaries at confluences, OR if they can not be found
**  by the existing algorithm (the like reason is the wrong finding of the 
**  node for search of the confluence)
**  
*****************************************************************************/
int bFlowNet::FindConfluenceDisjoints(std::list<bFlowNode*> &NodesLst)
{
  int cnt = 1;
  int flag = 0;
  int niterations = 0;
  bFlowNode *cn, *cmove;

  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> niter = nodeList->begin();

  // Loop through the nodelist until all the possibilities are checked, i.e.
  // loop until you find stream head that may lead us to a confluence, if it
  // does - get out and work with this tributary, if not - continue searching 

  while (cnt > 0 && flag == 0) {

    cnt  = 0;
    flag = 0;
    
    cn = (*niter);
    while ( niter != nodeList->getLastActive() && cnt == 0 ) {
      // Consider only nodes with (tracer == 0)
      if (cn->getBoundaryFlag() == kStream && cn->NoMoreTracers()) {
	if ( IsStreamHead(cn) ) {
	  cnt++;
	}
      }
      cmove = cn;
      niter++;
      cn = (*niter);
    }

    // If an unsettled stream head has been found --> we need 
    // to find where this tributary joins the main network of 
    // settled nodes. For that, go down stream of 'cmove' until   
    // a settled node is reached. Figure out how to join them. 

    if ( cnt ) {
      niterations = 0;
      do { 
	cn = cmove;
	cmove->setTracer(-1); // <-- Assign tracer to '-1'
	cmove = cmove->getDownstrmNbr();

	// Get out of the loop if you have found a settled node 
	// OR have reached the basin outlet 
	if ((cmove->getTracer() == SETTLED && 
	     cmove->getBoundaryFlag() == kStream) || cmove == OutletNode)
	  niterations = -9999;
        niterations++;
        assert( niterations < gridPtr->getNodeList()->getActiveSize() );
      } while ( niterations > 0 );

      // We now have found the outlet node for the tributary. In the 
      // vicinity of this node we should look for the node that is a 
      // true tributary outlet (the first one may be fictitious outlet) 
      // BUT: we will assume that 'cmove' IS the true outlet!
      // It is a strong assumption but in most cases would be true. 
      // If it is not - search for the true outlet among 'cmove' neighbors 
      // (most likely it would be upstream of 'cmove')

      // NOTE: Please note, that if no stream disjoints have been 
      // found in the vicinity of 'cmove', the algorithm will return 
      // (flag = 0) and therefore the tributary will not be ever found

      flag = FindStreamDisjoints(cmove, 1, NodesLst);
    }
  }
  return flag;
}

/*****************************************************************************
**  
**  FindStreamDisjoints()
**  
**  The function is intended to define if a current stream head node 'cn' is
**  actually a disjoint node in stream network.  For that, the function looks
**  at the hillslope neighbors of 'cn': if some of these nodes have connection
**  to a stream node that has "unsettled" value of tracer (< INSTACK) then the 
**  node 'cn' is a disjoint node. Correspondingly, an optimum path to the found
**  "unsettled" stream node is defined through a hillslope node (which becomes
**  'stream') and it is stored in the general stack of nodes 'NodesLst'.
**  
**  NOTE: - So far, it is assumed that there might be only ONE hillslope node 
**  that disjoints stream network. The function might be designed for recursive
**  calls to look for a stream node through several disjoining hillslope nodes.
**  Variable 'times' is intended for that.
** 
*****************************************************************************/
int bFlowNet::FindStreamDisjoints(bFlowNode *cn, int times, 
				  std::list<bFlowNode*> &NodesLst)
{
  int cnt = 0;
  double tempo = 0;
  double min = 1.0E+9;
  double edgeWeight1, edgeWeight2;
  bFlowNode *cnn;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge

  std::list<bEdge*> EdgeLst1; 
  std::_List_iterator<bEdge*> EdgeIter1;
  std::list<bEdge*> EdgeLst2; 
  std::_List_iterator<bEdge*> EdgeIter2;

  if (times > 0) {
    times--;
    if ((cn->getBoundaryFlag() == kStream &&
	 cn->getTracer() == SETTLED) || cn == OutletNode) {

      if (cn == OutletNode)
	firstedg = cn->getEdg();
      else
	firstedg = cn->getFlowEdg();
      cnn = (bFlowNode*)firstedg->getDestinationPtrNC();

      // Check the node
      cnt += IsStreamDisjoint(firstedg, cnn, EdgeLst1, EdgeLst2);

      curedg = firstedg->getCCWEdg();
      while (curedg != firstedg) {
	cnn = (bFlowNode*)curedg->getDestinationPtrNC();
	cnt += IsStreamDisjoint(curedg, cnn, EdgeLst1, EdgeLst2);
	curedg = curedg->getCCWEdg();
      }

      // Now, if (cnt > 0), then some nodes have been found
      if (cnt > 0) {
        for (EdgeIter1 = EdgeLst1.begin(), EdgeIter2 = EdgeLst2.begin();
             EdgeIter1 != EdgeLst1.end(); EdgeIter1++, EdgeIter2++) {

	  // Compute path weight
          edgeWeight1 = pow(((*EdgeIter1)->getLength() / 100.0), 1);
          edgeWeight2 = pow(((*EdgeIter2)->getLength() / 100.0), 1);
	  tempo  = edgeWeight1 + edgeWeight2;

	  // Choose, if it provides min path
	  // Store the value and the edges
	  if (tempo < min) {
	    min = tempo;
	    firstedg = (*EdgeIter1);
	    curedg   = (*EdgeIter2);
	  }
	}
	// By now, the min path has been chosen
	// Do appropriate assignments:

	cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
	// 1.) Change its boundary flag
	cnn->setBoundaryFlag( kStream );

	// 2.) Assign the flow edge - complimentary to 'firstedg'
	cnn->setFlowEdg( firstedg->FindComplement() );

	// 3.) Set tracer to "settled"
	cnn->setTracer(SETTLED);

	// 4.) Assign the weight
        edgeWeight1 = pow((firstedg->getLength() / 100.0), STEDGWEIGHT);
	tempo = cn->getContrArea() + edgeWeight1;
	cnn->setContrArea( tempo );

	// 5.) Put the found stream node in stack
	cnn = (bFlowNode*)curedg->getDestinationPtrNC();
	NodesLst.push_back( cnn );

	// 6.) Clean up temporary stacks
	EdgeLst1.erase(EdgeLst1.begin(), EdgeLst1.end());
	EdgeLst2.erase(EdgeLst2.begin(), EdgeLst2.end());
      }
    }
  }
  return cnt;
}

/*****************************************************************************
**  
**  IsStreamDisjoint()
**  
**  From a current stream head (in the calling function), analyze a hillslope
**  node 'cn', a destination node of edge 'edg'. If 'cn' has connection 
**  to "unsettled" stream nodes, put two edges in the stacks, return cnt > 0
**  
*****************************************************************************/
int bFlowNet::IsStreamDisjoint(bEdge *edg, bFlowNode *cn, 
			       std::list<bEdge*> &EdgeLst1,
			       std::list<bEdge*> &EdgeLst2)
{
  int cnt = 0;
  int ntimes;
  double tempo;
  bFlowNode *cnn;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge

  // If 'cn' is not a hillslope node, we do not consider it at all
  if (cn->getBoundaryFlag() == kNonBoundary) {

    // Now, consider edge by edge if there is 
    // any stream node that has not been assigned 
    // the flow direction in the weighted shortest path algorithm 
    firstedg = cn->getEdg();
    cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
    if (cnn->getTracer() < INSTACK && cnn->getBoundaryFlag() == kStream) {

      // Now, check if this node has positive elevation tendency,
      // i.e. minimum elevation of stream nodes increases
      ntimes = 1;
      tempo = CompElevationTendency(cn, cnn, &ntimes, RECURS);
      tempo /= ntimes;
         
      // IF the tendency is higher than the node elevation - fine
      // ELSE - we've found a node on another tributary 
      if (tempo > cnn->getZ()) {
        EdgeLst1.push_back( edg );        // to Hillslope node
	EdgeLst2.push_back( firstedg );   // to Stream node 
	cnt++;
      }
    }

    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg) {
      cnn = (bFlowNode*)curedg->getDestinationPtrNC();
      if (cnn->getTracer() < INSTACK && cnn->getBoundaryFlag() == kStream) {
        ntimes = 1;
	tempo = CompElevationTendency(cn, cnn, &ntimes, RECURS);
        tempo /= ntimes;
            
	if (tempo > cnn->getZ()) {
	  EdgeLst1.push_back( edg );
	  EdgeLst2.push_back( curedg );
	  cnt++;
	}
      }
      curedg = curedg->getCCWEdg();
    }
  }
  return cnt;
}

/*****************************************************************************
**  
**  IsConfluence
**
**  Checks if stream node 'cn' is a confluence node. Node 'cup' defines the 
**  upstream node for 'cn' of a stream link for which junction node is sought 
**
**  - Returns '1' if 'cn' is a confluence or outlet
**  - Returns '0' otherwise
**
*****************************************************************************/
int bFlowNet::IsConfluence(bFlowNode *cn, bFlowNode *cup) 
{
  int flag = 0;
  int niterations = 0;
  bFlowNode *cnn, *citer;
  bEdge  *firstedg; 
  bEdge  *curedg;

  // This flowedge has been checked before
  // It points to a downstream STREAM node
  if ( cn->getBoundaryFlag() == kOpenBoundary ) //Outlet!
    flag = 1;  // Outlet Node -- Exit!
  else {
    firstedg = cn->getFlowEdg();
    curedg = firstedg->getCCWEdg();
    while ((curedg != firstedg) && (!flag)) {

      // CHECK ## 1 : First, we look only at stream nodes, != upper node
      cnn = (bFlowNode*)curedg->getDestinationPtrNC();
      if (cnn != cup &&
	  cnn->getBoundaryFlag() == kStream) {

	// CHECK ## 2 : Slope condition may not work (ideally, 
	// upstream node from the tributary should be upslope); 
	// Check the contributing area condition
	if (cn->getContrArea() > cnn->getContrArea()) {

	  // CHECK ## 3 : The stream node which is the destination node 
	  // of 'curedg' must have FLOWEDGE pointed to 'cn' node (test node)
	  if (cnn->getFlowEdg()->getDestinationPtrNC() == (bNode*)cn ) {

	    // CHECK ## 4 : Last check. We have to make sure that 'cnn' 
	    // is not one of the downstream nodes which streamflow from 
	    // 'cn' would pass by on its way to Outlet (this would 
	    // actually duplicate Check ## 2)
	    flag = 1;
	    citer = cn;
	    do { 
	      citer = citer->getDownstrmNbr();
	      if (citer == cnn)
		flag = 0;
	      niterations++;
	      assert( niterations < gridPtr->getNodeList()->getActiveSize() );
	    } while ( niterations < 30 && citer != OutletNode && flag == 1 );
	  }
	}
      }
      curedg = curedg->getCCWEdg();
    }
  }
  return flag;
}

/*****************************************************************************
**  
**  CompElevationTendency()
**  
**  Computes average elevation of neighboring to 'cn_check' nodes which
**  are not shared with the 'cn' node ('cn' is usually upstream of 'cn_check'
**  and we therefore try to figure out if 'cn_check' leads us downslope)
**
*****************************************************************************/
double bFlowNet::CompElevationTendency(bFlowNode *cn, bFlowNode *cn_check, 
				       int *cnt, int flag)
{
  double Elev = 1.0E+6;
  double tempo = 0;
  bFlowNode *cnn, *cmm;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge
  std::list<bFlowNode*> DnLst; 
  std::_List_iterator<bFlowNode*> DnIter;

  std::list<int> Tracer;   // Tracer value in a stream node
  std::_List_iterator<int> TracerIter = Tracer.begin();

  if (flag > 0) { 
    flag--;

    firstedg = cn_check->getEdg();
    cnn = (bFlowNode*)firstedg->getDestinationPtrNC();

    // If 'cnn' is not shared with 'cn' 
    // use its elevation for computation
    if (cnn->getBoundaryFlag() == kOpenBoundary) {
       Elev = flag-10;
       flag = 0;
    }

    // 1.) If it is not an upstream node
    // 2.) If it is not connected to the upstream node
    // 3.) If it is a _Stream_ node
    // 4.) We have not previously stored it in the stack
    if (cnn != cn && !IsConnected(cn, cnn) && 
        cnn->getBoundaryFlag() == kStream && 
	cnn->getTracer() < 1) {
        Elev = cnn->getZ();
        DnLst.push_back( cnn );
        Tracer.push_back( cnn->getTracer() ); 
	cnn->setTracer(1);  // Set tracer to a "stored" value
    }
    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg) {
      cnn = (bFlowNode*)curedg->getDestinationPtrNC();
      // If it's an Outlet - make it the lowest
      if (cnn->getBoundaryFlag() == kOpenBoundary) {
	Elev = flag-10;
	flag = 0;
      }
      // If 'cnn' is not connected to the upper node 'cn'
      // --> check its elevation and add to the list
      // of nodes that are to be checked
      if (cnn != cn && !IsConnected(cn, cnn) && 
	  cnn->getBoundaryFlag() == kStream &&
	  cnn->getTracer() < 1) {
	if (cnn->getZ() < Elev)
	  Elev = cnn->getZ();
	DnLst.push_back( cnn );
	Tracer.push_back( cnn->getTracer() ); 

	// Set tracer to a "stored" value
	cnn->setTracer(1);
      }
      curedg = curedg->getCCWEdg();
    }

    // Now, for each of the nodes and compute its tendency
    // Estimate the average minimum tendency and add it to Elev

    double minElv = 0;
   
    if (DnLst.size() > 0) {
      for (DnIter = DnLst.begin(); DnIter != DnLst.end(); DnIter++) {
        cmm = (*DnIter);
	tempo = CompElevationTendency(cn_check, cmm, cnt, flag);
	if (tempo != 1.0E+6) {
	  minElv += tempo;  // if (tempo < minElv)
	  *cnt = *cnt + 1;
	}
      }
      // Use "average minimum" value 
      if ( minElv > 0 && minElv < 1.0E+6 ) {
	// 1.) minElv /= cnt;  // Elev = (Elev + minElv)/2.;
	// 2.) Elev = (Elev + minElv)/(cnt+1);
	Elev += minElv;
      }
      // Re-set tracer values to what they were...
      for (DnIter = DnLst.begin(), TracerIter = Tracer.begin(); 
	   DnIter != DnLst.end();
	   DnIter++, TracerIter++)
	(*DnIter)->setTracer( (*TracerIter) );
    }
    else
      flag = 0;

    // Release memory
    DnLst.erase(DnLst.begin(), DnLst.end());
    Tracer.erase(Tracer.begin(), Tracer.end());
  }
  return Elev;
}

/*****************************************************************************
**  
**  IsConnected()
**  
**  The routine checks if the node 'cn' has a common edge with 'cn_check'
**  Returns '1' if so, '0' otherwise
**
*****************************************************************************/
int bFlowNet::IsConnected(bFlowNode *cn, bFlowNode *cn_check) 
{
  int    flaggss = 0;
  bFlowNode *cnn;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge

  firstedg = cn->getFlowEdg();
  cnn = (bFlowNode*)firstedg->getDestinationPtrNC();
  if (cnn == cn_check)
    return 1;
  else {
    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg && flaggss == 0) {
      cnn = (bFlowNode*)curedg->getDestinationPtrNC();
      if (cnn == cn_check)
	flaggss = 1;    
      curedg = curedg->getCCWEdg();
    }
  }
  return flaggss;
}

/*****************************************************************************
**  
**  FindAngle()
**  
**  The routine finds angle for the nodes 'cn', 'cn1', and 'cn2'
**  It is assumed that they form vectors:
**           cn ---------- cn1
**             \ ) <--  alpha
**              \
**               \
**                \
**                cn2
**  Returned value is in degrees
**
*****************************************************************************/
double bFlowNet::FindAngle(bFlowNode *cn, bFlowNode *cn1, bFlowNode *cn2)
{
   double d1, d2, x1, y1, x2, y2, alpha;
   std::vector<double> xy(2), xy1(2), xy2(2);

   xy2 = cn2->get2DCoords();
   xy1 = cn1->get2DCoords();
   xy = cn->get2DCoords();

   x1  = xy1[0]-xy[0]; //<= THESE ARE COORDINATES
   y1  = xy1[1]-xy[1]; //<= OF THE VECTORS
   x2  = xy2[0]-xy[0]; //<= REQUIRED FOR FURTHER
   y2  = xy2[1]-xy[1]; //<= CALCULATIONS

   d1 = FindDistance(xy[0], xy[1], xy1[0], xy1[1]);
   d2 = FindDistance(xy[0], xy[1], xy2[0], xy2[1]);

   // Compute the angle between vectors cn-cn1 and cn-cn2
   // Angle is in degrees
   alpha = acos((x1*x2 + y1*y2)/(d1*d2))*180/(4*atan(1.0));  //Added decimal

   return alpha;
}

/***************************************************************************
**
**  FindDistance( )
**
**  Finds distance between two nodes located at (x1,y1) and (x2,y2)
**
***************************************************************************/
double bFlowNet::FindDistance(double x1, double y1, double x2, double y2) 
{
  return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
}

/***************************************************************************
**
**  TellAboutNode()
**
**  Prints node info
**
***************************************************************************/
void bFlowNet::TellAboutNode(bFlowNode *cn)
{
    if (cn->getBoundaryFlag() == kNonBoundary ||
	cn->getBoundaryFlag() == kStream) {
      cout<<"ID: " << cn->getID()
	  <<"\tPoint ("<<cn->getX()
	  <<"\t"<<cn->getY()
	  <<"\t"<<cn->getZ() << ")" 
          << "\tBoundary: " << cn->getBoundaryFlag()<< endl;
/*
      if (cn->getFlowEdg() != 0) {
	  cout <<"EdgeLen:\t"<<cn->getFlowEdg()->getLength()
	  <<"EdgeSlope:\t"<<cn->getFlowEdg()->getSlope() << endl;
          bFlowNode* o = (bFlowNode*) cn->getFlowEdg()->getOriginPtrNC();
          bFlowNode* d = (bFlowNode*) cn->getFlowEdg()->getDestinationPtrNC();
          if (o == 0) 
              cout << "Missing flow edge origin pointer" << endl;
          else
              cout << "Flow edge origin " << o->getID() << endl;
          if (d == 0) 
              cout << "Missing flow edge dest pointer" << endl;
          else
              cout << "Flow edge dest " << d->getID() << endl;
      
      } else {
          cout << "Missing flow edge" << endl;
      }
	  cout <<"\t"<<cn->getBoundaryFlag()<<endl<<flush;
*/
    } else 
      cout<<cn->getID()
	  <<"\t"<<cn->getX()
	  <<"\t"<<cn->getY()
	  <<"\t"<<cn->getZ()
          <<"\t-999\t-999"
	  <<"\t"<<cn->getBoundaryFlag()<<endl<<flush; 
  return;
}

/*****************************************************************************
**  
**  SortStreamNodes
**
**  Similar to SortNodesByStreamOrder but deals only with stream nodes.
**  Reason:  for some cases we have to re-adjust flow directions for stream
**           nodes and therefore it is necessary to ensure that the order 
**           of computations is correct
**  Algorithm: - Send all nodes to the back of the list such all hillslope 
**               nodes are preceding
**             - Sort stream nodes according to the relationship "drains to"
**             - Re-enumerate indices of the nodes
** 
*****************************************************************************/
void bFlowNet::SortStreamNodes() 
{
  int i;
  int nThisPass, nPassed; // Number moved in current iteration & in total
  int nStreams = 0;
  int done = false;

  bFlowNode * cn;
  bMeshList<bFlowNode*> *nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> listIter = nodeList->begin();
  std::_List_iterator<bFlowNode*> nodeToMove;

  int nUnsortedNodes = nodeList->getActiveSize(); // Number not yet sorted

  // Down stream node list for nodes extracted from mesh node list
  std::list<bFlowNode*> streamList;
  std::list<bFlowNode*> tempList;

  // First, move all the downstream nodes to the temporary list
  // following FIFO principle. Assign initial tracers
  cn = (*listIter);
  for ( i = 0; i < nUnsortedNodes; i++ ) {

    // If stream, move to temporary list
    if (cn->getBoundaryFlag() == kStream) {
      nStreams++;
      nodeToMove = listIter;
      listIter++;
      cn = (*listIter);

      // Move to temporary list
      streamList.push_back((*nodeToMove));
 
      // Remove from current node list
      nodeList->removeActive( nodeToMove );
    }
    else {
      listIter++;
      cn = (*listIter);
    }
  }

  // Activate tracers for the stream nodes
  std::list<bFlowNode*>::iterator streamIter;
  for (streamIter = streamList.begin(); 
       streamIter != streamList.end(); streamIter++)
        (*streamIter)->ActivateSortTracer();

  // Sort the temporary stream list using the water drop method
  // Move tracers downstream and sort until no nodes with tracers are left
  nPassed = 0;
  do {
    // --- Send tracers downstream ---
    streamIter = streamList.begin();
    cn = (*streamIter);
    for ( i = 0; i < nStreams; i++ ) {
        assert( cn!=0 );
	cn->MoveSortTracerDownstream();
        streamIter++;
        cn = (*streamIter);
    }

    // Scan for any nodes that have no tracers, and move them 
    // to the bottom of the stream list. This simply means that there 
    // were no any nodes for which they would lie downstream! 
    nThisPass = 0;
    done = true;
    streamIter = streamList.begin();
    cn = (*streamIter);

    for ( i = 0; i < nStreams; i++ ) {
      if (cn->NoMoreTracers() ) {

	// --- If no tracers, move to bottom of list ---
	nodeToMove = streamIter;
        streamIter++;
	cn = (*streamIter);

        // Add the node without the tracer to temporary list
        tempList.push_back((*nodeToMove));

        // Remove it from the current node list
	streamList.erase( nodeToMove );
	nThisPass++;
      }
      else {
        streamIter++;
	cn = (*streamIter);
        done = false;
      }
    }
    nStreams -= nThisPass;

    // Reset all active node tracers to value 1
    for (streamIter = streamList.begin();
         streamIter != streamList.end(); streamIter++) {
         (*streamIter)->ActivateSortTracer();
    }

  } while ( (done == false) );

  // Move temporary flow list in order to the active part of the node list
  // which should have only inactive nodes at this point
  std::list<bFlowNode*>::iterator tempIter;
  for (tempIter = tempList.begin(); tempIter != tempList.end(); tempIter++)
    nodeList->insertAtActiveBack((*tempIter));

  // Changed to make all IDs consistent
  bMeshList<bEdge*>* edgeList = gridPtr->getEdgeList();
  std::list<bTriangle*>* triList = gridPtr->getTriList();
                                                                                
  std::_List_iterator<bFlowNode*> niter;
  std::_List_iterator<bEdge*> eiter;
  std::_List_iterator<bTriangle*> titer;
                                                                                
  cout <<"\nRenumbering nodes, edges, triangles..."<<endl;
  int id;
  for (niter = nodeList->begin(), id = 0; 
       niter != nodeList->end(); niter++, id++)
    (*niter)->setID(id);
  for (eiter = edgeList->begin(), id = 0; 
       eiter != edgeList->end(); eiter++, id++)
    (*eiter)->setID(id);
  for (titer = triList->begin(), id = 0; 
       titer != triList->end(); titer++, id++)
    (*titer)->setID(id);
  return;
}

/*****************************************************************************
**  
**  UpdatePathVariable()
**  
**  The function updates the value of path variable for any neighbors of
**  'cn' that are "settled" and has flowedge to 'cn'
**  
*****************************************************************************/
void bFlowNet::UpdatePathVariable(bFlowNode *cn)
{
  double tempo, edgeWeight;
  bFlowNode *cnn;
  bEdge  *firstedg;   // pointer to first edge
  bEdge  *curedg;     // pointer to current edge

  if (cn->getBoundaryFlag() != kStream) 
    return;

  // If for any of the _settled_ neighbors 'cn'
  // is a better choice to flow to ->
  // then we need to update its weight 
  firstedg = cn->getFlowEdg();
  curedg = firstedg->getCCWEdg();
  while (curedg != firstedg) {
    cnn = (bFlowNode*)curedg->getDestinationPtrNC();

    if (cnn->getBoundaryFlag() == kStream && cnn->getTracer() == SETTLED) {
      edgeWeight = pow((curedg->FindComplement()->getLength() / 100.0), 
                        STEDGWEIGHT);
      tempo = cn->getContrArea() + edgeWeight;
      if ( tempo < cnn->getContrArea() ) {
	cnn->setContrArea( tempo );
	cnn->setFlowEdg( curedg->FindComplement() );
	UpdatePathVariable( cnn );
      }
    }
    curedg = curedg->getCCWEdg();
  }
  return;
}
       
#undef STEDGWEIGHT
#undef SETTLED
#undef INSTACK
#undef RECURS

/*************************************************************************
**
**  CheckVDrainageWidths
**
**  The function loops through the active node list in search of the 
**  nodes that have zero (or close to zero) voronoi edge length in the
**  steepest drainage direction. It calls a function that computes a 
**  pseudo width instead defined as the area of voronoi sector
**  located between the CW & CCW neighbors of the flow edge divided by 
**  the flow edge length
**
*************************************************************************/
void bFlowNet::CheckVDrainageWidths() 
{
  bFlowNode *cn;
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> nodIter;

  cout<<"bFlowNet: Checking Voronoi drainage widths..."<<endl;

  for (nodIter = nodeList->begin();
       nodIter != nodeList->getLastActive(); nodIter++) {
      cn = (*nodIter);
      FixVoronoiEdgeWidth( cn );
  }

  return;
}

/*************************************************************************
**
**  FixVoronoiEdgeWidth
**
**  The function computes a pseudo width for a node 'cn'. The pseudo width
**  is defined as the area of voronoi sector located between the CW & CCW 
**  neighbors of the flow edge divided by the flow edge length
**
*************************************************************************/
void bFlowNet::FixVoronoiEdgeWidth(bFlowNode *cn) 
{
  int i, j, flag, flagg, NV;
  int changeWidth;
  double sArea, wWidth;
  double **polyg;

  bEdge *flowedg, *ccwedg, *cwedg;
  bEdge  *firstedg, *curedg;
  std::vector<double> centroid(2);
  std::vector<double> vv_flow(2);
  std::vector<double> vv_ccw1(2), vv_ccw3(2), vv_cw2(2);
  std::vector<double> orig_ccw(2), dest_ccw(2);
  std::vector<double> orig_cw(2), dest_cw(2);
  std::vector<double> orig_flow(2), dest_flow(2);
  std::vector<double> intersect_ccw(2), intersect_cw(2);
  std::vector<double> vv_cur(2);

  cout<<setprecision(6);

  flag = flagg = 0;
  wWidth = sArea = 0;
  changeWidth = 0;

  NV = 5;
  polyg = new double*[2];
  for (int i = 0; i < 2; i++)
    polyg[i] = new double[NV];

  centroid = cn->get2DCoords();

  // Flow edge information
  flowedg = cn->getFlowEdg();
  vv_flow = flowedg->getRVtx();
  orig_flow = flowedg->getOriginPtrNC()->get2DCoords();
  dest_flow = flowedg->getDestinationPtrNC()->get2DCoords();

  // CCW edge information
  ccwedg = flowedg->getCCWEdg();		// CCW neighbor of flow edge
  orig_ccw = ccwedg->getOriginPtrNC()->get2DCoords();
  dest_ccw = ccwedg->getDestinationPtrNC()->get2DCoords();

  // Find the CW neighbor of flow edge by iterating around
  // Save voronoi vertices around the flow edge for later
  curedg = flowedg->getCCWEdg();
  while ( curedg != flowedg )  {

    vv_cur = curedg->getRVtx();		// Voronoi vertex of current edge

    if (vv_cur[0] != vv_flow[0] || vv_cur[1] != vv_flow[1]) {
      // Eventually this will be first CW voronoi vertex
      vv_cw2[0] = vv_cur[0]; 			// Current edge voronoi vertex
      vv_cw2[1] = vv_cur[1]; 

      // First time through loop
      if (flag == 0) {
        // Save first distinct CCW voronoi vertex
        vv_ccw1[0] = vv_cur[0]; 
        vv_ccw1[1] = vv_cur[1]; 
        flag++;

      // Second time through loop
      } else if (flagg == 0) {
        // Save the second distinct CCW voronoi vertex
        if (vv_cur[0] != vv_ccw1[0] || vv_cur[1] != vv_ccw1[1]) {
          vv_ccw3[0] = vv_cur[0];
          vv_ccw3[1] = vv_cur[1]; 
	  flagg++;
	}
      }
    }
    // When loop exits this will be the CW edge
    cwedg = curedg;
    curedg = curedg->getCCWEdg();
  }

  // CW edge information
  orig_cw = cwedg->getOriginPtrNC()->get2DCoords();
  dest_cw = cwedg->getDestinationPtrNC()->get2DCoords();

  flag = 0;
  flagg = 0;

// ===========================================================================
//
// Fix flow widths of zero
//
// ===========================================================================

  if (flowedg->getVEdgLen() < THRESH) {

    // Check to see if the node vv_flow falls in the 
    // sector bounded by 'ccw_edg' & 'cw_edg'. If it 
    // does not -> the Voronoi geometry is complicated
    // and we just ignore correction of the Voronoi width

    if (IsInTriangle(centroid, dest_cw,   dest_flow, vv_flow[0], vv_flow[1]) ||
	IsInTriangle(centroid, dest_flow, dest_ccw,  vv_flow[0], vv_flow[1]) ||
	IsInTriangle(centroid, dest_cw,   dest_flow, vv_ccw1[0], vv_ccw1[1]) ||
 	IsInTriangle(centroid, dest_flow, dest_ccw,  vv_ccw1[0], vv_ccw1[1]) ) {

      // First vertex in the sector area is the flow voronoi vertex
      j = 0;
      polyg[0][j] = vv_flow[0];
      polyg[1][j] = vv_flow[1];
      j++;

      // Now we have two lines that are defined by vertices 
      // vv_flow & vv_ccw1, vv_flow & vv_cw2, and two edges ccwedg & cwedg 
      // - neighbors of flowedg. We need to define where the
      // lines cross the edges and then find an area of a 
      // corresponding polygon. 

      // Find intersection of the first and second voronois with the CCW edge
      // If the intersection is in the triangle add one point to section area
      // otherwise add two points to section area
      intersect_ccw = FindIntersectionCoords(vv_flow, vv_ccw1,
                                             orig_ccw, dest_ccw);

      if (IsInTriangle(centroid, dest_flow, dest_ccw, 
                       intersect_ccw[0], intersect_ccw[1]) ) {

        // If that intersection is inside the triangle of centroid to dest node
        // of the flow edge to the dest node of the CCW edge
        if (intersect_ccw != vv_flow) {
          // And the intersection is not the flow voronoi vertex
          // then the intersection is the second point in the section area
 	  polyg[0][j] = intersect_ccw[0];
 	  polyg[1][j] = intersect_ccw[1];
	  j++;
	}
      } else {
        // Otherwise the second voronoi vertex is the second point in the
        // section area, and we look for a third point
 	polyg[0][j] = vv_ccw1[0];
 	polyg[1][j] = vv_ccw1[1];
	j++;

        // Find intersection of second and third voronoi vertices with CCW edge
	intersect_ccw = FindIntersectionCoords(vv_ccw1, vv_ccw3,
                                               orig_ccw, dest_ccw);

	if (IsInTriangle(centroid, dest_flow, dest_ccw, 
                         intersect_ccw[0], intersect_ccw[1])) {
	  // If intersection falls in the triangle it is the third sector point
 	  polyg[0][j] = intersect_ccw[0];
 	  polyg[1][j] = intersect_ccw[1];
	  j++;
	} else if (IsInTriangle(centroid, dest_flow, dest_ccw, 
                                vv_ccw3[0], vv_ccw3[1])) {
          // Otherwise if the second voronoi vertex falls in triangle it is
          // the third sector point
 	  polyg[0][j] = vv_ccw3[0];
 	  polyg[1][j] = vv_ccw3[1];
	  j++;
	}
      }
      // Now add the centroid as the fourth sector point
      polyg[0][j] = centroid[0];
      polyg[1][j] = centroid[1];
      j++;

      // Wrap around and deal with the far side of the CW triangle
      // Find intersection of the first voronoi and first CW voronoi with
      // the CW edge
      intersect_cw = FindIntersectionCoords(vv_flow, vv_cw2,
                                            orig_cw, dest_cw);

      if (IsInTriangle(centroid, dest_cw, dest_flow, 
                       intersect_cw[0], intersect_cw[1]) ) {
        // If that intersection is inside the CW triangle
        if (intersect_cw != vv_flow) {
          // And the intersection is not the flow voronoi vertex
          // then the intersection is the fifth sector point
 	  polyg[0][j] = intersect_cw[0];
 	  polyg[1][j] = intersect_cw[1];
	  j++;
	}
      } else {
        // If the intersection is outside the CW triangle
        if (vv_cw2 != vv_flow) {
          // And the last voronoi vertex is not the first voronoi vertex
          // then the last voronoi vertex is the fourth sector point
 	  polyg[0][j] = vv_cw2[0];
 	  polyg[1][j] = vv_cw2[1];
	  j++;
	}
      }
      sArea = polygonArea(polyg, j);
      changeWidth++; // <- It MUST be here, never eliminate!
    } else {
        cout<<"\tbFlowNet Cell: " << cn->getID() 
            << " Can not re-define zero flow width! "<<endl;
        cout<<"\tvv_flow node is not within the ccwedg - cwedg sector"<<endl;
        TellAboutNode(cn);
        cout<<endl;
    }
  }

// ===========================================================================
//
// Fix flow widths size greater than 0 (all others)
//
// ===========================================================================

  else {
    // Check to see if the voronoi vertex associated with the flow edge or the
    // voronoi vertex ccw from that fall in the sector bounded by ccwedg - cwedg
    // If they do not the Voronoi geometry is complicated and we skip the fix
    if (IsInTriangle(centroid, dest_cw,   dest_flow, vv_flow[0], vv_flow[1]) ||
	IsInTriangle(centroid, dest_flow, dest_ccw,  vv_flow[0], vv_flow[1]) ||
	IsInTriangle(centroid, dest_cw,   dest_flow, vv_ccw1[0], vv_ccw1[1]) ||
 	IsInTriangle(centroid, dest_flow, dest_ccw,  vv_ccw1[0], vv_ccw1[1]) ) {
	
      // Collect the 3, 4 or 5 points which will define the sector area
      // Always the centroid and either 1 or 2 points on CW and CCW sides
      j = 0;

      // CW edge intersects line from flow voronoi vertex to next ccw vertex
      if (!AreSegmentsParallel(vv_flow, vv_ccw1, orig_cw, dest_cw)) {

        // Find intersection of vv_flow to vv_ccw1 with CW edge
        intersect_cw = FindIntersectionCoords(vv_flow, vv_ccw1,
                                              orig_cw, dest_cw);
    
        if (IsBetweenEndPnts(vv_flow, vv_ccw1, orig_cw, dest_cw,
                             intersect_cw[0], intersect_cw[1])) {
          // Intersection is within minimum bounding box of 4 endpoints
 	  flag = 1;
          polyg[0][j] = intersect_cw[0];
          polyg[1][j] = intersect_cw[1];
          j++;
        }
      }

      if (flag == 0) {
        // Intersection is outside minimum bounding box or segments parallel
        // Use the flow voronoi vertex instead
        polyg[0][j] = vv_flow[0];
        polyg[1][j] = vv_flow[1];
        j++;
      }

      // CCW edge intersects line from flow voronoi vertex to next ccw vertex
      if (!AreSegmentsParallel(vv_flow, vv_ccw1, orig_ccw, dest_ccw)) {

        // Find intersection of vv_flow to vv_ccw1 with CCW edge
        intersect_ccw = FindIntersectionCoords(vv_flow, vv_ccw1, 
                                               orig_ccw, dest_ccw);
        if (IsBetweenEndPnts(vv_flow, vv_ccw1, orig_ccw, dest_ccw,
                             intersect_ccw[0], intersect_ccw[1])) {
          // Intersection is within minimum bounding box of 4 endpoints
          // Keep the intersection point
          flagg = 1;
          polyg[0][j] = intersect_ccw[0];
          polyg[1][j] = intersect_ccw[1];
          j++;
        }
      }

      // We didn't find a good intersection points on the CCW side
      // We'll add two points instead of one to sector polygon
      if (flagg == 0) {
        // Intersection is outside minimum bounding box or segments parallel
        // Use voronoi vertex CCW to flow voronoi vertex instead
        polyg[0][j] = vv_ccw1[0];
        polyg[1][j] = vv_ccw1[1];
        j++;
 
        // Find intersection of CCW edge segment with next pair of voronoi verts
        intersect_ccw = FindIntersectionCoords(vv_ccw1, vv_ccw3, 
                                               orig_ccw, dest_ccw);

        // Should there be a minimum bounding box test again ???
        polyg[0][j] = intersect_ccw[0];
        polyg[1][j] = intersect_ccw[1];
        j++;
      }

      // Add the centroid to the sector polygon
      polyg[0][j] = centroid[0];
      polyg[1][j] = centroid[1];
      j++;

      // We've wrapped around to the CW side again so pick up the extra point
      // if the initial intersection point was not usable
      if (flag == 0) {
        // Find intersection ov CW edge segment with flow vv and cw vv
        intersect_cw = FindIntersectionCoords(vv_flow, vv_cw2, 
                                              orig_cw, dest_cw);

        // Make sure the intersection is on the CW edge ??
        polyg[0][j] = intersect_cw[0];
        polyg[1][j] = intersect_cw[1];
        j++;
      }

      // Sector polygon was defined with 3 points if both intersections were
      // good or with 5 points if both were bad
      sArea = polygonArea(polyg, j);
      changeWidth++; //It MUST be here, never eliminate!
    }
    // If neither vv_flow or xy1 is within the sector of CCW to CW then we
    // aren't changing the width (flag will still == 0) and we won't
    // execute the next block of code
  }

// ===========================================================================
//
// If we DID fix the flow width for a current Voronoi cell,
// then proceed with the following:
//
// ===========================================================================

  // Just check of a situation when Voronoi cell is somewhat strange...
  // If we enter this code we will not change the voronoi width
  if (cn->getID() == -999 ||
     ((changeWidth > 0) && (sArea > 1.0E+6)) ) {
    
    firstedg = cn->getFlowEdg();
    curedg = firstedg->getCCWEdg();

    cout << endl;
    cout << "ID: " << cn->getID()
         << "\twidth1 = " << (firstedg->getVEdgLen()*1000.0)
         << "\twidth2 = " << (curedg->getVEdgLen()*1000.0) << endl;

    vv_cur = firstedg->getRVtx();
    cout<<"Voronoi Vertex = [ "<<vv_cur[0]<<" " << vv_cur[1] << "]" << endl;
         
    while ( curedg != firstedg )  {
      vv_cur = curedg->getRVtx();
      cout<<"Voronoi Vertex = [ "<<vv_cur[0]<<" " << vv_cur[1] << "]" << endl;
      curedg = curedg->getCCWEdg();
    } 

    curedg = firstedg;
    cout<<"Slope: " << curedg->getSlope()<<"\t";
    TellAboutNode((bFlowNode *)curedg->getDestinationPtrNC());
    cout<<"       -9999999\t";
    TellAboutNode(cn);
    curedg = firstedg->getCCWEdg();
    while ( curedg != firstedg )  {
      cout<<"Slope: " << curedg->getSlope()<<"\t";
      TellAboutNode((bFlowNode *)curedg->getDestinationPtrNC());
      cout<<"       -9999999\t";
      TellAboutNode(cn);
      curedg = curedg->getCCWEdg();
    }

    cout<<"\nvv_flow = "<<vv_flow[0]<<", "<<vv_flow[1]<<endl;
    cout<<"vv_ccw1 = "<<vv_ccw1[0]<<", "<<vv_ccw1[1]<<endl;
    cout<<"vv_cw2 = "<<vv_cw2[0]<<", "<<vv_cw2[1]<<endl;
    cout<<"vv_ccw3 = "<<vv_ccw3[0]<<", "<<vv_ccw3[1]<<endl;
 
    cout<<"intersect_ccw = "<<intersect_ccw[0]<<", "<<intersect_ccw[1]<<endl;
    cout<<"Centroid = "<<centroid[0]<<", "<<centroid[1]<<endl;
    cout<<"intersect_cw = "<<intersect_cw[0]<<", "<<intersect_cw[1]<<endl;

    if (changeWidth > 0) {
      NV = j;
      cout<<"changeWidth = "<<changeWidth<<"; flagg = "<<flagg<<"; j = "<<j<<endl;
      for (i=0; i < 2; i++) {
        for (j=0; j < NV; j++)
	  cout<<polyg[i][j]<<" ";
	  cout<<endl;
      }
    }
    changeWidth = 0;
  }

  if ( changeWidth ) {
    if (flowedg->getLength() > 0.0)
      wWidth = sArea/flowedg->getLength();
    else 
      wWidth = -99999; // To ensure we dont get value

    if (flowedg->getVEdgLen() < wWidth) {
/*
      cout<<"\tsArea = "<<sArea
          <<"\tWIDTH BEFORE =  "<<flowedg->getVEdgLen()
          <<"\twWIDTH = "<<wWidth<<endl;
*/
      flowedg->setVEdgLen(wWidth);
    }
  }
  return;
}

/***************************************************************************
**
** Function:  polygonArea
** Arguments: *poly[2] - coordinates of NOT closed polygon
**            - poly[0][] - X coordinates
**            - poly[1][] - Y coordinates
** Objective: calculates polygon area: polygon is 
**            considered to be NON-closed
** Return value: - area of a polygon
** Algorithm: [O'Rourke], page 21 
**
***************************************************************************/
double bFlowNet::polygonArea(double **poly, int L) 
{
  double sum = 0.0;
  for (int i=0; i < L; i++) {
    if (i == (L-1))
      sum += (poly[0][i]*poly[1][0] - poly[1][i]*poly[0][0]);
    else 
      sum += (poly[0][i]*poly[1][i+1] - poly[1][i]*poly[0][i+1]);
  }
  if (sum == 0.0)
    return 0.0;
  else
    return(fabs(sum/2.));
}

/*************************************************************************
**
**  IsBetweenEndPnts
**
**  The function checks if the point (x,y) falls between the endpoints 
**  of two line segments defined by xy1(x1,y1) - xy2(x2,y2) and 
**  xy3(xx1,yy1) - xy4(xx2,yy2).  Returns '1' if yes, '0' otherwise
**
*************************************************************************/
int bFlowNet::IsBetweenEndPnts(
		std::vector<double>& xy1, std::vector<double>& xy2,
		std::vector<double>& xy3, std::vector<double>& xy4,
		double x, double y)
{
  int result;
  double minx1, minx2, maxx1, maxx2;
  double miny1, miny2, maxy1, maxy2;
  double x1,x2,xx1,xx2,y1,y2,yy1,yy2;

  x1  = xy1[0];
  x2  = xy2[0];
  xx1 = xy3[0];
  xx2 = xy4[0];
  y1  = xy1[1];
  y2  = xy2[1];
  yy1 = xy3[1];
  yy2 = xy4[1];

  if (x1 < x2) {
    minx1 = x1;
    maxx1 = x2;
  }
  else { 
    minx1 = x2;  // <-- Can be a point if x1 = x2;
    maxx1 = x1;
  }
  if (xx1 < xx2) {
    minx2 = xx1;
    maxx2 = xx2;
  }
  else { 
    minx2 = xx2;  // <-- Can be a point if xx1 = xx2;
    maxx2 = xx1;
  }
  if (y1 < y2) {
    miny1 = y1;
    maxy1 = y2;
  }
  else { 
    miny1 = y2;  // <-- Can be a point if y1 = y2;
    maxy1 = y1;
  }
  if (yy1 < yy2) {
    miny2 = yy1;
    maxy2 = yy2;
  }
  else { 
    miny2 = yy2;  // <-- Can be a point if yy1 = yy2;
    maxy2 = yy1;
  }

  // In order to deal with the numerical issues...
  minx1 -= THRESH;
  maxx1 += THRESH;
  minx2 -= THRESH;
  maxx2 += THRESH;
  miny1 -= THRESH;
  maxy1 += THRESH;
  miny2 -= THRESH;
  maxy2 += THRESH;

  if ((x >= minx1 && x <= maxx1) &&
      (x >= minx2 && x <= maxx2) &&
      (y >= miny1 && y <= maxy1) &&
      (y >= miny2 && y <= maxy2))
    result = 1;
  else
    result = 0;
  return result;
}

/*************************************************************************
**
**  AreSegmentsParallel
**
**  The function checks if the two segments defined by xy1(x1,y1)-xy2(x2,y2) 
**  and xy3(xx1,yy1)-xy4(xx2,yy2) are parallel. Returns '1' if yes, '0' 
**  otherwise
**
*************************************************************************/
int bFlowNet::AreSegmentsParallel(
		std::vector<double>& xy1, std::vector<double>& xy2,
		std::vector<double>& xy3, std::vector<double>& xy4)
{
  int result;
  double dxa, dxb, dya, dyb, det;
  double x1,x2,xx1,xx2,y1,y2,yy1,yy2;
  x1  = xy1[0];
  x2  = xy2[0];
  xx1 = xy3[0];
  xx2 = xy4[0];
  y1  = xy1[1];
  y2  = xy2[1];
  yy1 = xy3[1];
  yy2 = xy4[1];

  // Check if  segments exist
  if ( ((x2 == x1) && (y1 == y2)) || 
     ( (xx1 == xx2) && (yy1 == yy2)) ) {
    cout<<"Segment doesn't exist:"<<endl;
    cout<<"X1 = "<<x1<<"; Y1 = "<<y1<<"; X2 = "<<x2<<"; Y2 = "<<y2
        <<"XX1 = "<<xx1<<"; YY1 = "<<yy1
        <<"; XX2 = "<<xx2<<"; YY2 = "<<yy2<<endl;
    result = 0;
  }
  else {
    dxa = x2 - x1;
    dxb = xx2 - xx1;
    dya = y2 - y1;
    dyb = yy2 - yy1;
    // Check if the lines are not parallel
    det = dya*dxb-dyb*dxa;
    if (det !=0 )
      result = 0;
    else
      result = 1;
  }
  return result;
}

/***************************************************************************
**
**  IsInTriangle()
**
**  Defines if a point (x,y) falls in the triangle. The algorithm exploits
**  the fact that the 3 triangle points are always in counter-clockwise
**  order, so that the point is contained within a given triangle (p0,p1,p2)
**  if and only if the point lies to the left of vectors p0->p1, p1->p2,
**  and p2->p0. Here's how it works:
** 
***************************************************************************/
int bFlowNet::IsInTriangle(std::vector<double>& xyp1,
                           std::vector<double>& xyp2,
                           std::vector<double>& xyp3,
			   double x, double y) 
{
  int k;
  double a, b, c;

  k = 1;
  for (int i=0; (i<3)&&(k>0) ; i++) {

    if (i == 0) {
      a = (xyp1[1] - y) * (xyp2[0] - x);
      b = (xyp1[0] - x) * (xyp2[1] - y);
      c = a - b;
    }
    else if (i == 1) {
      a = (xyp2[1] - y) * (xyp3[0] - x);
      b = (xyp2[0] - x) * (xyp3[1] - y);
      c = a - b;
    }
    else if (i == 2) {
      a = (xyp3[1] - y) * (xyp1[0] - x);
      b = (xyp3[0] - x) * (xyp1[1] - y);
      c = a - b;
    }

    if ( c > 0.0 )     // <--- Not to the LEFT
      k = -1;
    else { 
      if ( c == 0.0 )  // <--- on the BND
	;
    }
  }
  if (k > 0)
    return 1;
  else 
    return 0;
}

/***************************************************************************
**  
** bFlowNet::SetReachInformation() Function
**  
** Set the reach ID for each node
**
***************************************************************************/

void bFlowNet::SetReachInformation()
{
  // Assign reaches for stream nodes
  bFlowNode *chead, *coutlet, *cn;
  std::_List_iterator<bFlowNode*> HeadIter;
  std::_List_iterator<bFlowNode*> OutletIter;

  int reach;
  for (HeadIter = NodesLstH.begin(), OutletIter = NodesLstO.begin(), reach = 0;
       HeadIter != NodesLstH.end();
       HeadIter++, OutletIter++, reach++) {

    chead = (*HeadIter);
    coutlet = (*OutletIter);
    cn = chead;

    while (cn != coutlet) {
      cn->setReach(reach);
      cn = cn->getDownstrmNbr();
    }
  }

  // Set the final outlet node to the final reach for the next part
  OutletNode->setReach(NodesLstH.size() - 1);

  // Set all the other node reaches using the stream nodes
  bMeshList<bFlowNode*> *nlist = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> niter;
  for (niter = nlist->begin(); niter != nlist->getLastActive(); niter++) {
    cn = (*niter);

    if (cn->getBoundaryFlag() != kStream) {
      reach = cn->getStreamNode()->getReach();
      cn->setReach(reach);
    }
  }

  // Set the final outlet node back to "no reach" for the rest
  OutletNode->setReach(-1);
}

/***************************************************************************
**
** bFlowNet::WriteFlowNet() Function
**
** Write information needed to build tFlowNet in the simulator
**
***************************************************************************/
                                                                                
void bFlowNet::WriteFlowNet() const
{
  fstream flowStr("flow.meshb", ios::out | ios::binary);

  BinaryWrite(flowStr, hillvel);
  BinaryWrite(flowStr, streamvel);
  BinaryWrite(flowStr, velratio);
  BinaryWrite(flowStr, velcoef);
  BinaryWrite(flowStr, flowexp);
  BinaryWrite(flowStr, baseflow);
  BinaryWrite(flowStr, flowout);
  BinaryWrite(flowStr, maxttime);
  BinaryWrite(flowStr, dist_hill_max);
  BinaryWrite(flowStr, dist_stream_max);
  BinaryWrite(flowStr, BasArea);

  BinaryWrite(flowStr, OutletNode->getID());

  std::list<bFlowNode*>::const_iterator niter;
  std::list<int>::const_iterator iter;

  BinaryWrite(flowStr, (int) HeadsLst.size());
  for (niter = HeadsLst.begin(); niter != HeadsLst.end(); niter++)
     BinaryWrite(flowStr, (*niter)->getID());

  BinaryWrite(flowStr, (int) NodesLstH.size());
  for (niter = NodesLstH.begin(); niter != NodesLstH.end(); niter++)
     BinaryWrite(flowStr, (*niter)->getID());

  BinaryWrite(flowStr, (int) NodesLstO.size());
  for (niter = NodesLstO.begin(); niter != NodesLstO.end(); niter++)
     BinaryWrite(flowStr, (*niter)->getID());

  BinaryWrite(flowStr, (int) NNodes.size());
  for (iter = NNodes.begin(); iter != NNodes.end(); iter++)
     BinaryWrite(flowStr, (*iter));

  flowStr.close();
}

/***************************************************************************
**   
** bFlowNet::WriteGeometry() Function
**   
** Write information needed to for visualizer
**
***************************************************************************/
                                                                          
void bFlowNet::WriteGeometry() const
{    
  bMeshList<bFlowNode*>* nodeList = gridPtr->getNodeList();
  std::_List_iterator<bFlowNode*> niter;
  std::vector<double> voronoiVertex;
  bFlowNode *cn;
  bEdge *firstedg, *curedg;

  fstream geomStr("voronoi.meshb", ios::out | ios::binary);
     
  int nActiveNodes = nodeList->getActiveSize();
  int nPoints = 0;
  int* edgesPerCell = new int[nActiveNodes];

  // Collect the number of voronoi vertex points per cell
  int index = 0;
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++) {
    cn = (*niter);
    firstedg = cn->getFlowEdg();
    curedg = firstedg->getCCWEdg();
    edgesPerCell[index] = 1;

    while (curedg != firstedg) {
      edgesPerCell[index]++;
      curedg = curedg->getCCWEdg();
    }
    nPoints += edgesPerCell[index];
    index++;
  }

  // Write the counts of nodes and points
  BinaryWrite(geomStr, nActiveNodes);
  BinaryWrite(geomStr, nPoints);

  // Write all the voronoi vertices (points in unstructured grid)
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++) {

    firstedg = (*niter)->getFlowEdg();
    voronoiVertex = firstedg->getRVtx();
    BinaryWrite(geomStr, (float) voronoiVertex[0]);
    BinaryWrite(geomStr, (float) voronoiVertex[1]);

    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg) {
      voronoiVertex = curedg->getRVtx();
      BinaryWrite(geomStr, (float) voronoiVertex[0]);
      BinaryWrite(geomStr, (float) voronoiVertex[1]);
      curedg = curedg->getCCWEdg();
    }
  }

  // Write the edges per used cell
  index = 0;
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, edgesPerCell[index++]);

  // Write the centroid points per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getX());
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getY());

  // Centroid z is the elevation
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getZ());

  // Write the boundary per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getBoundaryFlag());

  // Write the reach number per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getReach());

  // Write the flood status per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getFloodStatus());

  // Write the flow width per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getFlowEdg()->getVEdgLen());

  // Write the flow length per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getFlowEdg()->getLength());

  // Write the contributing area per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getContrArea());

  // Write the slope per cell
  for (niter = nodeList->begin(); niter != nodeList->getLastActive(); niter++)
    BinaryWrite(geomStr, (float) (*niter)->getFlowEdg()->getSlope());

  delete [] edgesPerCell;
  geomStr.close();
} 

//=========================================================================
//
//
//                           End of bFlowNet 
//
//
//=========================================================================
