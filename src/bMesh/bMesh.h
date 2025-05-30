// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bMesh.h: Header file for class bMesh
**
**  bMesh is the master class that handles the implementation of the
**  triangulated mesh. The class includes lists of the mesh elements
**  (nodes, triangles, and edges; see meshElements.h/.cpp), and provides
**  functionality to:
**    - read in or create meshes, either from scratch, from a list of
**      points, from a pre-existing set of triangulation files (e.g., a
**      previous run), or an Arc/Info files
**    - move, add, and/or delete nodes
**    - update Delaunay and Voronoi geometry
**
\***************************************************************************/

#ifndef BMESH_H
#define BMESH_H

#include "src/Headers/Inclusions.h"
#include "src/bMesh/bTriangulator.h"
#include "src/bMeshElements/bMeshElements.h"
#include "src/bInOut/bInputFile.h"
#include "src/bListInputData/bListInputData.h"
//#include "src/bList/bList.h"

#ifdef ALPHA_64
  #include <stdlib.h>
  #include <strings.h>
  #include <list.h>
#elif defined LINUX_32
  #include <cstdlib>
  #include <string>
  #include <list>
  #include <vector>
#else 
  #include <stdlib.h>
  #include <strings.h>
  #include <list.h>
  #include <vector>
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 1: bMesh Class Declarations
//
//
//=========================================================================

template< class bSubNode >
class bMesh
{
  bMesh(const bMesh&);
  bMesh& operator=(const bMesh&);
public:
   bMesh();
   ~bMesh();
 
   // Creates mesh using Tipper triangulation and writes to files
   void MakeBasicMesh(bInputFile& infile);

   // Reads files from triangulator, builds mesh, enhances with voronoi
   void RefineBasicMesh();

   // Writes bNode mesh information
   void WriteBasicMesh();

   // Reads bNode mesh information but populates bFlowNode mesh
   void ReadBasicMesh();
   
   // Calculate counter clockwise edges so that spokes can be eliminated
   void MakeCCWEdges();

   // Creating Mesh from Existing Mesh 
   void MakeMeshFromInputData(bInputFile& infile);

   // Calculate voronoi vertex and areas and save in nodes and edges
   void setVoronoiVertices();
   void CalcVoronoiEdgeLengths();
   void CalcVAreas();
   void CheckMeshConsistency( int boundaryCheckFlag=1 );
   void TellAboutNode(bSubNode *);

   bTriangle *LocateTriangle( double, double );
   bTriangle *TriWithEdgePtr( bEdge * );

   bMeshList<bEdge*>*	  getEdgeList()		{ return &edgeList; }
   bMeshList<bSubNode*>*  getNodeList()		{ return &nodeList; }
   std::list<bTriangle*>* getTriList()		{ return &triList; }

#ifndef NDEBUG
   void DumpEdges();
   void DumpSpokes( bSubNode * );
   void DumpTriangles();
   void DumpNodes();
#endif
   
protected:
   int nnodes;			       	// # of nodes
   int nedges;			       	// # of edges
   int ntri;  			     	// # of triangles

   bMeshList<bSubNode*>  nodeList; 	// list of nodes
   bMeshList<bEdge*>	 edgeList;      // list of directed edges
   std::list<bTriangle*> triList;     	// list of triangles
   bTriangle* mSearchOriginTriPtr; 	// triangle to start searches from
};

// DEFAULT CONSTRUCTOR
template< class bSubNode >
bMesh< bSubNode >:: bMesh()
{
    nnodes = nedges = ntri = 0;
    mSearchOriginTriPtr=0;
}

// DESTRUCTOR
template< class bSubNode >
bMesh< bSubNode >::~bMesh()
{
    cout << "bMesh delete nodes" << endl;
    std::_List_iterator<bSubNode*> niter;
    for (niter = nodeList.begin(); niter != nodeList.end(); niter++)
        delete (*niter);
    nodeList.erase(nodeList.begin(), nodeList.end());

    cout << "bMesh delete edges" << endl;
    std::_List_iterator<bEdge*> eiter;
    for (eiter = edgeList.begin(); eiter != edgeList.end(); eiter++)
        delete (*eiter);
    edgeList.erase(edgeList.begin(), edgeList.end());

    cout << "bMesh delete triangles" << endl;
    std::_List_iterator<bTriangle*> titer;
    for (titer = triList.begin(); titer != triList.end(); titer++)
        delete (*titer);
    triList.erase(triList.begin(), triList.end());

    cout << "bMesh Object has been destroyed..." << endl;
}

//=========================================================================
//
//
//  bMesh::MakeMeshFromInputData()
//  
//   Orignal code from tRIBS, tMesh.cpp but major changes to the orignal code 
//   code were made to get the code working in MeshBuilder. CJC 2025
//   Constructs tListInputData object and makes mesh from data in that object.
//                    
//   Calls: tListInputData( infile ), CheckMeshConsistency()
//   Inputs: infile -- main input file from which various items are read
//=========================================================================

template< class bSubNode >
void bMesh< bSubNode >::MakeMeshFromInputData(bInputFile& infile)
{
    bListInputData< bSubNode > input( infile );

    int i;
    nnodes = input.x.size();
    nedges = input.orgid.size();
    ntri = input.p0.size();
        
    assert( nnodes > 0 );
    assert( nedges > 0 );
    assert( ntri > 0 );
        
    // Create the node list by creating a temporary node and iteratively
    // (1) assigning it values from the input data and (2) inserting it
    // onto the back of the node list.
        
    // --- NODES ---
    cout << "\nCreating node list..." << endl;
    // Populate nodeList 
    for( int i = 0; i < nnodes; i++ ) {
        bSubNode* tempnode_ptr = new bSubNode(); // <<< USE DEFAULT CONSTRUCTOR
        // Now set all properties for tempnode_ptr
        tempnode_ptr->set3DCoords( input.x[i], input.y[i], input.z[i] );
        tempnode_ptr->setID( i ); // Crucial: 0 to nnodes-1
        int bound = input.boundflag[i]; // boundflag comes from tListInputData
        tempnode_ptr->setBoundaryFlag( bound );

        // Use bMeshList's insertion methods
        if( (bound == kNonBoundary) || (bound == kStream) ) // Use defined constants if available
            nodeList.insertAtActiveBack( tempnode_ptr );
        else if( bound == kOpenBoundary )
            nodeList.insertAtBoundFront( tempnode_ptr );
        else // kClosedBoundary
            nodeList.insertAtBack( tempnode_ptr );
    }

    // Create NodeTable 
    std::vector<bSubNode*> NodeTable_vec(nnodes, nullptr);
    for (std::_List_iterator<bSubNode*> it = nodeList.begin(); it != nodeList.end(); ++it) {
        bSubNode* node_ptr = *it;
        if (node_ptr && node_ptr->getID() >= 0 && node_ptr->getID() < nnodes) {
            NodeTable_vec[node_ptr->getID()] = node_ptr;
        }
    }

    // --- EDGES ---
    cout << "\nCreating edge list..." << endl;
    // Populate edgeList (bMeshList<bEdge*>, i.e., std::list<bEdge*>)
    int miNextEdgID;
    for( miNextEdgID = 0; miNextEdgID < nedges; miNextEdgID++ ){ // Iterate for each edge from input
        bEdge* tempedge_ptr = new bEdge();
        tempedge_ptr->setID( miNextEdgID ); // ID is 0 to nedges-1, matching input.orgid's implicit index

        bSubNode *nodPtr1 = NodeTable_vec[ input.orgid[miNextEdgID] ]; // Use NodeTable_vec
        tempedge_ptr->setOriginPtr( nodPtr1 );
        int obnd = nodPtr1->getBoundaryFlag();

        bSubNode *nodPtr2 = NodeTable_vec[ input.destid[miNextEdgID] ]; // Use NodeTable_vec
        tempedge_ptr->setDestinationPtr( nodPtr2 );
        int dbnd = nodPtr2->getBoundaryFlag();
        
        // Set flow allowed and insert into edgeList
        if( obnd == kClosedBoundary || dbnd == kClosedBoundary
            || (obnd==kOpenBoundary && dbnd==kOpenBoundary) ) {
            tempedge_ptr->setFlowAllowed( 0 );
            edgeList.insertAtBack( tempedge_ptr );
        } else {
            tempedge_ptr->setFlowAllowed( 1 );
            edgeList.insertAtActiveBack( tempedge_ptr );
        }
    }
     // The original code incremented miNextEdgID by 2 because it created pairs.
     // Your tListInputData provides orgid[i], destid[i], nextid[i] for each individual edge 'i'.
     // So the loop for miNextEdgID should go from 0 to nedges-1.

    // Create EdgeTable (std::vector based)
    std::vector<bEdge*> EdgeTable_vec(nedges, nullptr);
    for (std::_List_iterator<bEdge*> it = edgeList.begin(); it != edgeList.end(); ++it) {
        bEdge* edge_ptr = *it;
        if (edge_ptr && edge_ptr->getID() >= 0 && edge_ptr->getID() < nedges) {
            EdgeTable_vec[edge_ptr->getID()] = edge_ptr;
        }
    }

    // --- SPOKE LISTS & CCW EDGES ---
    cout << "\nSetting up spoke lists and CCW edges..." << endl;
    for (std::_List_iterator<bSubNode*> it = nodeList.begin(); it != nodeList.end(); ++it) {
        bSubNode* curnode = *it;
        int node_id = curnode->getID(); // This is 0 to nnodes-1

        // input.edgid[node_id] is the ID of the first spoke for this node
        // input.nextid[edge_id] is the next CCW edge for 'edge_id'
        
        const int first_spoke_id = input.edgid[node_id];
        if (first_spoke_id < 0 || first_spoke_id >= nedges) { /* error */ continue; }

        bEdge *first_edge_ptr = EdgeTable_vec[first_spoke_id];
        if (!first_edge_ptr) { /* error */ continue; }

        curnode->setEdg( first_edge_ptr ); // Set first edge
        // The spoke list in bNode seems to build itself if setEdg and ccw links are correct.
        // Or, we might need to populate it manually if that's how tMesh did it.
        // The original MakeMeshFromInputData did:
        // curnode->insertBackSpokeList( edgPtr );
        // int ne;
        // for( ne = input.nextid[e1]; ne != e1; ne = input.nextid[ne] ) {
        //    bEdge *nextEdgPtr = EdgeTable_vec[ne];
        //    curnode->insertBackSpokeList( nextEdgPtr );
        // }
        // Let's assume bNode::makeCCWEdges() or similar handles full spoke list from first edge + CCW links.

        // If `bNode::insertBackSpokeList` is necessary:
        curnode->getSpokeListNC().clear(); // Clear any existing spokes
        bEdge* current_spoke_ptr = first_edge_ptr;
        int current_spoke_id = first_spoke_id;
        do {
            curnode->insertBackSpokeList(current_spoke_ptr);
            current_spoke_id = input.nextid[current_spoke_id]; // Get next CCW edge ID from .edges file data
            if (current_spoke_id < 0 || current_spoke_id >= nedges) { /* error, break */ break;}
            current_spoke_ptr = EdgeTable_vec[current_spoke_id];
            if (!current_spoke_ptr) { /* error, break */ break; }
        } while (current_spoke_id != first_spoke_id && curnode->getSpokeListNC().size() < 100); // Safety break
    }

    // Set CCW edges for each edge
    for (std::_List_iterator<bEdge*> it = edgeList.begin(); it != edgeList.end(); ++it) {
        bEdge* curedg = *it;
        int edge_id = curedg->getID(); // 0 to nedges-1
        int ccw_neighbor_id = input.nextid[edge_id]; // From .edges file data

        if (ccw_neighbor_id >= 0 && ccw_neighbor_id < nedges) {
            bEdge* ccwedg_ptr = EdgeTable_vec[ccw_neighbor_id];
            curedg->setCCWEdg( ccwedg_ptr );
        } else { /* error: invalid CCW edge ID */ }
    }


    // --- TRIANGLES ---
    cout << "\nSetting up triangle connectivity..." << endl;
    // Change triList to std::list<bTriangle*> or bMeshList<bTriangle*>
    // For now, assuming it's std::list<bTriangle*> as in the original bMesh.h code
    for (int i = 0; i < ntri; i++ ) {
        bTriangle* newtri_ptr = new bTriangle();
        newtri_ptr->setID( i ); // 0 to ntri-1
        newtri_ptr->setPPtr( 0, NodeTable_vec[ input.p0[i] ] );
        newtri_ptr->setPPtr( 1, NodeTable_vec[ input.p1[i] ] );
        newtri_ptr->setPPtr( 2, NodeTable_vec[ input.p2[i] ] );
        newtri_ptr->setEPtr( 0, EdgeTable_vec[ input.e0[i] ] );
        newtri_ptr->setEPtr( 1, EdgeTable_vec[ input.e1[i] ] );
        newtri_ptr->setEPtr( 2, EdgeTable_vec[ input.e2[i] ] );
        triList.push_back( newtri_ptr ); // If std::list
        // If bMeshList: triList.insertAtBack(newtri_ptr); or similar
    }

    // Create TriTable (std::vector based)
    std::vector<bTriangle*> TriTable_vec(ntri, nullptr);
    for (std::_List_iterator<bTriangle*> it = triList.begin(); it != triList.end(); ++it) {
        bTriangle* tri_ptr = *it;
        if (tri_ptr && tri_ptr->getID() >= 0 && tri_ptr->getID() < ntri) {
            TriTable_vec[tri_ptr->getID()] = tri_ptr;
        }
    }
    
    // Set triangle neighbors
    for (std::_List_iterator<bTriangle*> it = triList.begin(); it != triList.end(); ++it) {
        bTriangle* ct = *it;
        if (!ct) continue;
        int current_tri_id = ct->getID(); // 0 to ntri-1

        bTriangle* nbrtri;
        nbrtri = (input.t0[current_tri_id] >= 0 && input.t0[current_tri_id] < ntri) ? TriTable_vec[ input.t0[current_tri_id] ] : nullptr;
        ct->setTPtr( 0, nbrtri );
        nbrtri = (input.t1[current_tri_id] >= 0 && input.t1[current_tri_id] < ntri) ? TriTable_vec[ input.t1[current_tri_id] ] : nullptr;
        ct->setTPtr( 1, nbrtri );
        nbrtri = (input.t2[current_tri_id] >= 0 && input.t2[current_tri_id] < ntri) ? TriTable_vec[ input.t2[current_tri_id] ] : nullptr;
        ct->setTPtr( 2, nbrtri );
    }

    cout<<"\nTesting Mesh..."<<endl;
    CheckMeshConsistency(); 
}

//=========================================================================
//
//
//                  Section 5a: bMesh:: MakeMeshFromTriangulator()
//
//
//=========================================================================


/**************************************************************************
**
**   bMesh::MakeBasicMesh( infile )
**
**   Makes a mesh from a points file using Tipper's triangulation algorithm.
**
**   Created: 07/2002, Arnaud Desitter, Greg Tucker, Oxford
**   Modified: 08/2002, MIT
**
**************************************************************************/

// edge numbering translation
static int e_t2c(int ei, bool o){ // Tipper to child
    return o? 2*ei : 2*ei+1;
}
static inline int e_t2c(const oriented_edge &oe){
    return e_t2c(oe.e(), oe.o());
}

template< class bSubNode >
void bMesh< bSubNode >::MakeBasicMesh(bInputFile& infile)
{
    int i, numpts;                      // no. of points in mesh
    double x, y, z;
    int b;
    char pointFilenm[80];            // name of file containing (x,y,z,b) data
    ifstream pointfile;              // the file (stream) itself

    //Read Points
    infile.ReadItem( pointFilenm, "POINTFILENAME" );
    pointfile.open( pointFilenm );
    if ( !pointfile.good() ){
        cout << "\nPoint file name: '" << pointFilenm << "' not found\n";
        exit(1);
    }

    cout<<"\nReading in '"<<pointFilenm<<"' points file..."<<endl;
    pointfile >> numpts;
    cout << "Number of points to read: " << numpts << endl;
    bSubNode** pointArray = new bSubNode*[numpts];
    bSubNode* curnode;

    //Read point file, make Nodelist
    int streamNodeCount = 0;
    int openNodeCount = 0;
    int closedNodeCount = 0;
    int regularNodeCount = 0;

    for ( i=0; i<numpts; i++ ){

        if ( pointfile.eof() )
            cout << "\nReached end-of-file while reading points.\n" ;
        pointfile >> x >> y >> z >> b;

        curnode = new bSubNode();
        curnode->set3DCoords( x, y, z);
        // Check for legal boundary flagM
        if (b < kNonBoundary || b > kStream) {
            cout << "Boundary error: " << " Point " << i << " "
                 << x << " " << y << " "
                 << z << " " << b << endl;
            cout << "\nInvalid boundary code.\n"<<endl;
            cout << "\n\nExiting Program..."<<endl;
            exit(2);
        }

        if (b == kStream)
            streamNodeCount++;
        if (b == kOpenBoundary)
            openNodeCount++;
        if (b == kClosedBoundary)
            closedNodeCount++;
        if (b == kNonBoundary)
            regularNodeCount++;

        curnode->setBoundaryFlag( b );
        curnode->setID( i );

        if (b == kNonBoundary || b == kStream)
            nodeList.insertAtActiveBack(curnode);
        else if (b == kOpenBoundary)
            nodeList.insertAtBoundFront(curnode);
        else
            nodeList.insertAtBack(curnode);

        pointArray[i] = curnode;
    }
    pointfile.close();

    nnodes = numpts;
    cout << "Read in " << numpts << " points " << endl;
    cout << "Number of regular points " << regularNodeCount << endl;
    cout << "Number of stream points " << streamNodeCount << endl;
    cout << "Number of closed boundary points " << closedNodeCount << endl;
    cout << "Number of open boundary points " << openNodeCount << endl;

    if (openNodeCount == 0) {
        cout << "One outlet node must be present" << endl;
        exit(1);
    }

    // Create the lookup table of nodes indexed by node id
    cout << "Creating Node Lookup Table" << endl;
    bSubNode** NodeTable = new bSubNode*[nodeList.size()];
    std::_List_iterator<bSubNode*> node;
    for (node = nodeList.begin(); node != nodeList.end(); node++)
        NodeTable[(*node)->getID()] = (*node);

    // call mesh generator based on Tipper's method
    int nedgesl;
    int nelem;
    edge* edges(0);
    elem* elems(0);

    // Calculate the triangulation giving edges
    cout << "\nComputing triangulation..." << flush;
    tt_sort_triangulate(nnodes, pointArray, &nedgesl, &edges);

    // Every edge is bidirectional (is this necessary)
    nedges = 2*nedgesl;

    // Create and initialize the edge list by creating two edges
    // (which are complementary, ie share the same endpoints) and then
    // iteratively assigning values to the pair and inserting them onto the
    // back of the edgeList

    cout << "\nCreating edge list..." << endl<<flush;
    bEdge *edge1, *edge2;
    int origBoundary, destBoundary;

    for ( int iedge = 0; iedge < nedgesl; ++iedge ) {

        // Assign values: ID, origin and destination pointers
        edge1 = new bEdge();
        edge2 = new bEdge();
        edge1->setID( e_t2c(iedge,true) );
        edge2->setID( e_t2c(iedge,false) );

        bSubNode *nodPtr1 = NodeTable[pointArray[edges[iedge].from]->getID()];
        edge1->setOriginPtr( nodPtr1 );
        edge2->setDestinationPtr( nodPtr1 );
        origBoundary = nodPtr1->getBoundaryFlag();

        bSubNode *nodPtr2 = NodeTable[pointArray[edges[iedge].to]->getID()];
        edge1->setDestinationPtr( nodPtr2 );
        edge2->setOriginPtr( nodPtr2 );
        destBoundary = nodPtr2->getBoundaryFlag();

        // set the "flowallowed" status (false if either endpoint is a
        // closed boundary, or both are open boundaries)
        // and insert edge pair onto the list --- active
        // part of list if flow is allowed, inactive if not
        if (origBoundary == kClosedBoundary ||
            destBoundary == kClosedBoundary ||
            (origBoundary == kOpenBoundary && destBoundary == kOpenBoundary)) {
            edge1->setFlowAllowed( 0 );
            edge2->setFlowAllowed( 0 );
            edgeList.insertAtBack( edge1 );
            edgeList.insertAtBack( edge2 );
        } else {
            edge1->setFlowAllowed( 1 );
            edge2->setFlowAllowed( 1 );
            edgeList.insertAtActiveBack( edge1 );
            edgeList.insertAtActiveBack( edge2 );
        }
    }

    // set up the lists of edges (spokes) connected to each node
    // connectivity point - sorted point
    int* p2sp = new int[nnodes];
    for (int inodes = 0; inodes != nnodes; ++inodes)
        p2sp[pointArray[inodes]->getID()] = inodes;

    // Calculate the nodes oriented edges
    oriented_edge *oedge;
    cout << "Triangulation building spokes" << endl;
    tt_build_spoke(nnodes, nedgesl, edges, &oedge);

    // Create the lookup table of edges indexed by edge id
    cout << "Creating Edge Lookup Table" << endl;
    bEdge** EdgeTable = new bEdge*[edgeList.size()];
    std::_List_iterator<bEdge*> edge;
    for (edge = edgeList.begin(); edge != edgeList.end(); edge++)
        EdgeTable[(*edge)->getID()] = (*edge);

    // Write the temporary node file to create room to triangulate
    // Contains the spokelist with edge ids and not edge pointers
    //
    ofstream nodeStr;
    nodeStr.open ("nodes.temp", ios::out | ios::binary);
    nodeStr.write((char*) &nnodes, sizeof(int));

    std::vector<int> spokelist;
    int edgid1, edgid, ne;

    cout << "Adding spokes to nodes and writing temporary nodes file" << endl;
    std::_List_iterator<bSubNode*> nodIter = nodeList.begin();
    for (int ii=0; ii<nnodes; ii++) {
        curnode = (*nodIter);
        edgid1 = e_t2c(oedge[p2sp[curnode->getID()]]);
        // Check for wrong edge id
        if (edgid1 < 0 || edgid1 >= (2*nedgesl)) {
            cout << "First Accessing out of range edge id " << edgid1
                 << " for node " << curnode->getID() << endl;
            cout << "No edge pointer" << endl;
            cout << "ID = " << curnode->getID()
                 << " p2sp = " << p2sp[curnode->getID()]
                 << " oedge = " << oedge[p2sp[curnode->getID()]].e()
                 << endl;

            // Good first edge id, collect rest of spokes
        } else {
            bEdge *edgPtr = EdgeTable[edgid1];
            spokelist.push_back(edgPtr->getID());
            edgid = edgPtr->getID();

            // Build rest of spoke list
            oriented_edge& oe_ref = oedge[p2sp[curnode->getID()]];
            oriented_edge ccw_from = oe_ref.ccw_edge_around_from(edges);
            while ( ccw_from.e() != oe_ref.e()) {
                assert(ccw_from.e() < nedgesl);
                ne = e_t2c(ccw_from);
                if (ne < 0 || ne >= (2*nedgesl)) {
                    cout << "Second Accessing out of range edge id " << edgid1
                         << " for node " << curnode->getID() << endl;
                }
                bEdge *edgPtr = EdgeTable[ne];
                spokelist.push_back(edgPtr->getID());
                ccw_from = ccw_from.ccw_edge_around_from(edges);
            }
        }

        // Write out node with spoke list edge ids
        int temp = curnode->getID();
        x = curnode->getX();
        y = curnode->getY();
        z = curnode->getZ();
        b = curnode->getBoundaryFlag();

        nodeStr.write((char*) &temp, sizeof(int));
        nodeStr.write((char*) &x, sizeof(double));
        nodeStr.write((char*) &y, sizeof(double));
        nodeStr.write((char*) &z, sizeof(double));
        nodeStr.write((char*) &b, sizeof(int));

        temp = spokelist.size();
        nodeStr.write((char*)&temp, sizeof(int));
        for (int i = 0; i < spokelist.size(); i++) {
            temp = spokelist[i];
            nodeStr.write((char*)&temp, sizeof(int));
        }
        spokelist.clear();
        nodIter++;
    }
    delete [] oedge;
    delete [] p2sp;

    nodeStr.close();
    delete [] NodeTable;

    // Assign ccwedg connectivity (that is, tell each edge about its neighbor
    // immediately counterclockwise)
    int iedge;
    bEdge *curedge, *ccwedge;
    int ecnt = 0;
    int ccwedgid;

    cout << "Building counter clockwise edge pointer in edges" << endl;
    std::_List_iterator<bEdge*> edgIter = edgeList.begin();
    for (iedge = 0; iedge < nedgesl; ++iedge) {
        curedge = (*edgIter);

        oriented_edge ee1(iedge,true);
        oriented_edge ccw_from = ee1.ccw_edge_around_from(edges);
        ccwedgid = e_t2c(ccw_from);
        ccwedge = EdgeTable[ccwedgid];
        curedge->setCCWEdg( ccwedge );
        edgIter++;
        curedge = (*edgIter);

        oriented_edge ee2(iedge,false);
        oriented_edge ccw_to = ee2.ccw_edge_around_from(edges);
        ccwedgid = e_t2c(ccw_to);
        ccwedge = EdgeTable[ccwedgid];
        curedge->setCCWEdg( ccwedge );
        edgIter++;
    }

    // Write edges to file
    cout << "Write temporary edge file" << endl;
    ofstream edgeStr("edges.temp", ios::out | ios::binary);
    edgeStr.write((char*) &nedges, sizeof(int));
    int id, origID, destID, ccwID;

    for (edgIter = edgeList.begin(); edgIter != edgeList.end(); edgIter++) {
        curedge = (*edgIter);
        id = curedge->getID();
        origID = curedge->getOriginPtrNC()->getID();
        destID = curedge->getDestinationPtrNC()->getID();
        ccwID = curedge->getCCWEdg()->getID();

        edgeStr.write((char*) &id, sizeof(int));
        edgeStr.write((char*) &origID, sizeof(int));
        edgeStr.write((char*) &destID, sizeof(int));
        edgeStr.write((char*) &ccwID, sizeof(int));
    }
    edgeStr.close();

    // Can flush nodeList but can't delete bNodes which are used in pointarray
    nodeList.Flush();

    // Can flush and delete bEdge which is not used again
    std::_List_iterator<bEdge*> eiter;
    for (eiter = edgeList.begin(); eiter != edgeList.end(); eiter++)
        delete (*eiter);
    edgeList.Flush();

    delete [] EdgeTable;

    // Build the element table of triangles
    cout << "Triangulation building element table" << endl;
    tt_build_elem_table(nnodes, pointArray, nedgesl, edges, &nelem, &elems);

    delete [] edges;
    ntri = nelem;

    // Write out triangles
    cout << "Write temporary triangles file" << endl;
    ofstream triangleStr("triangles.temp", ios::out | ios::binary);
    triangleStr.write((char*) &nelem, sizeof(int));
    int nodeID[3], edgeID[3], triID[3];

    for (int ielem = 0; ielem < nelem; ++ielem) {
        nodeID[0] = pointArray[elems[ielem].p1]->getID();
        nodeID[1] = pointArray[elems[ielem].p2]->getID();
        nodeID[2] = pointArray[elems[ielem].p3]->getID();

        edgeID[0] = e_t2c(elems[ielem].e1, elems[ielem].eo1);
        edgeID[1] = e_t2c(elems[ielem].e2, elems[ielem].eo2);
        edgeID[2] = e_t2c(elems[ielem].e3, elems[ielem].eo3);

        triID[0] = ( elems[ielem].t1>=0 ) ? elems[ielem].t1 : -1;
        triID[1] = ( elems[ielem].t2>=0 ) ? elems[ielem].t2 : -1;
        triID[2] = ( elems[ielem].t3>=0 ) ? elems[ielem].t3 : -1;

        triangleStr.write((char*) &ielem, sizeof(int));
        triangleStr.write((char*) &nodeID, 3*sizeof(int));
        triangleStr.write((char*) &triID, 3*sizeof(int));
        triangleStr.write((char*) &edgeID, 3*sizeof(int));
    }
    triangleStr.close();

    // deallocation of Tipper triangulator data structures
    delete [] elems;

    // Have to delete individual bNodes also
    for (int i = 0; i < numpts; i++)
        delete pointArray[i];
    delete [] pointArray;
    cout << "Finished building basic mesh and writing to files" << endl;
}

/**************************************************************************
**
**   bMesh::RefineBasicMesh
**
**   Starting with the basic mesh structure from the triangulator,
**   check the consistency of the mesh and add some voronoi data
**
**************************************************************************/

template< class bSubNode >
void bMesh< bSubNode >::RefineBasicMesh()
{
    double x, y, z;
    int id, bound, nspokes, spoke;
    ifstream nodeStr, edgeStr, triangleStr;

    // Read the temporary node file which includes spoke ids
    cout << "Refining basic mesh, reading temporary node file" << endl;
    nodeStr.open("nodes.temp", ios::in | ios::binary);
    nodeStr.read((char*) &nnodes, sizeof(int));
    assert( nnodes > 0 );
    bSubNode* curnode;

    for (int i = 0; i < nnodes; i++ ) {
        nodeStr.read((char*) &id, sizeof(int));
        nodeStr.read((char*) &x, sizeof(double));
        nodeStr.read((char*) &y, sizeof(double));
        nodeStr.read((char*) &z, sizeof(double));
        nodeStr.read((char*) &bound, sizeof(int));
        nodeStr.read((char*) &nspokes, sizeof(int));

        // Edge ids of spokes can't be used yet because we have no edges
        for (int j = 0; j < nspokes; j++)
            nodeStr.read((char*) &spoke, sizeof(int));

        curnode = new bSubNode();
        curnode->set3DCoords( x, y, z );
        curnode->setID( id );
        assert( bound >= kNonBoundary && bound <= kStream );

        curnode->setBoundaryFlag( bound );
        if ((bound == kNonBoundary) || (bound == kStream) ) {
            nodeList.insertAtActiveBack( curnode );
        } else if ( bound == kOpenBoundary ) {
            nodeList.insertAtBoundFront( curnode );
        } else {
            nodeList.insertAtBack( curnode );       //kClosedBoundary
        }
    }
    nodeStr.close();

    // Create the lookup table of nodes indexed by node id
    bSubNode** NodeTable = new bSubNode*[nnodes];
    std::_List_iterator<bSubNode*> node;
    for (node = nodeList.begin(); node != nodeList.end(); node++)
        NodeTable[(*node)->getID()] = (*node);

    for (int i = 0; i < nnodes; i++)
        if (NodeTable[i] == 0)
            cout << "node " << i << " is missing" << endl;

    // Read the edge file and create all edges using NodeTable to locate nodes
    cout << "Refining basic mesh, reading temporary edge file" << endl;
    edgeStr.open("edges.temp", ios::in | ios::binary);
    edgeStr.read((char*) &nedges, sizeof(int));
    assert( nedges > 0 );
    int origID, destID, ccwID;
    int origBoundary, destBoundary;
    bEdge* curedge;

    for (int i = 0; i < nedges; i++) {
        edgeStr.read((char*) &id, sizeof(int));
        edgeStr.read((char*) &origID, sizeof(int));
        edgeStr.read((char*) &destID, sizeof(int));
        edgeStr.read((char*) &ccwID, sizeof(int));

        // Assign values: ID, origin and destination pointers
        curedge = new bEdge();
        curedge->setID( id );

        bSubNode *nodPtr1 = NodeTable[ origID ];
        curedge->setOriginPtr( nodPtr1 );
        origBoundary = (*nodPtr1).getBoundaryFlag();

        bSubNode *nodPtr2 = NodeTable[ destID ];
        curedge->setDestinationPtr( nodPtr2 );
        destBoundary = (*nodPtr2).getBoundaryFlag();

        // set the "flowallowed" status (false if either endpoint is a
        // closed boundary, or both are open boundaries)
        // and insert edge pair onto the list --- active
        // part of list if flow is allowed, inactive if not

        if (origBoundary == kClosedBoundary ||
            destBoundary == kClosedBoundary ||
            (origBoundary == kOpenBoundary && destBoundary == kOpenBoundary) ) {
            curedge->setFlowAllowed( 0 );
            edgeList.insertAtBack( curedge );
        } else {
            curedge->setFlowAllowed( 1 );
            edgeList.insertAtActiveBack( curedge );
        }
    }
    edgeStr.close();

    // Create the lookup table of edges indexed by edge id
    bEdge** EdgeTable = new bEdge*[nedges];
    std::_List_iterator<bEdge*> edge;
    for (edge = edgeList.begin(); edge != edgeList.end(); edge++)
        EdgeTable[(*edge)->getID()] = (*edge);

    // Read in the temporary node file and set up list of edges (spokes)
    cout << "Refining basic mesh, reading temporary node file" << endl;
    nodeStr.open("nodes.temp", ios::in | ios::binary);
    nodeStr.read((char*) &nnodes, sizeof(int));
    std::vector<int> spokeList;
    int first1, next1;
    bEdge* edgPtr;

    std::_List_iterator<bSubNode*> nodIter = nodeList.begin();
    for (nodIter = nodeList.begin(); nodIter != nodeList.end(); nodIter++) {
        bSubNode * curnode;
        curnode = (*nodIter);
        nodeStr.read((char*) &id, sizeof(int));
        nodeStr.read((char*) &x, sizeof(double));
        nodeStr.read((char*) &y, sizeof(double));
        nodeStr.read((char*) &z, sizeof(double));
        nodeStr.read((char*) &bound, sizeof(int));
        nodeStr.read((char*) &nspokes, sizeof(int));

        // Now we can process the spokes because we have the edges
        for (int j = 0; j < nspokes; j++) {
            nodeStr.read((char*)&spoke, sizeof(int));
            spokeList.push_back(spoke);
        }

        // First spoke coming off this node
        bEdge *edge = EdgeTable[spokeList[0]];
        curnode->insertBackSpokeList( edge );
        curnode->setEdg( edge );

        for (int k = 1; k < spokeList.size(); k++) {
            next1 = spokeList[k];
            if (next1 >= nedges) {
                cerr << "Warning: edge " << spokeList[0]
                     << " has non-existant ccw edge " << next1 << endl;
                cerr << "This is likely to be a problem in the edge input file"
                     << endl;
            }
            edge = EdgeTable[next1];
            curnode->insertBackSpokeList( edge );
        }
        spokeList.clear();
    }
    nodeStr.close();

    // Reread edge file because we can assign ccw connectivity
    // because we know all the edges now
    cout << "Refining basic mesh, reading temporary edge file" << endl;
    edgeStr.open("edges.temp", ios::in | ios::binary);
    edgeStr.read((char*) &nedges, sizeof(int));
    assert( nedges > 0 );

    for (int i = 0; i < nedges; i++) {
        edgeStr.read((char*) &id, sizeof(int));
        edgeStr.read((char*) &origID, sizeof(int));
        edgeStr.read((char*) &destID, sizeof(int));
        edgeStr.read((char*) &ccwID, sizeof(int));

        if (ccwID >= 0 && ccwID < nedges)
            EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
        else
            cout << "Edge " << id << " has bad ccwedge "
                 << ccwID << " should be 0 to " << nedges << endl;
    }
    edgeStr.close();

    // Calculate the counter clockwise edges
    cout << "Calculating counterclockwise edges" << endl;
    MakeCCWEdges();

    // Read the temporary triangle file and fill in nodes and edges
    int nodeID[3], edgeID[3], triID[3];
    triangleStr.open("triangles.temp", ios::in | ios::binary);
    triangleStr.read((char*) &ntri, sizeof(int));
    assert( ntri > 0 );

    // Set up the triangle connectivity
    bTriangle* curtriangle;
    for (int i = 0; i < ntri; i++ ) {
        triangleStr.read((char*) &id, sizeof(int));
        triangleStr.read((char*) &nodeID, 3*sizeof(int));
        triangleStr.read((char*) &triID, 3*sizeof(int));
        triangleStr.read((char*) &edgeID, 3*sizeof(int));

        curtriangle = new bTriangle();
        curtriangle->setID(id);
        curtriangle->setPPtr( 0, NodeTable[ nodeID[0] ] );
        curtriangle->setPPtr( 1, NodeTable[ nodeID[1] ] );
        curtriangle->setPPtr( 2, NodeTable[ nodeID[2] ] );
        curtriangle->setEPtr( 0, EdgeTable[ edgeID[0] ] );
        curtriangle->setEPtr( 1, EdgeTable[ edgeID[1] ] );
        curtriangle->setEPtr( 2, EdgeTable[ edgeID[2] ] );
        triList.push_back( curtriangle );
    }
    triangleStr.close();

    // Don't need the indexed node or edge lookup tables any more
    delete [] NodeTable;
    delete [] EdgeTable;

    // Create lookup table for triangles
    std::_List_iterator<bTriangle*> triIter;
    bTriangle** TriTable = new bTriangle*[ntri];
    for (triIter = triList.begin(); triIter != triList.end(); triIter++)
        TriTable[(*triIter)->getID()] = (*triIter);

    // Read in and set adjacent triangles for each triangle
    triangleStr.open("triangles.temp", ios::in | ios::binary);
    triangleStr.read((char*) &ntri, sizeof(int));
    assert( ntri > 0 );

    for (triIter = triList.begin(); triIter != triList.end(); triIter++) {
        triangleStr.read((char*) &id, sizeof(int));
        triangleStr.read((char*) &nodeID, 3*sizeof(int));
        triangleStr.read((char*) &triID, 3*sizeof(int));
        triangleStr.read((char*) &edgeID, 3*sizeof(int));

        curtriangle = ( triID[0] >= 0 ) ? TriTable[ triID[0] ] : 0;
        TriTable[id]->setTPtr( 0, curtriangle );
        curtriangle = ( triID[1] >= 0 ) ? TriTable[ triID[1] ] : 0;
        TriTable[id]->setTPtr( 1, curtriangle );
        curtriangle = ( triID[2] >= 0 ) ? TriTable[ triID[2] ] : 0;
        TriTable[id]->setTPtr( 2, curtriangle );
    }
    triangleStr.close();

    system("rm nodes.temp");
    system("rm edges.temp");
    system("rm triangles.temp");

    // Don't need indexed triangle lookup table any more
    delete [] TriTable;

    // assertions
    assert( nodeList.size() == nnodes );
    assert( edgeList.size() == nedges );
    assert( triList.size() == ntri );
    cout << nnodes << " nodes " << nedges << " edges "
         << ntri << " elements" << endl;

    // Calculate voronoi vertices, edge lengths and areas
    cout << "Calculate voronoi vertices" << endl;
    setVoronoiVertices();
    cout << "Calculate voronoi edge lengths" << endl;
    CalcVoronoiEdgeLengths();
    cout << "Calculate voronoi areas" << endl;
    CalcVAreas();

    // Check the consistency of the mesh
    cout << "Check consistency" << endl;
    CheckMeshConsistency();
}


//=========================================================================
//
//
//                  Section 9: bMesh:: CheckMeshConsistency( )
// 					ChangePointOrder( )
//
//=========================================================================

/*****************************************************************************
**
**  bMesh::CheckMeshConsistency
**
**  Performs a series of tests to make sure the mesh connectivity is correct.
**  Should be called immediately after reading in a user-defined mesh
**
**  The consistency checks include the following:
**
**  1) Each edge:
**     - Has valid origin and destination pointers
**     - Has a valid counter-clockwise edge, which shares the same origin but
**       not the same destination
**     - Is paired with its complement in the list
**
**  2) Each node:
**     - Points to a valid edge which has the node as its origin
**     - If the node is not a boundary, it has at least one neighbor that
**       is not a closed boundary (unless boundaryCheckFlag is false).
**     - Has a consistent spoke list (ie, you can go around the spokes and
**       get back to where you started)
**
**  3) Each triangle:
**     - Has 3 valid points and edges
**     - Each edge Ei has Pi as its origin and P((i+2)%3) as its
**       destination
**     - If an opposite triangle Ti exists, points P((i+1)%3) and
**       P((i+2)%3) are the same as points PO((n+2)%3) and PO((n+1)%3) in
**       the opposite triangle, where PO denotes a point in the opposite
**       triangle and n is the vertex ID (0, 1, or 2) of the point in the
**       opposite triangle that is opposite from the shared face.
**     - If an opposite triange Ti does not exist, points P((i+1)%3) and
**       and P((i+2)%3) should both be boundary points.
**
**      Parameters:  boundaryCheckFlag -- defaults to true; if false,
**                                        node connection to open node or
**                                        open boundary isn't tested
**
*****************************************************************************/

template<class bSubNode>
void bMesh< bSubNode >::CheckMeshConsistency( int boundaryCheckFlag )
{
    std::_List_iterator<bSubNode*> nodIter;
    std::_List_iterator<bEdge*> edgIter;
    std::_List_iterator<bTriangle*> triIter;

    std::_List_iterator<bEdge*> sIter;
    bNode * cn, * org, * dest;
    bEdge * ce, * cne, * ccwedg;
    bTriangle * ct, * optr;
    int boundary_check_ok, i, nvop;
    int kMaxSpokes = 100;

    // Edges: make sure complementary pairs are together in the list
    // (each pair Ei and Ei+1, for i=0,2,4,...nedges-1, should have the same
    // endpoints but the opposite orientation)

    for (edgIter = edgeList.begin(); edgIter != edgeList.end(); edgIter++) {
        ce = (*edgIter);
        edgIter++;
        cne = (*edgIter);
        if( ce->getOriginPtrNC() != cne->getDestinationPtrNC()
            || ce->getDestinationPtrNC() != cne->getOriginPtrNC() ){
            cerr << "EDGE #" << ce->getID()
                 << " must be followed by its complement in the list\n";
        }
    }

    // Edges: check for valid origin, destination, and ccwedg
    for (edgIter = edgeList.begin(); edgIter != edgeList.end(); edgIter++) {
        ce = (*edgIter);
        if( !(org=ce->getOriginPtrNC() ) ){
            cerr << "EDGE #" << ce->getID()
                 << " does not have a valid origin point\n";
        }
        if( !(dest=ce->getDestinationPtrNC() ) ){
            cerr << "EDGE #" << ce->getID()
                 << " does not have a valid destination point\n";
        }
        if( !(ccwedg=ce->getCCWEdg() ) ){
            cerr << "EDGE #" << ce->getID()
                 << " does not point to a valid counter-clockwise edge\n";
        }
        if( ccwedg->getOriginPtrNC()!=org ){
            cerr << "EDGE #" << ce->getID()
                 << " points to a CCW edge with a different origin\n";
        }
        if( ccwedg->getDestinationPtrNC()==dest ){
            cerr << "EDGE #" << ce->getID()
                 << " points to a CCW edge with the same destination\n";
        }
        if( org==dest ){
            cerr << "EDGE #" << ce->getID()
                 << " has the same origin and destination nodes\n";
        }
    }

    // Nodes: check for valid edg pointer, spoke connectivity, and connection
    // to at least one non-boundary or open boundary node

    for (nodIter = nodeList.begin(); nodIter != nodeList.end(); nodIter++) {
        cn = (*nodIter);
        // edg pointer
        if( !(ce = cn->getEdg()) ){
            cerr << "NODE #" << cn->getID()
                 << " does not point to a valid edge\n";
        }
        if( ce->getOriginPtrNC()!=cn ){
            cerr << "NODE #" << cn->getID()
                 << " points to an edge that has a different origin\n";
        }

        boundary_check_ok = ( cn->getBoundaryFlag()==kNonBoundary &&
                              boundaryCheckFlag ) ? 0 : 1;
        i = 0;
        // Loop around the spokes until we're back at the beginning
        do{

            if( ce->getDestinationPtrNC()->getBoundaryFlag()!=kClosedBoundary )
                boundary_check_ok = 1;  // OK, there's at least one open nbr
            i++;
            if( i>kMaxSpokes ){
                cerr << "NODE #" << cn->getID()
                     << ": infinite loop in spoke connectivity\n";
            }

            // Make sure node is the origin --- and not the destination
            if( ce->getOriginPtrNC()!=cn ){
                cerr << "EDGE #" << ce->getID()
                     << " is in the spoke chain of NODE " << cn->getID()
                     << " but does not have the node as an origin\n";
            }
            if( ce->getDestinationPtrNC()==cn ){
                cerr << "EDGE #" << ce->getID()
                     << " is in the spoke chain of NODE " << cn->getID()
                     << " but has the node as its destination\n";
            }

        } while( (ce=ce->getCCWEdg())!=cn->getEdg() );

        if( !boundary_check_ok ){
            std::vector<double> x(2);
            x = cn->get2DCoords();
            cerr << "NODE #" << cn->getID()
                 <<" ( "<<x[0]<< " , "<<x[1]<<" )"
                 << " is surrounded by closed boundary nodes\n";
        }

        //make sure node coords are consistent with edge endpoint coords:
        std::list<bEdge*>& spokList = cn->getSpokeListNC();
        for (sIter = spokList.begin(); sIter != spokList.end(); sIter++) {
            ce = (*sIter);
            if( ce->getOriginPtrNC()->getX() != cn->getX() ||
                ce->getOriginPtrNC()->getY() != cn->getY() ){
                cerr << "NODE #" << cn->getID()
                     << " coords don't match spoke origin coords\n";
            }
        }

    }

    // Triangles: check for valid points and connectivity

    for (triIter = triList.begin(); triIter != triList.end(); triIter++) {
        ct = (*triIter);
        for( i=0; i<=2; i++ ){
            // Valid point i?
            if( !(cn=ct->pPtr(i)) ){
                cerr << "TRIANGLE #" << ct->getID()
                     << " has an invalid point " << i << endl;
            }
            // Valid edge i?
            if( !(ce=ct->ePtr(i)) ){
                cerr << "TRIANGLE #" << ct->getID()
                     << " has an invalid edge " << i << endl;
                goto error;
            }
            // Edge and point consistency
            if( ce->getOriginPtrNC()!=cn ){
                cerr << "TRIANGLE #" << ct->getID()
                     << ": edge " << i << " does not have point " << i
                     << " as origin\n";
            }
            // changed from (i+1) to (i+2) for "right-hand" format
            if( ce->getDestinationPtrNC()!=ct->pPtr((i+2)%3) ){
                cerr << "TRIANGLE #" << ct->getID()
                     << ": edge " << i << " does not have point " << (i+1)%3
                     << " as destination\n";
            }
            // Opposite triangle: if it exists, check common points
            if( (optr = ct->tPtr(i)) ){
                nvop = optr->nVOp(ct); // Num (0,1,2) of opposite vertex in optr
                if( nvop < 3 ){
                    if( ct->pPtr((i+1)%3) != optr->pPtr((nvop+2)%3)
                        || ct->pPtr((i+2)%3) != optr->pPtr((nvop+1)%3) ){
                        cerr << "TRIANGLE #" << ct->getID()
                             << ": opposite triangle " << i << " does not share nodes "
                             << (ct->pPtr((i+1)%3))->getID() << " and "
                             << (ct->pPtr((i+2)%3))->getID() << endl;
                    }
                }
                else{
                    cerr << "TRIANGLE #" << ct->getID()
                         << ": opposite triangle " << i << ", triangle #"
                         << optr->getID() << ",does not have current tri as neighbor\n";
                }
            }
                // If no opposite triangle, make sure it really is a boundary
            else{
                if( (ct->pPtr((i+1)%3))->getBoundaryFlag()==kNonBoundary
                    || (ct->pPtr((i+2)%3))->getBoundaryFlag()==kNonBoundary )
                {
                    cerr << "TRIANGLE #" << ct->getID()
                         << ": there is no neighboring triangle opposite node "
                         << cn->getID() << " but one (or both) of the other nodes "
                         << "is a non-boundary point."
                         <<"\nX = "<<cn->getX()<<"\tY = "<<cn->getY()
                         <<"\tZ = "<<cn->getZ()<<endl;
                }
            }
        }
    }
    return;

    error:
    cerr<<"Error in mesh consistency." << endl;
}

//=========================================================================
//
//
//             Section 11: bMesh Functions using MeshElements
//
//
//=========================================================================

template< class bSubNode >
void bMesh< bSubNode >::MakeCCWEdges()
{
    std::_List_iterator<bSubNode*> nodIter;
    for(nodIter = nodeList.begin(); nodIter != nodeList.end(); nodIter++)
        (*nodIter)->makeCCWEdges();
}

/*****************************************************************************
**
**  bMesh::setVoronoiVertices
**
**  Each Delaunay triangle is associated with an intersection between
**  three Voronoi cells, called a Voronoi vertex. These Voronoi vertices
**  are used in computing the area of each Voronoi cell. The Voronoi
**  vertex associated with each triangle is the circumcenter of the
**  triangle. This routine finds the Voronoi vertex associated with
**  each triangle by finding the triangle's circumcenter.
**
**    Assumes: correct triangulation with valid edge pointers in each tri.
**    Data members modified: none
**    Other objects modified: Voronoi vertices set for each bEdge
**    Modifications:
**     - reverted to earlier triangle-based computation, from an edge-based
**       computation that takes 3x as long because NE = 3NT. In so doing,
**       the definition of the Voronoi vertex stored in a bEdge is changed
**       to "left-hand", meaning the V. vertex associated with the edge's
**       lefthand triangle (the vertex itself may or may not lie to the left
**       of the edge). 1/98 GT
**     - also moved circumcenter computation into a bTriangle mbr fn.
**     - copied function to bMesh member from tStreamNet member, gt 3/98.
**       Other fns now use "right-hand" definition; this fn may have to
**       be changed.
**
*****************************************************************************/

template <class bSubNode>
void bMesh<bSubNode>::setVoronoiVertices()
{
    std::vector<double> xy(2);

    // Find the Voronoi vertex associated with each Delaunay triangle
    std::_List_iterator<bTriangle*> ct;
    for( ct = triList.begin(); ct != triList.end(); ct++) {
        xy = (*ct)->FindCircumcenter();

        // Assign the Voronoi point as the left-hand point of the three edges
        // associated with the current triangle

        (*ct)->ePtr(0)->setRVtx( xy );
        (*ct)->ePtr(1)->setRVtx( xy );
        (*ct)->ePtr(2)->setRVtx( xy );
    }

}

/**************************************************************************
**
**  bMesh::CalcVoronoiEdgeLengths
**
**  Updates the length of the Voronoi cell edge associated with each
**  triangle edge. Because complementary edges are stored pairwise on
**  the edge list, we can save computation time by only computing the
**  vedglen once for the first of the pair, then assigning it to the
**  second. For boundary triangle edges, the corresponding Voronoi edge
**  is infinitely long, so the calculation is only done for interior
**  (active) edges.
**
**************************************************************************/

template <class bSubNode>
void bMesh<bSubNode>::CalcVoronoiEdgeLengths()
{
    bEdge *ce;
    double vedglen;
    std::_List_iterator<bEdge*> edgIter;

    for (edgIter = edgeList.begin();
         edgIter != edgeList.getLastActive(); edgIter++) {
        ce = (*edgIter);
        vedglen = ce->CalcVEdgLen();     // Compute Voronoi edge length
        ce->setVEdgLen( vedglen );
        edgIter++;
        ce = (*edgIter);
        ce->setVEdgLen( vedglen );       // assign the same edge length.
    }
}

/**************************************************************************
**
**  bMesh::CalcVAreas
**
**  Computes Voronoi area for each active (non-boundary) node in the
**  mesh (Voronoi area is only defined for interior nodes). Accomplishes
**  this by calling ComputeVoronoiArea for each node.
**
**************************************************************************/

template <class bSubNode>
void bMesh<bSubNode>::CalcVAreas()
{
    std::_List_iterator<bSubNode*> nodIter;
    for (nodIter = nodeList.begin();
         nodIter != nodeList.getLastActive(); nodIter++)
        (*nodIter)->ComputeVoronoiArea();
}

//=========================================================================
//
//
//                  Section 15: Mesh Data Debugging Functions
//
//
//=========================================================================

//#ifndef NDEBUG

/*****************************************************************************
**
**      DumpEdges(), DumpSpokes(), DumpTriangles(), DumpNodes(): debugging
**         routines which simply write out information pertaining to the mesh;
**      DumpNodes() calls DumpSpokes for each node;
**      DumpSpokes() takes a pointer to a node as an argument.
**
*****************************************************************************/

template<class bSubNode>
void bMesh<bSubNode>::DumpEdges()
{
    std::_List_iterator<bEdge*> edgIter;
    bEdge *ce;
    bTriangle *ct;
    int tid;
    for (edgIter = edgeList.begin(); edgIter != edgeList.end(); edgIter++) {
        ce = (*edgIter);
        ct = TriWithEdgePtr( ce );
        tid = ( ct != 0 ) ? ct->getID() : -1;
        cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
             << " to " << ce->getDestinationPtrNC()->getID() << "; in tri "
             << tid << " (flw " << ce->getBoundaryFlag() << ")" << endl;
    }
}

template<class bSubNode>
void bMesh<bSubNode>::DumpSpokes( bSubNode *cn )
{
    bEdge *ce;
    std::list<bEdge*>* spokList = cn->getSpokeListNC();
    std::_List_iterator<bEdge*> spokIter;

    cout << "node " << cn->getID() << " with spoke edges " << endl;
    for (spokIter = spokList->begin(); spokIter != spokList->end(); spokIter++) {
        ce = (*spokIter);
        cout << "   " << ce->getID()
             << " from node " << ce->getOriginPtrNC()->getID()
             << " to " << ce->getDestinationPtrNC()->getID() << endl;
    }
}

template<class bSubNode>
void bMesh<bSubNode>::DumpTriangles()
{
    std::_List_iterator<bTriangle*> triIter;
    bTriangle *ct, *nt;
    int tid0, tid1, tid2;
    cout << "triangles:" << endl;
    for (triIter = triList.begin(); triIter != triList.end(); triIter++) {
        ct = (*triIter);
        nt = ct->tPtr(0);
        tid0 = ( nt != 0 ) ? nt->getID() : -1;
        nt = ct->tPtr(1);
        tid1 = ( nt != 0 ) ? nt->getID() : -1;
        nt = ct->tPtr(2);
        tid2 = ( nt != 0 ) ? nt->getID() : -1;
        cout << ct->getID() << " with vertex nodes "
             << ct->pPtr(0)->getID() << ", "
             << ct->pPtr(1)->getID() << ", and "
             << ct->pPtr(2)->getID() << "; edges "
             << ct->ePtr(0)->getID() << ", "
             << ct->ePtr(1)->getID() << ", and "
             << ct->ePtr(2)->getID() << "; nbr triangles "
             << tid0 << ", "
             << tid1 << ", and "
             << tid2 << endl;
    }
}

template<class bSubNode>
void bMesh<bSubNode>::DumpNodes()
{
    std::_List_iterator<bSubNode*> nodIter;
    bSubNode *cn;
    cout << "nodes: " << endl;
    for (nodIter = nodeList.begin(); nodIter != nodeList.end(); nodIter++) {
        cn = (*nodIter);
        cout << " at " << cn->getX() << ", " << cn->getY() << ", " << cn->getZ()
             << "; bndy: " << cn->getBoundaryFlag() << "; ";
        DumpSpokes( cn );
    }
}

template<class bSubNode>
void bMesh<bSubNode>::TellAboutNode(bSubNode *cn)
{
    cout<<cn->getID()
        <<"\t"<<cn->getX()
        <<"\t"<<cn->getY()
        <<"\t"<<cn->getZ()
        <<"\t"<<cn->getBoundaryFlag()<<endl<<flush;
    return;
}

/***************************************************************************
**
** bMesh::WriteBasicMesh() Function
**
** Write the nodes and edges belonging to mesh to a file, structure only
**
***************************************************************************/

template<class bSubNode>
void bMesh<bSubNode>::WriteBasicMesh()
{
    // Write the nodes list
    fstream nodeStr("nodes.basic", ios::out | ios::binary);
    BinaryWrite(nodeStr, (int) nodeList.size());

    bSubNode * curnode;
    std::_List_iterator<bSubNode*> nodIter;

    cout << "Write basic mesh " << nnodes << " nodes" << endl;
    for (nodIter = nodeList.begin(); nodIter != nodeList.end(); nodIter++) {
        curnode = (*nodIter);
        BinaryWrite(nodeStr, curnode->getID());
        BinaryWrite(nodeStr, curnode->getBoundaryFlag());
        BinaryWrite(nodeStr, curnode->getX());
        BinaryWrite(nodeStr, curnode->getY());
        BinaryWrite(nodeStr, curnode->getZ());
        BinaryWrite(nodeStr, curnode->getVArea());
        BinaryWrite(nodeStr, curnode->getEdg()->getID());
    }
    nodeStr.close();

    // Write the edge list
    fstream edgeStr ("edges.basic", ios::out | ios::binary);
    BinaryWrite(edgeStr, (int) edgeList.size());

    std::vector<double> rvtx(2);
    bEdge* curedge;
    std::_List_iterator<bEdge*> edgIter;

    cout << "Write basic mesh " << nedges << " edges" << endl;
    for (edgIter = edgeList.begin(); edgIter != edgeList.end(); edgIter++) {
        curedge = (*edgIter);
        rvtx = curedge->getRVtx();
        BinaryWrite(edgeStr, curedge->getID());
        BinaryWrite(edgeStr, rvtx[0]);
        BinaryWrite(edgeStr, rvtx[1]);
        BinaryWrite(edgeStr, curedge->getVEdgLen());
        BinaryWrite(edgeStr, curedge->getOriginPtrNC()->getID());
        BinaryWrite(edgeStr, curedge->getDestinationPtrNC()->getID());
        BinaryWrite(edgeStr, curedge->getCCWEdg()->getID());
    }
    edgeStr.close();
}

/**************************************************************************
**
**   bMesh::ReadBasicMesh
**
**   Reads mesh which was already checked and enhanced with voronoi
**   information from data files.  Reads in a bNode but fills in a
**   bFlowNode with that information so that bFlowNet can be called.
**
**************************************************************************/

template< class bSubNode >
void bMesh< bSubNode >::ReadBasicMesh()
{
    int id, boundary, firstEdgeID;
    int origID, destID, ccwID, origBoundary, destBoundary;
    double x, y, z, varea, vedglen;
    std::vector<double> rvtx(2);
    fstream nodeStr, edgeStr;

    // Read the node file setting all variables except the pointer to first edge
    // Since that can't be found until edges are read
    //
    nodeStr.open("nodes.basic", ios::in | ios::binary);
    BinaryRead(nodeStr, nnodes);
    assert( nnodes > 0 );
    bSubNode* curnode;

    cout << "Read basic mesh " << nnodes << " nodes (Pass 1)" << endl;
    for (int i = 0; i < nnodes; i++) {
        curnode = new bSubNode();
        BinaryRead(nodeStr, id);
        BinaryRead(nodeStr, boundary);
        BinaryRead(nodeStr, x);
        BinaryRead(nodeStr, y);
        BinaryRead(nodeStr, z);
        BinaryRead(nodeStr, varea);
        BinaryRead(nodeStr, firstEdgeID);
        curnode->setID(id);
        curnode->setBoundaryFlag(boundary);
        curnode->setX(x);
        curnode->setY(y);
        curnode->setZ(z);
        if (varea < 0.0) {
            cout << "Problem setting varea of node " << id
                 << " varea=" << varea << endl;
            curnode->setVArea(0.0);
        } else {
            curnode->setVArea(varea);
        }

        if ((boundary == 0) || (boundary==3))
            nodeList.insertAtActiveBack(curnode);
        else if (boundary == kOpenBoundary)
            nodeList.insertAtBoundFront(curnode);
        else
            nodeList.insertAtBack(curnode);       //kClosedBoundary
    }
    nodeStr.close();

    // Create the lookup table of nodes indexed by node id
    // Used when assigning node pointers to edge origin and destination
    //
    bSubNode** NodeTable = new bSubNode*[nnodes];
    std::_List_iterator<bSubNode*> node;
    for (node = nodeList.begin(); node != nodeList.end(); node++)
        NodeTable[(*node)->getID()] = (*node);

    // Read the edge file setting all variables except ccw edge which can't
    // be found until all edges are read once
    //
    edgeStr.open("edges.basic", ios::in | ios::binary);
    edgeStr.read((char*) &nedges, sizeof(int));
    assert( nedges > 0 );
    bEdge* curedge;

    cout << "Read basic mesh " << nedges << " edges (Pass 1)" << endl;
    for (int i = 0; i < nedges; i++) {
        BinaryRead(edgeStr, id);
        BinaryRead(edgeStr, rvtx[0]);
        BinaryRead(edgeStr, rvtx[1]);
        BinaryRead(edgeStr, vedglen);
        BinaryRead(edgeStr, origID);
        BinaryRead(edgeStr, destID);
        BinaryRead(edgeStr, ccwID);

        curedge = new bEdge();
        curedge->setID(id);
        curedge->setRVtx(rvtx);
        curedge->setVEdgLen(vedglen);

        // Look up the node origin and destination for the given id
        bSubNode *origNode = NodeTable[origID];
        curedge->setOriginPtr(origNode);
        origBoundary = origNode->getBoundaryFlag();

        bSubNode *destNode = NodeTable[destID];
        curedge->setDestinationPtr(destNode);
        destBoundary = destNode->getBoundaryFlag();

        // set the "flowallowed" status (false if either endpoint is a
        // closed boundary, or both are open boundaries)
        // and insert edge pair onto the list --- active
        // part of list if flow is allowed, inactive if not

        if (origBoundary == kClosedBoundary ||
            destBoundary == kClosedBoundary ||
            (origBoundary == kOpenBoundary && destBoundary == kOpenBoundary) ) {
            curedge->setFlowAllowed(0);
            edgeList.insertAtBack(curedge);
        }
        else {
            curedge->setFlowAllowed(1);
            edgeList.insertAtActiveBack(curedge);
        }
    }
    edgeStr.close();

    // Create the lookup table of edges indexed by edge id
    bEdge** EdgeTable = new bEdge*[nedges];
    std::_List_iterator<bEdge*> edge;
    for (edge = edgeList.begin(); edge != edgeList.end(); edge++)
        EdgeTable[(*edge)->getID()] = (*edge);

    // Reread nodes file to get the id of the first edge which we now
    // can get a pointer to because of EdgeTable
    //
    nodeStr.open("nodes.basic", ios::in | ios::binary);
    BinaryRead(nodeStr, nnodes);

    cout << "Read basic mesh " << nnodes << " nodes (Pass 2)" << endl;
    for (int i = 0; i < nnodes; i++) {
        BinaryRead(nodeStr, id);
        BinaryRead(nodeStr, boundary);
        BinaryRead(nodeStr, x);
        BinaryRead(nodeStr, y);
        BinaryRead(nodeStr, z);
        BinaryRead(nodeStr, varea);
        BinaryRead(nodeStr, firstEdgeID);

        // Index of the first spoke coming off this node
        NodeTable[id]->setEdg(EdgeTable[firstEdgeID]);
    }
    nodeStr.close();

    // Reread edge file to get the pointer to the counter clockwise edge
    // which we can not get a pointer to because of EdgeTable
    edgeStr.open("edges.basic", ios::in | ios::binary);
    BinaryRead(edgeStr, nedges);
    assert( nedges > 0 );

    bEdge *curedg, *ccwedg;

    cout << "Read basic mesh " << nedges << " edges (Pass 2)" << endl;
    for (int i = 0; i < nedges; i++) {
        BinaryRead(edgeStr, id);
        BinaryRead(edgeStr, rvtx[0]);
        BinaryRead(edgeStr, rvtx[1]);
        BinaryRead(edgeStr, vedglen);
        BinaryRead(edgeStr, origID);
        BinaryRead(edgeStr, destID);
        BinaryRead(edgeStr, ccwID);

        if (ccwID >= 0 && ccwID < nedges)
            EdgeTable[id]->setCCWEdg(EdgeTable[ccwID]);
        else
            cout << "Edge " << id << " has bad ccwedge "
                 << ccwID << " should be 0 to " << nedges << endl;
    }
    edgeStr.close();

    system("/bin/rm edges.basic");
    system("/bin/rm nodes.basic");

    delete [] NodeTable;
    delete [] EdgeTable;
}
#endif

//=========================================================================
//
//
//                          End of bMesh.h
//
//
//=========================================================================
