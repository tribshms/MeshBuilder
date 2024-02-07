// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  meshElements.h: Header file for mesh elements bNode, bEdge, and bTriangle. 
**                  Each of these mesh elements is implemented as an object.
**
**  This file contains declarations of the three classes that collectively
**  make up the triangulated mesh. These classes are:
**   - bNode: nodes (ie, the points in the triangulation)
**   - bEdge: directed edges, described by a starting node and an
**            ending node
**   - bTriangle: triangles in the mesh, with each triangle maintaining
**                pointers to its 3 vertex nodes, its 3 neighboring
**                triangles, and its 3 clockwise-oriented edges
**
**  Lists of each of these 3 types of mesh element are maintained by the
**  bMesh class, which implements the mesh and its routines. Connectivity
**  between mesh elements is managed using pointers, as follows:
**   - Each bNode object points to one of its "spokes" (the bEdges that
**     originate at the node). In the current implementation, each
**     bNode also maintains a list of all its spokes. In a future
**     version such lists will only be created when needed by mesh
**     modification routines (to reduce memory overhead).
**   - Each bEdge points to its origin and destination nodes, and to
**     the bEdge that lies counterclockwise relative to the origin node.
**     bEdge objects also contain the coordinates of the the Voronoi
**     vertex that lies on the righthand side of the bEdge. (A Voronoi
**     vertex is the intersection between 3 Voronoi cells, and in a
**     Delaunay triangulation is found at the circumcenter of triangle).
**   - Each bTriangle object points to its 3 vertex nodes, its 3
**     neighboring triangles (or null if no neighboring triangle exists
**     across a given face), and the 3 clockwise-oriented bEdges. The
**     data structure uses the "opposite" numbering scheme, so that
**     triangle node 1 represents the vertex that is opposite to
**     neighboring triangle 1, and so on. Node 1 is also the origin for
**     edge 1, etc.
**
**************************************************************************/

#ifndef BMESHELEMENTS_H
#define BMESHELEMENTS_H

#include "src/Headers/Definitions.h"
#include "src/Mathutil/geometry.h"
#include "src/bMeshList/bMeshList.h"
#include "src/Headers/globalFns.h"

#ifdef ALPHA_64
  #include <iostream.h>
  #include <math.h> 
  #include <vector.h>
#elif defined LINUX_32
  #include <iostream>
  #include <cmath> 
  #include <vector>
#elif defined WIN
  #include <iostream.h>
  #include <math.h> 
  #include <vector.h>
#else 
  #include <iostream.h>
  #include <math.h> 
  #include <vector.h>
#endif

using namespace std;

class bEdge;

//=========================================================================
//
//
//                  Section 1: bNode() Class Declaration
//
//
//=========================================================================

/**************************************************************************
**  
**  bNode() 
**
**  bNodes are the points in a Delaunay triangulation, and their data
**  include x and y coordinates, a z value (which could be elevation or
**  some other variable), an ID, the point's Voronoi area and its reciprocal,
**  and a pointer to one of its "spokes" (edges that originate at the node).
**  Because each spoke points to its counter-clockwise neighbor, bNode
**  objects only need to point to one of their spokes. However, for
**  convenience, the current implementation of bNode also contains a list
**  of pointers to all spokes (a "spoke list"). In the
**  future, to conserve memory usage, these spoke lists will only be
**  allocated when needed by various mesh modification routines.
**
**  Other bNode variables include the area of the corresponding Voronoi
**  cell, an ID number, and a flag indicating the node's boundary status.
**  Possible boundary status codes are non-boundary (mesh interior), closed
**  (outer boundary, not counted as part of the solution domain), and
**  open (boundary node not part of the solution domain but representing
**  a valid exit point for mass or energy flows). Note that these boundary
**  codes are used by bMesh to segregate nodes according to whether they
**  are boundary or non-boundary points. Note also that while all hull
**  points must be boundaries, interior points do not necessarily have to
**  be flagged as non-boundaries (e.g., one could include an "island" of
**  boundary points in the interior of a mesh if needed).
**
**************************************************************************/

class bNode
{
  friend ostream &operator<<( ostream &, bNode & );
  friend istream &operator>>( istream &, bNode & );

public:
  bNode();                                   // default constructor
  bNode( const bNode & );                    // copy constructor
  virtual ~bNode();                          // destructor

  const bNode &operator=( const bNode & );   // assignment operator
  int operator<( const bNode & );            // comparison operator for sort

  std::vector<double> get3DCoords() const; // returns x,y,z
  std::vector<double> get2DCoords() const; // returns x,y
  int getID() const;                  // returns ID number
  double getX() const;                // returns x coord
  double getY() const;                // returns y coord
  double getZ() const;                // returns z value
  double getVArea() const;            // returns Voronoi area
  double getVArea_Rcp() const;        // returns 1/Voronoi area
  int getBoundaryFlag() const;        // returns boundary code
  bEdge * getEdg();                   // returns ptr to one spoke
  void getVoronoiVertexList( std::list<Point2D> * );    
                                      // Returns list of V vertices
  void getVoronoiVertexXYZList( std::list<Point3D> * ); 
                                      // As above plus interp z

  const std::list<bEdge*> &getSpokeList() const;  // returns ref to spokelist
  std::list<bEdge*> &getSpokeListNC();            // returns non-const ref "
  void insertFrontSpokeList( bEdge * ); // adds edge ptr to front of spokelist
  void insertBackSpokeList( bEdge * );  // adds edge ptr to back of spokelist
  void makeWheel();                	// makes spokelist circular
  void makeCCWEdges();             	// sets up CCWEdg connect from spokelst
  
  void setID( int );              // sets ID number
  void setX( double );            // sets x coord
  void setY( double );            // sets y coord
  void setZ( double );            // sets z value
  void ChangeZ( double );         // adds or sub from the z value
  void setVArea( double );        // sets Voronoi area
  void set2DCoords( double, double );    // sets x and y values
  void set3DCoords( double, double, double ); // sets x, y, and z values
  void setBoundaryFlag( int );   // sets boundary status flag
  void setEdg( bEdge * );         // sets ptr to one spoke

  double Dist( bNode *, bNode * ); // distance from node to line (node1,node2)
  bEdge *EdgToNod( bNode * );      // finds spoke connected to given node
  double ComputeVoronoiArea();     // calculates node's Voronoi area
  void ConvertToClosedBoundary();  // makes node a closed bdy & updates edges
  void AttachFirstSpoke( bEdge * ); // welcomes first spoke (sets edg)
  void FlushSpokeList();           // Empty spokes
  virtual void WarnSpokeLeaving( bEdge *); 
                                   // signals node that spoke is being deleted
  virtual void InitializeNode();   // used when new nodes are created

#ifdef PARALLEL_TRIBS
  void setResIndex(int rind);          //tResample (Parallel only)
  int getResIndex();
#endif

#ifndef NDEBUG
   void TellAll();  			// Debugging routine 
#endif
 
protected:
  int id;           			// ID number
  int boundary;     			// Boundary status code
  double x;				// x coordinate
  double y;				// y coordinate
  double z;				// z value
  double varea;     			// Voronoi cell area
 
  bEdge * edg;      			// Ptr to one edge
  std::list<bEdge*> spokeList; 		// list of connected edges (spokes)

#ifdef PARALLEL_TRIBS
  int resIndex;               // tResample index (parallel only)
#endif
};


class bFlowNode : public bNode
{
public:
  bFlowNode();                                   // default constructor
  bFlowNode( const bFlowNode & );                  // copy constructor
  virtual ~bFlowNode() {}                        // destructor

  bEdge* getFlowEdg();
  int getFloodStatus();
  bFlowNode* getDownstrmNbr();
  bFlowNode* getStreamNode();
  double getHillPath();
  double getTTime();
  double getSrf();
  double getHsrf();
  double getSbsrf();
  double getPsrf();
  double getSatsrf();
  double getRain();
  double getSoilMoistureSC();
  double getRootMoistureSC();
  double getSoilMoistureUNSC();
  double getNwtNew();
  double getEvapoTrans();
  double getContrArea();
  double getCurvature();
  double getStreamPath();
  int getTracer();
  int getReach();

  double getCentroidX();                                // Geometric Methods
  double getCentroidY();
  int  polyCentroid(double *, double *, int, double *, double *, double *);

  void setFlowEdg(bEdge *);
  void setFloodStatus(int status);
  void setStreamNode(bFlowNode*);
  void setStreamPath(double);
  void setHillPath(double);
  void setTTime(double);
  void setNwtNew(double);
  void setContrArea(double);
  void setCurvature(double);
  void setTracer(int);
  void setReach(int);

  void addContrArea(double value);

  void ActivateSortTracer();
  void DeactivateTracer();
  void MoveSortTracerDownstream();
  void AddTracer();
  int NoMoreTracers();

private:
  bEdge* flowedge;
  bFlowNode* StreamPtr;

  int flood;
  double hillpath;
  double traveltime;
  double srf;
  double hsrf;
  double psrf;
  double satsrf;
  double sbsrf;
  double EvapoTranspiration;
  double SoilMoistureSC;
  double SoilMoistureUNSC;
  double RootMoistureSC;
  double ContrArea;
  double NwtNew;
  double Rain;
  double Curvature;
  double streampath;
  int tracer;
  int reach;
  double xC;
  double yC;
};



//=========================================================================
//
//
//                  Section 2: bEdge() Class Declaration
//
//
//=========================================================================

/***************************************************************************
**  
**  bEdge()
**
**  bEdge objects represent the directed edges in a Delaunay triangulation
**  of bNode objects. "Directed" means that the edge has directionality
**  from one point to the other; one is the origin and the other the  
**  the destination. In addition to pointing to its origin and destination
**  nodes, each bEdge points to the bEdge that shares the same origin and
**  lies immediately counter-clockwise. This makes it possible to obtain,
**  given one bEdge, all of the bEdges connected to a given origin node.
**  Other data maintained by a bEdge include its length, slope (if
**  applicable), a boundary flag, the coordinates of the Voronoi vertex
**  vertex associated with the right-hand triangle, and the length of
**  the corresponding Voronoi cell edge.
**
**  The boundary status of a bEdge object depends on the boundary status
**  of the two nodes to which it is connected: if either node is a closed
**  boundary, the edge's is a "no flow" (boundary) edge; otherwise it is a 
**  "flow allowed" (non-boundary) edge.
**
**  Note that an edge's slope is defined as the (Zo - Zd)/L, where Zo and
**  Zd are the z values of the origin and destination nodes, respectively,
**  and L is the edge's (projected) length.
**
***************************************************************************/

class bEdge
{
  friend ostream &operator<<( ostream &, const bEdge & );

public:

  bEdge();                		// default constructor
  bEdge( const bEdge & ); 		// copy constructor

  const bEdge &operator=( const bEdge & );  // assignment operator
  int getID() const;            	    // returns ID number
  int getBoundaryFlag() const;  	    // returns boundary status
  double getLength();          // returns edge's length 
  double getSlope();           // slope from org to dest nodes 
  double getOrgZ();                  // returns origin's z value
  double getDestZ();                 // returns destination's z value
  const bNode *getOriginPtr() const; // returns ptr to origin node 
  const bNode *getDestinationPtr() const;   // returns ptr to dest node 
  bNode *getOriginPtrNC();      // returns ptr to origin node (non-const)
  bNode *getDestinationPtrNC(); // returns ptr to destination node (non-const)
  bEdge * getCCWEdg();          // returns ptr to counter-clockwise neighbor
  std::vector<double> getRVtx() const; // returns Voronoi vertex for RH triangle
  double getVEdgLen();    // returns length of assoc'd Voronoi cell edge
  int FlowAllowed();            // returns boundary status ("flow allowed")

  void setID( int );                 // sets ID number
  void setOriginPtr( bNode * );      // sets origin ptr
  void setDestinationPtr( bNode * ); // sets destination ptr
  void setFlowAllowed( int );        // sets boundary code
  void setCCWEdg( bEdge * edg );     // sets ptr to counter-clockwise neighbor
  void setRVtx(std::vector<double>); // sets coords of Voronoi vertex RH tri
  void setVEdgLen( double ); // sets length of corresponding Voronoi edge

  double CalcLength();               // computes & sets length
  double CalcSlope();                // computes & sets slope
  double CalcVEdgLen();      // computes, sets & returns length of V cell edg

  bEdge * FindComplement();  // returns ptr to edge's complement
  void WelcomeCCWNeighbor( bEdge * );  // Adds another edge ccw from this edge

#ifndef NDEBUG
  void TellCoords();  // debug routine that reports edge coordinates
#endif
  
private:
  int id;          	 // ID number
  bool flowAllowed; 	 // boundary flag, false when org & dest = closed bds 
  std::vector<double> rvtx; // (x,y) coords of Voronoi vertex in RH triangle
  double vedglen;        // length of Voronoi edge shared by org & dest cells

  bNode *org;		// ptrs to origin node
  bNode *dest;		// ptrs to destination node
  bEdge *ccwedg;         // ptr to counter-clockwise edge w/ same origin 
};

//=========================================================================
//
//
//                  Section 3: bTriangle() Class Declaration
//
//
//=========================================================================

/**************************************************************************
**  
**  bTriangle()
**
**  bTriangles are the Delaunay triangles that form the triangulated mesh.
**  Each bTriangle maintains pointers to its three nodes (vertices), three
**  adjacent triangles (if they exist), and to the three counter-clockwise
**  directed edges. For example, a bTriangle containing points a, b, and c
**  (where the order abc is counter-clockwise) would point to bNodes a,b,c
**  and bEdges a->b, b->c, and c->a, as well as to the bTriangles that share
**  sides ab, bc, and ac.
** 
**  Numbering convention: 
**   - points p0,p1,p2 are in counter-clockwise order
**   - adjacent triangle t0 lies opposite point p0
**   - directed edge e0 has points p0->p2
**
**************************************************************************/

class bTriangle
{
  friend ostream &operator<<( ostream &, const bTriangle & );
  friend istream &operator>>( istream &, bTriangle & );

public:
  bTriangle();                    	    // default constructor
  bTriangle( const bTriangle & ); 	    // copy constructor
  bTriangle( int, bNode*, bNode*, bNode* ); // copy constructor

  const bTriangle &operator=( const bTriangle & ); // assignment operator
  int getID() const;                 // returns ID number
  bNode *pPtr( int );                // returns ptr to given vertex (0,1, or 2)
  bEdge *ePtr( int );                // returns ptr to given clockwise edge
  bTriangle *tPtr( int );            // returns ptr to given neighboring tri
  void setID( int );                 // sets ID number
  void setPPtr( int, bNode * );      // sets ptr to given vertex
  void setEPtr( int, bEdge * );      // sets ptr to given clockwise edge
  void setTPtr( int, bTriangle * );  // sets ptr to given neighboring tri
  int nVOp( bTriangle * );    // returns side # (0,1 or 2) of nbr triangle
  int nVtx( bNode * );        // returns vertex # (0,1 or 2) of given node
  std::vector<double> FindCircumcenter(); // computes & returns tri's circumcenter

#ifndef NDEBUG
  void TellAll();  // debugging routine
#endif

private:
  int id;          // triangle ID number

  bNode *p[3];     // ptrs to 3 nodes (vertices)
  bEdge *e[3];     // ptrs to 3 clockwise-oriented edges
  bTriangle *t[3]; // ptrs to 3 neighboring triangles (or 0 if no nbr exists)
};


//=========================================================================
//
//
//                  Section 4. Inline Functions for bNode
//
//
//=========================================================================

/***********************************************************************
**
**  Constructors & destructors:
**
**  Default:  initializes values to zero.
**  Copy:  copies all values and makes duplicate spoke list
**
***********************************************************************/

//Default constructor
inline bNode::bNode() :  
  id(0), x(0.), y(0.), z(0.),
  varea(0.),
  boundary(0), edg(0), 
  spokeList()
{
}

//Copy constructor
inline bNode::bNode( const bNode &original ) :
  id(original.id),
  x(original.x), y(original.y), z(original.z),
  varea(original.varea),
  boundary(original.boundary), edg(original.edg),
  spokeList()
{
  if( &(original.spokeList) != 0 ){
    std::list<bEdge*>::const_iterator siter = original.spokeList.begin();
    for( int i=0; i<original.spokeList.size(); i++ ) {
      insertBackSpokeList((*siter));
      siter++;
    }
  }    
}


//Copy constructor
inline bFlowNode::bFlowNode( const bFlowNode &original ) : bNode(original)
{
}


//Destructor
inline bNode::~bNode()
{
}

/***********************************************************************
**
**  Overloaded operators:
**
**    assignment: copies all values (spokelist's assignment operator
**                creates duplicate copy of list)
**    right shift: takes input for x, y and z values from input stream
**    left shift: sends the following data to the output stream:
**                node ID, x, y, z values, and IDs of neighboring nodes,
**                which are obtained through the spokelist.
**
***********************************************************************/

inline const bNode &bNode::operator=( const bNode &right )
{
   if ( &right != this ) {
      id = right.id;
      x = right.x;
      y = right.y;
      z = right.z;
      boundary = right.boundary;
      varea = right.varea;
      edg = right.edg;
      spokeList = right.spokeList;
   }
   return *this;
}

inline int bNode::operator<(const bNode& node)
{
   if (x == node.x)
      return y < node.y;
   return x < node.x;
}

//right shift
inline istream &operator>>( istream &input, bNode &node )
{
   cout << "x y z:" << endl;
   input >> node.x >> node.y >> node.z;
   return input;
}

//left shift
inline ostream &operator<<( ostream &output, bNode &node )
{
   output << node.id 
          << " " << node.x << " " << node.y << " " << node.z
          << " " << node.boundary << " " << node.varea 
          << " " << node.spokeList.size() << endl;

   std::list<bEdge*>::iterator spokIter;
   for (spokIter = node.spokeList.begin();
        spokIter != node.spokeList.end(); spokIter++) {
      output << (*spokIter)->getID() << endl;
   }
   return output;
}


/***********************************************************************
**
**  bNode "get" functions:
**
**  get3DCoords - returns x, y, z as a 3-element array
**  get2DCoords - returns x & y coords as a 2-element array
**  getID - returns ID #
**  getX - returns node's x coord
**  getY - returns node's y coord
**  getZ - returns node's z value
**  getVeg - returns node's vegetation cover
**  getVArea - returns Voronoi area
**  getVArea_Rcp - returns 1 / Voronoi area
**  getBoundaryFlag - returns boundary code
**  getEdg - returns pointer to one spoke
**  getSpokeList - returns const reference to spoke list
**  getSpokeListNC - returns non-const reference to spoke list
**  getFirstSpokeNode - returns ptr to 1st spokelist item 
**                      (return type tPtrLisbNode *)
**  getFirstSpokeNodeNC - non-const version of the above
**  getNextSpokeNode - returns ptr to next spokelist item that follows
**                       prevedg
**  getNextSpokeNodeNC - non-const version of the above
**
***********************************************************************/

inline std::vector<double> bNode::get3DCoords() const
{
   std::vector<double> xyz(3);
   xyz[0] = x;
   xyz[1] = y;
   xyz[2] = z;
   return xyz;
}

inline std::vector<double> bNode::get2DCoords() const
{
   std::vector<double> xy(2);
   xy[0] = x;
   xy[1] = y;
   return xy;
}

inline double bNode::getX() const {return x;}
inline double bNode::getY() const {return y;}
inline double bNode::getZ() const {return z;}
inline int bNode::getID() const {return id;}                  
inline double bNode::getVArea() const {return varea;}            
inline double bNode::getVArea_Rcp() const {return 1.0 / varea;}    
inline int bNode::getBoundaryFlag() const {return boundary;}     
inline bEdge * bNode::getEdg() {return edg;}

inline const std::list<bEdge*> &                                    
bNode::getSpokeList() const {assert( &spokeList != 0 ); return spokeList;}

inline std::list<bEdge*> &                                          
bNode::getSpokeListNC() {assert( &spokeList != 0 ); return spokeList;}

/***********************************************************************
**
**  bNode "set" functions:
**
**  setID - sets ID number to val
**  setX - sets x coord to val
**  setY - sets y coord to val
**  setZ - sets z value to val
**  setVeg - sets vegetation cover to val
**  setVArea - sets Voronoi area to val
**  setVArea_Rcp - sets 1/Voronoi area to val
**  setBoundaryFlag - returns boundary code
**  set3DCoords - sets x, y, z to val1, val2, val3
**  set2DCoords - sets x & y to val1 and val2
**  setEdg - sets edge ptr to theEdg
**
***********************************************************************/

inline void bNode::setID( int val ) {id = val;}    
inline void bNode::setX( double val ) {x = val;}   
inline void bNode::setY( double val ) {y = val;}
inline void bNode::setZ( double val ) {z = val;}

inline void bNode::setVArea( double val ){
  assert( val>=0.0 );
  varea = val;
}

inline void bNode::setBoundaryFlag( int val ){
  assert( val>=0 && val<=3 );
  boundary = val;
}

inline void bNode::set2DCoords( double val1, double val2 ){
   setX( val1 );
   setY( val2 );
}

inline void bNode::set3DCoords( double val1, double val2, double val3 ){
   setX( val1 );
   setY( val2 );
   setZ( val3 );
}

inline void bNode::setEdg( bEdge * theEdg ){
   //assert( theEdg > 0 );
   edg = theEdg;
}

/***********************************************************************
**
**  bNode::ChangeZ:  Adds delz to current z value
**
***********************************************************************/

inline void bNode::ChangeZ( double delz ) { z += delz; }   

/*******************************************************************
**
**  bNode::WarnSpokeLeaving( bEdge * edglvingptr )
**
**  This function is called when an edge is being removed from the list.
**  If edg (the edge pointer member of bNode) is pointing to the edge
**  which will be removed, this edg must be updated.
**
**  edglvingptr is as it says, a pointer to the edge which will be
**  removed.
**
*******************************************************************/

inline void bNode::WarnSpokeLeaving( bEdge * edglvingptr ){
   if( edglvingptr == edg )
       edg = edg->getCCWEdg();
}

/**********************************************************************
**
**  bNode::InitializeNode()
**
**  A virtual function.
**  This functions doesn't do anything here, only in inherited classes.
**  Used for initializing things in newly created nodes that are set up
**  for the rest of the nodes when the mesh is created.
**
**********************************************************************/

inline void bNode::InitializeNode(){}

//=========================================================================
//
//
//                Section 5. Inline Functions for bEdge
//
//
//=========================================================================


/***********************************************************************
**
**  Constructors & destructors:
**
**  Default:  initializes values to zero and makes rvtx a 2-elem array
**  Copy:  copies all values
**
***********************************************************************/

//default constructor
inline bEdge::bEdge() : 
  id(0), flowAllowed(false),
  vedglen(0), rvtx(2), org(0), dest(0), ccwedg(0)
{
}

//copy constructor
inline bEdge::bEdge( const bEdge &original ) :
  id(original.id), flowAllowed(original.flowAllowed),
  rvtx(original.rvtx), vedglen(original.vedglen),
  org(original.org), dest(original.dest), ccwedg(original.ccwedg)
{
}

/***********************************************************************
**
**  Overloaded operators:
**
**    assignment: copies all values (spokelist's assignment operator
**                creates duplicate copy of list)
**    left shift: sends the following data to the output stream:
**                edge ID, length, slope, and origin and destination IDs
**
***********************************************************************/

inline const bEdge &bEdge::operator=( const bEdge &original )
{
   if( &original != this ){
      id = original.id;
      rvtx[0] = original.rvtx[0];
      rvtx[1] = original.rvtx[1];
      vedglen = original.vedglen;
      org = original.org;
      dest = original.dest;
      ccwedg = original.ccwedg;
      flowAllowed = original.flowAllowed;
   }
   return *this;
}

//left shift
inline ostream &operator<<( ostream &output, const bEdge &edge ){
   output << edge.id << " "
          << " " << edge.org->getID()
          << " " << edge.dest->getID() 
          << " " << (edge.flowAllowed ? 1 : 0)
          << " " << edge.ccwedg->getID()
          << " " << edge.rvtx[0] << " " << edge.rvtx[1] << endl;
   return output;
}

/***********************************************************************
**
**  bEdge "get" functions:
**
**  getID - returns ID #
**  getBoundaryFlag - returns boundary code
**  getLength - returns projectd length
**  getSlope - returns slope
**  getOriginPtr - returns const ptr to origin node
**  getDestinationPtr - returns const ptr to destination node
**  getOriginPtrNC - returns non-const ptr to origin node
**  getDestinationPtrNC - returns non-const ptr to destination node
**  getOrgZ- returns z value of origin node
**  getDestZ - returns z value of destination node
**  getCCWEdg - returns ptr to counterclockwise neighboring edge
**  FlowAllowed - returns the boundary flag, which indicates whether
**                or not the edge is an active flow conduit (which is
**                true as long as neither endpoint is a closed bdy node)
**  getRVtx - returns coordinates of right-hand Voronoi vertex as a
**            2-element array
**  getVEdgLen - returns the length of the corresponding Voronoi edge
**
***********************************************************************/

inline int bEdge::getID() const {return id;}                               

inline int bEdge::getBoundaryFlag() const             
{
   return !( flowAllowed == kFlowAllowed );
} 

inline double bEdge::getLength() {return CalcLength();}  

inline double bEdge::getSlope() {return CalcSlope();} 

inline const bNode *bEdge::getOriginPtr() const {return org;} 

inline const bNode *bEdge::getDestinationPtr() const {return dest;}

inline bNode *bEdge::getOriginPtrNC() {return org;}                

inline bNode *bEdge::getDestinationPtrNC() {return dest;}          

inline double bEdge::getOrgZ() {
   assert( org!=0 );
   return( org->getZ() );
}

inline double bEdge::getDestZ(){  
   assert( dest!=0 );
   return( dest->getZ() );
}

inline bEdge * bEdge::getCCWEdg() {
   return ccwedg;
}

inline int bEdge::FlowAllowed() {
   if (flowAllowed) return 1;
   else return 0;
}

inline std::vector<double> bEdge::getRVtx() const
{
   std::vector<double> x(2);
   x[0] = rvtx[0];
   x[1] = rvtx[1];
   return x;
}

inline double bEdge::getVEdgLen()
{
   return vedglen;
}


/***********************************************************************
**
**  bEdge "set" functions:
**
**  setID - sets ID # to val
**  setLength - sets length to val
**  setSlope - sets slope to slp
**  setOriginPtr - sets origin pointer to ptr (if ptr is nonzero)
**  setDestinationPtr - sets destination pointer to ptr (if nonzero)
**  setFlowAllowed - sets flowAllowed status to val
**  setCCWEdg - sets ptr to counter-clockwise neighbor to edg
**  setRVtx - sets the coordinates of the right-hand Voronoi vertex
**            (ie, the Voronoi vertex at the circumcenter of the RH
**            triangle) to the 1st two elements in arr, which is
**            assumed to be a 2-element array
**  setVEdgLen - sets vedglen to val (vedglen is the length of the
**               corresponding Voronoi cell edge)
**
**  Note: unless otherwise noted, no checking of range or validity is
**        performed in these routines (aside from assert statements)
**
***********************************************************************/

inline void bEdge::setID( int val ) {
   assert( val>=0 );
   id = val;
}           

inline void bEdge::setOriginPtr( bNode * ptr ) {if( ptr != 0 ) org = ptr;}

inline void bEdge::setDestinationPtr( bNode * ptr )
{if( ptr != 0 ) dest = ptr;}

inline void bEdge::setFlowAllowed( int val ){
   assert( val==0 || val==1 );
   flowAllowed = val;
}

inline void bEdge::setCCWEdg( bEdge * edg ){
   //assert( edg > 0 );
   ccwedg = edg;
}

inline void bEdge::setRVtx( std::vector<double> arr )
{
   assert( &arr != 0 );
   rvtx[0] = arr[0];
   rvtx[1] = arr[1];
}

inline void bEdge::setVEdgLen( double val )
{
   if (val < 0.0) cout << "setVEdgLen edge " << getID() << " orig " << getOriginPtr()->getID() << " dest " << getDestinationPtr()->getID() << " val " << val << endl;
   //assert( val>=0.0 );
   vedglen = val;
}

/**************************************************************************
**
**  bEdge::CalcSlope
**
**  Computes the slope of the edge as ( Zorg - Zdest ) / length.
**
**  Returns: the slope
**  Modifies: slope (data mbr)
**  Assumes: length >0; org and dest valid.
**
**************************************************************************/

inline double bEdge::CalcSlope(){ 
   assert( org!=0 );  // Failure = edge has no origin and/or destination node
   assert( dest!=0 );
   double slope;
   slope = ( org->getZ() - dest->getZ() ) / CalcLength();
   return slope;
}


/**************************************************************************
**
**  bEdge::CalcVEdgLen
**
**  Calculates the length of the Voronoi cell edge associated with the
**  current triangle edge. The Voronoi cell edge length is equal to the
**  distance between the Voronoi vertex of the right-hand triangle and
**  the Voronoi vertex of the left-hand triangle. The vertex for the
**  right-hand triangle is stored with rvtx[] (and is assumed to be up to
**  date), and the vertex for the left-hand triangle is stored in the
**  edge's counter-clockwise (left-hand) neighbor.
**
**************************************************************************/

inline double bEdge::CalcVEdgLen(){
   assert( ccwedg!=0 );	
   double dx, dy;
	
   dx = rvtx[0] - ccwedg->rvtx[0];
   dy = rvtx[1] - ccwedg->rvtx[1];
   double vedglen = sqrt( dx*dx + dy*dy );
   return( vedglen );
}


/**************************************************************************
**
**  bEdge::WelcomeCCWNeighbor
**
**  Welcomes a new spoke to the neighborhood! neighbor is a new edge to
**  be inserted counter-clockwise from this one. We point neighbor at
**  the edge we're currently pointing to, and then point ourself to
**  neighbor, thus maintaining the edge connectivity.
**
**************************************************************************/

inline void bEdge::WelcomeCCWNeighbor( bEdge * neighbor ){
   assert( neighbor!=0 );
   assert( neighbor->org == org );
   neighbor->ccwedg = ccwedg;
   ccwedg = neighbor;
}

//=========================================================================
//
//
//                   Section 6. Inline Functions for bTriangle
//
//
//=========================================================================

/***********************************************************************
**
**  Constructors & destructors:
**
**  Default:  initializes node, edge, and triangle ptrs to zero.
**  Copy:  copies all values
**  ID & Vertices: creates a triangle w/ pointers to three vertices,
**                 and sets up edge pointers as well. Does not set
**                 triangle pointers however (these are zero'd).
**
***********************************************************************/

//default
inline bTriangle::bTriangle() : 
  id(-1)
{
   assert( p != 0 && e != 0 && t != 0 );
   for( int i=0; i<3; i++ ){
      p[i] = 0;
      e[i] = 0;
      t[i] = 0;
   }
}

//copy constructor
inline bTriangle::bTriangle( const bTriangle &init ) :
   id(init.id)
{
   assert( p != 0 && e != 0 && t != 0 );
   if( &init != 0 ){
      id = init.id;
      for( int i=0; i<3; i++ ){
         p[i] = init.p[i];
         e[i] = init.e[i];
         t[i] = init.t[i];
      }
   }
}

// construct with id and 3 vertices
inline bTriangle::bTriangle( int num, bNode* n0, bNode* n1, bNode* n2 ) :
  id(num)
{
   //assert( n0 > 0 && n1 > 0 && n2 > 0 );
   p[0] = n0;
   p[1] = n1;
   p[2] = n2;
   setEPtr( 0, n0->EdgToNod( n2 ) );
   setEPtr( 1, n1->EdgToNod( n0 ) );
   setEPtr( 2, n2->EdgToNod( n1 ) );
   t[0] = t[1] = t[2] = 0;
}


/***********************************************************************
**
**  Overloaded operators:
**
**    assignment: copies all values 
**    left shift: sends the following data to the output stream:
**                triangle ID and the IDs of its 3 nodes, clockwise
**                edges, ad neighboring triangles (or -1 if no
**                neighboring triangle exists across a given face)
**    right shift: reads triangle ID and 3 other unspecified IDs from
**                 the input stream (the latter are not currently used
**                 for anything)
**
***********************************************************************/

//overloaded assignment operator
inline const bTriangle &bTriangle::operator=( const bTriangle &init ){
   if( &init != this ){
      id = init.id;
      for( int i=0; i<3; i++ ){
         p[i] = init.p[i];
         e[i] = init.e[i];
         t[i] = init.t[i];
      }
   }
   return *this;
}

//left shift
inline ostream &operator<<( ostream &output, const bTriangle &tri ){
   int i;
   output << tri.id << ":";
   for( i=0; i<3; i++ )
       output << " " << tri.p[i]->getID();
   output << ";";
   for( i=0; i<3; i++ )
       output << " " << tri.e[i]->getID();
   output << ";";
   for( i=0; i<3; i++ ){
      if( tri.t[i] != 0 ) output << " " << tri.t[i]->getID();
      else  output << " -1";
   }
   output << endl;
   return output;
}

inline istream &operator>>( istream &input, bTriangle &tri ){
   int id1, id2, id3;
   cout << "triangle id, origin id, dest id:";
   input >> tri.id >> id1 >> id2 >> id3; 
   return input;
}

/***********************************************************************
**
**  bTriangle "get" functions:
**
**  getID - returns ID #
**  pPtr - returns ptr to one of the 3 vertex nodes, as specified by
**         _index_ (index is 0, 1, or 2)
**  ePtr - returns ptr to one of the 3 clockwise edges, as specified by
**         _index_ (index is 0, 1, or 2)
**  tPtr - returns ptr to one of the 3 adjacent triangles, as specified
**         by _index_ (index is 0, 1, or 2)
**
***********************************************************************/

inline int bTriangle::getID() const {return id;}

inline bNode *bTriangle::pPtr( int index ){
   assert( index >= 0 && index <= 3 );
   return p[index];
}

inline bEdge *bTriangle::ePtr( int index ){
   assert( index >= 0 && index <= 3 );
   return e[index];
}

inline bTriangle *bTriangle::tPtr( int index ){
   assert( index >= 0 && index <= 3 );
   return t[index];
}

/***********************************************************************
**
**  bTriangle "set" functions:
**
**  setID - sets ID #
**  setPPtr - sets pointer to one of the 3 vertex nodes, as specified
**            by _index_ (index is 0, 1, or 2)
**  setEPtr - sets pointer to one of the 3 clockwise edges, as specified
**            by _index_ (index is 0, 1, or 2)
**  setTPtr - sets pointer to one of the 3 adjacent triangles, as
**            specified by _index_ (index is 0, 1, or 2)
**
***********************************************************************/

inline void bTriangle::setID( int val ) {id = ( val >= 0 ) ? val : 0;}

inline void bTriangle::setPPtr( int index, bNode * ndptr ){
   assert( index >= 0 && index <= 3 );
   p[index] = ndptr;
}

inline void bTriangle::setEPtr( int index, bEdge * egptr ){
   assert( index >= 0 && index <= 3 );
   e[index] = egptr;
}

inline void bTriangle::setTPtr( int index, bTriangle * trptr ){
   assert( index >= 0 && index <= 3 );
   t[index] = trptr;
}

/**************************************************************************
**
**  bTriangle::nVOp
**
**  Returns the side number (0, 1, or 2) of the neighboring triangle ct.
**  Assumes that ct _is_ one of the neighboring triangles.
**
**************************************************************************/

inline int bTriangle::nVOp( bTriangle *ct ){
   int i;
   for( i=0; i<4; i++ ){
      assert( i<3 );
      if( t[i] == ct ) return i;
   }
   return i;
}


/**************************************************************************
**
**  bTriangle::nVtx
**
**  Returns the vertex number (0, 1, or 2) associated with node cn.
**  (In other words, it says whether cn is vertex 0, 1, or 2 in the 
**  triangle).
**  Assumes that cn _is_ one of the triangle's vertices.
**
**************************************************************************/

inline int bTriangle::nVtx( bNode *cn ){
   int i;
   for( i=0; i<4; i++ ){
      assert( i<3 );
      if( p[i] == cn ) return i;
   }
   return i;
}

#endif

//=========================================================================
//
//
//                        End of meshElements.h
//
//
//=========================================================================
