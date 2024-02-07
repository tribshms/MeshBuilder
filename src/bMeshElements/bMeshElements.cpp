// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  meshElements.cpp: Functions for mesh element classes bNode, bEdge,
**                    and bTriangle. 
**
***************************************************************************/

#include "src/bMeshElements/bMeshElements.h"

//====================================================================
//
//
//                  Section 1: Mesh Elements Global Functions
//
//
//====================================================================

int PointsCCW( double*, double*, double* );

/*************************************************************************
**
**  FindIntersectionCoords
**
**  Finds and returns intersection of line segments
**  defined by endpoints given as arguments; 1st seg. endpts are xy1, xy2;
**  2nd seg. endpts are xy3, xy4. (SL) 
**
*************************************************************************/

std::vector<double> FindIntersectionCoords(
		std::vector<double> xy1, std::vector<double> xy2, 
                std::vector<double> xy3, std::vector<double> xy4)
{
   double dxa, dxb, dya, dyb, a, b, c, f, g, h;
   std::vector<double> intxy(2);
  
   dxa = xy2[0] - xy1[0];
   dxb = xy4[0] - xy3[0];
   dya = xy2[1] - xy1[1];
   dyb = xy4[1] - xy3[1];
   a = dya;
   b = -dxa;
   c = dxa * xy1[1] - dya * xy1[0];
   f = dyb;
   g = -dxb;
   h = dxb * xy3[1] - dyb * xy3[0];  

   if( fabs(dxa) > THRESH && fabs(dxb) > THRESH ){
      if( fabs(f - g * a / b) > 0 ){
         intxy[0] = (g * c / b - h) / (f - g * a / b);
         intxy[1] = (-c - a * intxy[0]) / b;
      }
   }
   else{
      if( fabs(dya) > THRESH && fabs(dyb) > THRESH ){
         if( fabs(g - f * b / a) > 0 ){
            intxy[1] = (f * c / a - h) / (g - f * b / a);
            intxy[0] = (-c - b * intxy[1]) / a;
         }
      }
      else {
         if( fabs(dya) <= THRESH ){
            intxy[0] = xy3[0];
            intxy[1] = xy1[1];
         }
         else{
            intxy[0] = xy1[0];
            intxy[1] = xy3[1];
         }
      }
   }
   return intxy;
}


//=========================================================================
//
//
//                  Section 2: bNode Class Functions
//
//
//=========================================================================            


/***********************************************************************
**
**  bNode::insertFrontSpokeList
**  bNode::insertBackSpokeList
**
**  Places eptr at the front or back of the spoke list (respectively)
**  and makes the list circular.
**
***********************************************************************/

void bNode::insertFrontSpokeList( bEdge *eptr )
{
   spokeList.insert(spokeList.begin(), eptr );
}

void bNode::insertBackSpokeList( bEdge *eptr )
{
   spokeList.push_back( eptr );
}

/*****************************************************************************
**
**  bNode::AttachFirstSpoke
**
**  Attaches the first spoke to the node by pointing edg to that spoke,
**  and then telling the spoke to point to itself. thespoke is the edge
**  being added.
**
*****************************************************************************/

void bNode::AttachFirstSpoke( bEdge *thespoke ){
   assert( thespoke!=0 );
   assert( thespoke->getOriginPtr()==this );
   edg = thespoke;
   thespoke->setCCWEdg( thespoke );
}

/*****************************************************************************
**
**  bNode::AttachFirstSpoke
**
**  Clear spoke list.
**
*****************************************************************************/

void bNode::FlushSpokeList() {
   spokeList.erase(spokeList.begin(), spokeList.end());
}

/*****************************************************************************
**
**  bNode::Dist
**
**  Computes the distance of the node from the line formed
**  by points p0 p1 using x y.
**
*****************************************************************************/

double bNode::Dist( bNode * p0, bNode * p1 ){
  double a,b,c,res;

  a=(p1->y)-(p0->y);  
  b=-((p1->x)-(p0->x));
  c=-((a*(p0->x))+(b*(p0->y)));
  res=(a*x + b*y + c) / sqrt(a*a + b*b);
  if (res<0) res=-res;
  return(res);
}

/*****************************************************************************
**
**  bNode::EdgToNod
**
**  Finds and returns the spoke (edge) that connects the current node to _nod_,
**  or zero if no such spoke is found.
**
*****************************************************************************/

bEdge *bNode::EdgToNod( bNode * nod )
{
   std::_List_iterator<bEdge*> spokIter;
   bEdge * ce;

   for (spokIter = spokeList.begin(); spokIter != spokeList.end(); spokIter++) {
      ce = (*spokIter);
      if( ce->getDestinationPtr()->getID() == nod->getID() ) return ce;
   }
   return 0;
}

/*****************************************************************************
**
**  bNode::ComputeVoronoiArea
**
**  Computes the node's Voronoi area by summing the area of embedded triangles.
**
**  The basic Voronoi polygon is described by the set of "right-hand
**  Voronoi vertices" associated with each spoke (edge). These vertices
**  are computed by setVoronoiVertices() as the intersection of the
**  perpendicular bisectors of consecutive spokes. However, in some cases
**  the basic polygon will be distorted, with consecutive vertices NOT
**  being counter-clockwise from one another. (This seems to be the result
**  of numerical errors (possibly in estimating the intersection of two
**  nearly parallel bisectors); according to Sugihara and Iri
**  [J. Comp. Geom. & Appl., 1994, v. 4, p. 179], each Delaunay triangle
**  should be associated with one Voronoi vertex). This is handled by
**  detecting "loops" in the Voronoi polygon and cutting them off by
**  taking the area of the closest (counterclockwise) intersection of 
**  perpendicular bisectors.
**
*****************************************************************************/

double bNode::ComputeVoronoiArea()
{
   int cw;
   double area = 0;
   double a, b, c, dx, dy, dx0, dx1, dy0, dy1, dx2, dy2;
   double vx, vy, x0, y0, x1, y1, x2, y2, m1, m2;
   double temparea; 

   bEdge *edgptr;
   std::list<bEdge*> vedgList;

   // Iterators to current edge, plus next 3 edges
   std::_List_iterator<bEdge*> iter1, iter2, iter3, iter4;
   bEdge *ce, *edge2, *edge3, *edge4;

   std::list<double> vcX, vcY;  
   std::_List_iterator<double> vcIx, vcIy, vcIxnext, vcIynext;

   int i;
   std::vector<double> xy(2), xyn(2), xynn(2), xynnn(2);
   std::vector<double> xy1(2), xy2(2), xy3(2), xy4(2);

   // Create a duplicate list of edges; we will modify this list to obtain
   // the correct vertices. In some cases, we may need to delete an edge
   // to get the correct list of vertices; we don't want to delete the
   // spoke ptr, so we make a duplicate list.

   ce = edg;
   do{
      //assert( ce>0 );
      vedgList.push_back( ce );
      ce = ce->getCCWEdg();
   } while( ce != edg );

   // Check boundary status: Voronoi area only defined for non-boundary nodes
   // edgList is not circular so we have to mimic that behavior
   if(( boundary == 0 ) || (boundary == 3)) {
      cw = true;
      do {
        cw = false;
        // For each edge, look at the next two to see if clockwise
        for (iter1 = vedgList.begin(); iter1 != vedgList.end(); iter1++) {
            // First edge
            ce = (*iter1);
            xy = ce->getRVtx();
            
            // Are any of the initial voronoi vertices equal to the centroid
            if (x == xy[0] && y == xy[1])
               cout << "Cell " << id << ": voronoi vertex = centroid" << endl;

            // Second edge
            iter2 = iter1;
            iter2++;
            if (iter2 == vedgList.end())
               iter2 = vedgList.begin();
            xyn = (*iter2)->getRVtx();

            // Third edge
            iter3 = iter2;
            iter3++;
            if (iter3 == vedgList.end())
               iter3 = vedgList.begin();
            xynn = (*iter3)->getRVtx();

            // Are vertices 1,2,3 clockwise
            // Differences between voronoi vertex 2 and 3
            dx0 = xynn[0] - xyn[0];
            dy0 = xynn[1] - xyn[1];

            // Differences between voronoi vertex 1 and 2
            dx1 = xy[0] - xyn[0];
            dy1 = xy[1] - xyn[1];

            // If the differences are very small set them to be 0
            if (fabs(dx0) < THRESH) dx0 = 0.0;
            if (fabs(dy0) < THRESH) dy0 = 0.0;
            if (fabs(dx1) < THRESH) dx1 = 0.0;
            if (fabs(dy1) < THRESH) dy1 = 0.0;

            // Voronoi vertices 1,2,3 clockwise
            if ( dy0 * dx1 > dx0 * dy1 ) {

               // Fourth edge
               iter4 = iter3;
               iter4++;
               if (iter4 == vedgList.end())
                  iter4 = vedgList.begin();
               xynnn = (*iter4)->getRVtx();

               // Are vertices 2,3,4 clockwise
               dx0 = xynnn[0] - xynn[0];
               dy0 = xynnn[1] - xynn[1];
               dx1 = xyn[0] - xynn[0];
               dy1 = xyn[1] - xynn[1];

               // If the differences are very small set them to be 0
               if (fabs(dx0) < THRESH) dx0 = 0.0;
               if (fabs(dy0) < THRESH) dy0 = 0.0;
               if (fabs(dx1) < THRESH) dx1 = 0.0;
               if (fabs(dy1) < THRESH) dy1 = 0.0;

               // Voronoi vertices 2,3,4 clockwise also (must repair)
               if( dy0 * dx1 > dx0 * dy1 ) {
                  //two consecutive clockwise vertices=>want intersection of
                  //bisectors of edges 2 and 4
                  cw = true;

                  // Centroid
                  x0 = x;
                  y0 = y;

                  // Location of first edge destination node
                  xy1 = ce->getDestinationPtr()->get2DCoords();

                  // Location of third edge destination node
                  xy2 = (*iter3)->getDestinationPtr()->get2DCoords();

                  // Midpoint between centroid and first voronoi vertex
                  x1 = ( x0 + xy1[0] ) / 2;
                  y1 = ( y0 + xy1[1] ) / 2;

                  // Midpoint between centroid and second voronoi vertex
                  x2 = ( x0 + xy2[0] ) / 2;
                  y2 = ( y0 + xy2[1] ) / 2;

                  // Differences between centroid and two midpoints
                  dx1 = x1 - x0;
                  dy1 = y1 - y0;
                  dx2 = x2 - x0;
                  dy2 = y2 - y0;

                  if( fabs(dy1)>0 && fabs(dy2) > 0 ){
                     m1 = -dx1/dy1;
                     m2 = -dx2/dy2;

		     if (m1 == m2){
			   cout<<"\n\n1: Point Correction"<<setprecision(14);
			   cout<<"\nOrigin ID = "<<ce->getOriginPtr()->getID();
			   cout<<"\nOrigin X = "<<ce->getOriginPtr()->getX();
			   cout<<"\nOrigin Y = "<<ce->getOriginPtr()->getY();
			   cout<<"\nOrigin Z = "<<ce->getOriginPtr()->getZ();
			   cout<<"\nOrigin B = "
                               <<ce->getOriginPtr()->getBoundaryFlag();
			   cout<<"\n\nDestination ID = "
                               <<ce->getDestinationPtr()->getID();
			   cout<<"\n\nDestination X = "
                               <<ce->getDestinationPtr()->getX();
			   cout<<"\nDestination Y = "
                               <<ce->getDestinationPtr()->getY();
			   cout<<"\nDestination Z = "
                               <<ce->getDestinationPtr()->getZ();
			   cout<<"\nDestination B = "
                               <<ce->getDestinationPtr()->getBoundaryFlag();
                           cout << "m1 = " << m1 << "  m2 = " << m2 << endl;
			   cout<<"\n\n";
      		     }
                     assert ( m1 != m2 );
                     vx = (y2-m2*x2-y1+m1*x1) / (m1-m2);
                     vy = m1*(vx-x1)+y1;
                  }
                  else{
                     if( fabs(dx1) > 0 && fabs(dx2) > 0 ){
                        m1 = dy1/dx1;
                        m2 = dy2/dx2;

			if (m1 == m2){
 			   cout<<"\n\n2: Point Correction"<<setprecision(14);
			   cout<<"\nOrigin ID = "<<ce->getOriginPtr()->getID();
			   cout<<"\nOrigin X = "<<ce->getOriginPtr()->getX();
			   cout<<"\nOrigin Y = "<<ce->getOriginPtr()->getY();
			   cout<<"\nOrigin Z = "<<ce->getOriginPtr()->getZ();
			   cout<<"\nOrigin B = "
                               <<ce->getOriginPtr()->getBoundaryFlag();
			   cout<<"\n\nDestination ID = "
                               <<ce->getDestinationPtr()->getID();
			   cout<<"\n\nDestination X = "
                               <<ce->getDestinationPtr()->getX();
			   cout<<"\nDestination Y = "
                               <<ce->getDestinationPtr()->getY();
			   cout<<"\nDestination Z = "
                               <<ce->getDestinationPtr()->getZ();
			   cout<<"\nDestination B = "
                               <<ce->getDestinationPtr()->getBoundaryFlag();
                           cout << "m1 = " << m1 << "  m2 = " << m2 << endl;
			   cout<<"\n\n";
      			}
                        assert ( m1 != m2 );
                        vy=(m1*y1+x1-m2*y2-x2)/(m1-m2);
                        vx= -vy*m1+m1*y1+x1;
                     }
                     else if( fabs(dx1) > 0 ){
                        vx = x1;
                        vy = y2;
                     }
                     else{
                        vx = x2;
                        vy = y1;
                     }
                  }
                  // Second edge
                  edgptr = (*iter2);
                  xyn[0] = vx;
                  xyn[1] = vy;
                  dx = xy[0] - vx;
                  dy = xy[1] - vy;

                  // Is new voronoi vertex equal to the centroid
                  if (x == vx && y == vy)
                    cout << "Cell " << id << ": new vertex = centroid" << endl;

                  edgptr->setVEdgLen( sqrt( dx*dx + dy*dy ) );
                  edgptr->setRVtx( xyn );

                  // Third edge
                  edgptr = (*iter3);
                  edgptr->setVEdgLen(0.0);
                  edgptr->setRVtx( xynnn );
                  edgptr = 0;
                  
                  // Remove third edge
                  vedgList.erase(iter3);
               }
            }
         }
      } while( cw ); //while we're still finding loops in the polygon

      //Before the next step, make a list of V. vertex coord. arrays.
      //In doing so, check for parts of the V. area lying outside the
      //mesh domain and cut them off by substituting coordinates of
      //two intersections of V. edges with boundary edge for the V.
      //vertex lying outside the boundary. This should take care of any
      //outlying area as long as all boundaries are convex.
      // Go through spokes and put RVtx of ccw edge in coord list, but
      // first check that the vtx lies within the bndies
      bEdge *ne, *nne;
      bNode *bn0, *bn1;
      for (iter1 = vedgList.begin(); iter1 != vedgList.end(); iter1++) {
         ce = (*iter1);
         ne = ce->getCCWEdg();
         xy1 = ne->getRVtx();

         //checking polygon edge is on boundary and ccw edge's RVtx is on
         //wrong side of bndy edge...
         if( ce->getBoundaryFlag() && ne->getBoundaryFlag() ) {
            bn0 = ce->getDestinationPtrNC();
            bn1 = ne->getDestinationPtrNC();
            xy2 = bn0->get2DCoords();
            xy3 = bn1->get2DCoords();

            if( !PointsCCW( xy1, xy2, xy3 ) ) {					

               xy = FindIntersectionCoords( ce->getRVtx(), xy1, xy2, xy3 );
               vcX.push_back( xy[0] );
               vcY.push_back( xy[1] );

               nne = ne->getCCWEdg();
               xy = FindIntersectionCoords( xy1, nne->getRVtx(), xy2, xy3);
               vcX.push_back( xy[0] );
               vcY.push_back( xy[1] );
            }
            else {
               vcX.push_back( xy1[0] );
               vcY.push_back( xy1[1] );
            }
         }
         else {
            vcX.push_back( xy1[0] );
            vcY.push_back( xy1[1] );
         }
      }
      
      // Now that we've found the correct vertices, make triangles to
      // fill the polygon; the sum of the tri areas is the v. area.
      // For a convex polygon, we can compute the total area as the
      // sum of the area of triangles [P(1) P(i) P(i+1)] for i=2,3...N-1.
  
      // coords of first vertex:
      vcIx = vcX.begin();
      vcIy = vcY.begin();
      xy[0] = (*vcIx);
      xy[1] = (*vcIy);
			
      // Find out # of vertices in polygon:
      int nverts = vcX.size(); 
      for( i=2; i<=nverts-1; i++ )
      {
         vcIx++; vcIy++;
         xyn[0] = (*vcIx);
         xyn[1] = (*vcIy);
         vcIxnext = vcIx;  vcIynext = vcIy;
         vcIxnext++;  vcIynext++;
         xynn[0] = (*vcIxnext);
         xynn[1] = (*vcIynext);

         dx = xyn[0] - xy[0];
         dy = xyn[1] - xy[1];
         a = sqrt( dx*dx + dy*dy );

         dx = xynn[0] - xyn[0];
         dy = xynn[1] - xyn[1];
         b = sqrt( dx*dx + dy*dy );

         dx = xynn[0] - xy[0];
         dy = xynn[1] - xy[1];
         c = sqrt( dx*dx + dy*dy );

         //Added by smr to take care of rounding error, leading to negative value
         temparea = 4*a*a*b*b - (c*c - (b*b + a*a))*(c*c - (b*b + a*a));
	 if (temparea < 0) temparea = 0;

         area += 0.25*sqrt( temparea);
      }
   }
   if (area < 0.0) {
      cout << "CalculateVoronoiArea: varea = " << area << endl;
      area = 0.0;
   }
   varea = area;
   return area;
}

/*******************************************************************
**
**  bNode::getVoronoiVertexList()
**
**  Creates and returns a list of (x,y) coordinates for the
**  Voronoi vertices associated with the node. The list is 
**  created by moving around the spokes and adding each spoke's
**  right-hand Voronoi vertex to the back of the list.
**  A pointer to the vertex list is passed as a parameter; any
**  prior contents are flushed before the list of points is created.
**
*******************************************************************/

void bNode::getVoronoiVertexList( std::list<Point2D> * vertexList )
{
   std::vector<double> vtxarr(2);
   Point2D vtx;
   assert( !boundary );
   vertexList->erase(vertexList->begin(), vertexList->end());

   // Loop around spokes, adding the right-hand Voronoi vertex of
   // each to the list, until we've gone all the way around
   bEdge *ce = edg;
   do{
      vtxarr = ce->getRVtx();
      vtx.x = vtxarr[0];
      vtx.y = vtxarr[1];
      vertexList->push_back( vtx );
      ce = ce->getCCWEdg();
   }
   while( ce!=edg );
   
   assert( vertexList->size()!=0 );
}


/*******************************************************************
**
**  bNode::getVoronoiVertexXYZList()
**
**  Creates and returns a list of (x,y,z) coordinates for the
**  Voronoi vertices associated with the node. The list is 
**  created by moving around the spokes and adding each spoke's
**  right-hand Voronoi vertex to the back of the list. The z
**  coordinate is obtained by linear interpolation from the 3
**  points of the triangle of which the vertex is the circumcenter.
**  A pointer to the vertex list is passed as a parameter; any
**  prior contents are flushed before the list of points is created.
**
*******************************************************************/

void bNode::getVoronoiVertexXYZList( std::list<Point3D> * vertexList )
{
  std::vector<double> vtxarr(2), zvals(3);
  Point3D vtx;
  bNode *n1, *n2;
  assert( !boundary );
  vertexList->erase(vertexList->begin(), vertexList->end());

  // Loop around spokes, adding the right-hand Voronoi vertex of
  // each to the list, until we've gone all the way around
  bEdge *ce = edg;
  n2 = ce->getDestinationPtrNC();
  do{
    ce = ce->getCCWEdg();
    n1 = n2;
    n2 = ce->getDestinationPtrNC();
    vtxarr = ce->getRVtx();
    vtx.x = vtxarr[0];
    vtx.y = vtxarr[1];
    zvals[0] = this->z;
    zvals[1] = n1->getZ();
    zvals[2] = n2->getZ();
    vtx.z = PlaneFit(vtx.x, vtx.y, 
	this->get2DCoords(), n1->get2DCoords(), n2->get2DCoords(), zvals);
    vertexList->push_back( vtx );

  }
  while( ce!=edg );
   
  assert( vertexList->size()!=0 );
}


/*******************************************************************
**
**  bNode::makeCCWEdges
**
**  This function provides for compatibility between the CCW Edge
**  data structure and the Spoke List data structure. It sets up
**  CCW edge connectivity from the spoke list data (which is 
**  assumed to be up to date) by: (1) setting the node's edg 
**  pointer to the first spoke on the list, and (2) setting the
**  ccwedg pointer for each spoke.
**
*******************************************************************/

void bNode::makeCCWEdges()
{
   bEdge *ce, *ccwe;
   std::_List_iterator<bEdge*> tempIter;
   std::_List_iterator<bEdge*> spokIter = spokeList.begin();
   
   ce = (*spokIter);
   assert( ce != 0 );
   setEdg( ce );

   for (; spokIter != spokeList.end(); spokIter++) {
      ce = (*spokIter);
      tempIter = spokIter;
      tempIter++;
      if (tempIter == spokeList.end())
         tempIter = spokeList.begin();
      ccwe = (*tempIter);
      assert( ccwe != 0 );
      ce->setCCWEdg( ccwe );
   }
}

/*******************************************************************
**
**  bNode::ConvertToClosedBoundary
**
**  Makes the node into a closed boundary by setting its boundary
**  status flag. The function also updates the boundary ("flow
**  allowed") status of the node's spokes and their complements.
**
******************************************************************/

void bNode::ConvertToClosedBoundary()
{
   bEdge *ce;   // an edge and its complement

   // Reset boundary flag
   boundary = kClosedBoundary;

   // Signal all connected edges and their complements to become no-flow
   ce = edg;
   do{
     assert( ce!=0 );
     if( ce->getBoundaryFlag()==kFlowAllowed ){
       ce->setFlowAllowed( 0 );
     }  
   } while( (ce=ce->getCCWEdg()) != edg );
}

#ifdef PARALLEL_TRIBS
/**************************************************************************
**
**  get or set ResIndex
**
**  The resample index may not be the same as the node ID
**  when running in parallel. So an additional index is
**  available for resampling.
**
**************************************************************************/
                                                                                             
void bNode::setResIndex(int rind) {
  resIndex = rind;
}
                                                                                             
int bNode::getResIndex() {
    return resIndex;
}

#endif

/**************************************************************************
**
**  TellAll
**
**  Debugging routine.
**
**************************************************************************/

#ifndef NDEBUG
void bNode::TellAll()
{
   cout << " NODE " << id << ":\n";
   cout << "  x=" << x << " y=" << y << " z=" << z;
   cout << "  boundary: " << boundary
        << "\n  varea: " << varea << endl;
   if( edg )
       cout << "  points to edg #" << edg->getID() << endl;
   else cout << "  edg is undefined!\n";
   
}
#endif

//=========================================================================
//
//
//                  Section 2: bFlowNode Class Functions
//
//
//=========================================================================
bFlowNode::bFlowNode() : bNode()
{
  NwtNew = 0.0;
  Rain = 0.0;
  srf=hsrf=psrf=satsrf=sbsrf=0.0;
  flowedge = 0;
  traveltime = hillpath = streampath = 0.0;
  EvapoTranspiration = 0.0;
  SoilMoistureSC = SoilMoistureUNSC = 0.0;
  RootMoistureSC = 0.0;
  ContrArea = Curvature = 0.0;
  StreamPtr = 0;
  tracer = flood = 0;
  reach = -1;
  xC = -1;
  yC = -1;
}

double bFlowNode::getNwtNew() { return NwtNew; }
double bFlowNode::getRain()   { return Rain;  }
double bFlowNode::getSrf()    { return srf;  }
double bFlowNode::getHsrf()   { return hsrf; }
double bFlowNode::getPsrf()   { return psrf; }
double bFlowNode::getSatsrf() { return satsrf; }
double bFlowNode::getSbsrf()  { return sbsrf;  }
double bFlowNode::getTTime()  { return traveltime; }
double bFlowNode::getHillPath()   { return hillpath; }
double bFlowNode::getStreamPath() { return streampath; }
double bFlowNode::getEvapoTrans()    { return EvapoTranspiration; }
double bFlowNode::getSoilMoistureSC()   { return SoilMoistureSC; }
double bFlowNode::getSoilMoistureUNSC() { return SoilMoistureUNSC; }
double bFlowNode::getRootMoistureSC()   { return RootMoistureSC; }
double bFlowNode::getContrArea() { return ContrArea; }
double bFlowNode::getCurvature() { return Curvature; }
int bFlowNode::getTracer() { return tracer; }
int bFlowNode::getReach() { return reach; }
int bFlowNode::getFloodStatus()  { return flood; }
bEdge * bFlowNode::getFlowEdg()     { return flowedge; }
bFlowNode * bFlowNode::getStreamNode()	{ return StreamPtr; }

void bFlowNode::setStreamPath(double l) { streampath = l; }
void bFlowNode::setFlowEdg(bEdge * edgs) { flowedge = edgs; }
void bFlowNode::setContrArea(double value)    { ContrArea = value; }
void bFlowNode::setCurvature(double value)    { Curvature = value; }
void bFlowNode::setTracer(int cnt) { tracer = cnt; }
void bFlowNode::setReach(int value) { reach = value; }
void bFlowNode::setFloodStatus( int status )  { flood = status; }
void bFlowNode::setHillPath(double l)   { hillpath = l; }
void bFlowNode::setTTime(double value)  { traveltime = value; }
void bFlowNode::setStreamNode(bFlowNode *cn) { StreamPtr = cn; }

void bFlowNode::addContrArea(double value)   { ContrArea += value; }

void bFlowNode::ActivateSortTracer() { tracer = 1; }
void bFlowNode::DeactivateTracer() { tracer = 0; }

void bFlowNode::MoveSortTracerDownstream()
{
  tracer--;
  getDownstrmNbr()->AddTracer();
}

void bFlowNode::AddTracer()
{
  if((boundary==0)||(boundary==3)) tracer++;
}

int bFlowNode::NoMoreTracers()
{
  return( tracer==0 );
}

bFlowNode * bFlowNode::getDownstrmNbr()
{
  if( flowedge == 0 ) return 0;
  return (bFlowNode*) (flowedge->getDestinationPtrNC());
}

/**************************************************************************
**
**  tCNode:: Get Centroid of associated voronoi polygon 
**
**************************************************************************/
double bFlowNode::getCentroidX()
{
  if (xC == -1){
    double areaT = 0.0;
        
    std::vector<double> xy(2);
    int nPoints;    
    int cnt = 0;  
    bEdge *firstedg;
    bEdge *curedg;
                
    firstedg = getFlowEdg();
    if (!firstedg) 
      firstedg = getEdg();
    curedg = firstedg->getCCWEdg();
    cnt++;      
    while (curedg != firstedg) { 
      curedg = curedg->getCCWEdg();
      cnt++;            
    }           

    nPoints = cnt;
                
    double *vXs = new double [cnt+1];
    double *vYs = new double [cnt+1];

    int iv = 0; 
    firstedg = getFlowEdg(); 

    xy = firstedg->getRVtx();
    vXs[iv] = xy[0];
    vYs[iv] = xy[1];
    iv++;       
    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg) { 
      xy = curedg->getRVtx();
      if (xy[0] == vXs[iv-1] && xy[1] == vYs[iv-1]) {
        iv--;         //If points coincide    
        nPoints--; //just skip it...
      }                 
      else {            
        vXs[iv] = xy[0];        
        vYs[iv] = xy[1];        
      }                 
      iv++;             
      curedg = curedg->getCCWEdg();
    }           
    vXs[iv] = vXs[0];
    vYs[iv] = vYs[0];

    cnt = polyCentroid(vXs, vYs, nPoints,
              &xC, &yC, &areaT);

    delete [] vXs;
    delete [] vYs;

  }
  return xC;
}

double bFlowNode::getCentroidY()
{
  if (yC == -1) {
    double areaT = 0.0;

    std::vector<double> xy(2);
    int nPoints;
    int cnt = 0;
    bEdge *firstedg;
    bEdge *curedg;

    firstedg = getFlowEdg();
    if (!firstedg)
      firstedg = getEdg();
    curedg = firstedg->getCCWEdg();
    cnt++;
    while (curedg != firstedg) {
      curedg = curedg->getCCWEdg();
      cnt++;
    }

    nPoints = cnt;

    double *vXs = new double [cnt+1];
    double *vYs = new double [cnt+1];

    int iv = 0;
    firstedg = getFlowEdg();

    xy = firstedg->getRVtx();
    vXs[iv] = xy[0];
    vYs[iv] = xy[1];
    iv++;
    curedg = firstedg->getCCWEdg();
    while (curedg != firstedg) {
      xy = curedg->getRVtx();
      if (xy[0] == vXs[iv-1] && xy[1] == vYs[iv-1]) {
        iv--;         //If points coincide
        nPoints--; //just skip it...
      }
      else {
        vXs[iv] = xy[0];
        vYs[iv] = xy[1];
      }
      iv++;
      curedg = curedg->getCCWEdg();
    }
    vXs[iv] = vXs[0];
    vYs[iv] = vYs[0];

    cnt = polyCentroid(vXs, vYs, nPoints,
              &xC, &yC, &areaT);

    delete [] vXs;
    delete [] vYs;

  }
  return yC;
}

int bFlowNode::polyCentroid(double x[], double y[], int n,
            double *xCentroid, double *yCentroid, double *area) {
    int i, j;
  double ai, atmp = 0.0, xtmp = 0.0, ytmp = 0.0;
  if (n < 3) return 1;
  for (i = n-1, j = 0; j < n; i = j, j++)
    {   
    ai = x[i] * y[j] - x[j] * y[i];
    atmp += ai; 
    xtmp += (x[j] + x[i]) * ai;
    ytmp += (y[j] + y[i]) * ai; 
    }           
  *area = atmp / 2.0;
  if (atmp != 0.0)
    {   
    *xCentroid =  xtmp / (3.0 * atmp);
    *yCentroid =  ytmp / (3.0 * atmp);
    return 0;   
    }           
  return 2;
}       


//=========================================================================
//
//
//                  Section 3: bEdge Class Functions
//                                                    
//                                                    
//
//=========================================================================

/**************************************************************************
**
**  bEdge::CalcLength
**
**  Computes the edge length and returns it. (Length is the projected
**  on the x,y plane). Assumes org and dest are valid.
**
**************************************************************************/

double bEdge::CalcLength()
{
   assert( org!=0 );  
   assert( dest!=0 );
   
   double dx = org->getX() - dest->getX();
   double dy = org->getY() - dest->getY();
   double len = sqrt( dx*dx + dy*dy );
   return len;
}

/**************************************************************************
**
**  bEdge::TellCoords
**
**  Debugging routine that reports edge ID and coords of endpoints.
**
**************************************************************************/

void bEdge::TellCoords(){
   cout << "EDGE " << id << ":\n";
   cout << "  " << org->getID() << " (" << org->getX() << ","
        << org->getY() << ") -> " << dest->getID() << " ("
        << dest->getX() << "," << dest->getY() << ")" << endl;
}


/**************************************************************************
**
**  bEdge::FindComplement
**
**  Finds and returns the edge's complement edge. Does this by checking
**  the spokes connected to its destination node and returning the one
**  that connects back to its origin.
**
**  Returns:  ptr to the complement edge
**  Assumes:  valid connectivity (one of destination's spokes connects
**            back to origin)
**
**************************************************************************/

bEdge * bEdge::FindComplement(){
   assert( org!=0 && dest!=0 && dest->getEdg()!=0 );
   
   bEdge * ce = dest->getEdg();
   while( ce->getDestinationPtrNC() != org ){
      ce = ce->getCCWEdg();
      assert( ce!=0 && ce!=dest->getEdg() );
   }
   return ce;
}

//=========================================================================
//
//
//                  Section 4: bTriangle Class Functions
//                                                    
//                                                    
//
//=========================================================================


/*****************************************************************************
**
**  bTriangle::FindCircumcenter
**
**  Finds the circumcenter of the triangle by finding the intersection of
**  the perpendicular bisectors of sides (p0,p1) and (p0,p2). Returns the
**  coordinates of the circumcenter as a 2-element array. Note that the
**  circumcenter is also the Voronoi cell vertex associated with the
**  triangle's three nodes (that's the point of computing it).
**
*****************************************************************************/

std::vector<double> bTriangle::FindCircumcenter()
{
   double x1, y1, x2, y2, dx1, dy1, dx2, dy2, m1, m2;
   std::vector<double> xyo(2), xyd1(2), xyd2(2), xy(2);

   assert( pPtr(0) && pPtr(1) && pPtr(2) );
   
   // Coordinates of triangle's nodes p0, p1, and p2
   xyo = pPtr(0)->get2DCoords();
   xyd1 = pPtr(1)->get2DCoords();
   xyd2 = pPtr(2)->get2DCoords();

   // Find the midpoints of the two sides (p0,p1) and (p0,p2) and store them 
   // in (x1,y1) & (x2,y2). Then get the distance between p0 and the 
   // midpoints of each side

   x1 = (xyo[0] + xyd1[0]) / 2;
   y1 = (xyo[1] + xyd1[1]) / 2;
   x2 = (xyo[0] + xyd2[0]) / 2;
   y2 = (xyo[1] + xyd2[1]) / 2;
   dx1 = x1-xyo[0];
   dy1 = y1-xyo[1];
   dx2 = x2-xyo[0];
   dy2 = y2-xyo[1];

   // Compute the intercept of the bisectors of the two sides:
   // Case: neither spoke is horizontal (ok to divide by dy1 and dy2)

   if( fabs(dy1)>0 && fabs(dy2)>0 ){
      assert( dy1!=0 && dy2!=0 );
      m1= -dx1/dy1;
      m2 = -dx2/dy2;
      if (m1 == m2){
        cout<<"\nEdges are parallel " << id;
	cout<<"\nPoint 0 X = "<<pPtr(0)->getX();
	cout<<"\nPoint 0 Y = "<<pPtr(0)->getY();
	cout<<"\nPoint 0 Z = "<<pPtr(0)->getZ();
	cout<<"\nPoint 0 B = "<<pPtr(0)->getBoundaryFlag();
	cout<<"\n\nPoint 1 X = "<<pPtr(1)->getX();
	cout<<"\nPoint 1 Y = "<<pPtr(1)->getY();
	cout<<"\nPoint 1 Z = "<<pPtr(1)->getZ();
	cout<<"\nPoint 1 B = "<<pPtr(0)->getBoundaryFlag();
	cout<<"\n\nPoint 2 X = "<<pPtr(2)->getX();
	cout<<"\nPoint 2 Y = "<<pPtr(2)->getY();
	cout<<"\nPoint 2 Z = "<<pPtr(2)->getZ();
	cout<<"\nPoint 2 B = "<<pPtr(0)->getBoundaryFlag();
	cout<<"\n\n";
      }
      assert( m1!=m2 );
      xy[0] = (y2 - m2 * x2 - y1 + m1 * x1) / (m1 - m2);
      xy[1] = m1 * (xy[0] - x1) + y1;
   }

   // Case: one spoke is horizontal, but neither are vertical

   else if( dx1!=0 && dx2!=0 ){
      assert( dx1!=0 && dx2!=0 );
      m1 = dy1/dx1;
      m2 = dy2/dx2;
      if (m1 == m2){
        cout<<"One edge is horizontal, neither vertical " << id;
	cout<<"\nPoint 0 X = "<<pPtr(0)->getX();
	cout<<"\nPoint 0 Y = "<<pPtr(0)->getY();
	cout<<"\nPoint 0 Z = "<<pPtr(0)->getZ();
	cout<<"\nPoint 0 B = "<<pPtr(0)->getBoundaryFlag();
	cout<<"\n\nPoint 1 X = "<<pPtr(1)->getX();
	cout<<"\nPoint 1 Y = "<<pPtr(1)->getY();
	cout<<"\nPoint 1 Z = "<<pPtr(1)->getZ();
	cout<<"\nPoint 1 B = "<<pPtr(0)->getBoundaryFlag();
	cout<<"\n\nPoint 2 X = "<<pPtr(2)->getX();
	cout<<"\nPoint 2 Y = "<<pPtr(2)->getY();
	cout<<"\nPoint 2 Z = "<<pPtr(2)->getZ();
	cout<<"\nPoint 2 B = "<<pPtr(0)->getBoundaryFlag();
	cout<<"\n\n";
      }
      assert( m1!=m2 );
      xy[1] = (m1 * y1 + x1 - m2 * y2 - x2) / (m1 - m2);
      xy[0] = -xy[1] * m1 + m1 * y1 + x1;
   }

   // Special case: one is vertical, the other horizontal

   else{
      if( dx1!=0 ){
         xy[0] = xyo[0] + dx1;
         xy[1] = xyo[1] + dy2;
         assert( dx2==0 && dy1==0 );
      }
      else{
         xy[0] = xyo[0] + dx2;
         xy[1] = xyo[1] + dy1;
         assert( dx1==0 && dy2==0 );
      }
   }
   assert( &xy != 0 );
   return xy;
}

/*****************************************************************************
**
**  bTriangle::TellAll()
** 
**  Debugging Routing
**
*****************************************************************************/

#ifndef NDEBUG
void bTriangle::TellAll()
{
   int i;
   
   assert( this!=0 );
   cout << "TRIANGLE #" << id << ":\n";
   for( i=0; i<3; i++ ){
      cout << "  P" << i << " ";
      if( p[i]!=0 ) cout << p[i]->getID() << " (" << p[i]->getX() << ","
                         << p[i]->getY() << ")";
      else cout << "(ndef)";
      cout << "  E" << i << " ";
      if( e[i]!=0 ) cout << e[i]->getID();
      else cout << "(ndef)";
      cout << "  T" << i << " ";
      if( t[i]!=0 ) cout << t[i]->getID();
      else cout << "(ndef)";
      cout << endl;
   }
}
#endif

//=========================================================================
//
//
//              End of MeshElements.cpp
//
//
//=========================================================================
