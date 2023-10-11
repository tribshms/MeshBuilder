/***************************************************************************
**
**  		      tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  globalFns.h: Global tRIBS Class Header
**
***************************************************************************/

#ifndef GLOBALFNS_H
#define GLOBALFNS_H

//=========================================================================
//
//
//                  Section 1: globalFns Include Statements
//
//
//=========================================================================
#include "src/bMeshElements/bMeshElements.h"
#include "src/Mathutil/predicates.h"

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
class bTriangle;

//=========================================================================
//
//
//                  Section 2: Function Declarations
//
//
//=========================================================================

extern Predicates predicate; 

double ran3( long * ); 

std::vector< double > UnitVector( bEdge* );

double FindCosineAngle0_2_1( std::vector< double > &, std::vector< double > &,
                             std::vector< double > & );

int TriPasses( std::vector< double > &, std::vector< double > &,
               std::vector< double > &, std::vector< double > & );

int PointsCCW( std::vector< double > &, std::vector< double > &, std::vector< double > & );

int NewTriCCW( bTriangle * );

int InNewTri( std::vector< double > &, bTriangle * );

int Intersect( bEdge *, bEdge * );

bEdge* IntersectsAnyEdgeInList( bEdge*, std::list<bEdge*>& );

double InterpSquareGrid( double, double, double**, int, int, int );

std::vector<double> FindIntersectionCoords(
		std::vector<double>, std::vector<double>, 
		std::vector<double>, std::vector<double>);

double PlaneFit(double x, double y, 
                std::vector<double> p0,
                std::vector<double> p1, 
	        std::vector<double> p2, 
		std::vector<double> zs);

double LineFit(double x1, double y1, double x2, double y2, double nx);

double DistanceBW2Points(double x1, double y1, double x2, double y2 );

//=========================================================================
//
//
//                  Section 3: Templated Function Declarations
//
//
//=========================================================================

template< class tSubNode >
int Next3Delaunay(std::list<tSubNode*> &, std::_List_iterator<tSubNode*>);

template< class tSubNode >
int PointAndNext2Delaunay( tSubNode &, std::list<tSubNode*> &,
                           std::_List_iterator<tSubNode*>);
template< class T > 
ostream &operator<<( ostream &, const std::vector< T > & );

/**********************************************************************
**
**  Templated function BinaryWrite
**
**********************************************************************/
                                                                                
template< class outDataType >
inline void BinaryWrite(ostream& outStream, const outDataType& outData)
{
  outStream.write(
    reinterpret_cast<const char*>(&outData), sizeof(outDataType));
}
                                                                                
/**********************************************************************
**
**  Templated function BinaryRead
**
**********************************************************************/
                                                                                
template< class inHolderType >
inline istream& BinaryRead(istream& inStream, inHolderType& inHolder)
{
   return inStream.read(
      reinterpret_cast<char*>(&inHolder), sizeof(inHolderType));
}

#endif

//=========================================================================
//
//
//                          End of globalFns.h 
//
//
//=========================================================================
