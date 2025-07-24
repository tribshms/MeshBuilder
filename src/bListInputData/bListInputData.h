// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bListInputData.h: 	Header file for class bListInputData
**
**  This class is used to read in lists of Delaunay-triangulated
**  mesh elements from three user-provided input files, which contain
**  the nodes, directed edges, and triangles in the mesh, respectively.
**
**************************************************************************/

#ifndef BLISTINPUTDATA_H
#define BLISTINPUTDATA_H
#define kTimeLineMark ' '

#include "src/Headers/Inclusions.h"
#include "src/Headers/Definitions.h"
#include "src/bInOut/bInputFile.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
using namespace std;

//=========================================================================
//
//
//                  Section 1: bListInputData Class Declaration
//
//
//=========================================================================

/**************************************************************************
**  
**  class bListInputData()
**
**  bListInputData reads an established triangulation from a set of four
**  files, and stores the data in a series of arrays. The files are:
**    <name>.nodes  --  node (point) data
**    <name>.edges  --  directed edge data
**    <name>.tri    --  triangle data
**    <name>.z      --  "z" value data (elevation or other)
**
**************************************************************************/

template< class bSubNode >
class bListInputData{
    //friend class bMesh< bSubNode >; 

public:
    bListInputData( bInputFile & ); 
    ~bListInputData();
    void GetKeyEntry();         // not currently supported
    void GetFileEntry();        // read data from files
    int nnodes, nedges, ntri;  	// # nodes, edges, & triangles
    ifstream nodeinfile;   	// node input file
    ifstream edgeinfile;   	// edge input file
    ifstream triinfile;    	// triangle input file
    ifstream zinfile;      	// "z" input file std::vector
    std::vector< double > x;      	// node x coords
    std::vector< double > y;      	// node y coords
    std::vector< double > z;      	// node z values
    std::vector< int > edgid;     	// node edge ID #s
    std::vector< int > boundflag; 	// node boundary codes
    std::vector< int > orgid;     	// directed edge origin node ID #s
    std::vector< int > destid;    	// directed edge destination node ID #s
    std::vector< int > nextid;    	// ID #s of next counter-clockwise edges
    std::vector< int > p0;     	// IDs of triangle node 0
    std::vector< int > p1;     	// IDs of triangle node 1
    std::vector< int > p2;     	// IDs of triangle node 2
    std::vector< int > e0;     	// IDs triangle clockwise-oriented edge 0
    std::vector< int > e1;     	// IDs triangle clockwise-oriented edge 1
    std::vector< int > e2;   	// IDs triangle clockwise-oriented edge 2
    std::vector< int > t0;     	// IDs of neighboring tri's opposite node 0
    std::vector< int > t1;     	// IDs of neighboring tri's opposite node 1
    std::vector< int > t2;     	// IDs of neighboring tri's opposite node 2 

};

#endif

//=========================================================================
//
//
//                  Section 1: bListInputData Constructor and Destructor
//
//
//=========================================================================

/**************************************************************************
**
**  bListInputData constructor
**
**  The constructor takes an input file and reads from it the base name 
**  containing the triangulation. It then opens <basename>.nodes,
**  <basename>.z, <basename>.edges, and <basename>.tri. Assuming the
**  files are valid, the desired time-slice is read from infile, and
**  the start of data for that time-slice is sought in each of the four
**  triangulation files. The arrays are dimensioned as needed, and
**  GetFileEntry() is called to read the data into the arrays. Note that
**  the time in each file is identified by a space character preceding it
**  on the same line.
**
**  Modifications for MeshBuilder simplifies routine by not having to search
**  through multiple time slices. MeshBuilder has only one time slice at 
**  time zero. Option for CHILD not maintained, input time assumed zero. 
**
**************************************************************************/

template< class bSubNode >
bListInputData< bSubNode >::
bListInputData( bInputFile &infile )                  
{
	int righttime;                   	//flag: found the right time slice
	double time, intime;             	//current & desired time
	char basename[80],               	//base name of input files
		inname[80];                  	//full name of an input file
	char nodeHeader[kMaxNameLength]; 	//header line read from input file
	char zHeader[kMaxNameLength]; 	//header line read from input file
	char edgHeader[kMaxNameLength]; 	//header line read from input file
	char triHeader[kMaxNameLength]; 	//header line read from input file
	
	cout<<"\nCreating bListInputData for reading in existing mesh..."<<endl;
	
	// Read base name for triangulation files from infile
	infile.ReadItem( basename, "INPUTDATAFILE" );

	// Open each of the four files
	strcpy( inname, basename );
	strcat( inname, ".nodes" );
	nodeinfile.open(inname);    	//Node input file pointer
	
	strcpy( inname, basename );
	strcat( inname, ".edges" );
	edgeinfile.open(inname);   	//Edge input file pointer
	
	strcpy( inname, basename );
	strcat( inname, ".tri" );
	triinfile.open( inname );   	//Triangle input file pointer
	
	strcpy( inname, basename );
	strcat( inname, ".z" );
	zinfile.open( inname );     	//Elevations input file pointer
	
	// Make sure we found them
	if( !nodeinfile.good() || !edgeinfile.good() || !triinfile.good()
		|| !zinfile.good() ){
		cerr << "Error: I can't find one or more of the following files:\n"
		<< "\t" << basename << ".nodes\n"
		<< "\t" << basename << ".edges\n"
		<< "\t" << basename << ".tri\n"
		<< "\t" << basename << ".z\n";
		cerr <<"Unable to open triangulation input file(s).";
	}
	
	// Timeslice hardcoded to 0 as the time funcationaility is not used.
	intime = 0;
	
	cout<<"\nReading in existing mesh files:"<<endl;
	cout<<"'"<< basename << ".nodes'"<<endl;
	cout<<"'"<< basename << ".edges'"<<endl;
	cout<<"'"<< basename << ".tri'"<<endl;
	cout<<"'"<< basename << ".z'"<<endl;
	
	
	// For tRIBS, use INPUTTIME = 0 as flag implying only one mesh
	// list (nodes, edges, tri, z) per file opened. Will read the
	// the number of elements directly, followed by the values in
	// the GetFileEntry() Function.
	
	
	if(intime == 0){
		nodeinfile >> time;
		nodeinfile >> nnodes;
		
		zinfile >> time;
		zinfile >> nnodes;
		
		edgeinfile >> time;
		edgeinfile >> nedges;
		
		triinfile >> time;
		triinfile >> ntri;
	}
	
	
	// Keep the following for CHILD compatibility
	
	else{    // Find specified input times in input data files and read # items.
		
		//Nodes:
		righttime = 0;
		time = 0;
		while( !( nodeinfile.eof() ) && !righttime ){
			nodeinfile.getline( nodeHeader, kMaxNameLength );
			cout<<"nodeheader[0] = "<<nodeHeader[0]<<endl;
			if( nodeHeader[0] == kTimeLineMark ){
				nodeinfile.seekg( -nodeinfile.gcount(), ios::cur );
				nodeinfile >> time;
				cout << "from node file, time = " << time << endl;
				cout << "from node file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( nodeinfile.eof() ) ) nodeinfile >> nnodes;
		else{
			cerr << "\nCouldn't find the specified input time in the node file\n";
		}
		
		//"z" values:
		righttime = 0;
		time = 0;
		while( !( zinfile.eof() ) && !righttime ){
			zinfile.getline( zHeader, kMaxNameLength );
			cout<<"zheader[0] = "<<zHeader[0]<<endl;
			if( zHeader[0] == kTimeLineMark ){
				zinfile.seekg( -zinfile.gcount(), ios::cur );
				zinfile >> time;
				cout << "from z file, time = " << time << endl;
				cout << "from z file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( zinfile.eof() ) ) zinfile >> nnodes;
		else{
			cerr << "Couldn't find the specified input time in elevation file\n";
		}
		
		//Edges:
		righttime = 0;
		time = 0;
		while( !( edgeinfile.eof() ) && !righttime ){
			edgeinfile.getline( edgHeader, kMaxNameLength );
			cout<<"edgheader[0] = "<<edgHeader[0]<<endl;
			if( edgHeader[0] == kTimeLineMark ){
				edgeinfile.seekg( -edgeinfile.gcount(), ios::cur );
				edgeinfile >> time;
				cout << "from edg file, time = " << time << endl;
				cout << "from edg file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( edgeinfile.eof() ) ) edgeinfile >> nedges;
		else{
			cerr << "Couldn't find the specified input time in the edge file\n";
		}
		
		//Triangles:
		righttime = 0;
		time = 0;
		while( !( triinfile.eof() ) && !righttime ){
			triinfile.getline( triHeader, kMaxNameLength );
			cout<<"triheader[0] = "<<triHeader[0]<<endl;
			if( triHeader[0] == kTimeLineMark ){
				triinfile.seekg( -triinfile.gcount(), ios::cur );
				triinfile >> time;
				cout << "from tri file, time = " << time << endl;
				cout << "from tri file, intime = " <<intime <<endl;
				if( time == intime ) righttime = 1;
			}
		}
		if( !( triinfile.eof() ) ) triinfile >> ntri;
		else{
			cerr << "Couldn't find the specified input time in the tri file\n";
		}
	}
	
	// Dimension the arrays 
	x.resize( nnodes );
	y.resize( nnodes );
	z.resize( nnodes );
	edgid.resize( nnodes );
	boundflag.resize( nnodes );
	orgid.resize( nedges );
	destid.resize( nedges );
	nextid.resize( nedges );
	p0.resize( ntri );
	p1.resize( ntri );
	p2.resize( ntri );
	e0.resize( ntri );
	e1.resize( ntri );
	e2.resize( ntri );
	t0.resize( ntri );
	t1.resize( ntri );
	t2.resize( ntri );
	
	// Read in data from file
	GetFileEntry();
	
	// Close the files
	nodeinfile.close();
	edgeinfile.close();
	triinfile.close();
	zinfile.close();
}

template< class bSubNode >
bListInputData< bSubNode >:: ~bListInputData(){ }

//=========================================================================
//
//
//                  Section 2: bListInputData Functions
//
//
//=========================================================================

/**************************************************************************
**
**  bListInputData::GetFileEntry()
**
**  Reads node, edge, and triangle data from the four triangulation input
**  files. Assumes that each files is open and valid and that the current
**  reading point in each corresponds the start of data for the desired
**  time-slice.
**
**************************************************************************/

template< class bSubNode >
void bListInputData< bSubNode >::
GetFileEntry()                 
{
	int i;
	for( i=0; i< nnodes; i++ ){
		nodeinfile >> x[i] >> y[i] >> edgid[i] >> boundflag[i];
		zinfile >> z[i];
	}
	
	for( i=0; i<nedges; i++ )
		edgeinfile >> orgid[i] >> destid[i] >> nextid[i];
	
	for( i=0; i< ntri; i++ )
		triinfile >> p0[i] >> p1[i] >> p2[i] >> t0[i] >> t1[i] >> t2[i] 
			>> e0[i] >> e1[i] >> e2[i];
	
}

//=========================================================================
//
//
//                      End of bListInputData.h
//
//
//=========================================================================

