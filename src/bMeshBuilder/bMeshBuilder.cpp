// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  bMeshBuilder.cpp: Functions for class MeshBuilder
**
***************************************************************************/

#include "src/bMeshBuilder/bMeshBuilder.h"
#include <ctime>

//=========================================================================
//
//
//                  Section 1: bMeshBuilder Constructors/Destructors
//
//
//=========================================================================

bMeshBuilder::bMeshBuilder(bInputFile& infile)
{
  time_t clock;
  clock = time(NULL);
  cout << "bMeshBuilder: Start ... " << ctime(&clock) << endl;

  // Check for basic mesh files and skip creation if they exist
  ifstream nodeFile("nodes.basic");
  ifstream edgeFile("edges.basic");

  if (nodeFile.fail() || edgeFile.fail()) {
     nodeFile.close();
     edgeFile.close();

    // Read in MESHINPUT Option, previously meshbuilder assumed only poitns files
    int read = infile.ReadItem( read, "OPTMESHINPUT" );  
    
    if( read < 1 || read > 9 ){
      cerr << "\nInvalid mesh input option requested.";
      cerr << "\nValid options for reading mesh input are:\n"
        << "  1 -- read mesh from input data files\n"
        << "  8 -- create mesh from points file using Triangulator\n\n";
      exit(1);
    }

    // Create mesh using only structural nodes
    baseMesh = new bMesh<bNode>();
     
    if( read == 1 ) {
      cout<<"\n\nCreating Mesh from Existing Mesh (Option 1)"<<endl;
      cout<<"---------------------------------------------------"<<endl;
      baseMesh = new bMesh<bNode>();      
      baseMesh->MakeMeshFromInputData( infile );
      cout<<"\nMakeMeshFromInputData Successful Using Option 1"<<endl<<flush;

      // NOW, the mesh is in 'baseMesh'. We need to calculate Voronoi properties
      // and then write it to .basic files.
      // This mirrors what RefineBasicMesh does after loading temporary files.
      cout << "bMesh: Calculating Voronoi for pre-built mesh ... " << ctime(&clock) << endl;
      baseMesh->setVoronoiVertices();
      baseMesh->CalcVoronoiEdgeLengths();
      baseMesh->CalcVAreas();

    } else if ( read == 8 ) { // "triangulate from points" option
      cout<<"\n\nCreating Mesh from Points File (Option 8)"<<endl;
      cout<<"---------------------------------------------------"<<endl;
      // Using a points file, triangulate to create nodes, edge, triangles
      // Write this basic structure to temporary files
      clock = time(NULL);
      cout << "bMesh: MakeBasicMesh ... " << ctime(&clock) << endl;
      baseMesh->MakeBasicMesh(infile);

      // Read in the temporary files building just nodes, edges, triangles
      // Do the voronoi calculations and add the results to nodes and edges
      clock = time(NULL);
      cout << "bMesh: RefineBasicMesh ... " << ctime(&clock) << endl;
      baseMesh->RefineBasicMesh(); // This already does Voronoi calculations

    }

     // Write the complete node and edge files which will be read in later
     // *** Writing bNode
     clock = time(NULL);
     cout << "bMesh: WriteBasicMesh ... " << ctime(&clock) << endl;
     baseMesh->WriteBasicMesh();

     delete baseMesh;

  } else {
     nodeFile.close();
     edgeFile.close();
     cout << "bMesh: Basic mesh files exist" << endl;
  }

  // Create the flow mesh from the temporary files in data structures
  flowMesh = new bMesh<bFlowNode>();

  // Rebuild the node and edge files complete with voronoi information
  // so that the flow network can be calculated
  // *** Reading bNode but creating bFlowNode
  clock = time(NULL);
  cout << "bMesh: ReadBasicMesh ... " << ctime(&clock) << endl;
  flowMesh->ReadBasicMesh();

  // Create the flow network
  clock = time(NULL);
  cout << "bFlowNet: Constructor ... " << ctime(&clock) << endl;
  flowNet = new bFlowNet(flowMesh, infile);

  // Create the graph connectivity using the reaches of the flow network
  clock = time(NULL);
  cout << "bGraph: Constructor ... " << ctime(&clock) << endl;
  graph = new bGraph(flowMesh, flowNet, infile);

  // Write the geometry file for the visualizer 
  cout << "bFlowNet: WriteGeometry ... " << ctime(&clock) << endl;
  flowNet->WriteGeometry();

  // Write flow and partition information for the simulator
  cout << "bFlowNet: WriteFlowNet ... " << ctime(&clock) << endl;
  flowNet->WriteFlowNet();

  // Write the flow network to a file such that individual reaches can
  // be read in on parallel processors so that no mesh is too large
  // Deletes nodes from main node list and puts into lists by reach
  // So we can't access nodeList after this
  // *** Writing bFlowNode
  clock = time(NULL);
  cout << "bGraph: WriteFlowMesh ... " << ctime(&clock) << endl;
  graph->WriteFlowMesh();
  
  // Option 9 will read bFlowNode and will create tCNode
  clock = time(NULL);
  cout << "bMeshBuilder: Finished ... " << ctime(&clock) << endl;
}

bMeshBuilder::~bMeshBuilder() 
{  
  delete flowMesh;
  delete flowNet;
  delete graph;
}

//=========================================================================
//
//
//                          End of bMeshBuilder.cpp
//
//
//=========================================================================
