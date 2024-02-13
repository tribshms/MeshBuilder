// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

#include "src/bMeshBuilder/bMeshBuilder.h"
#include "src/bInOut/bInputFile.h"
#include <iostream>

using namespace std;

//=========================================================================
//
//
//                  Section 1: Main Program
//
//
//=========================================================================

int main( int argc, char **argv )
{  

  #ifdef LINUX_32
    using std::cout;
    using std::endl;
  #endif

    if (argc < 2) {
        // No arguments provided, print an error message
        cerr << "Error: No input file provided." << endl;
        return 1; // Return an error code
    }


  // Check command-line arguments
  bInputFile InputFile(argv[1]);
  cout << "Reading input file " << argv[1] << endl;

    try {
        // Call the function to process the input file
        bMeshBuilder MeshBuilder(InputFile);

        // Continue with the rest of your code
    } catch (const exception& e) {
        // Handle the exception (print an error message, etc.)
        std::cerr << "Error: " << e.what() << std::endl;
        return 1; // Return an error code
    }



  cout<<"\n\nPart 9: Deleting Objects and Exiting Program"<<endl;
  cout<<"------------------------------------------------"<<endl<<endl;
}

//=========================================================================
//
//
//                         End of main.cpp
//
//
//=========================================================================
