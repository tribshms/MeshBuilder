/***************************************************************************
**
**  		   tRIBS Distributed Hydrologic Model
**
**            TIN-based Real-time Integrated Basin Simulator
**		       Ralph M. Parsons Laboratory
**  		  Massachusetts Institute of Technology
**  
**                      VERSION 3.0, JULY 2006
**
**  Designed and created by the Hydrology group under supervision of
**      prof. Rafael L. Bras, Department of Civil and Environmental 
**      Engineering, Massachusetts Institute of Technology, 
**      Cambridge, MA, U.S.A.
**
**  MAIN.CPP -- This file contains the main() routine that handles
**              top-level initialization.
**
**  NOTE: This source code is copyrighted material. It is distributed
**        solely for noncommercial research and educational purposes.
**        Use in whole or in part for commercial purposes without a
**        written license from the copyright holder(s) is expressly
**        prohibited. Copies of this source code or any of its components 
**        may not be transferred to any other individuals or organizations
**        without written consent. Copyright (C) Massachusetts Institute
**        of Technology, 2001-2003. All rights reserved.
**
***************************************************************************/

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

  // Check command-line arguments
  bInputFile InputFile(argv[1]);
  cout << "Reading input file " << argv[1] << endl;

  bMeshBuilder MeshBuilder(InputFile);

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
