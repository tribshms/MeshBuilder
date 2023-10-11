/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  bInputFile.h: Header for bInputFile class and objects
**
**  bInputFile Class used in tRIBS for inputing parameters and pathnames
**  using keywords in a *.in file read within various classes.
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: bInputFile Include and Define Statements
//
//
//=========================================================================

#ifndef BINPUTFILE_H
#define BINPUTFILE_H

#include "src/Headers/Definitions.h"

#ifdef ALPHA_64
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#elif defined LINUX_32
  #include <iostream>
  #include <fstream>
  #include <cassert>
  #include <string> 
  #include <cstdlib>
#elif defined WIN
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#else 
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 2: bInputFile Class Definitions
//
//
//=========================================================================

class bInputFile
{
public:
  bInputFile( const char * );	
  ~bInputFile();

  int    IsItemIn( const char * );
  int    ReadItem( const int &, const char * );
  long   ReadItem( const long &, const char * );
  double ReadItem( const double &, const char * );
  void   ReadItem( char *, const char * );
  void   CloseOldAndOpenNew( const char * ); 
  char*  GebInFileName() { return InFileName; } 

private:
  ifstream infile;
  char InFileName[kMaxNameSize]; 
};

#endif

//=========================================================================
//
//
//                          End of bInputFile.h 
//
//
//=========================================================================
