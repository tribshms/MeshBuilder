/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  bInputFile.cpp: Functions for class bInputFile (see bInputFile.h)
**
***************************************************************************/

#include "src/bInOut/bInputFile.h"
#include <cstring>

//=========================================================================
//
//
//                  Section 1: bInputFile Constructors/Destructors
//
//
//=========================================================================

bInputFile::bInputFile( const char *filename )
{
  strcpy(InFileName, filename);
  infile.open( filename );
  if( !infile.good() ){
    cout << "\n\nbInputFile: Unable to open file '"
	 <<InFileName<<"'."<<endl;
    cout << "Exiting Program...\n\n";
    exit(1);
  }
}

bInputFile::~bInputFile()
{
  cout<<"bInputFile Object has been destroyed..."<<endl<<flush;
} 

void bInputFile::CloseOldAndOpenNew( const char * filename) {
  infile.close();
  strcpy(InFileName, filename);
  infile.open( filename );
  if( !infile.good() ){
    cout << "\nbInputFile: Unable to open '"<< InFileName << "'." << endl;
    cout << "Exiting Program...\n\n";
    exit(1);
  }
}

//=========================================================================
//
//
//                  Section 2: bInputFile ReadItem Functions
//
//
//=========================================================================


/***************************************************************************
**  
**  bInputFile::ReadItem( const int &datType, const char *itemCode )
**
**  Reads one parameter from the file. The format is assumed to be a line
**  of text that begins with the keyword, followed by a line containing
**  the parameter to be read. The function is overloaded according to the
**  type of data desired. Arbitrary order of the items in the infile allowed.
**  Routine searches through list until it finds itemCode. 
**
***************************************************************************/
int bInputFile::ReadItem( const int &datType, const char *itemCode )
{
  int item;
  char headerLine[kMaxNameLength];
   
  assert( infile.good() );

  infile.getline( headerLine, kMaxNameLength );
  while( !( infile.eof() ) &&
	 ( headerLine[0]==kCommentMark ||
	   strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
    infile.getline( headerLine, kMaxNameLength );

  if( !( infile.eof() ) ){
    infile >> item;
    infile.ignore( 1, '\n' );
  }
  else{
    cout<<"\nError: Expected to read the parameter '"<<itemCode
	<<"', but reached EOF first"<<endl;
    cout<<"\nMissing parameter in the input file..."<<endl;
    infile.clear(); 
    infile.seekg( 0, ios::beg );
    return -9999;
   }

  infile.seekg( 0, ios::beg );
  return item;
}

/***************************************************************************
**  
**  bInputFile::ReadItem( const long &datType, const char *itemCode )
** 
***************************************************************************/
long bInputFile::ReadItem( const long &datType, const char *itemCode )
{
  long item;
  char headerLine[kMaxNameLength];
  
  assert( infile.good() );

  infile.getline( headerLine, kMaxNameLength );
  while( !( infile.eof() ) &&
	 ( headerLine[0]==kCommentMark ||
	   strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
    infile.getline( headerLine, kMaxNameLength );

  if( !( infile.eof() ) ){
    infile >> item;
    infile.ignore( 1, '\n' );
  }
  else{
    cout<<"\nError: Expected to read the parameter '"<<itemCode
	<<"', but reached EOF first"<<endl;
    cout<<"\nMissing parameter in input file..."<<endl;
    infile.clear(); 
    infile.seekg( 0, ios::beg );
    item = -9999;
   }

  infile.seekg( 0, ios::beg );
  return item;
}

/***************************************************************************
**  
**  bInputFile::ReadItem( const double &datType, const char *itemCode )
** 
***************************************************************************/
double bInputFile::ReadItem( const double &datType, const char *itemCode )
{
  double item;
  char headerLine[kMaxNameLength];
 
  assert( infile.good() );

  infile.getline( headerLine, kMaxNameLength );
  while( !( infile.eof() ) &&
	 ( headerLine[0]==kCommentMark ||
	   strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
    infile.getline( headerLine, kMaxNameLength );

  if( !( infile.eof() ) ){
    infile >> item;
    infile.ignore( 1, '\n' );
  }
  else{
    cout<<"\nError: Expected to read the parameter '" << itemCode
	<<"', but reached EOF first" << endl;
    cout<<"\nMissing parameter in input file..."<<endl;
    infile.clear(); 
    infile.seekg( 0, ios::beg );
    item = -999999.;
  }
  infile.seekg( 0, ios::beg );
  return item;
}

/***************************************************************************
**  
**  bInputFile::ReadItem( char *theString, const char *itemCode )
** 
***************************************************************************/
void bInputFile::ReadItem( char * theString, const char *itemCode )
{
  char headerLine[kMaxNameLength];

  assert( infile.good() );
     
  infile.getline( headerLine, kMaxNameLength );
  while( !( infile.eof() ) &&
	 ( headerLine[0]==kCommentMark ||
	   strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
    infile.getline( headerLine, kMaxNameLength ); 

  if( !( infile.eof() ) ){
    infile.getline( theString, kMaxNameLength );
  }
  else {
    cout<<"\nError: Expected to read the parameter '" << itemCode
	<< "', but reached EOF first" << endl;
    cout<<"\nMissing parameter in input file..."<<endl;
    char errr[] = "-999";
    strcpy(theString, errr);
    infile.clear();
    infile.seekg( 0, ios::beg );
    return;
  }

  infile.seekg( 0, ios::beg );
}

/***************************************************************************
**  
**  bInputFile::IsItemIn( const char *itemCode )
** 
***************************************************************************/
int bInputFile::IsItemIn( const char *itemCode )
{
  int item;
  char headerLine[kMaxNameLength];

  assert( infile.good() );

  infile.getline( headerLine, kMaxNameLength );
  while( !( infile.eof() ) &&
	 ( headerLine[0]==kCommentMark ||
	   strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
    infile.getline( headerLine, kMaxNameLength );

  if ( !( infile.eof() ) ) {
    infile.seekg( 0, ios::beg );
    return 1;
  }
  else {
    infile.clear(); 
    infile.seekg( 0, ios::beg );
    return 0;
  }
}

//=========================================================================
//
//
//                        End of InputFile.cpp
//
//
//=========================================================================
