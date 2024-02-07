// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  Inclusions.h: General tRIBS Include Statements
**
***************************************************************************/

#ifndef INCLUSIONS_H
#define INCLUSIONS_H

// INCLUDED LIBRARY HEADER FILES

#ifdef ALPHA_64
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>

#elif defined LINUX_32
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <memory>

#elif defined WIN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory>

#else 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory>

#endif

#endif

//=========================================================================
//
//
//                   End of Inclusions.h  
//
//
//=========================================================================
