// MeshBuilder
//
// Copyright (c) 2024. tRIBS Developers
//
// See LICENSE file in the project root for full license information.

/***************************************************************************
**
**  Definitions.h: General tRIBS Definitions
**
***************************************************************************/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <string>

//=========================================================================
//
//
//                  Section 1: Definitions and Macros
//
//
//=========================================================================

// Definitions

const std::string  VERSION	= "MeshBuilder 1.0";

// File reading constants
const int kMaxNameLength	= 1024;
const int kMaxNameSize		= 1024;
const char kCommentMark		= '#';

// Types of nodes
const int kNonBoundary		= 0;
const int kClosedBoundary	= 1;
const int kOpenBoundary		= 2;
const int kStream		= 3;

// Flow network
int const kFlooded		= 1;  	
int const kNotFlooded		= 0;  	
int const kCurrentLake		= 2;
int const kSink			= 3;  	
int const kOutletFlag		= 4;  	
int const kOutletPreFlag	= 5;  	

int const kFlowNotAllowed	= 0;
int const kFlowAllowed		= 1;

// Numerical constants
const double PI			= 3.1415926;
const double EPS		= 2.2204e-16;
const double THRESH		= 1e-6;

const double kMaxElevation	= 100000.0;

#endif

//=========================================================================
//
//
//                         End of Definitions.h
//
//
//=========================================================================
