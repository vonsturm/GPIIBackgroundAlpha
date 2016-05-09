/*
 * DetectorPhaseII.cxx
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */


#include <iostream>

#include "DetectorPhaseII.h"

using namespace std;

// constructor
DetectorPhaseII::DetectorPhaseII()
{
	SetDetectorAttributes( "default", "default", 100, 100, "OFF" );
};

DetectorPhaseII::DetectorPhaseII( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, string DetectorAnalysisStatus )
{
	SetDetectorAttributes( DetectorName, DetectorType, DataChannel, MCChannel, DetectorAnalysisStatus );
};

// destructor
DetectorPhaseII::~DetectorPhaseII(){};

// set parameters
void DetectorPhaseII::SetDetectorAttributes( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, string DetectorAnalysisStatus )
{
	fDetectorName = DetectorName;
	fDetectorType = DetectorType;

	fDetectorAnalysisStatus = DetectorAnalysisStatus;

	fDataChannel = DataChannel;
	fMCChannel = MCChannel;
}

string DetectorPhaseII::PrintInfo()
{
	string detinfo = fDetectorName;
	detinfo += " "; detinfo += fDetectorType; detinfo += " "; detinfo += fDataChannel; detinfo += " ";
	detinfo += fMCChannel; detinfo	+= " "; detinfo += fDetectorAnalysisStatus;

	return detinfo;
}
