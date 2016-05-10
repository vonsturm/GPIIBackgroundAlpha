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
	SetDetectorAttributes( "default", "default", 100, 100, "OFF", 0, 0 );
};

DetectorPhaseII::DetectorPhaseII( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, string DetectorAnalysisStatus,
			double TotalMass, double ActiveVolumeFraction )
{
	SetDetectorAttributes( DetectorName, DetectorType, DataChannel, MCChannel, DetectorAnalysisStatus,
			TotalMass, ActiveVolumeFraction );
};

// destructor
DetectorPhaseII::~DetectorPhaseII(){};

// set parameters
void DetectorPhaseII::SetDetectorAttributes( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, string DetectorAnalysisStatus,
			double TotalMass, double ActiveVolumeFraction )
{
	fDetectorName = DetectorName;
	fDetectorType = DetectorType;

	fDetectorAnalysisStatus = DetectorAnalysisStatus;

	fDataChannel = DataChannel;
	fMCChannel = MCChannel;

	fMass = TotalMass;
	fAVfraction = ActiveVolumeFraction;
	fActiveMass = TotalMass * ActiveVolumeFraction;
}

string DetectorPhaseII::PrintInfo()
{
	string detinfo = fDetectorName;
	detinfo += " "; detinfo += fDetectorType; detinfo += " "; detinfo += fDataChannel; detinfo += " ";
	detinfo += fMCChannel; detinfo	+= " "; detinfo += fDetectorAnalysisStatus; detinfo	+= " ";
	detinfo += fMass; detinfo	+= " "; detinfo += fAVfraction; detinfo	+= " "; detinfo	+= fActiveMass;

	return detinfo;
}
