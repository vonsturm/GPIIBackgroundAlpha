/*
 * DetectorPhaseII.cxx
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */


#include <iostream>
#include <fstream>
#include <cmath> 
#include <ctgmath>

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
	ReadParametersFromDatabase();
};

// destructor
DetectorPhaseII::~DetectorPhaseII(){};

// set parameters that can change with each run
void DetectorPhaseII::SetDetectorAttributes( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, string DetectorAnalysisStatus )
{
	fDetectorName = DetectorName;
	fDetectorType = DetectorType;

	fDetectorAnalysisStatus = DetectorAnalysisStatus;

	fDataChannel = DataChannel;
	fMCChannel = MCChannel;

	return;
}

// Read detector parameters which cannot change like mass and active volume fraction
void DetectorPhaseII::ReadParametersFromDatabase()
{
	string DETECTOR_DATABASE;

	if( getenv("DETECTOR_DATABASE") )
		DETECTOR_DATABASE = getenv("DETECTOR_DATABASE");
	else
	{
		cout << "ERROR: DETECTOR_DATABASE environment variable not set." << endl;
		return;
	}

	ifstream file( DETECTOR_DATABASE );

	if( !file.is_open() )
	{
		cout << "ERROR: Detector Database File not found " << DETECTOR_DATABASE << endl;
		return;
	}

	string dummy = "";
	double ddummy;

	while( dummy != fDetectorName && !file.eof() ) file >> dummy;

	if( file.eof() && dummy != fDetectorName )
	{
		cout << "Detector not found in DB " << fDetectorName << endl;
		return;
	}

	file >> fMass;
	file >> fAVfraction;
	file >> ddummy;
	fEndDynamicRange = floor( ddummy );

	fActiveMass = fMass * fAVfraction/100.;

	file.close();

	return;
}


string DetectorPhaseII::PrintInfo()
{
	string detinfo = fDetectorName;
	detinfo += " "; detinfo += fDetectorType; detinfo += " "; detinfo += fDataChannel; detinfo += " ";
	detinfo += fMCChannel; detinfo	+= " "; detinfo += fDetectorAnalysisStatus; detinfo	+= " ";
	detinfo += fMass; detinfo	+= " "; detinfo += fAVfraction; detinfo	+= " "; detinfo	+= fActiveMass;

	return detinfo;
}
