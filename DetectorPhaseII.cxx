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
	SetDetectorAttributes( "default", "default", 100, 100, false, true );
};

DetectorPhaseII::DetectorPhaseII( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, bool OnlyACFlag,
			bool SwitchedOffFlag )
{
	SetDetectorAttributes( DetectorName, DetectorType, DataChannel, MCChannel,
			OnlyACFlag, SwitchedOffFlag );
};

// destructor
DetectorPhaseII::~DetectorPhaseII(){};

// set parameters
void DetectorPhaseII::SetDetectorAttributes( string DetectorName, string DetectorType,
			uint DataChannel, uint MCChannel, bool OnlyACFlag, bool SwitchedOffFlag )
{
	fDetectorName = DetectorName;
	fDetectorType = DetectorType;

	fSwitchedOffFlag = SwitchedOffFlag;
	fOnlyACFlag = OnlyACFlag && !SwitchedOffFlag;
	fAnalysisFlag = ! ( OnlyACFlag || SwitchedOffFlag );

	fDataChannel = DataChannel;
	fMCChannel = MCChannel;
}

// check the analysis flags !!!Only one flag has to be set!!!
int DetectorPhaseII::CheckFlags( int verbose )
{
	if( ( fAnalysisFlag && fOnlyACFlag ) || ( fAnalysisFlag && fSwitchedOffFlag ) ||
			( fSwitchedOffFlag && fOnlyACFlag ) )
	{
		if( verbose > 0 )
		{
			cout << "ERROR Analysis Flags: \n\t fAnalysisFlag = " << fAnalysisFlag;
			cout << "\n\t fOnlyACFlag = " << fOnlyACFlag;
			cout << "\n\t fSwitchedOffFlag = " << fSwitchedOffFlag << endl;
			cout << "Only one flag can be set at a time!" << endl;
		}
		return 1;
	}
	else if( ! ( fAnalysisFlag || fOnlyACFlag || fSwitchedOffFlag ) )
	{
		if( verbose > 0 )
		{
			cout << "ERROR Analysis Flags: \n\t fAnalysisFlag = " << fAnalysisFlag;
			cout << "\n\t fOnlyACFlag = " << fOnlyACFlag;
			cout << "\n\t fSwitchedOffFlag = " << fSwitchedOffFlag << endl;
			cout << "One flag can be set at a time!" << endl;
		}
		return 1;
	}

	if( verbose > 0 )
	{
		cout << "Analysis Flags are fine: \n\t fAnalysisFlag = " << fAnalysisFlag;
		cout << "\n\t fOnlyACFlag = " << fOnlyACFlag;
		cout << "\n\t fSwitchedOffFlag = " << fSwitchedOffFlag << endl;
	}

	return 0;
}


string DetectorPhaseII::PrintInfo()
{
	string detinfo = fDetectorName;
	detinfo += " "; detinfo += fDetectorType; detinfo += " "; detinfo += fDataChannel; detinfo += " ";
	detinfo += fMCChannel; detinfo	+= " "; detinfo += fAnalysisFlag; detinfo += fOnlyACFlag;
	detinfo += fSwitchedOffFlag;

	return detinfo;
}
