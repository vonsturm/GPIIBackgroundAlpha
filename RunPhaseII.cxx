/*
 * RunPhaseII.cxx
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "TString.h"

#include "RunPhaseII.h"


using namespace std;

RunPhaseII::RunPhaseII()
{
	SetRunSetup( 0 );
}

RunPhaseII::RunPhaseII( int RunNumber )
{
	SetRunSetup( RunNumber );
}

RunPhaseII::RunPhaseII( int RunNumber, string DetectorStatusFile,
		string DataKeysAnalysisFile, string DataKeysAllFile )
{
	SetRunSetup( RunNumber, DetectorStatusFile, DataKeysAnalysisFile, DataKeysAllFile );
	ParseDetectorStatusFile( 1 );
	SetExposures();
};

void RunPhaseII::SetRunSetup( int RunNumber, string DetectorStatusFile,
		string DataKeysAnalysisFile, string DataKeysAllFile )
{
	fRunNumber = RunNumber;

	SetGERDA_META_DATA();

	fDetectorStatusFile = DetectorStatusFile;
	fDataKeysAnalysisFile = DataKeysAnalysisFile;
	fDataKeysAllFile = DataKeysAllFile;

	fDetectors = vector<DetectorPhaseII*>(0);
}


RunPhaseII::~RunPhaseII()
{
	for( auto det : fDetectors ) delete det;
	fDetectors.clear();
};


void RunPhaseII::AddDetector( string DetectorName, string DetectorType,
		uint DataChannel, uint MCChannel, string DetectorAnalysisStatus, int verbose )
{
	for( auto det : fDetectors )
	{
		if( det->GetDetectorName() == DetectorName )
		{
			cout << "WARNING: Detector already added " << DetectorName << endl;
			return;
		}
	}

	if( verbose > 0 )
	{
		cout << "******************" << endl;
		cout << "Adding detector " << DetectorName << " to Run " << fRunNumber << endl;
		cout << "******************" << endl;
		cout << "Type: " << DetectorType << endl;
		cout << "DataCH: " << DataChannel << " MCCH: " << MCChannel << endl;
		cout << "Analysis Status: " << DetectorAnalysisStatus << endl;
	}

	DetectorPhaseII * det = new DetectorPhaseII( DetectorName, DetectorType, DataChannel, MCChannel,
			DetectorAnalysisStatus );

	fDetectors.push_back( det );

	if( verbose > 0 ) cout << "******************" << endl;

	return;
}


DetectorPhaseII* RunPhaseII::GetDetectorInDataChannel( int channel )
{
	if( fDetectors.size() <= 0 )
	{
		cout << "ERROR: No detectors found." << endl;
		return NULL;
	}

	for( auto det : fDetectors )
	{
		if( det->GetDataChannel() == channel )
			return det;
	}
}


int RunPhaseII::ParseDetectorStatusFile( int verbose )
{
	if( fGERDA_META_DATA.empty() )
	{
		cout << "WARNING: Environment variable GERDA_META_DATA not set." << endl;
		return 1;
	}

	string statusFileName = fGERDA_META_DATA; statusFileName += fDetectorStatusFile;
	ifstream detectorStatusFile( statusFileName );

	if( !detectorStatusFile.is_open() )
	{
		cout << "Could not open file " << statusFileName << endl;
		return 1;
	}
	else
	{
		if( verbose > 1 ) cout << "File is open " << statusFileName << endl;
	}

	string line, DetectorName, DetectorType, DetectorAnalysisStatus;
	uint DataChannel, MCChannel;

	getline(detectorStatusFile, line);

	detectorStatusFile >> line >> fLiveTime;
	for( int i = 0; i < 4; i++ ) detectorStatusFile >> line;

	getline(detectorStatusFile, line);

	int counterON = 0, counterOFF = 0, counterAConly = 0;

	detectorStatusFile >> DetectorName >> DetectorType >> DataChannel;
	detectorStatusFile >> MCChannel >> DetectorAnalysisStatus;

	while( !detectorStatusFile.eof() )
	{
		if( DetectorAnalysisStatus == "ON" ) counterON++;
		else if( DetectorAnalysisStatus == "OFF" ) counterOFF++;
		else if( DetectorAnalysisStatus == "AC-only" ) counterAConly++;

		AddDetector( DetectorName, DetectorType, DataChannel, MCChannel, DetectorAnalysisStatus );

		detectorStatusFile >> DetectorName >> DetectorType >> DataChannel;
		detectorStatusFile >> MCChannel >> DetectorAnalysisStatus;
	}

	detectorStatusFile.close();

	if( verbose > 0 )
	{
		cout << "******************" << endl;
		cout << "Found " << fDetectors.size() << " detectors in Run " << fRunNumber << endl;
		cout << "ON " << counterON << " OFF " << counterOFF << " AC only " << counterAConly << endl;
		cout << "******************" << endl;
	}

	return 0;
}


int RunPhaseII::PrintDetectorStatusFile( string filename )
{
	ifstream test( filename );
	if( test.is_open() )
	{
		cout << "ERROR: File exists " << filename << endl;
		return 1;
	}
	else
		test.close();

	ofstream detectorStatusFile( filename );
	cout << "Creating file in current directory." << endl;

	detectorStatusFile << "Detector Type DataCh MCCh AnalysisFlags Mass fAV ActiveMass" << endl;

	for( auto det : fDetectors ) detectorStatusFile << det->PrintInfo() << endl;

	detectorStatusFile.close();

	return 0;
}


void RunPhaseII::SetExposures()
{
	fExposure = 0;
	fExposureBEGE = 0;
	fExposureCOAX = 0;

	for( auto det : fDetectors )
	{
		string type = det->GetDetectorType();
		string status = det->GetDetectorAnalysisStatus();

		double mass = det->GetTotalMass();

		if( status == "ON" )
		{
			fExposure += mass;

			if( type == "enrBEGe" )
				fExposureBEGE += mass;
			else if( type == "enrCoax" )
				fExposureCOAX += mass;
		}
	}

	fExposure *= fLiveTime / 365.2425 / 1000.;
	fExposureBEGE *= fLiveTime / 365.2425 / 1000.;
	fExposureCOAX *= fLiveTime / 365.2425 / 1000.;
}
