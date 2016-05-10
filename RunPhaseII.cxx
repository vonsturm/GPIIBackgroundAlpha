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
	ParseDetectorStatusFile();
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
		uint DataChannel, uint MCChannel, string DetectorAnalysisStatus )
{
	for( auto det : fDetectors )
		if( det->GetDetectorName() == DetectorName )
		{
			cout << "WARNING: Detector already added " << DetectorName << endl;
			return;
		}

	DetectorPhaseII * det = new DetectorPhaseII( DetectorName, DetectorType, DataChannel, MCChannel,
			DetectorAnalysisStatus );

	cout << "******************" << endl;
	cout << "Added detector " << DetectorName << " to Run " << fRunNumber << endl;
	cout << "******************" << endl;
	cout << "Type: " << DetectorType << endl;
	cout << "DataCH: " << DataChannel << " MCCH: " << MCChannel << endl;
	cout << "Analysis Status: " << DetectorAnalysisStatus << endl;

	fDetectors.push_back( det );

	return;
}


int RunPhaseII::ParseDetectorStatusFile()
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
		cout << "File is open " << statusFileName << endl;
	}

	string line, DetectorName, DetectorType, DetectorAnalysisStatus;
	uint DataChannel, MCChannel;

	getline(detectorStatusFile, line);

	int counterON = 0, counterOFF = 0, counterAConly = 0;

	while( !detectorStatusFile.eof() )
	{
		detectorStatusFile >> DetectorName >> DetectorType >> DataChannel;
		detectorStatusFile >> MCChannel >> DetectorAnalysisStatus;

		if( DetectorAnalysisStatus == "ON" ) counterON++;
		else if( DetectorAnalysisStatus == "OFF" ) counterOFF++;
		else if( DetectorAnalysisStatus == "AC-only" ) counterAConly++;

		AddDetector( DetectorName, DetectorType, DataChannel, MCChannel, DetectorAnalysisStatus );
	}

	cout << "Found " << fDetectors.size() << " detectors in Run" << endl;
	cout << "ON " << counterON << " OFF " << counterOFF << " AC only " << counterAConly << endl;

	return 0;
}


int RunPhaseII::PrintDetectorStatusFile()
{
	string outputFilename = fGERDA_META_DATA;
	outputFilename += Form( "/run%04d-phy-detStatus.txt", fRunNumber );

	ofstream detectorStatusFile( fGERDA_META_DATA + "/" );
	cout << "Creating file in current directory." << endl;

	detectorStatusFile << "Detector Type DataCh MCCh AnalysisFlags" << endl;

	for( auto det : fDetectors ) detectorStatusFile << det->PrintInfo() << endl;

	return 0;
}
