/*
 * RunPhaseII.cxx
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */


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

	fDetectorStatusFile += DetectorStatusFile;
	fDataKeysAnalysisFile += DataKeysAnalysisFile;
	fDataKeysAllFile += DataKeysAllFile;

	fDetectors = vector<DetectorPhaseII*>(0);
}


RunPhaseII::~RunPhaseII()
{
	for( auto det : fDetectors ) delete det;
	fDetectors.clear();
};


void RunPhaseII::AddDetector( string DetectorName, string DetectorType,
		uint DataChannel, uint MCChannel, bool OnlyACFlag, bool SwitchedOffFlag )
{
	for( auto det : fDetectors )
		if( det->GetDetectorName() == DetectorName )
		{
			cout << "WARNING: Detector already added " << DetectorName << endl;
			return;
		}

	DetectorPhaseII * det = DetectorPhaseII( DetectorName, DetectorType, DataChannel, MCChannel,
			OnlyACFlag, SwitchedOffFlag );

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

	ifstream detectorStatusFile( fGERDA_META_DATA + "/" + fDetectorStatusFile );
	string line, DetectorName, DetectorType;
	uint DataChannel, MCChannel;
	bool OnlyACFlag, SwitchedOffFlag;

	getline(detectorStatusFile, line);

	detectorStatusFile >> DetectorName >> DetectorType >> DataChannel >> MCChannel;

	detectorStatusFile >> line;

	if( line.back() == '1' )
	{
		AddDetector( DetectorName, DetectorType, DataChannel, MCChannel, 0, 1 );
	}
	else if( line.at( (int)line.end() - 2 )  == '1' )
	{
		AddDetector( DetectorName, DetectorType, DataChannel, MCChannel, 1, 0 );
	}
	else AddDetector( DetectorName, DetectorType, DataChannel, MCChannel);

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
