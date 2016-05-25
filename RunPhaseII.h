/*
 * Run_PhaseII.h
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */

#ifndef RUNPHASEII_H_
#define RUNPHASEII_H_

#include <vector>
#include <cstdlib>

#include "DetectorPhaseII.h"

class RunPhaseII
{
public:

	RunPhaseII();
	RunPhaseII( int RunNumber );
	RunPhaseII( int RunNumber, std::string DetectorStatusFile, std::string DataKeysAnalysisFile,
			std::string DataKeysAllFile);
	~RunPhaseII();

	void SetRunSetup( int RunNumber, std::string DetectorStatusFile = "default.txt",
			std::string DataKeysAnalysisFile = "default.txt", std::string DataKeysAllFile  = "default.txt");

	void SetGERDA_DETECTOR_STATUS()
	{
		fGERDA_DETECTOR_STATUS = getenv("GERDA_DETECTOR_STATUS");
		fGERDA_DETECTOR_STATUS += "/";
	};

	void SetGERDA_DATA_SETS()
	{
		fGERDA_DATA_SETS = getenv("GERDA_DATA_SETS");
		fGERDA_DATA_SETS += "/";
	};

	// Getters and Setters
	int GetRunNumber(){ return fRunNumber; };

	double GetLiveTime(){ return fLiveTime; };		// in days
	double GetExposure(){ return fExposure; };		// in kg*yr
	double GetExposureBEGE(){ return fExposureBEGE; }; 	// in kg*yr
	double GetExposureCOAX(){ return fExposureCOAX; }; 	// in kg*yr

	std::string GetGERDA_DETECTOR_STATUS(){ return fGERDA_DETECTOR_STATUS; };
	std::string GetGERDA_DATA_SETS(){ return fGERDA_DATA_SETS; };

	std::string GetDetectorStatusFile(){ return fGERDA_DETECTOR_STATUS + "/" + fDetectorStatusFile; };
	std::string GetDataKeysAllFile(){ return fGERDA_DATA_SETS + "/" + fDataKeysAllFile; };
	std::string GetDataKeysAnalysisFile(){ return fGERDA_DATA_SETS + "/" + fDataKeysAnalysisFile; };

	std::vector<DetectorPhaseII*> GetDetectors(){ return fDetectors; };
	DetectorPhaseII* GetDetectorInDataChannel( int channel );

	//
	void AddDetector( std::string DetectorName, std::string DetectorType,
			uint DataChannel, uint MCChannel, std::string DetectorAnalysisStatus = "ON", int verbose = 0 );
	int ParseDetectorStatusFile( int verbose = 0 );
	int PrintDetectorStatusFile( std::string filename );

private:

	int fRunNumber;

	double fLiveTime;		// in days
	double fExposure; 		// in kg*yr
	double fExposureBEGE; 	// in kg*yr
	double fExposureCOAX; 	// in kg*yr

	std::string fGERDA_DETECTOR_STATUS;
	std::string fGERDA_DATA_SETS;

	std::string fDataKeysAllFile;
	std::string fDataKeysAnalysisFile;
	std::string fDetectorStatusFile;

	std::vector<DetectorPhaseII*> fDetectors;

	void SetExposures();
};

#endif /* RUNPHASEII_H_ */
