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

	void SetGERDA_META_DATA()
	{
		fGERDA_META_DATA = getenv("GERDA_META_DATA");
		fGERDA_META_DATA += "/";
	};

	// Getters and Setters
	int GetRunNumber(){ return fRunNumber; };
	std::string GetGERDA_META_DATA(){ return fGERDA_META_DATA; };
	std::string GetDataKeysAllFile(){ return fGERDA_META_DATA + "/" + fDataKeysAllFile; };
	std::string GetDataKeysAnalysisFile(){ return fGERDA_META_DATA + "/" + fDataKeysAnalysisFile; };
	std::string GetDetectorStatusFile(){ return fGERDA_META_DATA + "/" + fDetectorStatusFile; };
	std::vector<DetectorPhaseII*> GetDetectors(){ return fDetectors; };

	//
	void AddDetector( std::string DetectorName, std::string DetectorType,
			uint DataChannel, uint MCChannel, std::string DetectorAnalysisStatus = "ON",
			double TotalMass, double ActiveVolumeFraction );
	int ParseDetectorStatusFile();
	int PrintDetectorStatusFile();

private:

	int fRunNumber;

	double fLiveTime;		// in days
	double fExposure; 		// in kg*yr
	double fExposureBEGE; 	// in kg*yr
	double fExposureCOAX; 	// in kg*yr

	std::string fGERDA_META_DATA;

	std::string fDataKeysAllFile;
	std::string fDataKeysAnalysisFile;
	std::string fDetectorStatusFile;

	std::vector<DetectorPhaseII*> fDetectors;
};

#endif /* RUNPHASEII_H_ */
