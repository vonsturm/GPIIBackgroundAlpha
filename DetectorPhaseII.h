/*
 * DetectorPhaseII.h
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */

#ifndef DETECTORPHASEII_H_
#define DETECTORPHASEII_H_

#include <string>

class DetectorPhaseII
{
public:

	DetectorPhaseII();
	DetectorPhaseII( std::string DetectorName, std::string DetectorType,
			uint DataChannel, uint MCChannel, std::string DetectorAnalysisStatus = "ON" );
	~DetectorPhaseII();

	void SetDetectorAttributes( std::string DetectorName, std::string DetectorType,
			uint DataChannel, uint MCChannel, std::string DetectorAnalysisStatus = "ON" );

	// Getters and Setters
	std::string GetDetectorName(){ return fDetectorName; };
	std::string GetDetectorType(){ return fDetectorType; };

	std::string GetDetectorAnalysisStatus(){ return fDetectorAnalysisStatus; };

	uint GetDataChannel(){ return fDataChannel; };
	uint GetMCChannel(){ return fMCChannel; };

	double GetTotalMass(){ return fMass; };
	double GetAVfraction(){ return fAVfraction; };
	double GetActiveMass(){ return fActiveMass; };

	int GetDynamicRange(){ return fEndDynamicRange; };

	//
	void ReadParametersFromDatabase();

	std::string PrintInfo();


private:

	std::string fDetectorName;
	std::string fDetectorType;

	std::string fDetectorAnalysisStatus;

	uint fDataChannel;
	uint fMCChannel;

	double fMass;
	double fAVfraction;
	double fActiveMass;

	int fEndDynamicRange;
};



#endif /* DETECTORPHASEII_H_ */
