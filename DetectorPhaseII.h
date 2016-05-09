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
			uint DataChannel, uint MCChannel, bool OnlyACFlag = 0, bool SwitchedOffFlag = 0 );
	~DetectorPhaseII();

	void SetDetectorAttributes( std::string DetectorName, std::string DetectorType,
			uint DataChannel, uint MCChannel, bool OnlyACFlag = 0, bool SwitchedOffFlag = 0 );

	// Getters and Setters
	std::string GetDetectorName(){ return fDetectorName; };
	std::string GetDetectorType(){ return fDetectorType; };

	bool GetAnalysisFlag(){ return fAnalysisFlag; };
	bool GetOnlyACFlag(){ return fOnlyACFlag; };
	bool GetSwitchedOffFlag(){ return fSwitchedOffFlag; };

	uint GetDataChannel(){ return fDataChannel; };
	uint GetMCChannel(){ return fMCChannel; };

	//
	int CheckFlags( int verbose = 0 );

	std::string PrintInfo();


private:

	std::string fDetectorName;
	std::string fDetectorType;

	bool fAnalysisFlag;
	bool fOnlyACFlag;
	bool fSwitchedOffFlag;

	uint fDataChannel;
	uint fMCChannel;
};



#endif /* DETECTORPHASEII_H_ */
