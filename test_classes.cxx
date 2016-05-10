/*
 * test_classes.cxx
 *
 *  Created on: May 9, 2016
 *      Author: sturm
 */

#include <string>

#include "RunPhaseII.h"

using namespace std;

int main( int argc, char* argv[] )
{
	RunPhaseII * run60 = new RunPhaseII( 60, "run0060-phy-detStatus.txt",
			"run0060-phy-analysis.txt", "run0060-phy-allFiles.txt" );

	run60->ParseDetectorStatusFile();
}

