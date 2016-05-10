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
	RunPhaseII * run53 = new RunPhaseII( 53, "run0053-phy-detStatus.txt",
			"run0053-phy-analysis.txt", "run0053-phy-allFiles.txt" );

	RunPhaseII * run54 = new RunPhaseII( 54, "run0054-phy-detStatus.txt",
			"run0054-phy-analysis.txt", "run0054-phy-allFiles.txt" );

	RunPhaseII * run55 = new RunPhaseII( 55, "run0055-phy-detStatus.txt",
			"run0055-phy-analysis.txt", "run0055-phy-allFiles.txt" );

	RunPhaseII * run56 = new RunPhaseII( 56, "run0056-phy-detStatus.txt",
			"run0056-phy-analysis.txt", "run0056-phy-allFiles.txt" );

	RunPhaseII * run57 = new RunPhaseII( 57, "run0057-phy-detStatus.txt",
			"run0057-phy-analysis.txt", "run0057-phy-allFiles.txt" );

	RunPhaseII * run58 = new RunPhaseII( 58, "run0058-phy-detStatus.txt",
			"run0058-phy-analysis.txt", "run0058-phy-allFiles.txt" );

	RunPhaseII * run59 = new RunPhaseII( 59, "run0059-phy-detStatus.txt",
			"run0059-phy-analysis.txt", "run0059-phy-allFiles.txt" );

	RunPhaseII * run60 = new RunPhaseII( 60, "run0060-phy-detStatus.txt",
			"run0060-phy-analysis.txt", "run0060-phy-allFiles.txt" );

	RunPhaseII * run61 = new RunPhaseII( 61, "run0061-phy-detStatus.txt",
			"run0061-phy-analysis.txt", "run0061-phy-allFiles.txt" );

}

