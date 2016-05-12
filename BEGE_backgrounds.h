// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BEGE_BACKGROUNDS__H
#define __BAT__BEGE_BACKGROUNDS__H

// C/C++ includes
#include <vector>

// ROOT includes
#include "TH1D.h"
#include "TChain.h"

// BAT includes
#include <BAT/BCModel.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>


// This is a BEGE_backgrounds header file.
// Model source code is located in file BEGE_backgrounds/BEGE_backgrounds.cxx

// ---------------------------------------------------------
class BEGE_backgrounds : public BCModel
{
   public:

      // Constructors and destructor
      BEGE_backgrounds();
      BEGE_backgrounds(const char * name);
      ~BEGE_backgrounds();

      // Methods to overload, see file BEGE_backgrounds.cxx
      void DefineParameters();
      double LogAPrioriProbability(std::vector <double> parameters);
      double LogLikelihood(const std::vector <double> & parameters);
      // void MCMCIterationInterface();

      // my own methods
      void SetHistogramParameters(int hnumbins, double hemin, double hemax);
      int ReadData( std::string meta_filename = "" );
      int FillDataArray();
      int ReadMC();
      int AddMC(std::string name);
      int FillMCArrays();

      double getndets() {return f_ndets;};

      double EstimatePValue();
      void DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename);

 private:

      int f_ndets;
      int f_hnumbins;
      double f_hemin;
      double f_hemax;

      std::vector<TH1D*> f_hdata;
      TH1D* f_hdataSum;
      TH1D* f_hdataSum_fine;

      std::string f_MC_FileName;

      std::vector<TH1D*> f_MC;
      std::vector<TH1D*> f_MCfine;
      std::vector<TH1D*> f_MCall;
      std::vector<std::string> f_MCname;

      std::vector<double> f_vdata;
      std::vector<double> f_vMC;

      std::vector<double> f_lowerlimits;
      std::vector<double> f_upperlimits;
      
};
// ---------------------------------------------------------

#endif

