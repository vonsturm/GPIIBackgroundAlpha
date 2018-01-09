// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__GPIIBACKGROUNDALPHA__H
#define __BAT__GPIIBACKGROUNDALPHA__H

// C/C++ includes
#include <vector>

// ROOT includes
#include "TH1D.h"
#include "TChain.h"

// BAT includes
#include <BAT/BCModel.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>


// This is a GPIIBackgroundAlpha header file.
// Model source code is located in file GPIIBackgroundAlpha/GPIIBackgroundAlpha.cxx

// ---------------------------------------------------------
class GPIIBackgroundAlpha : public BCModel
{
   public:

      // Constructors and destructor
      GPIIBackgroundAlpha();
      GPIIBackgroundAlpha(const char * name);
      ~GPIIBackgroundAlpha();

      // Methods to overload, see file GPIIBackgroundAlpha.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector <double> & parameters);
      // void MCMCIterationInterface();

      // my own methods
      void SetHistogramParameters(int hnumbins, double hemin, double hemax);
      int ReadDataEnrBEGe( std::vector<int> runlist );
      int ReadDataEnrCoax( std::vector<int> runlist );
      int ReadDataNatCoax( std::vector<int> runlist );
      int FillDataArray();
      int ReadMCAlpha();
      int AddMC( std::string name );
      int AddMCSingle( std::string name, std::string histoname );
      int FillMCArrays();

      std::vector<std::string> GetMCParNames(){ return f_MCname; };

      double getndets() {return f_ndets;};

      double EstimatePValue();
      void DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename);

 private:

      int f_ndets;
      int f_hnumbins;
      double f_hemin;
      double f_hemax;
      std::vector<int> f_binsToSkip;

      std::vector<double> fDetectorLiveTime;
      std::vector<int> fDetectorDynamicRange;

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
