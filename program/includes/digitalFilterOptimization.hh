#ifndef digitalFilterOptimization_h
#define digitalFilterOptimization_h 1

#include "global.h"
#include "digitalFilters.h"
#include "trap_parameter.hh"
#include "resolutionStrip.hh"
#include "dssdData.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TList.h"
#include "TDirectory.h"

#include <stdio.h>//C library to perform Input/Output operations
#include <algorithm>    // std::sort
#include <exception> 
#include <fstream>
#include <iomanip> //setw,setprecision
//#include <Math/Interpolator.h>



class digitalFilterOptimization
{
	private:
		myGlobal *s1;
		digitalFilters *filter;
		dssdData *dData;
		//variables for the branches of the input tree
		ULong64_t ftime;
		UInt_t feventNo;
		UShort_t ftraceSize;
		UShort_t *ftrace;
		UShort_t fgain;
		UShort_t fboardID;
		UShort_t fchannelID;
		//variables for the branches of the trap tree
		Double_t trapAmplitude; 
		UShort_t fboardIdx;
		UShort_t k, m,l, kIdx, mIdx;
		//
		UShort_t kStart, kStop, kStep;
		UShort_t mStart, mStop, mStep;
		UShort_t ktimes;
		UShort_t mtimes;

		TFile* file_tree;
		TFile* file_spectra;
		TTree* trapTree; // store trapezoidal amplitudes
		TRandom3 * rand;

		//create histograms
		int nPeaks;
		Double_t *fPositionX;
		Double_t *fPositionY;
		//Fit functions
		TF1 *fLinear;
		TF1 *G0;
		TF1 *G1;
		TF1 *G2;
		// graph for calibration
		TGraph* gr_cal;
		Double_t Energy[2];
		TSpectrum *spectrum;
		TCanvas* canvas;
		TCanvas* canvas2;


		TList* list_res;

		TH1F ***** hRaw;
		TH1F ***** hCalib;
		TH2F *** hResolution;

		TList* list_calib;
		TH2F *** hCalibStrip;
		TH2F* h_m;
		TH2F* h_c;

		//calibration parameters
		Double_t ****gain;
		Double_t ****offset;
		void initialize_TObjects();
		void delete_TObjects();

		
	public:
		digitalFilterOptimization(UShort_t *);
		~digitalFilterOptimization();

		void optimize_trapezoidal_filter_parameters(const char *, double, double);
};

#endif
