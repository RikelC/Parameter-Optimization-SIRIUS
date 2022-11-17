#include "digitalFilterOptimization.hh"
#define PBSTR "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
#define PBWIDTH 60
bool compareFWHM(trap_parameter i,trap_parameter j) { return (i.FWHM < j.FWHM); }


bool compareFWHM2(resolution_strip i,resolution_strip j) { return (i.FWHM < j.FWHM); }


int extract_run_number_from_filename(const char * filename){
	std::string runN =  filename;
	std::size_t pos = runN.find("r-");
	if(pos!=std::string::npos)runN.erase(0, pos+2);
	pos = runN.length();
	if(pos!=std::string::npos)runN.erase(4, pos);
	return std::stoi(runN);
}

int extract_subrun_number_from_filename( const char * filename){
	std::string runN =  filename;
	std::size_t pos = runN.find("s-");
	if(pos!=std::string::npos)runN.erase(0, pos+2);
	pos = runN.length();
	size_t pos1 = runN.find_first_of("-Numexo2");
	if(pos1!=std::string::npos)runN.erase(pos1, pos);
	return std::stoi(runN);
}


	template<class T>
inline static void printProgress (ULong64_t i, T* t )
{
	static Int_t percentage = 0;

	static Int_t noEntries_dividedby100 = t->GetEntries()/100;
	if( percentage > 100) {
		percentage = 0;
		noEntries_dividedby100 = t->GetEntries()/100;
	}
	if(i == percentage * noEntries_dividedby100){
		Int_t val = static_cast<Int_t>(percentage);
		Int_t lpad = static_cast<Int_t>(percentage * PBWIDTH / 100);
		Int_t rpad = PBWIDTH - lpad;
		printf ("\rProgress:%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
		fflush (stdout);
		percentage++;
	}

}

digitalFilterOptimization::digitalFilterOptimization(UShort_t * par)
{
	s1 = myGlobal::getInstance();
	dData = new dssdData();
	filter = new digitalFilters();
	// Parameter range	
	kStart= par[0];
	kStop = par[1];
	kStep = par[2];
	mStart = par[3];
	mStop = par[4];
	mStep = par[5];

	ktimes = (kStop - kStart)/kStep + 1;
	mtimes = (mStop - mStart)/mStep + 1;
	nPeaks = 3;
	fPositionX = new Double_t[nPeaks];
	fPositionY = new Double_t[nPeaks];
	//3 alpha source 5156.6----241Am:  5485.8----244Cm: 5804.8
	ftrace =  new UShort_t[s1->TRACE_SIZE];
	Energy[0]= 5156.6;
	Energy[1]= 5804.8;
}
digitalFilterOptimization::~digitalFilterOptimization()
{
	delete dData;
	delete filter;
	delete[] fPositionX;
	delete[] fPositionY;
	delete [] ftrace;
}

void digitalFilterOptimization::initialize_TObjects()
{

	//calibration parameters
	gain = new Double_t***[s1->NBOARDS_DSSD];
	offset = new Double_t***[s1->NBOARDS_DSSD];
	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		gain[i] = new Double_t**[NCHANNELS];
		offset[i] = new Double_t**[NCHANNELS];
		for(UShort_t j = 0; j < NCHANNELS;j++){
			gain[i][j] = new Double_t*[ktimes];
			offset[i][j] = new Double_t*[ktimes];
			for(UShort_t k1 =0; k1 < ktimes; k1++){
				gain[i][j][k1] = new Double_t[mtimes];
				offset[i][j][k1] = new Double_t[mtimes];
				for(UShort_t m1 = 0; m1 < mtimes; m1++ ){
					gain[i][j][k1][m1] = 0.;
					offset[i][j][k1][m1] = 0.;
				}
			}
		}
	}

	//Fit functions
	fLinear = new TF1("fLinear","pol1",0,10000);
	G0    = new TF1("G0","gaus",0,10000);
	G1    = new TF1("G1","gaus",0,10000);
	G2    = new TF1("G2","gaus",0,10000);
	// graph for calibration
	gr_cal = new TGraph();
	//3 alpha source 5156.6----241Am:  5485.8----244Cm: 5804.8
	const Int_t nbins = 2500;
	Double_t xmin     = 0;
	Double_t xmax     = 2500.;

	spectrum = new TSpectrum(nPeaks);
	canvas = new TCanvas;
	canvas2 = new TCanvas;
	list_res = new TList;

	hRaw = new TH1F****[s1->NBOARDS_DSSD];
	hCalib = new TH1F****[s1->NBOARDS_DSSD];
	hResolution = new TH2F**[s1->NBOARDS_DSSD];

	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		UShort_t board =  s1->boardList_DSSD[i];
		hRaw[i] = new TH1F***[NCHANNELS];
		hCalib[i] = new TH1F***[NCHANNELS];
		hResolution[i] = new TH2F*[NCHANNELS];
		for(UShort_t j = 0; j < NCHANNELS;j++){
			hResolution[i][j] = new TH2F(Form("hRes_b%d_ch%d",board, j),Form("hRes_b%d_ch%d ; k; m",board, j), ktimes, kStart, kStop+kStep, mtimes, mStart, mStop+mStep);
			list_res->Add(hResolution[i][j]);
			hRaw[i][j]= new TH1F**[ktimes];
			hCalib[i][j]= new TH1F**[ktimes];
			for(UShort_t k1 =0; k1 < ktimes; k1++){
				hRaw[i][j][k1]= new TH1F*[mtimes];
				hCalib[i][j][k1]= new TH1F*[mtimes];
				for(UShort_t m1 = 0; m1 < mtimes; m1++ ){
					//convert to  m,k values
					UShort_t k_val = kStart + kStep * k1;
					UShort_t m_val = mStart + mStep * m1;
					// cout<<"i "<<i<<"  j "<<j<<"  k1 "<<k1<<" m1 "<<m1<<endl;
					//  cout<<"k val "<<k_val<<"  m val "<<m_val<<endl;
					hRaw[i][j][k1][m1]= new TH1F(Form("hRaw_b%d_c%d_k%d_m%d", board, j, k_val, m_val),Form("hRaw_boardId%d_channel%d_k%d_m%d", board,j,k_val,m_val),nbins,xmin,xmax);
					hCalib[i][j][k1][m1]= new TH1F(Form("hCalib_b%d_c%d_k%d_m%d", board,j, k_val, m_val),Form("hCalib_boardId%d_channel%d_k%d_m%d", board,j,k_val,m_val),1500,0,6000);
				}
			}
		}
	} 
	//calibrated spectra
	list_calib = new TList;
	hCalibStrip = new TH2F**[ktimes];
	for(UShort_t k1 =0; k1 < ktimes; k1++){
		hCalibStrip[k1]= new TH2F*[mtimes];
		for(UShort_t m1 = 0; m1 < mtimes; m1++ ){
			//convert to  m,k values
			UShort_t k_val = kStart + kStep * k1;
			UShort_t m_val = mStart + mStep * m1;
			hCalibStrip[k1][m1] = new TH2F(Form("hCalib_k%d_m%d", k_val, m_val), Form("hCalib_k%d_m%d", k_val, m_val),1500,0,6000,s1->NSTRIPS_DSSD,0,s1->NSTRIPS_DSSD);
			list_calib->Add( hCalibStrip[k1][m1]);
		}
	}

	if(ktimes > 1 && mtimes > 1){
		h_m = new TH2F("h_m", "calib par: m;k;m",ktimes, kStart, kStop, mtimes, mStart, mStop);
		h_c = new TH2F("h_c", "calib par: c;k;m",ktimes, kStart, kStop, mtimes, mStart, mStop);
	}
}

void digitalFilterOptimization::delete_TObjects()
{
	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		for(UShort_t j = 0; j < NCHANNELS;j++){
			delete [] gain[i][j];
			delete [] offset[i][j];

		}
	}
	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		delete [] gain[i];
		delete [] offset[i];
	}

	delete [] gain;
	delete [] offset;

	delete gr_cal;
	delete spectrum;
	delete canvas;
	delete canvas2;
	delete G0; delete G1; delete G2; delete fLinear;
	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		delete [] hResolution[i];
		for(UShort_t j = 0; j < NCHANNELS;j++){
			for(UShort_t k1 =0; k1 < ktimes; k1++){
				delete [] hRaw[i][j][k1];
				delete [] hCalib[i][j][k1];
			}
		}
	}
	delete [] hResolution;

	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		for(UShort_t j = 0; j < NCHANNELS;j++){
			delete [] hRaw[i][j];
			delete [] hCalib[i][j];        
		}
	}
	for(UShort_t i = 0; i < s1->NBOARDS_DSSD;i++){
		delete [] hRaw[i];
		delete [] hCalib[i];
	}

	delete [] hRaw;
	delete [] hCalib;
	for(UShort_t k1 =0; k1 < ktimes; k1++){
		delete[] hCalibStrip[k1];

	}
	delete[] hCalibStrip;
	list_res->Clear();
	list_calib->Clear();
	delete list_res;
	delete list_calib;
}
void digitalFilterOptimization::optimize_trapezoidal_filter_parameters(const char * file, double sr1, double sr2)
{
	//Extract run numbers

	std::string inDir = "/data/siriusX/test/acquisition/RootFiles/";
	std::string inFileName = inDir+std::string(file);

	int runNo = extract_run_number_from_filename( file);

	int subRunNo =extract_subrun_number_from_filename(file);
	// open a text file to save the minimization parameters
	TString minParFile;
	minParFile.Form("../results/min_trapezoidal_parameters_run_%d_r%d_%s.txt",  runNo, subRunNo, s1->filterAlgorithm.c_str());
	ofstream myfile (minParFile.Data());
	myfile << "boardID"<<setw(20)<<"channelID"<<setw(20)<<"k"<<setw(20)<<"m\n";

	//create an output file to save the spectra
	TString oFileName;
	oFileName.Form("../results/Spectra_run_%d_r%d_kStart%dkStop%dkStep%d_mStart%dmStop%dmStep%d_%s.root", runNo, subRunNo, kStart,kStop, kStep, mStart,mStop, mStep, s1->filterAlgorithm.c_str());

	file_spectra= new TFile(oFileName,"RECREATE");
	//create histograms, functions and lists
	initialize_TObjects();

	vector< trap_parameter> myvector[s1->NBOARDS_DSSD][NCHANNELS];
	//check if trapezoidal amplitude file is already created
	char FileName[200];
	sprintf(FileName,"../results/run_%d_r%d_kStart%dkStop%dkStep%d_mStart%dmStop%dmStep%d_%s.root", runNo, subRunNo, kStart,kStop, kStep, mStart,mStop, mStep, s1->filterAlgorithm.c_str());


	file_tree= new TFile(FileName,"RECREATE");
	trapTree = new TTree("trapData", "trapData");
	trapTree->Branch("boardId", &fboardID,"boardId/s");
	trapTree->Branch("boardIdx", &fboardIdx,"boardIdx/s");
	trapTree->Branch("channelId", &fchannelID,"channelId/s");
	trapTree->Branch("trapAmplitude", &trapAmplitude,"trapAmplitude/D");
	trapTree->Branch("k", &k,"k/s");
	trapTree->Branch("m", &m,"k/s");
	trapTree->Branch("kIdx", &kIdx,"kIdx/s");
	trapTree->Branch("mIdx", &mIdx,"mIdx/s");

	TFile* inFile = new TFile(inFileName.data(), "READ");

	cout<<inFileName<<endl;
	if(inFile->IsZombie()){
		//Close and clear memory
		file_spectra->cd();
		file_spectra->Delete("T;*");
		file_spectra->Close();
		remove(oFileName.Data());
		delete_TObjects();
		delete file_spectra;
		file_tree->cd();
		file_tree->Delete("T;*");
		file_tree->Close();
		delete trapTree;
		delete file_tree;
		std::cerr << "ERROR: could not open input file\n";
		std::terminate();

	}
	else{


		cout<< "file opened "<<inFileName<<endl;
		TTree * inTree = (TTree*) inFile->Get("rawDataTree");
		inTree->SetBranchAddress("Time", &ftime);
		inTree->SetBranchAddress("EventNo", &feventNo);
		inTree->SetBranchAddress("TraceSize", &ftraceSize);
		inTree->SetBranchAddress("Trace", ftrace);
		inTree->SetBranchAddress("Gain", &fgain);
		inTree->SetBranchAddress("BoardID", &fboardID);
		inTree->SetBranchAddress("ChannelID", &fchannelID);
		inTree->SetBranchStatus("*",1);
		// number of entries in the tree
		Int_t nEntries = inTree->GetEntries();
		//Start event loops
		// nEntries = 1;
		cout<<"Reading file.............."<<endl;
		cout<<"number of entries in the file: "<<nEntries<<endl;

		for(Int_t entry = 0; entry < nEntries; entry++){
			printProgress (entry,  inTree );
			inTree->GetEntry(entry);
			fboardIdx = s1->boardIndex_DSSD[fboardID];
			dData->set_channelID(fchannelID);
			dData->set_boardID(fboardID);
			dData->set_boardIdx(fboardIdx);
			dData->set_timestamp( ftime);
			dData->set_eventnumber( feventNo);
			dData->set_gain(fgain);
			for (int i = 0; i < s1->TRACE_SIZE; i++) {
				dData->set_trace_value(i, ftrace[i]);
			}
			dData->GetSignalInfo();

			//cout<<"baseline "<<baseline<<endl;
			for(k = kStart; k <= kStop; k += kStep){
				for(m= mStart; m <= mStop; m += mStep){


					if((2*k + m) >= s1->TRACE_SIZE - dData->get_Trigger()) break;


					kIdx = (k - kStart)/kStep;
					mIdx = (m - mStart)/mStep;
					trapAmplitude = filter->trapezoidal_filter_algorithm1(dData, k, m, NULL);
					hRaw[fboardIdx][fchannelID][kIdx][mIdx]->Fill(trapAmplitude);
					trapTree->Fill();
				}//m loop ends here
			}//k loop ends here

		}//event loop ends here

		//close input files
		inFile->cd();
		inFile->Close();
	}//read raw data files codition ends here




	cout<< "Calibration started.."<<endl;

	//----------------------
	// Calibrate
	//-----------------
	for (UShort_t bID = 0; bID < s1->NBOARDS_DSSD; bID++) {
		for(UShort_t chID = 0; chID < NCHANNELS; chID++) {
			for(UShort_t kId = 0; kId < ktimes; kId++){
				for(UShort_t mId = 0; mId < mtimes; mId++){
					if(hRaw[bID][chID][kId][mId]->GetEntries() < 100) continue;
					hRaw[bID][chID][kId][mId]->GetXaxis()->SetRangeUser(sr1,sr2);
					canvas->cd();
					hRaw[bID][chID][kId][mId]->Draw("");
					//Find the peaks
					Int_t nfound = spectrum->Search(hRaw[bID][chID][kId][mId],  2.5, "Markov", 0.2);
					Float_t *xpeaks = spectrum->GetPositionX();
					int idx[nfound];
					TMath::Sort(nfound,spectrum->GetPositionX(),idx,false);
					for (int i = 0; i < nfound; i++) {
						fPositionX[i] = xpeaks[idx[i]];
						Int_t bin = hRaw[bID][chID][kId][mId]->GetXaxis()->FindBin(fPositionX[i]);
						fPositionY[i] = hRaw[bID][chID][kId][mId]->GetBinContent(bin);
						//cout<<"x  "<<fPositionX[i]<<"  y "<<fPositionY[i]<<endl; 
					}

					//fit 1 peak first
					G0->SetRange(fPositionX[0]-10,fPositionX[0]+10);
					G0->SetParameters(fPositionY[0],fPositionX[0],3.);
					hRaw[bID][chID][kId][mId]->Fit(G0,"QLR");
					double sigma = G0->GetParameter(2);
					//Set Range for the fit functions
					G1->SetRange(fPositionX[0]-1*sigma,fPositionX[0]+2*sigma);
					G2->SetRange(fPositionX[2]-1*sigma,fPositionX[2]+2*sigma);
					//set parameters

					//G1->SetParameters(fPositionY[0],fPositionX[0],2.5,1,3, fPositionY[1],fPositionX[1],2.5,1,3);
					//G2->SetParameters(fPositionY[4],fPositionX[4],2.5,1,3, fPositionY[5],fPositionX[5],2.5,1,3);
					G1->SetParameters(fPositionY[0],fPositionX[0],sigma);
					G2->SetParameters(fPositionY[2],fPositionX[2],sigma);
					G1->SetParLimits(2, 1., 20.);
					G2->SetParLimits(2, 1., 20.);
					hRaw[bID][chID][kId][mId]->Fit(G1,"QLR+");
					hRaw[bID][chID][kId][mId]->Fit(G2,"QLR+");

					// Calibration
					fPositionX[0] = G1->GetParameter(1);
					fPositionX[1] = G2->GetParameter(1);
					for (int pt(0); pt<2; pt++) {
						gr_cal->SetPoint(pt,fPositionX[pt], Energy[pt]);
					}
					canvas2->cd();
					gr_cal->Draw();
					gr_cal->Fit(fLinear,"MQ","same");

					//Get Parameters
					gain[bID][chID][kId][mId] =  fLinear->GetParameter(1);
					offset[bID][chID][kId][mId] = fLinear->GetParameter(0);

				}//m loop ends here
			}//k loop ends here
		}//channel loop ends here
	}// board loop ends here

	// fill calibrated spectra
	for(Int_t i = 0; i < trapTree->GetEntries();i++){
		trapTree->GetEntry(i);
		Double_t calibData = gain[fboardIdx][fchannelID][kIdx][mIdx]*trapAmplitude + offset[fboardIdx][fchannelID][kIdx][mIdx];
		hCalib[fboardIdx][fchannelID][kIdx][mIdx]->Fill(calibData);
		// UShort_t sID = (boardIdx *4) + fchannelID;
		// hCalib_tot->Fill(calibData, sID);
	}

	//optimize k and m
	for(UShort_t bID = 0; bID < s1->NBOARDS_DSSD; bID++){
		for(UShort_t chID = 0; chID < NCHANNELS; chID++){
			for(UShort_t kId = 0; kId < ktimes; kId++){
				for(UShort_t mId = 0; mId < mtimes; mId++){
					if( hCalib[bID][chID][kId][mId]->GetEntries() < 100) continue;
					hCalib[bID][chID][kId][mId]->GetXaxis()->SetRangeUser(5700.,6000.);
					canvas->cd();
					hCalib[bID][chID][kId][mId]->Draw();
					double height= hCalib[bID][chID][kId][mId]->GetMaximum();
					G0->SetRange(5800-30.,5800+50.);
					G0->SetParameters(height, 5805,10.);
					hCalib[bID][chID][kId][mId]->Fit(G0,"QLR");
					double sigma = G0->GetParameter(2);

					G1->SetRange(5805-1*sigma,5805+3*sigma);
					//int bin =  hCalib[bID][chID][kIdx][mIdx]->GetXaxis()->FindBin(5770);
					//double height1= hCalib[bID][chID][kIdx][mIdx]->GetBinContent(bin);

					//G1->SetParameters(height1, 5770, 10,1, 3, height, 5805,10.,1.3);
					G1->SetParameters(height, 5805,sigma);
					hCalib[bID][chID][kId][mId]->Fit(G1,"QLR+");
					k = kId * kStep + kStart;
					m = mId * mStep + mStart;
					double FWHM = G1->GetParameter(2) * 2.35;
					trap_parameter r (k, m, FWHM);
					myvector[bID][chID].push_back(r);
					hResolution[bID][chID]->Fill(k, m,FWHM );
				}
			}   
		}
	}
	//----------------
	// Sort in decrasing order----------------
	//------------
	for(UShort_t bID = 0; bID < s1->NBOARDS_DSSD; bID++){
		for(UShort_t chID = 0; chID < NCHANNELS; chID++){
			if(myvector[bID][chID].size() != 0){

				std::sort (myvector[bID][chID].begin(), myvector[bID][chID].end(),compareFWHM);
				cout<<"bID "<< s1->boardList_DSSD[bID]<<"  chId "<<chID<< " FWHM "<<myvector[bID][chID][0].FWHM<<"  k "<<myvector[bID][chID][0].k <<" m "<<myvector[bID][chID][0].m<<endl;

			}
		}
	}
	// sort ends here
	cout<<"Saving plots......"<<endl;
	file_tree->cd();
	trapTree->Write();
	file_tree->Close();

	file_spectra->cd();
	//file_spectra->Write();
	list_res->Write();
	list_calib->Write();
	file_spectra->Close();
	cout<<"output file closed."<<endl;
	//clear memory
	cout<<"clearing memory..................";

	for(UShort_t bID = 0; bID < s1->NBOARDS_DSSD; bID++){
		for(UShort_t chID = 0; chID < NCHANNELS; chID++){
			myvector[bID][chID].clear();
		}
	}
	delete_TObjects();
	// delete inFile;
	delete file_tree;
	delete file_spectra;

	cout<<"done\n"<<endl;
	//close 
	myfile.close();

}

