#include <TMath.h>

const double Pi = TMath::Pi();
const double Pi2 = 2*Pi;

//PbPb centrality binning
const double centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,90,100};

const double ptbins[] = {
	        0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,
		        10.0, 12.0, 16.0, 20.0, 30.0, 40.0, 100.0}; // pPb pt binning nPtBins = 17;
const Int_t nPtBins = sizeof(ptbins)/sizeof(double)-1;

const double etabins[] = {
	        -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4
};
const Int_t nEtaBins = sizeof(etabins)/sizeof(double)-1;

const double fakpt[] = {
	        0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 4, 6, 8, 12
}; // 39 nbins=38

const Int_t nCentBins = 20;

const Int_t CentNoffCut[] = {100000, 350, 320, 300, 260, 240, 220, 185, 150, 120, 100, 80, 60, 50, 40, 30, 20, 10, 0};
//const Int_t CentNoffCut[] = {100000,         300, 260,      220, 185, 150, 120, 110, 90, 35, 0};
const Int_t nCentNoff = sizeof(CentNoffCut)/sizeof(Int_t);
