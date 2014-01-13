#define CORRELATIONS_CLOSED_ENABLE_U8 1
#include <QWAna/QWCumuV2/interface/correlations/Types.hh>
#include <QWAna/QWCumuV2/interface/correlations/Result.hh>
#include <QWAna/QWCumuV2/interface/correlations/QVector.hh>
#include <QWAna/QWCumuV2/interface/correlations/recursive/FromQVector.hh>
#include <QWAna/QWCumuV2/interface/correlations/recurrence/FromQVector.hh>
#include <QWAna/QWCumuV2/interface/correlations/closed/FromQVector.hh>
#include <TComplex.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtupleD.h>
#include <TRandom3.h>
#include <TFile.h>
#include "QWConstV2.h"
//
// constants, enums and typedefs
//

//#define QW_DEBUG 1
//#define QW_PEREVENT 1

#define PR(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << " = " << (x) << endl;
#define PRD(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << endl;
//
// class declaration
//

const int NMAX_TRK = 5000;
typedef struct QWEvent_ {
	int     Cent;
	int     Mult;
	double  vz;
	double  Pt[NMAX_TRK];
	double  Eta[NMAX_TRK];
	double  Phi[NMAX_TRK];
	int     Charge[NMAX_TRK];
	double	rEff[NMAX_TRK];
	double	rFak[NMAX_TRK];
	double	weight[NMAX_TRK];
} QWEvent;

///////////////// Class ////////////////////////////

class QWCumuV2 : public edm::EDAnalyzer {
	public:
		explicit QWCumuV2(const edm::ParameterSet&);
		~QWCumuV2();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	/////////////////////////////////////////////
		int getNoffCent(const edm::Event&, const edm::EventSetup&, int& Noff);
		//TRandom3 * gRandom;
		// ----------member data ---------------------------
		edm::InputTag tracks_; //used to select what tracks to read from configuration file
		edm::InputTag centrality_;	// centrality
		edm::InputTag vertexSrc_;
		edm::InputTag correctHist_;

		edm::InputTag fweight_;
		edm::InputTag facceptance_;
	/////////////////////////////////////////////
		double 	minvz_, maxvz_;
		double 	dzdzerror_;
		//double 	d0d0error_;
		double 	chi2_;
		double 	pterrorpt_;
		double 	rfpmineta_, rfpmaxeta_;
		double 	poimineta_, poimaxeta_;
		double 	rfpptmin_, rfpptmax_;
		double 	poiptmin_, poiptmax_;
		int 	charge_;

		bool	bFak;
		bool	bEff;
		bool	bacc;
		bool	bPhiEta;
		bool	bCentNoff;
		bool	bSim_;
		int 	Noffmin_;
		int 	Noffmax_;
		int	cmode_;
		int	cweight_;

		double	effCut_;
		QWEvent * t;
		TFile	* fEffFak;
		TFile	* facc;
	/////////////////////////////////////////////
		TH1D * hPt[nCentBins];
		TH2D * hPhiEta[nCentBins][nPtBins][2];

		TH2D	* hdNdPtdEta[nCentBins];
		TH2D	* hdNdPtdEtaPt[nCentBins];

		TH2D * hEff_cbin[nCentBins];
		TH2D * hFak_cbin[nCentBins];

		TH2D * hacc[nCentBins][nPtBins][2];

		TNtupleD * ntResult;

		correlations::HarmonicVector 	hc2;
		correlations::HarmonicVector 	hc3;

		correlations::QVector	q2;
		correlations::QVector	q3;
		correlations::FromQVector 	*cq2;
		correlations::FromQVector 	*cq3;

		void initQ();
		void doneQ();
		void Sim();
};



