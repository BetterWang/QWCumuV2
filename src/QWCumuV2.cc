// -*- C++ -*-
//
// Package:    QWCumuV2
// Class:      QWCumuV2
// 
/**\class QWCumuV2 QWCumuV2.cc QWAna/QWCumuV2/src/QWCumuV2.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Quan Wang
//         Created:  Tue Oct 16 16:33:30 EDT 2012
// $Id: QWCumuV2.cc,v 1.1 2013/01/15 15:56:58 qwang Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TComplex.h"


#include "QWAna/QWCumuV2/interface/QWCumuV2.h"


using namespace std;

//#ifdef QW_DEBUG
//
// constructors and destructor
//
QWCumuV2::QWCumuV2(const edm::ParameterSet& iConfig)
	:
		tracks_(iConfig.getUntrackedParameter<edm::InputTag>("tracks_"))
	,	centrality_(iConfig.getParameter<edm::InputTag>("centrality_"))
	,	vertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc_"))
	,	bacc(false)
	,	q2(0, 0, true)
	,	q3(0, 0, true)
{
	//now do what ever initialization is needed
	minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 3.);
	d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
	chi2_ = iConfig.getUntrackedParameter<double>("chi2_", 40);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt_", 0.1);
	rfpmineta_ = iConfig.getUntrackedParameter<double>("rfpmineta_", -2.);
	rfpmaxeta_ = iConfig.getUntrackedParameter<double>("rfpmaxeta_", 2.);
	poimineta_ = iConfig.getUntrackedParameter<double>("poimineta_", rfpmineta_);
	poimaxeta_ = iConfig.getUntrackedParameter<double>("poimaxeta_", rfpmaxeta_);
	fweight_ = iConfig.getUntrackedParameter<edm::InputTag>("fweight_", string("NA"));
	facceptance_ = iConfig.getUntrackedParameter<edm::InputTag>("facceptance_", string("NA"));
	rfpptmin_ = iConfig.getUntrackedParameter<double>("rfpptmin_", 0.1);
	rfpptmax_ = iConfig.getUntrackedParameter<double>("rfpptmax_", 100);
	poiptmin_ = iConfig.getUntrackedParameter<double>("poiptmin_", rfpptmin_);
	poiptmax_ = iConfig.getUntrackedParameter<double>("poiptmax_", rfpptmax_);
	charge_ = iConfig.getUntrackedParameter<int>("charge_", 0);
	bFak = iConfig.getUntrackedParameter<bool>("bFak_", false);
	bEff = iConfig.getUntrackedParameter<bool>("bEff_", false);
	bPhiEta = iConfig.getUntrackedParameter<bool>("bPhiEta_", false);
	bCentNoff = iConfig.getUntrackedParameter<bool>("bCentNoff_", false);
	bSim_ = iConfig.getUntrackedParameter<bool>("bSim_", false);
	Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
	Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 10000);
	effCut_ = iConfig.getUntrackedParameter<double>("effCut_", -1.0);
	cmode_ = iConfig.getUntrackedParameter<int>("cmode_", 1);
	cweight_ = iConfig.getUntrackedParameter<int>("cweight_", 1);
	bGen_ = iConfig.getUntrackedParameter<bool>("bGen_", false);
	nvtx_ = iConfig.getUntrackedParameter<int>("nvtx_", 100);

	if ( cweight_ == 0 ) {
		q2 = correlations::QVector(0, 0, false);
		q3 = correlations::QVector(0, 0, false);
	}

	string streff = fweight_.label();
	if ( streff == string("NA") ) {
		bFak = false;
		bEff = false;
		fEffFak = 0;
	} else {
		fEffFak = new TFile(streff.c_str());
		if ( !fEffFak->IsOpen() ) {
			bFak = false;
			bEff = false;
		} else {
			cout << "!!! Using particle weight " << streff << endl;
			if ( bFak ) cout << "!!! Apply Fak correction" << endl;
			if ( bEff ) cout << "!!! Apply Eff correction" << endl;
			for ( int i = 0; i < 20; i++ ) {
				if ( streff == string("TrackCorrections_HIJING_538_OFFICIAL_Mar24.root") || streff == string("trkEff_pp_all_42X_origin.root") ) {
					hEff_cbin[i] = (TH2D*) fEffFak->Get("rTotalEff3D");
					hFak_cbin[i] = (TH2D*) fEffFak->Get(Form("rFak_cbin%i", i));
				}
				if ( streff == string("trkEffNew2012_HI_hiGoodTightMerged_xsec_smoothv5true.root") ) {
					hEff_cbin[i] = (TH2D*) fEffFak->Get("Tot_4");
					hFak_cbin[i] = (TH2D*) fEffFak->Get("Fak_4");
				}
			}
			cout << "!!! eff histo done" << endl;
		}
	}
	string stracc = facceptance_.label();
	if ( stracc == string("NA") ) {
		cout << "!!! acc NA" << endl;
		bacc = false;
		facc = 0;
	} else {
		facc = new TFile(stracc.c_str());
		if ( !facc->IsOpen() ) {
			bacc = false;
		} else {
			cout << "!!! Using acceptance weight " << stracc << endl;
			bacc = true;
			for ( int cent = 0; cent < nCentBins; cent++ ) {
				for ( int ipt = 0; ipt < nPtBins; ipt++ ) {
					hacc[cent][ipt][0] = (TH2D*) facc->Get(Form("hPhiEta_%i_%i_0", cent, ipt));
					hacc[cent][ipt][1] = (TH2D*) facc->Get(Form("hPhiEta_%i_%i_1", cent, ipt));
				}
			}
		}
	}


	//
	//cout << __LINE__ << "\t" << tracks_.label().c_str() << "\t|" << tracks_.instance() << "\t|" << tracks_.process() << endl;
	//
//	cout << "!! before t" << endl;
	t = new QWEvent;
	memset(t, 0, sizeof(QWEvent));
//	cout << "!! after t" << endl;
	//
	for ( int cent = 0; cent < nCentBins; cent++ ) {
		hPt[cent]       = fs->make<TH1D>(Form("hPt_%i", cent), "", 20000, 0, 100);
		if ( bPhiEta ) {
			for ( int i = 0; i < nPtBins; i++ ) {
				hPhiEta[cent][i][0] = fs->make<TH2D>(Form("hPhiEta_%i_%i_0", cent, i), "", 512, -Pi, Pi, 480, -2.4, 2.4);
				hPhiEta[cent][i][1] = fs->make<TH2D>(Form("hPhiEta_%i_%i_1", cent, i), "", 512, -Pi, Pi, 480, -2.4, 2.4);
			}
		}
	}
	for ( int cent = 0; cent < 20; cent++ ) {
		hdNdPtdEta[cent] = fs->make<TH2D>(Form("hdNdPtdEta_%i", cent), Form("hdNdPtdEta_%i", cent), nEtaBins, etabins, 38, fakpt );
		hdNdPtdEtaPt[cent] = fs->make<TH2D>(Form("hdNdPtdEtaPt_%i", cent), Form("hdNdPtdEta_%i", cent), nEtaBins, etabins, 38, fakpt );
	}
//	cout << "!! before nt" << endl;

	initQ();
//	cout << cq2->name() << endl;
	ntResult = new TNtupleD("ntResult",cq2->name(),"Noff:Mult:Cent:C22:C24:C26:C28:iC22:iC24:iC26:iC28:wC22:wC24:wC26:wC28:C32:C34:C36:C38:iC32:iC34:iC36:iC38:wC32:wC34:wC36:wC38");
	ntResult->SetAutoFlush(-3000000);
	ntResult->SetAutoSave(-30000000);

//	cout << "!! initQ" << endl;
}


QWCumuV2::~QWCumuV2()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

int
QWCumuV2::getNoffCent(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& Noff)
{
	// very hard coded Noff track centrality cut
	using namespace edm;
	using namespace reco;
//	int Noff = 0;

	Handle<VertexCollection> vertexCollection;
	iEvent.getByLabel(vertexSrc_, vertexCollection);
	const VertexCollection * recoVertices = vertexCollection.product();

	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();


	Handle<TrackCollection> tracks;
	iEvent.getByLabel(tracks_,tracks);
	for(TrackCollection::const_iterator itTrack = tracks->begin();
		itTrack != tracks->end();                      
		++itTrack) {

		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
		if ( itTrack->charge() == 0 ) continue;
		if ( itTrack->pt() < 0.4 ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( fabs( dz/dzerror ) > 3. ) continue;
		if ( fabs( d0/derror ) > 3. ) continue;
		if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
//		bool b_pix = itTrack->numberOfValidHits() < 7;
//		if ( b_pix ) {
//			if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;
//			if ( itTrack->normalizedChi2() > chi2_ ) continue;
//		} else {
//			// full track
//			if ( fabs( dz/dzerror ) > 3. ) continue;
//			if ( fabs( d0/derror ) > 3. ) continue;
//			if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
//			if ( itTrack->numberOfValidHits() < 12 ) continue;
//		}

		Noff++;
	}

	int cent = nCentNoff-1;
	while ( CentNoffCut[cent] <= Noff ) cent--;
	return cent;
}


void
QWCumuV2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if ( bGen_ ) analyzeGen(iEvent, iSetup);
	else analyzeData(iEvent, iSetup);
	if ( t->Mult == 0 ) return;
	for ( int i = 0; i < t->Mult; i++ ) {
		q2.fill(t->Phi[i], t->weight[i]);
		q3.fill(t->Phi[i], t->weight[i]);
	}
	correlations::Result r22 = cq2->calculate(2, hc2);
	correlations::Result r24 = cq2->calculate(4, hc2);
	correlations::Result r26 = cq2->calculate(6, hc2);
	correlations::Result r28 = cq2->calculate(8, hc2);

	correlations::Result r32 = cq3->calculate(2, hc3);
	correlations::Result r34 = cq3->calculate(4, hc3);
	correlations::Result r36 = cq3->calculate(6, hc3);
	correlations::Result r38 = cq3->calculate(8, hc3);

	double C22,C24,C26,C28,iC22,iC24,iC26,iC28,wC22,wC24,wC26,wC28,C32,C34,C36,C38,iC32,iC34,iC36,iC38,wC32,wC34,wC36,wC38;

	C22 = r22._sum.real();
	C24 = r24._sum.real();
	C26 = r26._sum.real();
	C28 = r28._sum.real();

	iC22 = r22._sum.imag();
	iC24 = r24._sum.imag();
	iC26 = r26._sum.imag();
	iC28 = r28._sum.imag();

	wC22 = r22._weights;
	wC24 = r24._weights;
	wC26 = r26._weights;
	wC28 = r28._weights;

	C32 = r32._sum.real();
	C34 = r34._sum.real();
	C36 = r36._sum.real();
	C38 = r38._sum.real();

	iC32 = r32._sum.imag();
	iC34 = r34._sum.imag();
	iC36 = r36._sum.imag();
	iC38 = r38._sum.imag();

	wC32 = r32._weights;
	wC34 = r34._weights;
	wC36 = r36._weights;
	wC38 = r38._weights;

	double res[27] = { double(t->Noff), double(t->Mult), double(t->Cent),
		C22,C24,C26,C28,iC22,iC24,iC26,iC28,wC22,wC24,wC26,wC28,C32,C34,C36,C38,iC32,iC34,iC36,iC38,wC32,wC34,wC36,wC38};

	fs->cd();
	ntResult->Fill(res);

	doneQ();

}


void
QWCumuV2::analyzeGen(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	t->Mult = 0;
	// track
	Handle< std::vector<GenParticle> > tracks;
	iEvent.getByLabel(tracks_,tracks);
	t->Noff = 0;
	for ( std::vector<GenParticle>::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack
			)
	{
		if ( itTrack->status()!=1 ) continue;
		if ( itTrack->charge() == 0 ) continue;
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( itTrack->pt() < 0.4 ) continue;
		t->Noff++;
	}
	t->Cent = 0;
	for ( std::vector<GenParticle>::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack
			)
	{
		if ( itTrack->status()!=1 ) continue;
		if ( itTrack->charge() == 0 ) continue;
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		t->Pt[t->Mult] = itTrack->pt();
		if ( (t->Pt[t->Mult] < rfpptmin_) || (t->Pt[t->Mult] > rfpptmax_) || itTrack->eta() < rfpmineta_ || itTrack->eta() > rfpmaxeta_ ) continue;

		t->weight[t->Mult] = 1.;
		t->Phi[t->Mult] = itTrack->phi();

		t->Mult++;
	}
}


// ------------ method called for each event  ------------
	void
QWCumuV2::analyzeData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	t->Mult = 0;
	// vertex
	Handle<VertexCollection> vertexCollection;
	iEvent.getByLabel(vertexSrc_, vertexCollection);
	const VertexCollection * recoVertices = vertexCollection.product();
	if ( recoVertices->size() > nvtx_ ) return;

	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();

//	for ( unsigned int i = 0; i < recoVertices->size(); i++ ) {
//		size_t daughter = (*recoVertices)[i].tracksSize();
//		cout << "i = " << i << "\tnTracks = " << daughter <<"\t vz = " << (*recoVertices)[i].position().z() << endl;
//		//cout << "i = " << i << "\ttrkSize = " << "\t vz = " << (*recoVertices)[i].position().z() << endl;
//	}
	double vz = (*recoVertices)[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		return;
	}
	
	// centrality
	int bin = 0;
	int cbin = 0;
	t->Noff = 0;

	if ( bCentNoff ) {
		cbin = getNoffCent( iEvent, iSetup, t->Noff);
		if ( (t->Noff < Noffmin_) or (t->Noff >= Noffmax_) ) {
			return;
		}
	} else {
		edm::Handle<int> ch;
		iEvent.getByLabel(centrality_,ch);
		bin = *(ch.product());
		while ( centbins[cbin+1] < bin*2.5+0.1 ) cbin++;
	}
	bin = cbin;

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByLabel(tracks_,tracks);
	t->Cent = bin;
	t->vz = vz;
	//cout << __LINE__ << "\t" << bin << endl;

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();                      
			++itTrack) {
//		cout << "!!! " << __LINE__ << endl;
		if ( itTrack->charge() == 0 ) continue;
		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;

//		cout << "!!! " << __LINE__ << endl;
		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);

//		cout << "!!! " << __LINE__ << endl;
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;
		if ( fabs( d0/derror ) > d0d0error_ ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;

//		cout << "!!! " << __LINE__ << endl;
		t->Charge[t->Mult] = itTrack->charge();
		if ( (charge_ == 1) && (t->Charge[t->Mult]<0) ) continue;
		if ( (charge_ == -1) && (t->Charge[t->Mult]>0) ) continue;

//		cout << "!!! " << __LINE__ << endl;
		t->Pt[t->Mult] = itTrack->pt();
		if ( t->Pt[t->Mult] >= ptbins[nPtBins] || t->Pt[t->Mult] <= ptbins[0] ) continue;
		t->Eta[t->Mult] = itTrack->eta();

		if (effCut_>0.)  {
			double eff = hEff_cbin[bin]->GetBinContent( hEff_cbin[bin]->FindBin(t->Eta[t->Mult], t->Pt[t->Mult] ) ) ;
			if ( eff > effCut_ ) {
				if ( gRandom->Rndm() < (eff-effCut_)/eff ) continue;
			}
		}

		if ( bEff ) {
			t->rEff[t->Mult] = hEff_cbin[bin]->GetBinContent( hEff_cbin[bin]->FindBin(t->Eta[t->Mult], t->Pt[t->Mult] ) );
		} else {
			t->rEff[t->Mult] = 1.;
		}
		if ( bFak ) {
			t->rFak[t->Mult] = hFak_cbin[bin]->GetBinContent( hFak_cbin[bin]->FindBin(t->Eta[t->Mult], t->Pt[t->Mult] ) );
		} else {
			t->rFak[t->Mult] = 0.;
		}
		if ( t->rEff[t->Mult] <= 0.1 or TMath::IsNaN(t->rEff[t->Mult]) ) continue;
		double weight = (1.-t->rFak[t->Mult])/t->rEff[t->Mult];

		double phi = itTrack->phi();

		double wacc = 1.;
		int ipt=0;
		while ( t->Pt[t->Mult] > ptbins[ipt+1] ) ipt++;
		if ( bacc ) {
			wacc = 1./hacc[t->Cent][ipt][t->Charge[t->Mult]>0]->GetBinContent(hacc[t->Cent][ipt][t->Charge[t->Mult]>0]->FindBin(phi, t->Eta[t->Mult]));
		}
		if ( bPhiEta ) hPhiEta[t->Cent][ipt][t->Charge[t->Mult]>0]->Fill(phi, t->Eta[t->Mult], wacc);

		weight *= wacc;

		if ( (t->Pt[t->Mult] < rfpptmin_) || (t->Pt[t->Mult] > rfpptmax_) || itTrack->eta() < rfpmineta_ || itTrack->eta() > rfpmaxeta_ ) continue;

		t->weight[t->Mult] = weight;

		hdNdPtdEta[bin]->Fill(t->Eta[t->Mult], t->Pt[t->Mult]);
		hdNdPtdEtaPt[bin]->Fill(t->Eta[t->Mult], t->Pt[t->Mult], t->Pt[t->Mult]);

		t->Phi[t->Mult] = phi;
		hPt[t->Cent]->Fill(t->Pt[t->Mult]);

		t->Mult++;
	}
	if ( bSim_ ) Sim();
}


void
QWCumuV2::initQ()
{
	hc2 = correlations::HarmonicVector(8);
	hc2[0] = 2;
	hc2[1] = -2;
	hc2[2] = 2;
	hc2[3] = -2;
	hc2[4] = 2;
	hc2[5] = -2;
	hc2[6] = 2;
	hc2[7] = -2;

	hc3 = correlations::HarmonicVector(8);
	hc3[0] = 3;
	hc3[1] = -3;
	hc3[2] = 3;
	hc3[3] = -3;
	hc3[4] = 3;
	hc3[5] = -3;
	hc3[6] = 3;
	hc3[7] = -3;

	q2.resize(hc2);
	q3.resize(hc3);
	switch ( cmode_ ) {
		case 1:
			cq2 = new correlations::closed::FromQVector(q2);
			cq3 = new correlations::closed::FromQVector(q3);
			break;
		case 2:
			cq2 = new correlations::recurrence::FromQVector(q2);
			cq3 = new correlations::recurrence::FromQVector(q3);
			break;
		case 3:
			cq2 = new correlations::recursive::FromQVector(q2);
			cq3 = new correlations::recursive::FromQVector(q3);
			break;
	}
}

void
QWCumuV2::doneQ()
{
	q2.reset();
	q3.reset();
}

void
QWCumuV2::Sim()
{
	// QC{2} = -0.071428
	// QC{4} = -0.007142
	// QC{6} = 0.0357142
	// QC{8} = 0.0571429
	t->Cent = 100;
	t->Mult = 8;
	for ( int i = 0; i < 8; i++ ) {
		t->weight[i] = 1.;
	}
	t->Phi[0] = -Pi;
	t->Phi[1] = -Pi*5/6.;
	t->Phi[2] = -Pi/2.;
	t->Phi[3] = -Pi/6.;
	t->Phi[4] = 0;
	t->Phi[5] = Pi/6.;
	t->Phi[6] = Pi/2.;
	t->Phi[7] = Pi*5/6.;
}

// ------------ method called once each job just before starting event loop  ------------
	void 
QWCumuV2::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
QWCumuV2::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
QWCumuV2::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
QWCumuV2::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
QWCumuV2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
QWCumuV2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QWCumuV2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//////////////////////////////////////////


//define this as a plug-in
DEFINE_FWK_MODULE(QWCumuV2);
