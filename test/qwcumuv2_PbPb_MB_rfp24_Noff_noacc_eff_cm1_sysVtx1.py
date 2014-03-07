import FWCore.ParameterSet.Config as cms

process = cms.Process("CumuV2")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V27A::All'

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

#fN = cms.untracked.vstring();
#for line in open('flist').read().splitlines():
#	fN.append('file:'+line);
#
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring()
)

#import FWCore.PythonUtilities.LumiList as LumiList
#import FWCore.ParameterSet.Types as CfgTypes
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#JSONfile = 'Cert_210498-211631_HI_PromptReco_Collisions13_JSON_v2.txt'
#myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
#process.source.lumisToProcess.extend(myLumis)


process.cumulant = cms.EDAnalyzer('QWCumuV2'
	, tracks_ = cms.untracked.InputTag('generalTracks')
	, centrality_ = cms.InputTag("centralityBin")
	, chi2_ = cms.untracked.double(40.)
	, vertexSrc_ = cms.untracked.InputTag('offlinePrimaryVertices', "")
	, rfpptmin_ = cms.untracked.double(0.3)
	, rfpptmax_ = cms.untracked.double(3.0)
	, rfpmineta_ = cms.untracked.double(-2.4)
	, rfpmaxeta_ = cms.untracked.double(2.4)
	, bPhiEta_ = cms.untracked.bool(True)
	, bCentNoff_ = cms.untracked.bool(True)
	, poimineta_ = cms.untracked.double(-2.4)
	, poimaxeta_ = cms.untracked.double(2.4)
	, pterrorpt_ = cms.untracked.double(0.1)
	, Noffmin_ = cms.untracked.int32(0)
	, Noffmax_ = cms.untracked.int32(500)
	, fweight_ = cms.untracked.InputTag('TrackCorrections_HIJING_538_OFFICIAL_Mar24.root')
	, bEff_ = cms.untracked.bool(True)
	, cmode_ = cms.untracked.int32(1)
	, minvz_ = cms.untracked.double(3)
	, maxvz_ = cms.untracked.double(15)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cumu.root')
)


import sys

process.p = cms.Path(process.cumulant)

process.schedule = cms.Schedule(
	process.p
)
#process.schedule = cms.Schedule(
#	process.p100,
#	process.p130,
#	process.p160,
#	process.p190,
#	process.p220,
#)
