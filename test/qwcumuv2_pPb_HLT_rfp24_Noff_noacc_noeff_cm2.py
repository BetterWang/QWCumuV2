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
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V43::All'

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

fN = cms.untracked.vstring();
for line in open('flist').read().splitlines():
	fN.append('file:'+line);

process.source = cms.Source("PoolSource",
	fileNames = fN
)

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
JSONfile = 'Cert_210498-211631_HI_PromptReco_Collisions13_JSON_v2.txt'
myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
process.source.lumisToProcess.extend(myLumis)


import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltHM100 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM100.HLTPaths = [
	"HLT_PAPixelTracks_Multiplicity100_v*",
#	"HLT_PAPixelTracks_Multiplicity130_v*",
#	"HLT_PAPixelTracks_Multiplicity160_v*",
#	"HLT_PAPixelTracks_Multiplicity190_v*",
#	"HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM100.andOr = cms.bool(True)
process.hltHM100.throw = cms.bool(False)

process.hltHM130 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM130.HLTPaths = [
	"HLT_PAPixelTracks_Multiplicity100_v*",
	"HLT_PAPixelTracks_Multiplicity130_v*",
#	"HLT_PAPixelTracks_Multiplicity160_v*",
##	"HLT_PAPixelTracks_Multiplicity190_v*",
#	"HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM130.andOr = cms.bool(True)
process.hltHM130.throw = cms.bool(False)


process.hltHM160 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM160.HLTPaths = [
	"HLT_PAPixelTracks_Multiplicity100_v*",
	"HLT_PAPixelTracks_Multiplicity130_v*",
	"HLT_PAPixelTracks_Multiplicity160_v*",
#	"HLT_PAPixelTracks_Multiplicity190_v*",
#	"HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM160.andOr = cms.bool(True)
process.hltHM160.throw = cms.bool(False)




process.hltHM190 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM190.HLTPaths = [
	"HLT_PAPixelTracks_Multiplicity100_v*",
	"HLT_PAPixelTracks_Multiplicity130_v*",
	"HLT_PAPixelTracks_Multiplicity160_v*",
	"HLT_PAPixelTracks_Multiplicity190_v*",
#	"HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM190.andOr = cms.bool(True)
process.hltHM190.throw = cms.bool(False)



process.hltHM220 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM220.HLTPaths = [
	"HLT_PAPixelTracks_Multiplicity100_v*",
	"HLT_PAPixelTracks_Multiplicity130_v*",
	"HLT_PAPixelTracks_Multiplicity160_v*",
	"HLT_PAPixelTracks_Multiplicity190_v*",
	"HLT_PAPixelTracks_Multiplicity220_v*"
]
process.hltHM220.andOr = cms.bool(True)
process.hltHM220.throw = cms.bool(False)




process.cumulant100 = cms.EDAnalyzer('QWCumuV2'
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
	, Noffmin_ = cms.untracked.int32(120)
	, Noffmax_ = cms.untracked.int32(150)
#	, fweight_ = cms.untracked.InputTag('TrackCorrections_HIJING_538_OFFICIAL_Mar24.root')
	, bEff_ = cms.untracked.bool(True)
	, cmode_ = cms.untracked.int32(2)
	, cweight_ = cms.untracked.int32(0)
)

process.cumulant130 = process.cumulant100.clone()
process.cumulant160 = process.cumulant100.clone()
process.cumulant190 = process.cumulant100.clone()
process.cumulant220 = process.cumulant100.clone()

process.cumulant130.Noffmin_ = cms.untracked.int32(150)
process.cumulant130.Noffmax_ = cms.untracked.int32(185)
process.cumulant160.Noffmin_ = cms.untracked.int32(185)
process.cumulant160.Noffmax_ = cms.untracked.int32(220)
process.cumulant190.Noffmin_ = cms.untracked.int32(220)
process.cumulant190.Noffmax_ = cms.untracked.int32(260)
process.cumulant220.Noffmin_ = cms.untracked.int32(260)
process.cumulant220.Noffmax_ = cms.untracked.int32(1000)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cumu.root')
)


#process.p = cms.Path(process.PAcollisionEventSelection*process.pACentrality*process.cumulant)
#process.p100 = cms.Path(process.hltHM100*process.cumulant100)
#process.p130 = cms.Path(process.hltHM130*process.cumulant130)
#process.p160 = cms.Path(process.hltHM160*process.cumulant160)
#process.p190 = cms.Path(process.hltHM190*process.cumulant190)
#process.p220 = cms.Path(process.hltHM220*process.cumulant220)

import sys

switch = sys.argv[2];
if switch == '0':
	process.cumulant = process.cumulant100.clone()
	process.p = cms.Path(process.hltHM100*process.cumulant)
	process.TFileService.fileName = cms.string('cumu_100.root')
if switch == '1':
	process.cumulant = process.cumulant130.clone()
	process.p = cms.Path(process.hltHM130*process.cumulant)
	process.TFileService.fileName = cms.string('cumu_130.root')
if switch == '2':
	process.cumulant = process.cumulant160.clone()
	process.p = cms.Path(process.hltHM160*process.cumulant)
	process.TFileService.fileName = cms.string('cumu_160.root')
if switch == '3':
	process.cumulant = process.cumulant190.clone()
	process.p = cms.Path(process.hltHM190*process.cumulant)
	process.TFileService.fileName = cms.string('cumu_190.root')
if switch == '4':
	process.cumulant = process.cumulant220.clone()
	process.p = cms.Path(process.hltHM220*process.cumulant)
	process.TFileService.fileName = cms.string('cumu_220.root')

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
