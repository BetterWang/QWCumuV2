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

process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.HiGenCommon.VtxSmearedRealisticPPbBoost8TeVCollision_cff')
process.load('GeneratorInterface.HiGenCommon.AfterBurnerGenerator_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("EmptySource")


process.generator = cms.EDFilter("HijingGeneratorFilter",
                                     frame = cms.string('CMS     '),
                                     targ = cms.string('P       '),
                                     izp = cms.int32(82),
                                     bMin = cms.double(0),
                                     izt = cms.int32(1),
                                     proj = cms.string('A       '),
                                     comEnergy = cms.double(5020.0),
                                     iat = cms.int32(1),
                                     bMax = cms.double(10),
                                     iap = cms.int32(208),
                                     rotateEventPlane = cms.bool(True)
                                 )


process.cumulant = cms.EDAnalyzer('QWCumuV2'
	, tracks_ = cms.untracked.InputTag('genParticles')
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
#	, fweight_ = cms.untracked.InputTag('TrackCorrections_HIJING_538_OFFICIAL_Mar24.root')
	, bEff_ = cms.untracked.bool(False)
	, cmode_ = cms.untracked.int32(1)
	, bGen_ = cms.untracked.bool(True)
)

process.AftBurner.modv1 = cms.InputTag("0.0")
process.AftBurner.modv2 = cms.InputTag("0.05")
process.AftBurner.fluct_v1 = cms.double(0.0)
process.AftBurner.fluct_v1 = cms.double(0.0)
process.AftBurner.modmethod = cms.int32(1)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cumuHijing.root')
)

process.pgen_hijing = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")+process.VertexSmearing+process.AfterBurner+process.GeneInfo)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()


#process.generation_step = cms.Path(process.pgen)
process.generation_step = cms.Path(process.pgen_hijing+process.cumulant)

#process.p = cms.Path(process.generation_step, process.cumulant)
#	process.generation_step + process.cumulant

process.schedule = cms.Schedule(
	process.generation_step
)
for path in process.paths:
        getattr(process,path)._seq = process.generator * getattr(process,path)._seq
