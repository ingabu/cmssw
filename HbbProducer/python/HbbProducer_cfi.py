import FWCore.ParameterSet.Config as cms

HbbProducer = cms.EDProducer('HbbProducer',
                             rhoSource=cms.InputTag('fixedGridRhoFastjetAll'),

                             packedCandidateSource=cms.InputTag('chs'),

                             AK4Source =cms.InputTag('patJetsAK4PFCHS'),
                             AK8Source =cms.InputTag('patJetsAK8PFCHS'),
                             AK10Source=cms.InputTag('patJetsAK10PFCHS'),
                             AK12Source=cms.InputTag('patJetsAK12PFCHS'),
                             AK15Source=cms.InputTag('patJetsAK15PFCHS'),

                             muonSource=cms.InputTag('selectedMuons'),
)
