import FWCore.ParameterSet.Config as cms

HbbProducer = cms.EDProducer('HbbProducer',
                             rhoSource=cms.InputTag('fixedGridRhoFastjetAll'),

                             packedCandidateSource=cms.InputTag('chs'),

                             AK4Source =cms.InputTag('patJetsAK4PFCHS'),
                             AK8Source =cms.InputTag('patJetsAK8PFCHS'),
                             AK10Source=cms.InputTag('patJetsAK10PFCHS'),
                             AK12Source=cms.InputTag('patJetsAK12PFCHS'),
                             AK15Source=cms.InputTag('patJetsAK15PFCHS'),

                             AK8PackedSource =cms.InputTag('patJetsAK8PFCHSFilteredPacked'),
                             AK10PackedSource=cms.InputTag('patJetsAK10PFCHSFilteredPacked'),
                             AK12PackedSource=cms.InputTag('patJetsAK12PFCHSFilteredPacked'),
                             AK15PackedSource=cms.InputTag('patJetsAK15PFCHSFilteredPacked'),
                             
                             muonSource=cms.InputTag('selectedMuons'),
)
