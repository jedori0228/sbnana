from SelectionInfo import *

SampleSelectionInfos = []

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC",
    Cut = "cutIsNuMuCC",
    Comment = "NuMu-CC (+FV), 1mu1p"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC1mu1p",
    Cut = "cutIsNuMuCC && cutTFiducial && (cutHasTruthMuon && cutHasTruthProton)",
    Comment = "NuMu-CC (+FV), 1mu1p"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC1mu1pTruthCut",
    Cut = "cutIsNuMuCC && cutTFiducial && (cutHasTruthMuon && cutHasTruthProton) && (cutTruthMuonTCut && cutTruthProtonTCut)",
    Comment = "NuMu-CC (+FV), 1mu1p, truthcut ( == Signal)"
  )
)

