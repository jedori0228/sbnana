from SelectionInfo import *

SampleSelectionInfos = []

SampleSelectionInfos.append( 
  SelectionInfo(
    Name = "AllSamples",
    Cut = "kNoCut",
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC",
    Cut = "cutIsNuMuCC",
    Comment = "NuMu-CC (+FV), 1mu1p"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCCContained",
    Cut = "cutIsNuMuCC && cutTruthMuonContained",
    Comment = "NuMu-CC (+FV), 1mu1p, truth contained"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCCFailedTFV",
    Cut = "cutIsNuMuCC && !cutTFiducial",
    Comment = "All NuMu-CC but fails FV"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC1muTruthCut",
    Cut = "cutIsNuMuCC && cutTFiducial && cutHasTruthMuon && cutTruthMuonTCut",
    Comment = "NuMu-CC (+FV), 1mu, truthcut"
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

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC1mu1pTruthCutContained",
    Cut = "cutIsNuMuCC && cutTFiducial && (cutHasTruthMuon && cutHasTruthProton) && (cutTruthMuonTCut && cutTruthProtonTCut) && cutTruthMuonContained",
    Comment = "NuMu-CC (+FV), 1mu1p, truthcut ( == Signal), truth contained"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCCSignalDef",
    Cut = "cutIsNuMuCC && cutIsNuMuCCSignalDef",
    Comment = "NuMu-CC with 1mu1p (+truth cut) : NuMuCCSignalDef"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCC1mu1pFailedTruthCut",
    Cut = "cutIsNuMuCC && cutTFiducial && (cutHasTruthMuon && cutHasTruthProton) && !(cutTruthMuonTCut && cutTruthProtonTCut)",
    Comment = "NuMu-CC (+FV), 1mu1p, failed truthcut"
  )
)

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCCNot1mu1p",
    Cut = "cutIsNuMuCC && cutTFiducial && !(cutHasTruthMuon && cutHasTruthProton)",
    Comment = ""
  )
)

for name in ["NuMuNC", "NuE", "NuMuCCQE", "NuMuCCRes", "NuMuCCMEC", "NuMuCCDIS"]:

  SampleSelectionInfos.append(
    SelectionInfo(
      Name = name,
      Cut = "cutIs%s"%(name),
      Comment = name
    )
  )

SampleSelectionInfos.append(
  SelectionInfo(
    Name = "NuMuCCOtherInt",
    Cut = "cutIsNuMuCC && !(cutIsQE||cutIsRes||cutIsMEC||cutIsRes)",
    Comment = ""
  )
)
SampleSelectionInfos.append(
  SelectionInfo(
    Name = "Cosmic",
    Cut = "kIsCosmic",
    Comment = "Cosmic from the overliad sample"
  )
)

