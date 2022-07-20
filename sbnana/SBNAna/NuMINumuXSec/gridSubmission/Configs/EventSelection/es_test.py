from SelectionInfo import *

EventSelectionInfos = []

EventSelectionInfos.append( 
  SelectionInfo(
    Name = "NoCut",
    Cut = "kNoCut",
    SpillCut = "kNoSpillCut",
    Comment = "0. No cut",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutHasProton && cutMuonProtonCosineTheta",
    Comment = "7. Proton track found",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_FailMuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutHasProton && !cutMuonProtonCosineTheta",
    Comment = "7. Proton track found",
  )
)
