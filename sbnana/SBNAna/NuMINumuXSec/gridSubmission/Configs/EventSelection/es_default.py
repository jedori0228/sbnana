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
    Name = "RecoFiducial",
    Cut = "kNoCut && cutRFiducial",
    SpillCut = "kNoSpillCut",
    Comment = "1. Reco fiducial volume cut",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "FMScore",
    Cut = "kNoCut && cutRFiducial && cutFMScore",
    SpillCut = "kNoSpillCut",
    Comment = "2. Flash matching score",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "FMTime",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime",
    SpillCut = "kNoSpillCut",
    Comment = "3. Flash matching time",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "TopCRTVeto",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime",
    SpillCut = "spillcutFDTopCRTHitVeto",
    Comment = "4.CRT veto with top CRT",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "SideCRTVeto",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime",
    SpillCut = "spillcutFDSideCRTHitVeto",
    Comment = "4.CRT veto with side CRT",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "CRTVeto",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "5. CRT veto with top+side CRTs",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NuScore",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutNuScore",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "6. Nu score",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "CRLongestTrackDirY",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "6. Longest track y-direction under CR hypo.",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "7. Muon track found",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon_Contained",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutMuonContained",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "7-1. Muon track found, contained",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon_Exiting",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && !cutMuonContained",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "7-2. Muon track found, exiting",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasProton",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutHasProton",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "8. Proton track found",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutHasProton && cutMuonProtonCosineTheta",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "9. Reject back-to-back tracks",
  )
)

## Muon contained and exiting separately

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon_Contained_HasProton",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutMuonContained && cutHasProton",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "8-1. Muon track found, contained + proton",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon_Exiting_HasProton",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && !cutMuonContained && cutHasProton",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "8-2. Muon track found, exiting + proton",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon_Contained_HasProton_MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutMuonContained && cutHasProton && cutMuonProtonCosineTheta",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "9-1. Muon track found, contained + proton + reject back-to-back tracks",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "HasMuon_Exiting_HasProton_MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && !cutMuonContained && cutHasProton && cutMuonProtonCosineTheta",
    SpillCut = "spillcutFDTopCRTHitVeto && spillcutFDSideCRTHitVeto",
    Comment = "9-2. Muon track found, exiting + proton + reject back-to-back tracks",
  )
)



## Below is a no-crt version

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_NuScore",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutNuScore",
    Comment = "4. Nu score",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_CRLongestTrackDirY",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY",
    Comment = "4. Longest track y-direction under CR hypo.",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon",
    Comment = "5. Muon track found",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon_Contained",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutMuonContained",
    Comment = "5-1. Muon track found, contained",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon_Exiting",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && !cutMuonContained",
    Comment = "5-2. Muon track found, exiting",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasProton",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutHasProton",
    Comment = "6. Proton track found",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutHasProton && cutMuonProtonCosineTheta",
    Comment = "7. Proton track found",
  )
)

## Muon contained and exiting separately

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon_Contained_HasProton",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutMuonContained && cutHasProton",
    Comment = "6-1. Muon track found, contained + proton",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon_Exiting_HasProton",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && !cutMuonContained && cutHasProton",
    Comment = "6-2. Muon track found, exiting + proton",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon_Contained_HasProton_MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && cutMuonContained && cutHasProton && cutMuonProtonCosineTheta",
    Comment = "7-1. Muon track found, contained + proton + reject back-to-back tracks",
  )
)

EventSelectionInfos.append(
  SelectionInfo(
    Name = "NoCRT_HasMuon_Exiting_HasProton_MuonProtonCosineTheta",
    Cut = "kNoCut && cutRFiducial && cutFMScore && cutFMTime && cutSliceCRLongestTrackDirY && cutHasMuon && !cutMuonContained && cutHasProton && cutMuonProtonCosineTheta",
    Comment = "7-2. Muon track found, exiting + proton + reject back-to-back tracks",
  )
)
