from SelectionInfo import *

EventSelectionInfos = []

#------------------------------------------------------------------
# Truth muon is found from the primary particles

# Truth muon is contained
# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContainedHasRecoTrack",
    Cut = "cutTruthMuonContained && cutHasTruthMuonHasRecoTrack",
    Comment = "Has truth muon contained, recoed as a track",
  )
)

#------------------------------------------------------------------
# Truth proton is found from the primary particles

# Truth proton is contained
# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthProtonContainedHasRecoTrack",
    Cut = "cutTruthProtonContained && cutHasTruthProtonHasRecoTrack",
    Comment = "Has truth proton contained, recoed as a track",
  )
)

#------------------------------------------------------------------
# Truth charged pion is found from the primary particles

# Truth charged pion is contained
# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionContainedHasRecoTrack",
    Cut = "cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack",
    Comment = "Has truth charged pion contained, recoed as a track",
  )
)










