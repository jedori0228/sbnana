from SelectionInfo import *

EventSelectionInfos = []

#------------------------------------------------------------------
# No cut

EventSelectionInfos.append( 
  SelectionInfo(
    Name = "NoCut",
    Cut = "kNoCut",
    SpillCut = "kNoSpillCut",
    Comment = "No cut",
  )
)

#------------------------------------------------------------------
# Truth muon is found from the primary particles

EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuon",
    Cut = "cutHasTruthMuon",
    Comment = "Has truth muon",
  )
)

# Truth muon is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContained",
    Cut = "cutTruthMuonContained",
    Comment = "Has truth muon contained",
  )
)

# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContainedHasRecoTrack",
    Cut = "cutTruthMuonContained && cutHasTruthMuonHasRecoTrack",
    Comment = "Has truth muon contained, recoed as a track",
  )
)

# ---> that reco track is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContainedHasRecoTrackContained",
    Cut = "cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && cutTruthMuonMatchedTrackContained",
    Comment = "Has truth muon contained, recoed as a contained track",
  )
)

# -----> Applying nominal reco cut
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContainedHasRecoTrackContained_NominalContainedMuon",
    Cut = "cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && cutTruthMuonMatchedTrackContained && cutNominal_ContainedMuon",
    Comment = "Has truth muon contained, recoed as a contained track",
  )
)

# ---> that reco track is exiting
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContainedHasRecoTrackExiting",
    Cut = "cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && !cutTruthMuonMatchedTrackContained",
    Comment = "Has truth muon contained, recoed as a exiting track",
  )
)

# -----> Applying nominal reco cut
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonContainedHasRecoTrackExiting_NominalExitingMuon",
    Cut = "cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && !cutTruthMuonMatchedTrackContained && cutNominal_ExitingMuon",
    Comment = "Has truth muon contained, recoed as a exiting track",
  )
)

# Truth muon is exiting
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonExiting",
    Cut = "!cutTruthMuonContained",
    Comment = "Has truth muon exiting",
  )
)

# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonExitingHasRecoTrack",
    Cut = "!cutTruthMuonContained && cutHasTruthMuonHasRecoTrack",
    Comment = "Has truth muon exiting, recoed as a track",
  )
)

# ---> that reco track is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonExitingHasRecoTrackContained",
    Cut = "!cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && cutTruthMuonMatchedTrackContained",
    Comment = "Has truth muon exiting, recoed as a contained track",
  )
)

# -----> Applying nominal reco cut
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonExitingHasRecoTrackContained_NominalContainedMuon",
    Cut = "!cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && cutTruthMuonMatchedTrackContained && cutNominal_ContainedMuon",
    Comment = "Has truth muon exiting, recoed as a contained track",
  )
)

# ---> that reco track is exiting
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonExitingHasRecoTrackExiting",
    Cut = "!cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && !cutTruthMuonMatchedTrackContained",
    Comment = "Has truth muon exiting, recoed as a exiting track",
  )
)

# -----> Applying nominal reco cut
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthMuonExitingHasRecoTrackExiting_NominalExitingMuon",
    Cut = "!cutTruthMuonContained && cutHasTruthMuonHasRecoTrack && !cutTruthMuonMatchedTrackContained && cutNominal_ExitingMuon",
    Comment = "Has truth muon exiting, recoed as a exiting track",
  )
)

#------------------------------------------------------------------
# Truth proton is found from the primary particles

EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthProton",
    Cut = "cutHasTruthProton",
    Comment = "Has truth proton",
  )
)

# Truth proton is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthProtonContained",
    Cut = "cutTruthProtonContained",
    Comment = "Has truth proton contained",
  )
)

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

EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPion",
    Cut = "cutHasTruthChargedPion",
    Comment = "Has truth charged pion",
  )
)

# Truth charged pion is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionContained",
    Cut = "cutTruthChargedPionContained",
    Comment = "Has truth charged pion contained",
  )
)

# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionContainedHasRecoTrack",
    Cut = "cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack",
    Comment = "Has truth charged pion contained, recoed as a track",
  )
)

# ---> that reco track is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionContainedHasRecoTrackContained",
    Cut = "cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack && cutTruthChargedPionMatchedTrackContained",
    Comment = "Has truth charged pion contained, recoed as a contained track",
  )
)

# ---> that reco track is exiting
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionContainedHasRecoTrackExiting",
    Cut = "cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack && !cutTruthChargedPionMatchedTrackContained",
    Comment = "Has truth charged pion contained, recoed as a exiting track",
  )
)

# Truth charged pion is exiting
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionExiting",
    Cut = "!cutTruthChargedPionContained",
    Comment = "Has truth charged pion exiting",
  )
)

# -> Matched to a reco track
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionExitingHasRecoTrack",
    Cut = "!cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack",
    Comment = "Has truth charged pion exiting, recoed as a track",
  )
)

# ---> that reco track is contained
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionExitingHasRecoTrackContained",
    Cut = "!cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack && cutTruthChargedPionMatchedTrackContained",
    Comment = "Has truth charged pion exiting, recoed as a contained track",
  )
)

# ---> that reco track is exiting
EventSelectionInfos.append(
  SelectionInfo(
    Name = "TruthChargedPionExitingHasRecoTrackExiting",
    Cut = "!cutTruthChargedPionContained && cutHasTruthChargedPionHasRecoTrack && !cutTruthChargedPionMatchedTrackContained",
    Comment = "Has truth charged pion exiting, recoed as a exiting track",
  )
)










