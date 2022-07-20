from SelectionInfo import *

EventSelectionInfos = []


MyCutLists = [
["cutRFiducial","RecoFiducial"],
["cutFMScore","FMScore"],
["cutFMTime","FMTime"],
["cutSliceCRLongestTrackDirY","CRLongestTrackDirY"],
["cutHasMuon","HasMuon"],
["cutHasProton","HasProton"],
["cutMuonProtonCosineTheta","MuonProtonCosineTheta"],
]

FullSelection = "kNoCut"
for c in MyCutLists:
  cSTR = c[0]
  cName = c[1]
  FullSelection += " && "+cSTR

EventSelectionInfos.append(
  SelectionInfo(
    Name = "FullSelection",
    Cut = FullSelection,
    Comment = "Full selection",
  )
)

for i_c in range(0,len(MyCutLists)):
  thisCutSTR = MyCutLists[i_c][0]
  thisCutName = MyCutLists[i_c][1]

  thisSelection = "kNoCut"
  for j_c in range(0,len(MyCutLists)):
    tmpCutSTR = MyCutLists[j_c][0]
    tmpCutName = MyCutLists[j_c][1]
    if tmpCutSTR==thisCutSTR:
      continue
    thisSelection += " && "+tmpCutSTR

  EventSelectionInfos.append(
  SelectionInfo(
    Name = "NMinus_"+thisCutName,
    Cut = thisSelection,
    Comment = "Full excepct %s"%(MyCutLists[i_c][1]),
  )
)

'''
for EventSelectionInfo in EventSelectionInfos:
  print("============================================")
  print("Name = %s"%(EventSelectionInfo.Name))
  print("Cut = %s"%(EventSelectionInfo.Cut))
'''
