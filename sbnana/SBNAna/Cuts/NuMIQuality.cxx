#include "sbnana/SBNAna/Cuts/NuMIQuality.h"

namespace ana {

  //// ----------------------------------------------
  NuMIQuality::NuMIQuality()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIQuality: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    filename_bad_triggered_spill = std::string(sbndata) +
                   "beamData/NuMIdata/bad_triggered_spills.csv";

    std::cout << "[NuMIQuality::NuMIQuality] filename_bad_triggered_spill = " << filename_bad_triggered_spill << std::endl;

  }

  NuMIQuality::~NuMIQuality()
  {
  }

  bool NuMIQuality::IsBadTriggeredSpill(unsigned int run, unsigned int subrun, unsigned int event) const {

    std::ifstream csv_bad_triggered_spill(filename_bad_triggered_spill);
    std::string line;
    std::getline(csv_bad_triggered_spill, line); // skip first line

    while (std::getline(csv_bad_triggered_spill, line)) {
      std::istringstream iss(line);
      unsigned int entry, currentRun, currentSubrun, currentEvt;
      char comma;  // to read the commas in the CSV
      if (!(iss >> entry >> comma >> currentRun >> comma >> currentSubrun >> comma >> currentEvt)) {
        std::cerr << "Error reading line from file: " << filename_bad_triggered_spill << std::endl;
        abort();
      }

      // Check if the pair matches
      if (currentRun == run && currentSubrun == subrun && currentEvt == event) {
        printf("[BAD Trigger] (run, subrun, event) = (%d, %d, %d)\n", run, subrun, event);
        return true;  // Pair found
      }
    }
    return false;

  }


  const SpillCut kNuMINotBadTriggeredSpill( [](const caf::SRSpillProxy *sr) {

    if(sr->hdr.ismc) return true;

    if( kNuMIQuality.IsBadTriggeredSpill(sr->hdr.run, sr->hdr.subrun, sr->hdr.evt) ) return false;
    else return true;

  });


}
