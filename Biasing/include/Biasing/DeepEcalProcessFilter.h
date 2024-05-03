/**
 * @file DeepEcalProcessFilter.h
 * @class DeepEcalProcessFilter
 * @brief Class defining a UserActionPlugin that allows a user to filter out
 *        events where the interaction happened deep in the ECAL
 * @author Tamas Almos Vami (UCSB)
 * @author Anmol Sandhu, Tyler Horoho (UVA)
 */

#ifndef BIASING_DEEPECALFILTER_H
#define BIASING_DEEPECALFILTER_H

//----------------//
//   C++ StdLib   //
//----------------//
#include <algorithm>

/*~~~~~~~~~~~~~*/
/*   SimCore   */
/*~~~~~~~~~~~~~*/
#include "SimCore/UserAction.h"

/*~~~~~~~~~~~~~~~*/
/*   Framework   */
/*~~~~~~~~~~~~~~~*/
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h"

namespace biasing {

/**
 * User action that allows a user to filter out events where the interaction happened deep in the ECAL
 */
class DeepEcalProcessFilter : public simcore::UserAction {
 public:
  /// Constructor
  DeepEcalProcessFilter(const std::string& name,
                    framework::config::Parameters& parameters);

  /// Destructor
  virtual ~DeepEcalProcessFilter() = default;
  
  /// Method to set flags in the beginning of the event
  void BeginOfEventAction(const G4Event* event) final override;

  /**
   * Implement the stepping action which performs the target volume biasing.
   * @param step The Geant4 step.
   */
  void stepping(const G4Step* step) final override;

//  /**
//   * Method called at the end of every event.
//   * @param event Geant4 event object.
//   */
//  void EndOfEventAction(const G4Event*) final override;
  
  void NewStage() final override;

  /// Retrieve the type of actions this class defines
  std::vector<simcore::TYPE> getTypes() final override {
    return {simcore::TYPE::STACKING,
            simcore::TYPE::STEPPING,simcore::TYPE::EVENT};
  }

 private:
  /// Minimal energy the products  should have
  double bias_threshold_{1500.};
  ///  The allowed processes that can happen deep inside the ECAL, default is
  ///  conversion (conv) and photoelectron (photo)
  std::vector<std::string> processes_{"conv","phot"};
  /// Minimum Z location where the deep process should happen
  double ecal_min_Z_{400.};
  /// Enable logging
  enableLogging("DeepEcalProcessFilter")

};  // DeepEcalProcessFilter
}  // namespace biasing

#endif  // BIASING_DEEPECALFILTER_H