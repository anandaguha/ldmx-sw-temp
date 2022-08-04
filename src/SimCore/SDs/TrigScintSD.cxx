#include "SimCore/SDs/TrigScintSD.h"

/*~~~~~~~~~~~~~~*/
/*   DetDescr   */
/*~~~~~~~~~~~~~~*/
#include "DetDescr/TrigScintID.h"

/*~~~~~~~~~~~~*/
/*   Geant4   */
/*~~~~~~~~~~~~*/
#include "G4Step.hh"
#include "G4StepPoint.hh"

namespace simcore {

TrigScintSD::TrigScintSD(const std::string& name,
                         simcore::ConditionsInterface& ci,
                         const framework::config::Parameters& p)
    : SensitiveDetector(name, ci, p) {
      auto which = p.getParameter<std::string>("which");
      collection_name_ = which+"SimHits"; // default collection name
      if (which != "Target") {
        collection_name_ = "TriggerPad"+collection_name_;
      }
      collection_name_ = p.getParameter<std::string>("collection_name",collection_name_);


      static const std::map<std::string,std::string> which_to_vol = {
        {"Up", "trigger_pad_up_bar_volume"},
        {"Down", "trigger_pad_dn_bar_volume"},
        {"Tagger", "trigger_pad_tag_bar_volume"},
        {"Target", "target"}
      };
      try {
        vol_name_ = which_to_vol.at(which);
      } catch (const std::out_of_range&) {
        EXCEPTION_RAISE("SDConfig","Trigger Pad '"+which+"' not one of Up, Down, or Tagger");
      }
    }

TrigScintSD::~TrigScintSD() {}

G4bool TrigScintSD::ProcessHits(G4Step* step, G4TouchableHistory* history) {
  // Get the energy deposited by the particle during the step
  auto energy{step->GetTotalEnergyDeposit()};

  // If a non-Geantino particle doesn't deposit energy during the step,
  // skip processing it.
  if (energy == 0 and not isGeantino(step)) return false;

  // Create a new instance of a calorimeter hit
  //  emplace_back returns a *reference* to the hit that was constructed
  //  and we should keep that reference so that we are editing the correct hit
  ldmx::SimCalorimeterHit& hit = hits_.emplace_back();

  // Set the hit position
  auto position{0.5 * (step->GetPreStepPoint()->GetPosition() +
                       step->GetPostStepPoint()->GetPosition())};
  auto volumePosition{step->GetPreStepPoint()
                          ->GetTouchableHandle()
                          ->GetHistory()
                          ->GetTopTransform()
                          .Inverse()
                          .TransformPoint(G4ThreeVector())};
  hit.setPosition(position[0], position[1], volumePosition.z());

  // Get the track associated with this step
  auto track{step->GetTrack()};

  // Set the ID on the hit.
  int module = getModuleID(track->GetVolume()->GetLogicalVolume());
  auto bar{track->GetVolume()->GetCopyNo()};
  ldmx::TrigScintID id(module, bar);
  hit.setID(id.raw());

  // add single contrib to this calorimeter hit
  //  IncidentID - this track's ID
  //  Track ID
  //  PDG ID
  //  energy deposited
  //  global time of this hit
  hit.addContrib(track->GetTrackID(), track->GetTrackID(),
                 track->GetParticleDefinition()->GetPDGEncoding(), energy,
                 track->GetGlobalTime());

  return true;
}

int TrigScintSD::getModuleID(G4LogicalVolume* vol) const {
  /**
   * tag <-> 1
   * up  <-> 2
   * dn  <-> 3
   */
  if (vol->GetName().contains("tag")) return 1;
  if (vol->GetName().contains("up" )) return 2;
  if (vol->GetName().contains("dn" )) return 3;
  if (vol->GetName().contains("target")) return 4;
}

}  // namespace simcore

DECLARE_SENSITIVEDETECTOR(simcore::TrigScintSD)
