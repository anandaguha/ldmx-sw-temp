"""Configuration classes for sensitive detectors"""

from LDMX.SimCore import simcfg

class ScoringPlaneSD(simcfg.SensitiveDetector) :
    """Scoring plane SD

    Simply collecting tracker-hit equivalent for scoring planes that
    enclose different subsystems

    Parameters
    ----------
    subsystem : str
        Name of subsystem to store scoring plane hits for
        Names must match what is in gdml for <subsystem>_sp
    """
    def __init__(self,subsystem) :
        super().__init__(f'{subsystem}_sp','simcore::ScoringPlaneSD','SimCore_SDs')

        self.collection_name = f'{subsystem.capitalize()}ScoringPlaneHits'
        self.match_substr = f'sp_{subsystem.lower()}' #depends on gdml

    def ecal() :
        return ScoringPlaneSD('ecal')

    def hcal() :
        return ScoringPlaneSD('hcal')

    def target() :
        return ScoringPlaneSD('target')

    def magnet() :
        return ScoringPlaneSD('magnet')

    def tracker() :
        sp = ScoringPlaneSD('tracker')
        sp.match_substr = 'sp_recoil'
        return sp

class TrackerSD(simcfg.SensitiveDetector) :
    """SD for the recoil and tagging trackers

    Parameters
    ----------
    subsystem : str
        Recoil or Tagger
    subdet_id : int
        ID number for the subsystem
    """
    def __init__(self,subsystem,subdet_id) :
        super().__init__(f'{subsystem}_TrackerSD','simcore::TrackerSD','SimCore_SDs')

        self.subsystem = subsystem
        self.subdet_id = subdet_id

        self.collection_name = f'{subsystem}SimHits'

    def tagger() :
        return TrackerSD('Tagger',1)

    def recoil() :
        return TrackerSD('Recoil',4)

class HcalSD(simcfg.SensitiveDetector) :
    """SD for the HCal

    Separate from the other calorimeters since it includes a Birks law
    estimate. No other parameters
    
    """
    def __init__(self, gdml_identifiers = ['ScintBox', 'scint_box']) :
        super().__init__('hcal_sd', 'simcore::HcalSD','SimCore_SDs')
        self.gdml_identifiers = gdml_identifiers

class EcalSD(simcfg.SensitiveDetector) :
    """SD for the ECal

    The two configurable parameters are inherited from a legacy method of
    merging simulated hit contribs. We have plans to update this hit merging
    in the future.

    Parameters
    ----------
    enableHitContribs : bool, optional
        Should the simulation save contributions to Ecal sim hits?
    compressHitContribs : bool, optional
        Should the simulation compress contributions to Ecal sim hits by PDG ID?
    """
    def __init__(self) :
        super().__init__('ecal_sd', 'simcore::EcalSD','SimCore_SDs')
        self.enableHitContribs = True
        self.compressHitContribs = True

class TrigScintSD(simcfg.SensitiveDetector) :
    """Trigger Scintillaotr Sensitive Detector

    used for both the trigger pad modules as well as collecting hits
    within the target itself

    Parameters
    ----------
    module : int
        ID number for the module we are collecting hits from
    name : str
        Short name to be used in building collection name
    vol : str
        Name of logical volume(s) that this SD should be attached to
        DEPENDS ON GDML
    """
    def __init__(self, module, name, vol) :
        super().__init__(f'trig_scint_{name}_sd', 'simcore::TrigScintSD','SimCore_SDs')
        self.module_id = module
        self.volume_name = vol

        coll = name+'SimHits'
        if name != 'Target' :
            coll = 'TriggerPad'+coll

        self.collection_name = coll

    def up() :
        return TrigScintSD(2,'Up','trigger_pad_up_bar_volume')

    def tag() :
        return TrigScintSD(1,'Tagger','trigger_pad_tag_bar_volume')

    def down() :
        return TrigScintSD(3,'Down','trigger_pad_dn_bar_volume')

    def target() :
        return TrigScintSD(4,'Target','target')

