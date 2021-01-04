/**
 * @file EcalTrigPrimDigiProducer.h
 * @brief Class that performs emulation of the EcalTriggerPrimitives
 * @author Jeremiah Mans, University of Minnesota
 */

#ifndef ECAL_ECALTRIGPRIMDIGIPRODUCER_H_
#define ECAL_ECALTRIGPRIMDIGIPRODUCER_H_

//----------------//
//   LDMX Core    //
//----------------//
#include "Framework/EventProcessor.h"

namespace ldmx {

/**
 * @class EcalRecProducer
 * @brief Performs basic ECal reconstruction
 *
 * Reconstruction is done from the EcalDigi samples.
 * Some hard-coded parameters are used for position and energy calculation.
 */
class EcalTrigPrimDigiProducer : public Producer {
 public:
  /**
   * Constructor
   */
  EcalTrigPrimDigiProducer(const std::string& name, Process& process);

  /**
   * Grabs configure parameters from the python config file.
   *
   * Parameter        Default
   * inputDigiCollName     EcalDigis
   * inputDigiPassName     "" <-- blank means take any pass if only one
   * collection exists
   */
  virtual void configure(Parameters&);

  /**
   * Produce EcalHits and put them into the event bus using the
   * EcalDigis as input.
   */
  virtual void produce(Event& event);

 private:
  /** Digi Collection Name to use as input */
  std::string digiCollName_;

  /** Digi Pass Name to use as input */
  std::string digiPassName_;

  /** Conditions object for the calibration information */
  std::string condObjName_;
};
}  // namespace ldmx

#endif  // EVENTPROC_ECALTRIGPRIMDIGIPRODUCER_H_INC
