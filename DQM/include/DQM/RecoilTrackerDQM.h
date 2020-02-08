/**
 * @file RecoilTrackerDQM.h
 * @brief Analyzer used for DQM of the Recoil tracker. 
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _DQM_RECOIL_TRACKER_DQM_H_
#define _DQM_RECOIL_TRACKER_DQM_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <unordered_map>
#include <utility>

//----------//
//   LDMX   //
//----------//
#include "Tools/AnalysisUtils.h"
#include "Framework/EventProcessor.h"
#include "Framework/Event.h"
#include "Framework/HistogramPool.h"
#include "Event/EventDef.h"

//----------//
//   ROOT   //
//----------//
#include "TH1.h"
#include "TVector3.h"

namespace ldmx { 

    class RecoilTrackerDQM : public Analyzer { 


        public: 

            /** Constructor */
            RecoilTrackerDQM(const std::string &name, Process &process);

            /** Destructor */
            ~RecoilTrackerDQM(); 
            
            /** 
             * Configure the processor using the given user specified parameters.
             * 
             * @param parameters Set of parameters used to configure this processor.
             */
            void configure(std::map < std::string, std::any > parameters) final override;

            /**
             * Process the event and make histograms ro summaries.
             *
             * @param event The event to analyze.
             */
            void analyze(const Event& event);

            /** Method executed before processing of events begins. */
            void onProcessStart();


        private: 

            /** Singleton used to access histograms. */
            HistogramPool* histograms_{nullptr}; 

            /** Name of ECal veto collection. */
            std::string ecalVetoCollectionName_{"EcalVeto"}; 
    
    }; // RecoilTrackerDQM 
    
} // ldmx

#endif // _DQM_RECOIL_TRACKER_DQM_H_
