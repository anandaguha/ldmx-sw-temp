#ifndef DQM_SAMPLEVALIDATION_H
#define DQM_SAMPLEVALIDATION_H

//LDMX Framework
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h"

namespace dqm {

    /**
     * @class SampleValidation
     * @brief
     */

    class SampleValidation : public framework::Analyzer {
        public:

            SampleValidation(const std::string& name, framework::Process& process) : Analyzer(name, process) {}

            virtual void configure(framework::config::Parameters& ps);

            virtual void analyze(const framework::Event& event);
        private:

    };

}

#endif
