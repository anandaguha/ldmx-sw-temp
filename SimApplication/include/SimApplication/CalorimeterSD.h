#ifndef CALORIMETERSD_SIMAPPLICATION_H_
#define CALORIMETERSD_SIMAPPLICATION_H_ 1

// Geant4
#include "G4VSensitiveDetector.hh"

// LDMX
#include "DetDescr/DetectorID.h"
#include "Event/Event.h"
#include "SimApplication/G4CalorimeterHit.h"


class CalorimeterSD: public G4VSensitiveDetector {

    public:

        CalorimeterSD(G4String theName, G4String theCollectionName, int theSubdetId, DetectorID* theDetId);

        virtual ~CalorimeterSD();

        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

        void Initialize(G4HCofThisEvent* hcEvent);

        void EndOfEvent(G4HCofThisEvent* hcEvent);

    private:

        G4CalorimeterHitsCollection* hitsCollection;
        Event* currentEvent;
        int subdetId;
        DetectorID* detId;
};

#endif