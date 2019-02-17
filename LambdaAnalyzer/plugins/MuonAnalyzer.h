#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

class MuonAnalyzer {
    public:
        MuonAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~MuonAnalyzer();
        virtual std::vector<pat::Muon> FillMuonVector(const edm::Event&);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Muon> > MuonToken;

}; 


#endif
