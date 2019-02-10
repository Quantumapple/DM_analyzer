// system include files
#include <memory>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//ROOT required dependencies
#include "TH1D.h"
#include "TH1F.h"
#include "TTree.h"
#include "TLorentzVector.h"

//Helper
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

// Analysis-specified
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DM_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DM_Analyzer(const edm::ParameterSet&);
      ~DM_Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override; // Job starts
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override; // Loop every events
      virtual void endJob() override; // filling trees

      // ----------member data ---------------------------

  std::vector<reco::GenParticle> IndexByPt(std::vector<reco::GenParticle> vector);

  struct comp {
    bool operator() (reco::GenParticle& i,reco::GenParticle& j) { return ( (i.pt()) > (j.pt()) ); } // sort in descending order
  };

  //Token
  edm::EDGetTokenT<reco::GenParticleCollection> tok_gen; // object interacts with python <-> c++
  edm::EDGetTokenT<pat::METCollection> tok_met; // object interacts with python <-> c++
  edm::EDGetTokenT<pat::JetCollection> tok_jet;
  edm::EDGetTokenT<pat::GenJetCollection>  tok_genjet;

  //InputTag
  edm::InputTag genParticleTag_;
  edm::InputTag metTag_;
  edm::InputTag jetTag_;
  edm::InputTag genJetTag_;

  //LorentzVector
  TLorentzVector dmSystem;
  TLorentzVector dm;

  //TH1D
  TH1D *Njet;
  TH1D *pdg;
  TH1D *status;

  TH1F *ZpPt;
  TH1F *ZpPhi;
  TH1F *ZpMass;

  TH1F *DmPt;
  TH1F *DmMass;

  TH1F *hsPt;
  TH1F *hsMass;

  // Missing transverse energy
  TH1F *metPt;

  //unit cross section pb
  float unitxsec=1;

  //Tree for acceptance
  Double_t met[7];
  Float_t acc[7];
  TTree *tree;

  //variable
  double countEvent = 0;
  double countEventPass[7]={0.,0.,0.,0.,0.,0.,0.};

  int counter;

  // Check DM is working
  std::vector<bool> DM_flag;

};
