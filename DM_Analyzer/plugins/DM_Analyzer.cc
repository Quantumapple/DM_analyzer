// -*- C++ -*-
//
// Package:    EDMAnalyzer/DM_Analyzer
// Class:      DM_Analyzer
// 
/**\class DM_Analyzer DM_Analyzer.cc EDMAnalyzer/DM_Analyzer/plugins/DM_Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Siew Yan Hoh
//         Created:  Wed, 11 May 2016 15:18:40 GMT
//
//


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
  //edm::EDGetTokenT<pat::MET> tok_met; // object interacts with python <-> c++
  edm::EDGetTokenT<pat::METCollection> tok_met; // object interacts with python <-> c++

  //InputTag
  edm::InputTag genParticleTag_;
  edm::InputTag metTag_;
  edm::InputTag jetTag_;
  edm::InputTag genJetTag_;

  //vector
  //std::vector<reco::GenParticle> bjets;
  //std::vector<reco::GenParticle> thirdjet;

  //LorentzVector
  TLorentzVector dmSystem;
  TLorentzVector dm;

  //TH1D
  TH1D *Njet;
  TH1D *pdg;
  TH1D *status;

  //TH1F
  //TH1F *ptb1;
  //TH1F *ptb2;
  //TH1F *etab1;
  //TH1F *etab2;
  //TH1F *phib1;
  //TH1F *phib2;

  TH1F *phipt;
  TH1F *etaphi;
  TH1F *phimass;

  TH1F *chichipt;
  TH1F *chichimass;

  TH1F *hspt;

  // Missing transverse energy
  TH1F *metEt;
  TH1F *metpT;
  TH1F *diffAng;

  //light jet
  TH1F *phij1;
  TH1F *phij2;
  TH1F *etaj1;
  TH1F *etaj2;
  TH1F *ptj1;
  TH1F *ptj2;
  TH1F *deltaRj1j2;

  //bjet1 daughter
  TH1D *Ndoug;
  TH1D *dougstatus;
  TH1D *dougpdgid;

  //unit cross section pb
  float unitxsec=1;

  //Tree for acceptance
  Double_t met[7];
  Float_t acc[7];
  TTree *a;

  //variable
  double countEvent = 0;
  double countEventPass[7]={0.,0.,0.,0.,0.,0.,0.};

  int counter;

  // Check DM is working
  std::vector<bool> DM_flag;  

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DM_Analyzer::DM_Analyzer(const edm::ParameterSet& iConfig):
  genParticleTag_(iConfig.getUntrackedParameter<edm::InputTag>("genParticleTag")),
  metTag_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  jetTag_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  genJetTag__(iConfig.getUntrackedParameter<edm::InputTag>("genJetTag"))

{
   //now do what ever initialization is needed
   usesResource("TFileService"); // CMSSW service to save histograms, leave this!!

   //Token
   tok_gen = consumes<reco::GenParticleCollection>(genParticleTag_); // bridge python <-> c++
   tok_met = consumes<pat::METCollection>(metTag_); // bridge python <-> c++
   tok_jet = consumes<pat::JetCollection>(jetTag_);
   tok_genjet = consumes<pat::GenJetCollection>(genjetTag_);

   edm::Service<TFileService> fs;
   
   //global 
   //Njet       = fs->make<TH1D>( "Njet"   , "Jet Multiplicity"     , 10  , -0.5 , 9.5 );
   pdg        = fs->make<TH1D>( "pdg"    , "pdgid in event"       , 30  , -30. , 30. );
   status     = fs->make<TH1D>( "status" , "status code of Event" , 100 , 0.   , 100 );

   //process kinematics
   //ptb1       = fs->make<TH1F>( "ptb1"  , "P_{T} of the leading b jet; Pt_{bjet1} [GeV/c]; Events"    , 100 , 0.  , 1000. );
   //ptb2       = fs->make<TH1F>( "ptb2"  , "P_{T} of second leading b jet; Pt_{bjet2} [GeV/c]; Events" , 100 , 0.  , 1000. );
   //etab1      = fs->make<TH1F>( "etab1" , "#eta of the leading b jet; #eta_{bjet1}; Events"            , 100 , -8. , 8.   );
   //etab2      = fs->make<TH1F>( "etab2" , "#eta of the second leading b jet; #eta_{bjet2} ; Events"    , 100 , -8. , 8.   );
   //phib1      = fs->make<TH1F>( "phib1" , "#phi of the leading b jet; #phi_{bjet1}; Events"            , 100 , -3.5 , 3.5   );
   //phib2      = fs->make<TH1F>( "phib2" , "#phi of the second leading b jet; #phi_{bjet2} ; Events"    , 100 , -3.5 , 3.5   );

   phipt      = fs->make<TH1F>( "phipt"   , "P_{T} of Z'; Pt_{Z'} [GeV/c]; Events"      , 100 , 0.  , 500.  );
   etaphi     = fs->make<TH1F>( "etaphi"  , "Eta of Z'; #eta_{Z'} ; Events"             , 100 , -5. , 5.    );
   phimass    = fs->make<TH1F>( "phimass" , "Inv mass of Z'; M(Z') [GeV/c^{2}]; Events" , 100 , 0.  , 1500. );

   chichipt   = fs->make<TH1F>( "chichipt"   , "P_{T} of #chi#tilde{#chi}; Pt_{#chi#tilde{#chi}} [GeV/c]; Events" , 100 , 0. , 500.  );
   chichimass = fs->make<TH1F>( "chichimass" , "Inv mass of #chi#tilde{#chi}; M(#chi#chi) [GeV/c^{2}]; Events"    , 100 , 0. , 1500. );

   hspt      = fs->make<TH1F>( "hspt"   , "P_{T} of hs; Pt_{hs} [GeV/c]; Events"      , 100 , 0.  , 500.  );

   // light jet
   ptj1       = fs->make<TH1F>( "ptj1"  , "P_{T} of the first light jet; Pt_{jet1} [GeV/c]; Events"   , 100 , 0.  , 1000. );
   ptj2       = fs->make<TH1F>( "ptj2"  , "P_{T} of the second light jet; Pt_{jet2} [GeV/c]; Events"  , 100 , 0.  , 1000. );
   etaj1      = fs->make<TH1F>( "etaj1" , "Eta of the first light jet; #eta_{jet1}; Events"           , 100 , -8. , 8.   );
   etaj2      = fs->make<TH1F>( "etaj2" , "Eta of the second light jet; #eta_{jet2}; Events"          , 100 , -8. , 8.   );
   phij1      = fs->make<TH1F>( "phij1" , "Phi of the first light jet; #phi_{jet1}; Events"           , 100 , -3.5, 3.5  );
   phij2      = fs->make<TH1F>( "phij2" , "Phi of the second light jet; #phi_{jet2}; Events"          , 100 , -3.5, 3.5  );
   deltaRj1j2 = fs->make<TH1F>( "deltaRj1j2" , "#Delta R of the first and the second light jet; #Delta R {jet1, jet2}; Events" , 100 , 0. , 5.);

   //bjet1 daughter
   //Ndoug = fs->make<TH1D>( "Ndoug"   , "leading bjet number of daugther"     , 10  , -0.5 , 9.5 );
   //dougstatus = fs->make<TH1D>( "dougstatus"   , "leading bjet daughter status"     , 100  , -0.5 , 99.5 );
   //dougpdgid = fs->make<TH1D>( "dougpdgid"   , "leading bjet daughter pdgid"     , 10  , -0.5 , 9.5 );

   // MET
   metEt    = fs->make<TH1F>(  "metEt"    , "Missing transverse energy; Energy [GeV]; Events", 100 , 0., 1000. );
   metpT    = fs->make<TH1F>(  "metpT"    , "P_{T} of the missing transverse energy; Pt_{MET} [GeV/c]; Events", 100 , 0., 1000. );
   diffAng  = fs->make<TH1F>(  "diffAng"  , "Angle difference between MET and 1st leading jet; #Delta#phi; Events", 100, -3.5, 3.5 );

   //Tree
   a = fs->make<TTree>("a","Acceptance");
   a->Branch("met",met,"met[7]/D");
   a->Branch("acc",acc,"acc[7]/F");

   //variable initialization
   //bjets.clear();
   //thirdjet.clear();
}


DM_Analyzer::~DM_Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
std::vector<reco::GenParticle> DM_Analyzer::IndexByPt(std::vector<reco::GenParticle> vector)
{
  comp comparator;

  std::sort (vector.begin() , vector.end() , comparator);
  return vector;
}

// ------------ method called for each event  ------------
void
DM_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // PYTHIA convention -> check internet
  //status code:                                                                                                                                           
  // 1 = existing entry - not decayed or fragmented, represents the final state as given by the generator                                                  
  // 2 = decayed or fragmented entry (i.e. decayed particle or parton produced in shower.)                                                                 
  // 3 = identifes the "hard part" of the interaction, i.e. the partons that are used in the matrix element calculation, including immediate decays of resonances.                                                                   
   // CMSSW service name
   using namespace edm;
   using namespace reco;
   using namespace pat; 
   using namespace std;
   //using reco::GenParticleCollection;

   counter=0;

   //need to declare a vector storing reco::candidates
   //const reco::Candidate *bjet = 0;
   //std::vector<reco::GenParticle> bjets;
   std::vector<reco::GenParticle> thirdjet;
   //bjets.clear();
   thirdjet.clear();

   //TLorentzvector
   dmSystem.SetPxPyPzE(0.,0.,0.,0.);
   dm.SetPxPyPzE(0.,0.,0.,0.);

   // CMSSW service leave this!!
   edm::Handle<GenParticleCollection> genpart;
   edm::Handle<METCollection> met;
   edm::Handle<JetCollection> jet;
   edm::Handle<GenJetCollection> genjet;
   
  
   iEvent.getByToken(tok_gen, genpart);
   iEvent.getByToken(tok_met, met);
   iEvent.getByToken(tok_jet, jet);
   iEvent.getByToken(tok_genjet, genjet);

   for( size_t i=0 ; i < genpart->size() ; ++i )
     {
       const reco::GenParticle &p = (*genpart)[i];
       const reco::Candidate *mom = p.mother();

       dm.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());

       pdg->Fill(p.pdgId());
       status->Fill(p.status());  


       //bjets
       //if ( abs(p.pdgId()) == 5 && ( p.status() > 70 && p.status() < 80 ) ) { 
       //    bjets.push_back(p); 
       //}

       //if( abs(p.pdgId()) > 50 && abs(p.pdgId()) < 60 ) std::cout << "Particle ID: " << p.pdgId() << ", number of daughters: " << p.numberOfDaughters() << std::endl;

       //dark matter pair
       if ( /*p.fromHardProcessFinalState() && */ abs(p.pdgId())== 52 && p.numberOfDaughters() == 0 ) { 
           dmSystem+=dm;
           chichipt->Fill(dmSystem.Pt(),unitxsec);
           chichimass->Fill(dmSystem.M(),unitxsec);
       }
       //mediator
       if ( p.status() == 22 && abs(p.pdgId())==55 )
       {
           phipt->Fill(p.pt(),unitxsec);
           etaphi->Fill(p.rapidity(),unitxsec);
           phimass->Fill(p.mass(),unitxsec);
       }
       //dark higgs
       //if ( abs(p.pdgId())==54 ) {
       
       //}
       
     }


   //Njet->Fill(counter);
   
   //std::vector<reco::GenParticle> Sorted_bjets = IndexByPt(bjets);
   std::vector<reco::GenParticle> Sorted_thirdjet = IndexByPt(thirdjet);
   
   //
   for( size_t j=0 ; j < met->size() ; ++j )
   {
       const pat::MET &MeT = (*var2)[j];

       metEt->Fill(MeT.energy(),unitxsec);
       metpT->Fill(MeT.pt(),unitxsec);
       
       if( Sorted_thirdjet.size() > 0 ) {
	   float dphi = Sorted_thirdjet[0].phi() - MeT.phi();
	   if( dphi >= M_PI ) dphi -= 2.*M_PI;
	   if( dphi < -M_PI ) dphi += 2.*M_PI;
           diffAng->Fill(dphi,unitxsec);
       }
       //cout << "Missing transverse energy: " << MeT.energy() << endl;
   }

   /*
   if (Sorted_bjets.size() > 0){
     ptb1->Fill(Sorted_bjets[0].pt(),unitxsec);
     etab1->Fill(Sorted_bjets[0].rapidity(),unitxsec);
     phib1->Fill(Sorted_bjets[0].phi(),unitxsec);
     ptb2->Fill(Sorted_bjets[1].pt(),unitxsec);
     etab2->Fill(Sorted_bjets[1].rapidity(),unitxsec);
     phib2->Fill(Sorted_bjets[1].phi(),unitxsec);
   }
   */
   if (Sorted_thirdjet.size() > 0 ){
     ptj1->Fill(Sorted_thirdjet[0].pt(),unitxsec); 
     ptj2->Fill(Sorted_thirdjet[1].pt(),unitxsec); 
     etaj1->Fill(Sorted_thirdjet[0].rapidity(),unitxsec);
     etaj2->Fill(Sorted_thirdjet[1].rapidity(),unitxsec);
     phij1->Fill(Sorted_thirdjet[0].phi(),unitxsec);
     phij2->Fill(Sorted_thirdjet[1].phi(),unitxsec);

     if( Sorted_thirdjet[0].pt() != 0 && Sorted_thirdjet[1].pt() != 0 ) {
	float dphi = Sorted_thirdjet[0].phi() - Sorted_thirdjet[1].phi();
	if( dphi >= M_PI ) dphi -= 2.*M_PI;
	if( dphi < -M_PI ) dphi += 2.*M_PI;
	float deta = Sorted_thirdjet[0].rapidity() - Sorted_thirdjet[1].rapidity();
        float dR = sqrt(pow(dphi,2)+pow(deta,2));
        deltaRj1j2->Fill(dR,unitxsec);
     }
   }
   

   /*   Ndoug->Fill(Sorted_bjets[0].numberOfDaughters());

   size_t n = Sorted_bjets[0].numberOfDaughters();
   const reco::GenParticle *pp = &Sorted_bjets[0];

   for(size_t j = 0; j < n; ++ j) {
     const reco::Candidate &d = (*pp->daughter(j));
     
     dougstatus->Fill(d.status());
     dougpdgid->Fill(d.pdgId());
     }*/
   
   
   //chichipt->Fill(dmSystem.Pt(),unitxsec);
   //chichimass->Fill(dmSystem.M(),unitxsec);

   //acceptance
   countEvent++;
   for(Int_t i=0; i<7; i++)
   {
       if(dmSystem.Pt()<(double) 200+(i*50)) continue;
       countEventPass[i]+=1;
   }
     
   
   
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DM_Analyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DM_Analyzer::endJob() 
{

  for(Int_t i=0; i<7; i++)
    { 
      met[i] = (double) (200.+(i*50.));
      acc[i] = (100.*countEventPass[i]/countEvent); 
    }
  
  a->Fill();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DM_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DM_Analyzer);
