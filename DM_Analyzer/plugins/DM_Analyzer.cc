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
#include "DM_analyzer/DM_Analyzer/interface/DM_Analyzer.h"
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
   Njet       = fs->make<TH1D>( "Njet"   , "Jet Multiplicity"     , 10  , -0.5 , 9.5 );
   pdg        = fs->make<TH1D>( "pdg"    , "pdgid in event"       , 30  , -30. , 30. );
   status     = fs->make<TH1D>( "status" , "status code of Event" , 100 , 0.   , 100 );

   ZpPt      = fs->make<TH1F>( "ZpPt"   , "P_{T} of Z'; Pt_{Z'} [GeV/c]; Events"      , 100 , 0.  , 500.  );
   ZpPhi     = fs->make<TH1F>( "ZpPhi"  , "Eta of Z'; #eta_{Z'} ; Events"             , 100 , -5. , 5.    );
   ZpMass    = fs->make<TH1F>( "ZpMass" , "Inv mass of Z'; M(Z') [GeV/c^{2}]; Events" , 100 , 0.  , 1500. );

   DmPt   = fs->make<TH1F>( "DmPt"   , "P_{T} of #chi#tilde{#chi}; Pt_{#chi#tilde{#chi}} [GeV/c]; Events" , 100 , 0. , 500.  );
   DmMass = fs->make<TH1F>( "DmMass" , "Inv mass of #chi#tilde{#chi}; M(#chi#chi) [GeV/c^{2}]; Events"    , 100 , 0. , 1500. );

   hsPt   = fs->make<TH1F>( "DmPt"   , "darkhiggs pt" , 100 , 0. , 500.  );
   hsMass = fs->make<TH1F>( "DmMass" , "darkhiggs mass"    , 100 , 0. , 1500. );

   metPt = fs->make<TH1F>( "metPt"   , "MET pt"      , 100 , 0.  , 500.  );

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

   //GenParticle Collection
   for( size_t i=0 ; i < genpart->size() ; ++i )
     {
       const reco::GenParticle &p = (*genpart)[i];
       const reco::Candidate *mom = p.mother();

       dm.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());

       pdg->Fill(p.pdgId());
       status->Fill(p.status());

       //dark matter pair
       //stable dark matter --> not in the case of darkhiggs
       if ( /*p.fromHardProcessFinalState() && */ abs(p.pdgId())== 52 && p.numberOfDaughters() == 0 ) {
           dmSystem+=dm;
       }

       // Z' mediator
       if ( p.status() == 22 && abs(p.pdgId())==55 )
       {
           ZpPt->Fill(p.pt(),unitxsec);
           ZpEta->Fill(p.rapidity(),unitxsec);
           ZpMass->Fill(p.mass(),unitxsec);
       }

       //Dark Higgs
       if ( abs(p.pdgId())==54 ) {
         hsPt->Fill(p.pt(),unitxsec);
         hsMass->Fill(p.mass(),unitxsec);
       }

     }

   //MET object is a property of event, thus it does not exist more then one
   metPt->Fill(met.pt())

   //acceptance
   countEvent++;
   for(Int_t i=0; i<7; i++)
   {
       if(dmSystem.Pt()<(double) 200+(i*50)) continue;
       countEventPass[i]+=1;
   }

   //Filling njet
   Njet->Fill(counter);
   //Filling dmSystem
   DmPt->Fill(dmSystem.Pt(),unitxsec);
   DmMass->Fill(dmSystem.M(),unitxsec);

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
