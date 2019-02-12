// -*- C++ -*-
//
// Package:    edmAna/LambdaAnayzer
// Class:      LambdaAnayzer
// 
/**\class LambdaAnayzer LambdaAnayzer.cc edmAna/LambdaAnayzer/plugins/LambdaAnayzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Siew Yan Hoh
//         Created:  Tue, 12 Feb 2019 20:23:44 GMT
//
//

#include "LambdaAnalyzer.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
LambdaAnalyzer::LambdaAnalyzer(const edm::ParameterSet& iConfig):
  GenPSet(iConfig.getParameter<edm::ParameterSet>("genSet")),
  JetPSet(iConfig.getParameter<edm::ParameterSet>("jetSet")),
  FatJetPSet(iConfig.getParameter<edm::ParameterSet>("fastJetSet")),
  HistFile(iConfig.getParameter<std::string>("histFile"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   //TFileDirectory allDir=fs->mkdir("All/");

   // Initialize Objects
   theGenAnalyzer      = new GenAnalyzer(GenPSet, consumesCollector());
   theJetAnalyzer      = new JetAnalyzer(JetPSet, consumesCollector());

   //  ---------- Plots Initialization ----------
   TFileDirectory allDir=fs->mkdir("All/");
   TFileDirectory genDir=fs->mkdir("Gen/");
   TFileDirectory eleDir=fs->mkdir("Electrons/");
   TFileDirectory muoDir=fs->mkdir("Muons/");
   TFileDirectory tauDir=fs->mkdir("Taus/");
   TFileDirectory phoDir=fs->mkdir("Photons/");
   TFileDirectory jetDir=fs->mkdir("Jets/");
   TFileDirectory kinDir=fs->mkdir("Kin/");
   
   //make TH1F
   std::vector<std::string> nLabels={"All", "Trigger", "Iso Lep #geq 2", "Z cand ", "Jets #geq 2", "Z mass ", "h mass ", "Top veto", "bJets #geq 1", "bJets #geq 2"};

   int nbins;
   float min, max;
   std::string name, title, opt;

   std::ifstream histFile(HistFile);
   if(!histFile.is_open()) {
     throw cms::Exception("DMbb Analyzer", HistFile + " file not found");
   }
   while(histFile >> name >> title >> nbins >> min >> max >> opt) {
     if(name.find('#')==std::string::npos) {
       while(title.find("~")!=std::string::npos) title=title.replace(title.find("~"), 1, " "); // Remove ~                                                          
       if(name.substr(0, 2)=="a_") Hist[name] = allDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max); //.substr(2)
       if(name.substr(0, 2)=="g_") Hist[name] = genDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="e_") Hist[name] = eleDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="m_") Hist[name] = muoDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="t_") Hist[name] = tauDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="p_") Hist[name] = phoDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="j_") Hist[name] = jetDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="k_") Hist[name] = kinDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);

     Hist[name]->Sumw2();
     Hist[name]->SetOption(opt.c_str());
     if(name=="a_nEvents" || name=="e_nEvents" || name=="m_nEvents") for(unsigned int i=0; i<nLabels.size(); i++) Hist[name]->GetXaxis()->SetBinLabel(i+1, nLabels[i].c_str());
     }
   }
   histFile.close();

   std::cout << "---------- STARTING ----------" << std::endl;
}

LambdaAnalyzer::~LambdaAnalyzer(){
  
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  std::cout << "---------- ENDING  ----------" << std::endl;
  delete theGenAnalyzer;
  delete theJetAnalyzer;
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
LambdaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Initialize event variable
   nJets=0;
   
   //GEN
   // Gen Weight
   std::map<int, float> GenWeight = theGenAnalyzer->FillWeightsMap(iEvent);

   // Lhe Particles 
   std::map<std::string, float> LheMap = theGenAnalyzer->FillLheMap(iEvent);

   // Gen Particle
   std::vector<reco::GenParticle> GenPVect = theGenAnalyzer->FillGenVector(iEvent);

   // Gen candidates
   reco::Candidate* theGenZ = theGenAnalyzer->FindGenParticle(GenPVect, 23);
   reco::Candidate* theGenW = theGenAnalyzer->FindGenParticle(GenPVect, 24);
   reco::Candidate* theGenTop     = theGenAnalyzer->FindGenParticle(GenPVect, 6);
   reco::Candidate* theGenAntiTop = theGenAnalyzer->FindGenParticle(GenPVect, -6);

   std::vector<int> LepIds = {11,13,15,-11,-13,-15};
   std::vector<int> NeuIds = {12,14,16,-12,-14,-16};
   std::vector<int> HadIds = {1,2,3,4,5,-1,-2,-3,-4,-5};

   reco::GenParticle* theGenLep = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, LepIds);
   reco::GenParticle* theGenNeu = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, NeuIds);
   reco::GenParticle* theGenHad = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, HadIds);


   // Jet

   std::vector<pat::Jet> JetsVect = theJetAnalyzer->FillJetVector(iEvent);
   nJets = JetsVect.size();
   pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);


}


// ------------ method called once each job just before starting event loop  ------------
void 
LambdaAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LambdaAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LambdaAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LambdaAnalyzer);
