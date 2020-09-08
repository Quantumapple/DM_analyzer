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
  ElectronPSet(iConfig.getParameter<edm::ParameterSet>("electronSet")),
  MuonPSet(iConfig.getParameter<edm::ParameterSet>("muonSet")),
  HistFile(iConfig.getParameter<std::string>("histFile"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   //TFileDirectory allDir=fs->mkdir("All/");

   // Initialize Objects
   theGenAnalyzer      = new GenAnalyzer(GenPSet, consumesCollector());
   theJetAnalyzer      = new JetAnalyzer(JetPSet, consumesCollector());
   theElectronAnalyzer = new ElectronAnalyzer(ElectronPSet, consumesCollector());
   theMuonAnalyzer     = new MuonAnalyzer(MuonPSet, consumesCollector());

   //  ---------- Plots Initialization ----------
   TFileDirectory allDir=fs->mkdir("All/");
   TFileDirectory genDir=fs->mkdir("Gen/");
   TFileDirectory recoDir=fs->mkdir("Reco/");
   
   //make TH1F
   //std::vector<std::string> nLabels={"All", "Trigger", "Iso Lep #geq 2", "Z cand ", "Jets #geq 2", "Z mass ", "h mass ", "Top veto", "bJets #geq 1", "bJets #geq 2"};
   //std::vector<std::string> nLabels={"All", "MET > 200", "MET > 250", "Z cand ", "Jets #geq 2", "Z mass ", "h mass ", "Top veto", "bJets #geq 1", "bJets #geq 2"};
   std::vector<std::string> nLabels={"All", "MET > 200", "MET > 250", "Z' p_{T} > 250", "Jets #geq 2", "Z mass ", "h mass ", "Top veto", "bJets #geq 1", "bJets #geq 2"};

   int nbins;
   float min, max;
   std::string name, title, opt;

   std::ifstream histFile(HistFile);
   if(!histFile.is_open()) {
     throw cms::Exception("Analyzer", HistFile + " file not found");
   }
   while(histFile >> name >> title >> nbins >> min >> max >> opt) {
       //std::cout << name << ", " << title << ", " << nbins << ", " << min << ", " << ", " << max << ", " << opt << std::endl;

     if(name.find('#')==std::string::npos) {
       while(title.find("~")!=std::string::npos) title=title.replace(title.find("~"), 1, " "); // Remove ~
       if(name.substr(0, 2)=="a_") Hist[name] = allDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="g_") Hist[name] = genDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
       if(name.substr(0, 2)=="r_") Hist[name] = recoDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max); 

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
  delete theElectronAnalyzer;
  delete theMuonAnalyzer;
  
}


//
// member functions
//

float LambdaAnalyzer::deltaPhi(float phi1, float phi2){
  float PHI = fabs(phi1-phi2);
  if (PHI<=3.14159265)
    return PHI;
  else
    return 2*3.14159265-PHI; 
}

std::vector<const reco::GenParticle*> LambdaAnalyzer::IndexByPtGen(std::vector<const reco::GenParticle*> vector)
{
  compgen comparatorgen;

  std::sort (vector.begin() , vector.end() , comparatorgen);
  return vector;
}

std::vector<pat::Jet> LambdaAnalyzer::IndexByPtPat(std::vector<pat::Jet> vector)
{
  comppat comparatorpat;

  std::sort (vector.begin() , vector.end() , comparatorpat);
  return vector;
}

float LambdaAnalyzer::deltaR(float phi1, float eta1, float phi2, float eta2){
  return sqrt(pow((eta2-eta1),2)+pow(deltaPhi(phi1,phi2),2));
}

// ------------ method called for each event  ------------
void
LambdaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
 
   //Initialize event variable
   nJets=0; EventWeight=1.;
  
   Hist["a_nEvents"]->Fill(1.,EventWeight);
   //std::cout << "Event filled" << std::endl;
   //std::cout << std::endl;

   //GEN
   // Gen Weight
   //std::map<int, float> GenWeight = theGenAnalyzer->FillWeightsMap(iEvent);

   // Lhe Particles 
   //std::map<std::string, float> LheMap = theGenAnalyzer->FillLheMap(iEvent);

   // Gen Particle
   std::vector<reco::GenParticle> GenPVect = theGenAnalyzer->FillGenVector(iEvent);

   // Gen candidates
   //for(unsigned int i = 0; i < GenPVect.size(); i++)
   //{
   //    if(GenPVect[i].pdgId() > 50 ) std::cout << "Pdg ID: " << GenPVect[i].pdgId() << std::endl;
   //}
   //std::cout << std::endl;

   std::vector<reco::GenParticle> theDM = theGenAnalyzer->SelectGenVector(GenPVect, 52);
   std::vector<reco::GenParticle> thehs = theGenAnalyzer->SelectGenVector(GenPVect, 54);
   std::vector<reco::GenParticle> theZp = theGenAnalyzer->SelectGenVector(GenPVect, 55);
   //std::vector<reco::GenParticle> theGenW = theGenAnalyzer->SelectGenVector(GenPVect, 24);
   //std::vector<reco::GenParticle> theGenZ = theGenAnalyzer->SelectGenVector(GenPVect, 23);
   
   /*
   float Zp1pt = 0.;
   if( theZp.size() != 0 ) { 
       Zp1pt = theZp[0].pt();
       //Hist["g_Zpmass"]->Fill(theZp[0].mass(), EventWeight);
   }
   if( Zp1pt > 250. ) 
   {
       Hist["a_nEvents"]->Fill(4.,EventWeight);
       //Hist["g_hsmass"]->Fill(thehs[0].mass(), EventWeight);
       //Hist["g_hspt"]->Fill(thehs[0].pt(), EventWeight);
       //Hist["g_WPt"]->Fill(theGenW[0].pt());
       //Hist["g_ZPt"]->Fill(theGenZ[0].pt());
   }
   */

   //std::cout << theZp.size() << ", " << thehs.size() << ", " << theDM.size() << std::endl;
   if( theZp.size() > 0 && thehs.size() > 0 && theDM.size() > 0 ) std::cout << theZp[0].mass() << ", " << thehs[0].mass() << ", " << theDM[0].mass() << std::endl;
   //if( theZp.size() > 0 ) std::cout << theZp[0].mass() << std::endl;
   //if( thehs.size() > 0 ) std::cout << thehs[0].mass() << std::endl;
   //if( theDM.size() > 0 ) std::cout << theDM[0].mass() << std::endl;
   
   //for(unsigned int i = 0; i < theZp.size(); i++)
   //{
   //    if(thehs[0].mass() != 70) continue;
   //    Hist["g_Zpmass"]->Fill(theZp[0].mass(), EventWeight);
   //}


   //for(unsigned int i = 0; i < thehs.size(); i++)
   //{
   //    if(thehs[0].mass() != 70) continue;
   //    Hist["g_hsmass"]->Fill(thehs[0].mass(),EventWeight);
   //}
   //
   //for(unsigned int i = 0; i < theDM.size(); i++)
   //{
   //    if(thehs[0].mass() != 70) continue;
   //    Hist["g_DMmass"]->Fill(theDM[0].mass(),EventWeight);
   //    //Hist["g_DMmass"]->Fill(theDM[1].mass(),EventWeight);
   //}
   
   // Electron
   std::vector<pat::Electron> ElecVect = theElectronAnalyzer->FillElectronVector(iEvent);

   // Muon
   std::vector<pat::Muon> MuonVect = theMuonAnalyzer->FillMuonVector(iEvent);

   // Jet
   std::vector<pat::Jet> JetsVect = theJetAnalyzer->FillJetVector(iEvent);
   Hist["r_nJetb"]->Fill(JetsVect.size(),EventWeight);
   theJetAnalyzer->CleanJetsFromMuons(JetsVect, MuonVect, 0.4);
   theJetAnalyzer->CleanJetsFromElectrons(JetsVect, ElecVect, 0.4);

   // GenJet
   std::vector<reco::GenJet> GenJetsVect = theJetAnalyzer->FillGenJetVector(iEvent);
   
   //MET
   pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);   
   
   /*
   Hist["a_met"]->Fill(MET.pt(), EventWeight);
   // Save MET phi
   METphi = MET.phi();
   METeta = MET.eta();
   //std::cout << "MET filled" << std::endl;
   //std::cout << std::endl;
   
   // Fill number of events when MET > 200 GeV
   if ( MET.pt() > 200. ) Hist["a_nEvents"]->Fill(2.,EventWeight);
   if ( MET.pt() > 250. ) Hist["a_nEvents"]->Fill(3.,EventWeight);

   //JET MC Truth
   //The pruned genParticles are the ones pointed by the MC matching of the high level patObjectes (e.g. pat::Electron::genParticle()) 
   //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#MC_Truth

   //Matching by R may not work reliably in dense environments, such as jets. For studies needing high quality matching of reconstructed tracks with true tracks, it is possible to base the matching either on the number of hits that they share in common, or on a comparison of the 5 helix parameters describing the track. 

   std::vector<const reco::GenParticle*> HardGenJetsVect;
   std::vector<pat::Jet> HardJetsVect;

   for(unsigned int i = 0; i < JetsVect.size(); i++){
     //std::cout<<"jet pt = "<<JetsVect[i].pt()<<std::endl;
     pat::Jet j = JetsVect[i];
     if (j.genParton() == NULL)
       continue;
     
     //std::cout<<"genParton pdgid = "<< ( j.genParton() )->pdgId() <<std::endl;
     //std::cout<<"Robust Flags"<<std::endl;
     //std::cout<<"isPromptFinalState() = "<< ( j.genParton() )->isPromptFinalState() <<std::endl;
     //std::cout<<"Less robust Flags"<<std::endl;
     //std::cout<<"isHardProcess() = "<< ( j.genParton() )->isHardProcess() <<std::endl;
     //std::cout<<"fromHardProcessFinalState() = "<< ( j.genParton() )->fromHardProcessFinalState() <<std::endl;
     //std::cout<<"fromHardProcessDecayed() = "<< ( j.genParton() )->fromHardProcessDecayed() <<std::endl;
     //std::cout<<"isLastCopy() = "<< ( j.genParton() )->isLastCopy() <<std::endl;

     //std::cout<<"hadronFlavour = "<<(JetsVect[i].hadronFlavour())<<std::endl;
     //std::cout<<"partonFlavour = "<<(JetsVect[i].partonFlavour())<<std::endl;
     

     if (j.genParton()->isHardProcess()){
       HardGenJetsVect.push_back( j.genParton() );
       HardJetsVect.push_back(j);
     }
   }

   //std::cout << "Before Sorting" << std::endl; 
   //std::cout<<"Total number of match reco jet = "<< HardJetsVect.size() <<std::endl;
   for(unsigned int i = 0; i < HardJetsVect.size(); i++){
     pat::Jet j = HardJetsVect[i];
     //std::cout<<"Jet number "<<i<<" with j.pt() = "<<j.pt()<<std::endl;
   }
   //Sort vector
   JetsVect.clear();
   std::vector<const reco::GenParticle*> JetsMCmatch;
   JetsVect = IndexByPtPat(HardJetsVect);
   JetsMCmatch = IndexByPtGen(HardGenJetsVect);
   //std::cout << "After Sorting" << std::endl;
   //std::cout<<"Total number of match reco jet = "<< JetsVect.size() <<std::endl;
   for(unsigned int i = 0; i < JetsVect.size(); i++){
     pat::Jet j = JetsVect[i];
     //std::cout<<"Jet number "<<i<<" with j.pt() = "<<j.pt()<<std::endl;
   }

   //TLorentzVector v1;
   //TLorentzVector v2;
   //TLorentzVector v3;

   //Filling
  
   //GenJet
   //if( Zp1pt > 250. ) Hist["g_nJet"]->Fill(JetsMCmatch.size(), EventWeight);
   Hist["g_nJet"]->Fill(JetsMCmatch.size(), EventWeight);
   for(unsigned int i = 0; i < JetsMCmatch.size(); i++){
     if (i>2) break;
     //if( Zp1pt >  250. ) Hist[("g_Jet"+std::to_string(i+1)+"pt").c_str()]->Fill(JetsMCmatch[i]->pt(), EventWeight);
     //if( Zp1pt >  250. ) Hist[("g_Jet"+std::to_string(i+1)+"eta").c_str()]->Fill(JetsMCmatch[i]->eta(), EventWeight);
     Hist[("g_Jet"+std::to_string(i+1)+"pt").c_str()]->Fill(JetsMCmatch[i]->pt(), EventWeight);
     Hist[("g_Jet"+std::to_string(i+1)+"eta").c_str()]->Fill(JetsMCmatch[i]->eta(), EventWeight);
     if (JetsMCmatch.size() >= 2 )
     {
         //v1.SetPtEtaPhiE(JetsMCmatch[0]->pt(), JetsMCmatch[0]->eta(), JetsMCmatch[0]->phi(). JetsMCmatch[0]->energy());
         //v2.SetPtEtaPhiE(JetsMCmatch[1]->pt(), JetsMCmatch[1]->eta(), JetsMCmatch[1]->phi(). JetsMCmatch[1]->energy());
         Hist["g_J1J2dPhi"]->Fill(deltaPhi(JetsMCmatch[0]->phi(), JetsMCmatch[1]->phi()), EventWeight);
         Hist["g_J1J2dEta"]->Fill(JetsMCmatch[0]->eta() - JetsMCmatch[1]->eta(), EventWeight);
         Hist["g_J1J2dR"]->Fill(deltaR(JetsMCmatch[0]->phi(), JetsMCmatch[0]->eta(), JetsMCmatch[1]->phi(), JetsMCmatch[1]->eta()), EventWeight);
         //if( Zp1pt > 250. ) Hist["g_J1J2dPhi"]->Fill(deltaPhi(JetsMCmatch[0]->phi(), JetsMCmatch[1]->phi()), EventWeight);
         //if( Zp1pt > 250. ) Hist["g_J1J2dEta"]->Fill(JetsMCmatch[0]->eta() - JetsMCmatch[1]->eta(), EventWeight);
         //if( Zp1pt > 250. ) Hist["g_J1J2dR"]->Fill(deltaR(JetsMCmatch[0]->phi(), JetsMCmatch[0]->eta(), JetsMCmatch[1]->phi(), JetsMCmatch[1]->eta()), EventWeight);
     }
    }

   //v3 = v1 + v2;
   //if( Zp1pt > 250. ) Hist["g_dijetM"]->Fill(v3.M());
   
   //std::cout << "Genjet finished" << std::endl;
   //std::cout << std::endl;
  
   nJets=JetsVect.size();
   //RecoJet
   Hist["r_nJet"]->Fill(nJets, EventWeight);
   for(unsigned int i = 0; i < JetsVect.size(); i++){ 
     if (i>2) break;
     Hist[("r_Jet"+std::to_string(i+1)+"pt").c_str()]->Fill(JetsVect[i].pt(), EventWeight); 
     Hist[("r_Jet"+std::to_string(i+1)+"eta").c_str()]->Fill(JetsVect[i].eta(), EventWeight); 
     RecoJphi = JetsVect[0].phi();
     RecoJeta = JetsVect[0].eta();
     if (JetsVect.size() >= 2 )
     {
         Hist["r_J1J2dPhi"]->Fill(deltaPhi(JetsVect[0].phi(), JetsVect[1].phi()), EventWeight);
         Hist["r_J1J2dEta"]->Fill(JetsVect[0].eta() - JetsVect[1].eta(), EventWeight);
         Hist["r_J1J2dR"]->Fill(deltaR(JetsVect[0].phi(), JetsVect[0].eta(), JetsVect[1].phi(), JetsVect[1].eta()), EventWeight);
     }
   }
   //std::cout << "Recojet finished" << std::endl;
   //std::cout << std::endl;

   // Calculate angle difference (delta phi) of MET and 
   Hist["r_dPhi"]->Fill(deltaPhi(RecoJphi, METphi), EventWeight);
   Hist["r_dEta"]->Fill(RecoJeta - METeta, EventWeight);
   Hist["r_dR"]->Fill(deltaR(METphi, METeta, RecoJphi, RecoJeta), EventWeight);
   */

   tree->Fill();
   
   //std::cout << "Filling finished" << std::endl;
   //std::cout << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
LambdaAnalyzer::beginJob()
{
  tree=fs->make<TTree>("tree", "tree");
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
