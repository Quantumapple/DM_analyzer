// -*- C++ -*-
//
// Package:    DM_analyzer/LambdaAnalyzer
// Class:      SignalScanProducer
// 
/**\class SignalScanProducer SignalScanProducer.cc DM_analyzer/LambdaAnalyzer/plugins/SignalScanProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jongho Lee
//         Created:  Wed, 09 Sep 2020 17:13:29 GMT
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/GetterOfProducts.h"
#include "FWCore/Framework/interface/ProcessMatch.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
//#include "FWCore/Utilities/interface/StreamID.h"

// STL include files
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

enum class signal_type {
    None=0,
    darkhiggs=1
};

//
// class declaration
//

class SignalScanProducer : public edm::stream::EDProducer<> {
   public:
      explicit SignalScanProducer(const edm::ParameterSet&);
      ~SignalScanProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      void reset();
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      void getLHEComment(const LHEEventProduct& lhe);
      void getComment(const GenLumiInfoHeader& gen);

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::GetterOfProducts<LHEEventProduct> getterOfProducts_;
      edm::EDGetTokenT<GenLumiInfoHeader> genLumiHeaderToken_;
      bool shouldScan_, debug_, isLHE_;
      signal_type type_;
      std::vector<double> signalParameters_;
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
SignalScanProducer::SignalScanProducer(const edm::ParameterSet& iConfig) :
    getterOfProducts_(edm::ProcessMatch("*"), this),
    genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator"))),
    shouldScan_(true),
    debug_(iConfig.getParameter<bool>("debug")),
    isLHE_(iConfig.getParameter<bool>("isLHE"))
{
    std::string stype(iConfig.getParameter<std::string>("signalType"));
    if(stype=="None") { type_ = signal_type::None; shouldScan_ = false; }
    else if(stype=="darkhiggs") type_ = signal_type::darkhiggs;
    callWhenNewProductsRegistered(getterOfProducts_);
    produces<std::vector<double>>("SignalParameters");
}


SignalScanProducer::~SignalScanProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SignalScanProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    if(isLHE_ && shouldScan_){
        reset();
        if(debug_) edm::LogInfo("TreeMaker") << "SignalScanProducer: checking LHEEventProduct";
        std::vector<edm::Handle<LHEEventProduct> > handles;
        getterOfProducts_.fillHandles(iEvent, handles);

        if(!handles.empty()){
            edm::Handle<LHEEventProduct> product = handles[0];
            if(type_==signal_type::darkhiggs) getLHEComment(*product);
        }
    }

    auto signalParameters = std::make_unique<std::vector<double>>(signalParameters_);
    iEvent.put(std::move(signalParameters), "SignalParameters");
}

/*
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SignalScanProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SignalScanProducer::endStream() {
}
*/

// ------------ method called when starting to processes a run  ------------
/*
void
SignalScanProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SignalScanProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------

void
//SignalScanProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
SignalScanProducer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&)
{
    if(!isLHE_ && shouldScan_){
        reset();
        if(debug_) edm::LogInfo("TreeMaker") << "SignalScanProducer: checking GenLumiInfoHeader";
        edm::Handle<GenLumiInfoHeader> gen_header;
        iLumi.getByToken(genLumiHeaderToken_, gen_header);
        if(type_==signal_type::darkhiggs) getComment(*gen_header);
    }
}

 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SignalScanProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SignalScanProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//reset output vars
void SignalScanProducer::reset(){
    signalParameters_.clear();
}

void SignalScanProducer::getComment(const GenLumiInfoHeader& gen){
    std::string model = gen.configDescription();
    if(debug_) edm::LogInfo("TreeMaker") << model;
    std::cout << model << std::endl;
    //std::vector<std::string> fields;
}


//define this as a plug-in
DEFINE_FWK_MODULE(SignalScanProducer);
