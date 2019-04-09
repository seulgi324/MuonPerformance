// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include<map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// GEM
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
// CSC
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
// GEMCSC
#include "DataFormats/GEMRecHit/interface/GEMCSCSegmentCollection.h"
// RPC
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
// Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
// L1
#include "DataFormats/L1TMuon/interface/L1MuBMTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimMuon/MCTruth/interface/MuonToSimAssociatorByHits.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/Framework/interface/HCTypeTag.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

/// from MuonSimAnalyser
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h" 
#include "TDatabasePDG.h"
///

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//

using namespace std;
using namespace edm;

class MuonTrackAnalyser : public edm::EDAnalyzer {
public:
  explicit MuonTrackAnalyser(const edm::ParameterSet&);
  ~MuonTrackAnalyser();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  void initMuonValue();
  void initValue();

  // ----------member data ---------------------------
  edm::EDGetTokenT<CSCSegmentCollection>       cscSegments_;

  edm::EDGetTokenT<GEMSegmentCollection>       gemSegments_;
  
  edm::EDGetTokenT<GEMCSCSegmentCollection>       gemcscSegments_;

  edm::EDGetTokenT<TrackingParticleCollection> simToken_;

  edm::EDGetTokenT<edm::View<reco::Muon> >     muonToken_;
  
  edm::EDGetTokenT<reco::GenParticleCollection>     genParticlesToken_;

  const reco::VertexCollection* vertexes_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  
  edm::Service<TFileService> fs;

  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nCSCSegments, b_nGEMSegments, b_nGEMCSCSegments;
  float b_mu_pt, b_mu_eta;
  int b_mu_pdgId;

  /* CSC */
  TTree *t_CSC_seg;
  int   b_CSC_Seg_region, b_CSC_Seg_station, b_CSC_Seg_ring, b_CSC_Seg_chamber, b_CSC_Seg_layer, b_CSC_Seg_nRecHits;
  float b_CSC_Seg_eta,    b_CSC_Seg_phi;
  int b_CSC_Rec_cls;

  /* GEM */
  TTree *t_GEM_seg;
  int   b_GEM_Seg_region, b_GEM_Seg_station, b_GEM_Seg_ring, b_GEM_Seg_chamber, b_GEM_Seg_layer, b_GEM_Seg_nRecHits;
  float b_GEM_Seg_eta,    b_GEM_Seg_phi;
  int b_GEM_Rec_cls;
  
  /* GEMCSC */
  TTree *t_GEMCSC_seg;
  int   b_GEMCSC_Seg_region, b_GEMCSC_Seg_station, b_GEMCSC_Seg_ring, b_GEMCSC_Seg_chamber, b_GEMCSC_Seg_layer, b_GEMCSC_Seg_nRecHits, b_GEMCSC_Seg_nCSCRecHits, b_GEMCSC_Seg_nGEMRecHits;
  float b_GEMCSC_Seg_eta,    b_GEMCSC_Seg_phi;
  int b_GEMCSC_Rec_cls;
 
  /* Muon */
  TTree *t_genMuon;
  float b_genMuon_pt,      b_genMuon_eta,      b_genMuon_phi;
  int   b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose, b_genMuon_isTracker, b_genMuon_isGlobal, b_genMuon_isME0, b_genMuon_nMatchedStationLayer;

  TTree *t_Muon;
  float b_muon_pt,      b_muon_eta,      b_muon_phi;
  int   b_muon_isTight, b_muon_isMedium, b_muon_isLoose, b_muon_isTracker, b_muon_isGlobal, b_muon_isME0, b_muon_nMatchedStationLayer;
};

MuonTrackAnalyser::MuonTrackAnalyser(const edm::ParameterSet& iConfig)
{ 
  cscSegments_  = consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"));

  gemSegments_  = consumes<GEMSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemSegments"));
  
  gemcscSegments_  = consumes<GEMCSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemcscSegments"));

  simToken_     = consumes<TrackingParticleCollection>(iConfig.getParameter<InputTag>("simLabel"));
  
  muonToken_    = consumes<edm::View<reco::Muon> >(iConfig.getParameter<InputTag>("muonLabel"));
  
  genParticlesToken_  = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparts"));
  
  vtxToken_     = consumes<vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("primaryVertex"));
  
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nCSCSegments",  &b_nCSCSegments,  "nCSCSegments/I");
  t_event->Branch("nGEMSegments",  &b_nGEMSegments,  "nGEMSegments/I");
  t_event->Branch("nGEMCSCSegments",  &b_nGEMCSCSegments,  "nGEMCSCSegments/I");
  t_event->Branch("mu_pt",   &b_mu_pt,   "mu_pt/F");
  t_event->Branch("mu_eta",   &b_mu_eta,   "mu_eta/F");
  t_event->Branch("mu_pdgId",   &b_mu_pdgId,   "mu_pdgId/I");

  /*CSC*/
  t_CSC_seg = fs->make<TTree>("CSC_Segment", "CSC_Segment");
  t_CSC_seg->Branch("Seg_eta",      &b_CSC_Seg_eta,      "Seg_eta/F");
  t_CSC_seg->Branch("Seg_phi",      &b_CSC_Seg_phi,      "Seg_phi/F");
  t_CSC_seg->Branch("Seg_nRecHits", &b_CSC_Seg_nRecHits, "Seg_nRecHits/I");
  t_CSC_seg->Branch("Seg_region",   &b_CSC_Seg_region,   "Seg_region/I");
  t_CSC_seg->Branch("Seg_station",  &b_CSC_Seg_station,  "Seg_station/I");
  t_CSC_seg->Branch("Seg_ring",     &b_CSC_Seg_ring,     "Seg_ring/I");
  t_CSC_seg->Branch("Seg_chamber",  &b_CSC_Seg_chamber,  "Seg_chamber/I");
  t_CSC_seg->Branch("Seg_layer",    &b_CSC_Seg_layer,    "Seg_layer/I");
  t_CSC_seg->Branch("Rec_cls",    &b_CSC_Rec_cls,    "Rec_cls/I");


  /*GEM*/
  t_GEM_seg = fs->make<TTree>("GEM_Segment", "GEM_Segment");
  t_GEM_seg->Branch("Seg_eta",      &b_GEM_Seg_eta,      "Seg_eta/F");
  t_GEM_seg->Branch("Seg_phi",      &b_GEM_Seg_phi,      "Seg_phi/F");
  t_GEM_seg->Branch("Seg_nRecHits", &b_GEM_Seg_nRecHits, "Seg_nRecHits/I");
  t_GEM_seg->Branch("Seg_region",   &b_GEM_Seg_region,   "Seg_region/I");
  t_GEM_seg->Branch("Seg_station",  &b_GEM_Seg_station,  "Seg_station/I");
  t_GEM_seg->Branch("Seg_ring",     &b_GEM_Seg_ring,     "Seg_ring/I");
  t_GEM_seg->Branch("Seg_chamber",  &b_GEM_Seg_chamber,  "Seg_chamber/I");
  t_GEM_seg->Branch("Seg_layer",    &b_GEM_Seg_layer,    "Seg_layer/I");
  t_GEM_seg->Branch("Rec_cls",    &b_GEM_Rec_cls,    "Rec_cls/I");
  
  /*GEMCSC*/
  t_GEMCSC_seg = fs->make<TTree>("GEMCSC_Segment", "GEMCSC_Segment");
  t_GEMCSC_seg->Branch("Seg_eta",      &b_GEMCSC_Seg_eta,      "Seg_eta/F");
  t_GEMCSC_seg->Branch("Seg_phi",      &b_GEMCSC_Seg_phi,      "Seg_phi/F");
  t_GEMCSC_seg->Branch("Seg_nRecHits", &b_GEMCSC_Seg_nRecHits, "Seg_nRecHits/I");
  t_GEMCSC_seg->Branch("Seg_nCSCRecHits", &b_GEMCSC_Seg_nCSCRecHits, "Seg_nCSCRecHits/I");
  t_GEMCSC_seg->Branch("Seg_nGEMRecHits", &b_GEMCSC_Seg_nGEMRecHits, "Seg_nGEMRecHits/I");
  t_GEMCSC_seg->Branch("Seg_region",   &b_GEMCSC_Seg_region,   "Seg_region/I");
  t_GEMCSC_seg->Branch("Seg_station",  &b_GEMCSC_Seg_station,  "Seg_station/I");
  t_GEMCSC_seg->Branch("Seg_ring",     &b_GEMCSC_Seg_ring,     "Seg_ring/I");
  t_GEMCSC_seg->Branch("Seg_chamber",  &b_GEMCSC_Seg_chamber,  "Seg_chamber/I");
  t_GEMCSC_seg->Branch("Seg_layer",    &b_GEMCSC_Seg_layer,    "Seg_layer/I");
  t_GEMCSC_seg->Branch("Rec_cls",    &b_GEMCSC_Rec_cls,    "Rec_cls/I");

  /*Muon*/
  t_genMuon = fs->make<TTree>("gen_Muon", "gen_Muon");
  t_genMuon->Branch("pt",                   &b_genMuon_pt,                   "pt/F");
  t_genMuon->Branch("eta",                  &b_genMuon_eta,                  "eta/F");
  t_genMuon->Branch("phi",                  &b_genMuon_phi,                  "phi/F");
  t_genMuon->Branch("nMatchedStationLayer", &b_genMuon_nMatchedStationLayer, "nMatchedStationLayer/I");
  t_genMuon->Branch("isTracker",            &b_genMuon_isTracker,            "isTracker/I");
  t_genMuon->Branch("isGlobal",             &b_genMuon_isGlobal,             "isGlobal/I");
  t_genMuon->Branch("isME0",                &b_genMuon_isME0,                "isME0/I");
  t_genMuon->Branch("isTight",              &b_genMuon_isTight,              "isTight/I");
  t_genMuon->Branch("isMedium",             &b_genMuon_isMedium,             "isMedium/I");
  t_genMuon->Branch("isLoose",              &b_genMuon_isLoose,              "isLoose/I");

  t_Muon = fs->make<TTree>("Muon", "Muon");
  t_Muon->Branch("pt",                   &b_muon_pt,                   "pt/F");
  t_Muon->Branch("eta",                  &b_muon_eta,                  "eta/F");
  t_Muon->Branch("phi",                  &b_muon_phi,                  "phi/F");
  t_Muon->Branch("nMatchedStationLayer", &b_muon_nMatchedStationLayer, "nMatchedStationLayer/I");
  t_Muon->Branch("isTracker",            &b_muon_isTracker,            "isTracker/I");
  t_Muon->Branch("isGlobal",             &b_muon_isGlobal,             "isGlobal/I");
  t_Muon->Branch("isME0",                &b_muon_isME0,                "isME0/I");
  t_Muon->Branch("isTight",              &b_muon_isTight,              "isTight/I");
  t_Muon->Branch("isMedium",             &b_muon_isMedium,             "isMedium/I");
  t_Muon->Branch("isLoose",              &b_muon_isLoose,              "isLoose/I");
}

MuonTrackAnalyser::~MuonTrackAnalyser()
{
}

void
MuonTrackAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* CSC Geometry */
  edm::ESHandle<CSCGeometry> hCSCGeom;
  iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  const CSCGeometry* CSCGeometry_ = &*hCSCGeom;

  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegments_, cscSegments);

  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<GEMSegmentCollection> gemSegments;
  iEvent.getByToken(gemSegments_, gemSegments);

  /* GEMCSC Geometry */
  edm::Handle<GEMCSCSegmentCollection> gemcscSegments;
  iEvent.getByToken(gemcscSegments_, gemcscSegments);

  /* Muon */
  Handle<TrackingParticleCollection> simHandle;
  iEvent.getByToken(simToken_, simHandle);

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices); 

  vertexes_ = vertices.product();
  reco::Vertex pv0 = vertexes_->at(0);

  initValue();

  /* CSC */
  for (auto ch : CSCGeometry_->chambers()) {
    // CSC seg
    CSCDetId cId = ch->id();

    auto segsRange = cscSegments->get(cId);
    auto cscSeg = segsRange.first;
    if (cscSeg == segsRange.second) continue;
    //std::cout << cscSeg->isValid() << std::endl;
    for (auto seg = cscSeg; seg != segsRange.second; ++seg) { 
      auto segLd = seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_CSC_Seg_eta      = gp.eta();
      b_CSC_Seg_phi      = gp.phi();
      b_CSC_Seg_nRecHits = seg->nRecHits();
      b_CSC_Seg_region   = cId.endcap();
      b_CSC_Seg_station  = cId.station();
      b_CSC_Seg_ring     = cId.ring();
      b_CSC_Seg_chamber  = cId.chamber();
      b_CSC_Seg_layer    = cId.layer();
      
      b_nCSCSegments++;
      t_CSC_seg->Fill();
    }
  }
  /* GEM */
  for (auto ch : GEMGeometry_->chambers()) {
    // GEM
    GEMDetId cId = ch->id();
    auto segsRange = gemSegments->get(cId);
    auto gemSeg = segsRange.first;
    for (auto seg = gemSeg; seg != segsRange.second; ++seg) {
      auto segLd = seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_GEM_Seg_eta      = gp.eta();
      b_GEM_Seg_phi      = gp.phi();
      b_GEM_Seg_nRecHits = seg->nRecHits();
      b_GEM_Seg_region   = cId.region();
      b_GEM_Seg_station  = cId.station();
      b_GEM_Seg_ring     = cId.ring();
      b_GEM_Seg_chamber  = cId.chamber();
      b_GEM_Seg_layer    = cId.layer();
      
      b_nGEMSegments++;
      t_GEM_seg->Fill();
    }
  }

  /* GEMCSC */

  for (auto ch : CSCGeometry_->chambers()) {
    // GEMCSC
    CSCDetId cId = ch->id();
    if (!gemcscSegments.isValid()) continue;
    auto segsRange = gemcscSegments->get(cId);
    std::cout <<"he3e1"<<std::endl;
    std::cout << typeid(segsRange.first).name() << ", " << typeid(segsRange.second).name() << std::endl;
    auto gemcscSeg = segsRange.first;
    //if (gemcscSeg == NULL) continue;
    //std::cout << gemcscSeg->isValid() << std::endl;
    std::cout <<"heeq"<<std::endl;
    for (auto seg = gemcscSeg; seg != segsRange.second; ++seg) {
      std::cout <<"heeqc"<<std::endl;
      auto segLd = seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_GEMCSC_Seg_eta      = gp.eta();
      b_GEMCSC_Seg_phi      = gp.phi();
      b_GEMCSC_Seg_nRecHits = seg->nRecHits();
      b_GEMCSC_Seg_nCSCRecHits = seg->cscRecHits().size();
      b_GEMCSC_Seg_nGEMRecHits = seg->gemRecHits().size();
      b_GEMCSC_Seg_region   = cId.endcap();
      b_GEMCSC_Seg_station  = cId.station();
      b_GEMCSC_Seg_ring     = cId.ring();
      b_GEMCSC_Seg_chamber  = cId.chamber();
      b_GEMCSC_Seg_layer    = cId.layer();
      b_nGEMCSCSegments++;
      t_GEMCSC_seg->Fill();
      /*
      for (auto recHit : seg->cscRecHits()) {
        b_GEMCSC_Rec_cls = recHit.nStrips();
        //std::cout << b_GEMCSC_Rec_cls << endl; 
      }
      for (auto recHit : seg->gemRecHits()) {
        b_GEMCSC_Rec_cls = recHit.clusterSize();
        //std::cout << b_GEMCSC_Rec_cls << endl; 
      }*/
      if (b_GEMCSC_Seg_station==1 && b_GEMCSC_Seg_ring == 1){
          
        ++b_nGEMCSCSegments;
        t_GEMCSC_seg->Fill();
      }
    }
  }

  /* get muon information*/
  reco::GenParticleCollection::const_iterator mu = genParticles->begin();
  b_mu_pt = mu->pt(); 
  b_mu_eta = mu->eta();
  b_mu_pdgId = mu->pdgId();

  std::cout << "pdgId " << mu->pdgId() << std::endl;

  /* gen Muon */
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    b_genMuon_pt  = simTP->pt();
    b_genMuon_eta = simTP->eta();
    b_genMuon_phi = simTP->phi();
    //std::cout << simTP->genParticles().empty() << std::endl;
    bool isSignalMuon = abs(simTP->pdgId())==13 && !simTP->genParticles().empty() && (simTP->eventId().event() == 0) && (simTP->eventId().bunchCrossing() == 0);
    if (!isSignalMuon) continue;
    t_genMuon->Fill();
  }
  /* reco Muon */
  for(size_t i = 0; i < muonHandle->size(); ++i) {
    initMuonValue();
    edm::RefToBase<reco::Muon> muRef = muonHandle->refAt(i);
    const reco::Muon* mu = muRef.get();
    b_muon_isTight              = muon::isTightMuon(*mu, pv0);
    b_muon_isMedium             = muon::isMediumMuon(*mu);
    b_muon_isLoose              = muon::isLooseMuon(*mu);
    b_muon_pt                   = mu->pt();
    b_muon_eta                  = mu->eta();
    b_muon_phi                  = mu->phi();
    b_muon_isGlobal             = mu->isGlobalMuon();
    b_muon_isTracker            = mu->isTrackerMuon();
    b_muon_isME0                = mu->isME0Muon();
    b_muon_nMatchedStationLayer = mu->numberOfMatchedStations();
    t_Muon->Fill();
  } 

  t_event->Fill();
}

void MuonTrackAnalyser::beginJob(){}
void MuonTrackAnalyser::endJob(){}

void MuonTrackAnalyser::beginRun(const edm::Run& run, const edm::EventSetup& iSetup){//(Run const& run, EventSetup const&){
}
void MuonTrackAnalyser::endRun(Run const&, EventSetup const&){}

 void MuonTrackAnalyser::initMuonValue() {
  b_genMuon_pt  = -9; b_genMuon_eta = -9; b_genMuon_phi = -9;

  b_muon_isTight = -1; b_muon_isMedium = -1; b_muon_isLoose = -1; b_muon_isTracker = -1; b_muon_isGlobal = -1; b_muon_isME0 = -1;
  b_muon_pt = -9; b_muon_eta = -9; b_muon_phi = -9;
  b_muon_nMatchedStationLayer = -9;

} 
void MuonTrackAnalyser::initValue() {
  b_nCSCSegments  = 0;
  b_nGEMSegments  = 0;
  b_nGEMCSCSegments  = 0;
  b_mu_pt = 0;
  b_mu_eta = -9;
  b_mu_pdgId = 0;

  /*CSC seg*/
  b_CSC_Seg_region = -9; b_CSC_Seg_station = -9; b_CSC_Seg_ring = -9; b_CSC_Seg_chamber = -9; b_CSC_Seg_layer = -9; b_CSC_Seg_nRecHits = -9;
  b_CSC_Seg_eta    = -9; b_CSC_Seg_phi = -9;
  b_CSC_Rec_cls = -1;

  /*GEM seg*/
  b_GEM_Seg_region = -9; b_GEM_Seg_station = -9; b_GEM_Seg_ring = -9; b_GEM_Seg_chamber = -9; b_GEM_Seg_layer = -9; b_GEM_Seg_nRecHits = -9;
  b_GEM_Seg_eta    = -9; b_GEM_Seg_phi = -9;
  b_GEM_Rec_cls = -1;

  /*GEMCSC seg*/
  b_GEMCSC_Seg_region = -9; b_GEMCSC_Seg_station = -9; b_GEMCSC_Seg_ring = -9; b_GEMCSC_Seg_chamber = -9; b_GEMCSC_Seg_layer = -9; b_GEMCSC_Seg_nRecHits = -9;
  b_GEMCSC_Seg_eta    = -9; b_GEMCSC_Seg_phi = -9;
  b_GEMCSC_Rec_cls = -1;

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTrackAnalyser);
