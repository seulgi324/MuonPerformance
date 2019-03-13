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
// ME0
#include "DataFormats/GEMDigi/interface/ME0DigiCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
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
// DT
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
//#include "DataFormats/DTRecHit/interface/DTRecSegment2DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/MuonDetId/interface/DTBtiId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTSectCollId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/DTTracoId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
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

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

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
  bool isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut);

  // ----------member data ---------------------------
  edm::EDGetTokenT<ME0SegmentCollection>       me0Segments_;

  edm::EDGetTokenT<CSCSegmentCollection>       cscSegments_;

  edm::EDGetTokenT<GEMSegmentCollection>       gemSegments_;

  edm::EDGetTokenT<DTRecSegment4DCollection>   dt4DSegments_;

  edm::EDGetTokenT<TrackingParticleCollection> simToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> >     muonToken_;

//  edm::EDGetTokenT<edm::View<L1MuBMTrack> >     muonL1Token_;

  const reco::VertexCollection* vertexes_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  
  edm::Service<TFileService> fs;

  const ME0Geometry* ME0Geometry_;

  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nME0Segments, b_nCSCSegments, b_nGEMSegments, b_nDT4DSegments;

  /* ME0 */
  TTree *t_ME0_seg;
  int   b_ME0_Seg_region, b_ME0_Seg_chamber, b_ME0_Seg_layer, b_ME0_Seg_nRecHits;
  float b_ME0_Seg_eta,    b_ME0_Seg_phi;

  /* CSC */
  TTree *t_CSC_seg;
  int   b_CSC_Seg_region, b_CSC_Seg_station, b_CSC_Seg_ring, b_CSC_Seg_chamber, b_CSC_Seg_layer, b_CSC_Seg_nRecHits;
  float b_CSC_Seg_eta,    b_CSC_Seg_phi;

  /* GEM */
  TTree *t_GEM_seg;
  int   b_GEM_Seg_region, b_GEM_Seg_station, b_GEM_Seg_ring, b_GEM_Seg_chamber, b_GEM_Seg_layer, b_GEM_Seg_nRecHits;
  float b_GEM_Seg_eta,    b_GEM_Seg_phi;

  /* DT */
  TTree *t_DT_seg;
  int   b_DT_4DSeg_sector, b_DT_4DSeg_station, b_DT_4DSeg_wheel, b_DT_4DSeg_chamber, b_DT_4DSeg_nRecHits, b_DT_4DSeg_hasZed;
  float b_DT_4DSeg_eta,    b_DT_4DSeg_phi;

  /* Muon */
  TTree *t_genMuon;
  float b_genMuon_pt,      b_genMuon_eta,      b_genMuon_phi;
  int   b_genMuon_isTight, b_genMuon_isMedium, b_genMuon_isLoose, b_genMuon_isTracker, b_genMuon_isGlobal, b_genMuon_isME0, b_genMuon_isLooseME0, b_genMuon_nMatchedStationLayer;

  TTree *t_Muon;
  float b_muon_pt,      b_muon_eta,      b_muon_phi;
  int   b_muon_isTight, b_muon_isMedium, b_muon_isLoose, b_muon_isTracker, b_muon_isGlobal, b_muon_isME0, b_muon_isLooseME0, b_muon_nMatchedStationLayer;
};

MuonTrackAnalyser::MuonTrackAnalyser(const edm::ParameterSet& iConfig)
{ 
  me0Segments_  = consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0Segments"));

  cscSegments_  = consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"));

  gemSegments_  = consumes<GEMSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemSegments"));

  dt4DSegments_ = consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dt4DSegments"));

  simToken_     = consumes<TrackingParticleCollection>(iConfig.getParameter<InputTag>("simLabel"));
  muonToken_    = consumes<edm::View<reco::Muon> >(iConfig.getParameter<InputTag>("muonLabel"));
  
//  muonL1Token_  = consumes<edm::View<L1MuBMTrack> >(iConfig.getParameter<InputTag>("muonL1Label"));

  vtxToken_     = consumes<vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("primaryVertex"));
  
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nME0Segments",  &b_nME0Segments,  "nME0Segments/I");
  t_event->Branch("nCSCSegments",  &b_nCSCSegments,  "nCSCSegments/I");
  t_event->Branch("nGEMSegments",  &b_nGEMSegments,  "nGEMSegments/I");
  t_event->Branch("nDT4DSegments", &b_nDT4DSegments, "nDT4DSegments/I");

  /*ME0*/
  t_ME0_seg = fs->make<TTree>("ME0_Segment", "ME0_Segment");
  t_ME0_seg->Branch("Seg_eta",      &b_ME0_Seg_eta,      "Seg_eta/F");
  t_ME0_seg->Branch("Seg_phi",      &b_ME0_Seg_phi,      "Seg_phi/F");
  t_ME0_seg->Branch("Seg_nRecHits", &b_ME0_Seg_nRecHits, "Seg_nRecHits/I");
  t_ME0_seg->Branch("Seg_region",   &b_ME0_Seg_region,   "Seg_region/I");
  t_ME0_seg->Branch("Seg_chamber",  &b_ME0_Seg_chamber,  "Seg_chamber/I");
  t_ME0_seg->Branch("Seg_layer",    &b_ME0_Seg_layer,    "Seg_layer/I");

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

  /*DT*/
  t_DT_seg = fs->make<TTree>("DT_4DSegment", "DT_4DSegment");
  t_DT_seg->Branch("Seg_eta",      &b_DT_4DSeg_eta,      "Seg_eta/F");
  t_DT_seg->Branch("Seg_phi",      &b_DT_4DSeg_phi,      "Seg_phi/F");
  t_DT_seg->Branch("Seg_nRecHits", &b_DT_4DSeg_nRecHits, "Seg_nRecHits/I");
  t_DT_seg->Branch("Seg_sector",   &b_DT_4DSeg_sector,   "Seg_sector/I");
  t_DT_seg->Branch("Seg_station",  &b_DT_4DSeg_station,  "Seg_station/I");
  t_DT_seg->Branch("Seg_wheel",    &b_DT_4DSeg_wheel,    "Seg_wheel/I");
  t_DT_seg->Branch("Seg_chamber",  &b_DT_4DSeg_chamber,  "Seg_chamber/I");
  t_DT_seg->Branch("Seg_hasZed",   &b_DT_4DSeg_hasZed,   "Seg_hasZed/I");

  /*Muon*/
  t_genMuon = fs->make<TTree>("gen_Muon", "gen_Muon");
  t_genMuon->Branch("pt",                   &b_genMuon_pt,                   "pt/F");
  t_genMuon->Branch("eta",                  &b_genMuon_eta,                  "eta/F");
  t_genMuon->Branch("phi",                  &b_genMuon_phi,                  "phi/F");
  t_genMuon->Branch("nMatchedStationLayer", &b_genMuon_nMatchedStationLayer, "nMatchedStationLayer/I");
  t_genMuon->Branch("isTracker",            &b_genMuon_isTracker,            "isTracker/I");
  t_genMuon->Branch("isGlobal",             &b_genMuon_isGlobal,             "isGlobal/I");
  t_genMuon->Branch("isME0",                &b_genMuon_isME0,                "isME0/I");
  t_genMuon->Branch("isLooseME0",           &b_genMuon_isLooseME0,           "isLooseME0/I");
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
  t_Muon->Branch("isLooseME0",           &b_muon_isLooseME0,           "isLooseME0/I");
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
  edm::ESHandle<ME0Geometry> hME0Geom;
  iSetup.get<MuonGeometryRecord>().get(hME0Geom);
  ME0Geometry_ = &*hME0Geom;

  /* ME0 Geometry */
  edm::Handle<ME0SegmentCollection> me0Segments;
  iEvent.getByToken(me0Segments_, me0Segments);

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

  /* DT Geometry */
  edm::ESHandle<DTGeometry> hDTGeom;
  iSetup.get<MuonGeometryRecord>().get(hDTGeom);
  const DTGeometry* DTGeometry_ = &*hDTGeom;

  edm::Handle<DTRecSegment4DCollection> dt4DSegments;
  iEvent.getByToken(dt4DSegments_, dt4DSegments);

  /* Muon */
  Handle<TrackingParticleCollection> simHandle;
  iEvent.getByToken(simToken_, simHandle);

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

//  edm::Handle<edm::View<L1MuBMTrack> > muonL1Handle;
//  iEvent.getByToken(muonL1Token_, muonL1Handle);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices); 

  vertexes_ = vertices.product();
  reco::Vertex pv0 = vertexes_->at(0);

  initValue();

  for (auto ch : ME0Geometry_->chambers()) {
    /* ME0 seg */
    ME0DetId cId = ch->id();
    auto segsRange = me0Segments->get(cId);
    auto me0Seg = segsRange.first; 
    for (auto seg = me0Seg; seg != segsRange.second; ++seg) {
      auto segLd = seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_ME0_Seg_eta      = gp.eta();
      b_ME0_Seg_phi      = gp.phi();
      b_ME0_Seg_nRecHits = seg->nRecHits();
      b_ME0_Seg_region   = cId.region();
      b_ME0_Seg_chamber  = cId.chamber();
      b_ME0_Seg_layer    = cId.layer();
      //std::cout << b_ME0_Seg_chamber << " ==> seg x : " << segLd.x() << " seg y : " << segLd.y() << " seg eta : " << gp.eta() << " nRecHits : " << b_ME0_Seg_nRecHits << std::endl;
      b_nME0Segments++;
      t_ME0_seg->Fill();
    }
  }

  /* CSC */
  for (auto ch : CSCGeometry_->chambers()) {
    // CSC seg
    CSCDetId cId = ch->id();
    auto segsRange = cscSegments->get(cId);
    auto cscSeg = segsRange.first;
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

  /* DT */
  for (auto ch : DTGeometry_->chambers()) {
    // DT seg
    DTChamberId cId = ch->id();
    //b_DT_4DSeg_chamber       = ch->id();
    auto segsRange = dt4DSegments->get(cId);
    auto dt4DSeg = segsRange.first;
    for (auto seg = dt4DSeg; seg != segsRange.second; ++seg) {
      auto segLd = seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_DT_4DSeg_eta      = gp.eta();
      b_DT_4DSeg_phi      = gp.phi();
      b_DT_4DSeg_nRecHits = seg->recHits().size();
      b_DT_4DSeg_sector   = cId.sector();
      b_DT_4DSeg_station  = cId.station();
      b_DT_4DSeg_wheel    = cId.wheel();
      //b_DT_4DSeg_chamber  = cId.chamber();
      b_DT_4DSeg_hasZed   = seg->hasZed();
      b_nDT4DSegments++;
      t_DT_seg->Fill();
    }
  }

  /* gen Muon */
  for (TrackingParticleCollection::size_type i=0; i<simHandle->size(); i++) {
    TrackingParticleRef simRef(simHandle, i);
    const TrackingParticle* simTP = simRef.get();
    b_genMuon_pt  = simTP->pt();
    b_genMuon_eta = simTP->eta();
    b_genMuon_phi = simTP->phi();
    bool isSignalMuon = abs(simTP->pdgId())==13 && !simTP->genParticles().empty() && (simTP->eventId().event() == 0) && (simTP->eventId().bunchCrossing() == 0);
    if (!isSignalMuon) continue;
    t_genMuon->Fill();
  }
  /* reco Muon */
  for(size_t i = 0; i < muonHandle->size(); ++i) {
    initMuonValue();
    edm::RefToBase<reco::Muon> muRef = muonHandle->refAt(i);
    const reco::Muon* mu = muRef.get();

    double mom = mu->p();
    double dPhiCut_ = std::min(std::max(1.2/mom,1.2/100),0.056);
    double dPhiBendCut_ = std::min(std::max(0.2/mom,0.2/100),0.0096);
    b_muon_isLooseME0           = isME0MuonSelNew(*mu, 0.077, dPhiCut_, dPhiBendCut_);

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
  /*L1 trigger*/

  t_event->Fill();
}

void MuonTrackAnalyser::beginJob(){}
void MuonTrackAnalyser::endJob(){}

void MuonTrackAnalyser::beginRun(const edm::Run& run, const edm::EventSetup& iSetup){//(Run const& run, EventSetup const&){
}
void MuonTrackAnalyser::endRun(Run const&, EventSetup const&){}

void MuonTrackAnalyser::initMuonValue() {
  b_genMuon_pt  = -9; b_genMuon_eta = -9; b_genMuon_phi = -9;

  b_muon_isTight = -1; b_muon_isMedium = -1; b_muon_isLoose = -1; b_muon_isTracker = -1; b_muon_isGlobal = -1; b_muon_isME0 = -1; b_muon_isLooseME0 = -1;
  b_muon_pt = -9; b_muon_eta = -9; b_muon_phi = -9;
  b_muon_nMatchedStationLayer = -9;

}
void MuonTrackAnalyser::initValue() {
  b_nME0Segments  = 0;
  b_nCSCSegments  = 0;
  b_nGEMSegments  = 0;
  b_nDT4DSegments = 0;

  /*ME0 seg*/
  b_ME0_Seg_region = -9; b_ME0_Seg_chamber = -9; b_ME0_Seg_layer = -9; b_ME0_Seg_nRecHits = -9;
  b_ME0_Seg_eta    = -9; b_ME0_Seg_phi = -9;

  /*CSC seg*/
  b_CSC_Seg_region = -9; b_CSC_Seg_station = -9; b_CSC_Seg_ring = -9; b_CSC_Seg_chamber = -9; b_CSC_Seg_layer = -9; b_CSC_Seg_nRecHits = -9;
  b_CSC_Seg_eta    = -9; b_CSC_Seg_phi = -9;

  /*GEM seg*/
  b_GEM_Seg_region = -9; b_GEM_Seg_station = -9; b_GEM_Seg_ring = -9; b_GEM_Seg_chamber = -9; b_GEM_Seg_layer = -9; b_GEM_Seg_nRecHits = -9;
  b_GEM_Seg_eta    = -9; b_GEM_Seg_phi = -9;

  /* DT seg*/
  b_DT_4DSeg_sector = -9; b_DT_4DSeg_station = -9; b_DT_4DSeg_wheel = -9; b_DT_4DSeg_chamber = -9; b_DT_4DSeg_nRecHits = -9;
  b_DT_4DSeg_eta    = -9; b_DT_4DSeg_phi = -9;
  b_DT_4DSeg_hasZed = -1;
}

bool MuonTrackAnalyser::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
{
  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

      if (chamber->detector() == 5){

        for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

          LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
          LocalPoint seg_loc_coord(segment->x, segment->y, 0);
          LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
          LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

          const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);

          GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
          GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);

          //double segDPhi = segment->me0SegmentRef->deltaPhi();
          // need to check if this works
          double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
          double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);

          deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
          deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
          deltaPhiBend = std::abs(segDPhi - trackDPhi);

          if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;

        }
      }
    }

  }

  return result;

}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(MuonTrackAnalyser);
