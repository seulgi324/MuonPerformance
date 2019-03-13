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

class MuonDetHitAnalyser : public edm::EDAnalyzer {
public:
  explicit MuonDetHitAnalyser(const edm::ParameterSet&);
  ~MuonDetHitAnalyser();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  void initValue();

  TProfile* hME0Area_;
  TH1D* hME0Counts_;

  TProfile* hGEMArea_;
  TH1D* hGEMCounts_;

  TProfile* hCSCArea_;
  TH1D* hCSCCounts_;

  TProfile* hRPCArea_;
  TH1D* hRPCCounts_;

  TH1D* hEvents_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<ME0DigiCollection>          me0Digis_;
  edm::EDGetTokenT<ME0RecHitCollection>        me0RecHits_;

  edm::EDGetTokenT<CSCRecHit2DCollection>      csc2DRecHits_;

  edm::EDGetTokenT<GEMDigiCollection>          gemDigis_;
  edm::EDGetTokenT<GEMRecHitCollection>        gemRecHits_;

  edm::EDGetTokenT<DTDigiCollection>           dtDigis_;
  edm::EDGetTokenT<DTRecHitCollection>         dtRecHits_;

  edm::EDGetTokenT<RPCDigiCollection>          rpcDigis_;
  edm::EDGetTokenT<RPCRecHitCollection>        rpcRecHits_;
  
  edm::Service<TFileService> fs;

  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nME0Digis, b_nME0Segments, b_nME0RecHits, b_nCSCSegments, b_nCSC2DRecHits, b_nGEMDigis, b_nGEMSegments, b_nGEMRecHits, b_nDTDigis, b_nDT4DSegments, b_nDTRecHits, b_nRPCDigis, b_nRPCRecHits;

  /* ME0 */
  TTree *t_ME0_digi;
  int   b_ME0_Digi_region, b_ME0_Digi_chamber, b_ME0_Digi_layer, b_ME0_Digi_etaPartition, b_ME0_Digi_bx;
  float b_ME0_Digi_eta,    b_ME0_Digi_phi;
  TTree *t_ME0_rec;
  int   b_ME0_RecHit_region, b_ME0_RecHit_chamber, b_ME0_RecHit_layer, b_ME0_RecHit_etaPartition;
  float b_ME0_RecHit_eta,    b_ME0_RecHit_phi; 

  /* CSC */
  TTree *t_CSC_rec;
  int   b_CSC_2DRecHit_region, b_CSC_2DRecHit_station, b_CSC_2DRecHit_ring, b_CSC_2DRecHit_chamber, b_CSC_2DRecHit_layer;
  float b_CSC_2DRecHit_eta,    b_CSC_2DRecHit_phi; 

  /* GEM */
  TTree *t_GEM_digi;
  int   b_GEM_Digi_region, b_GEM_Digi_station, b_GEM_Digi_ring, b_GEM_Digi_chamber, b_GEM_Digi_layer, b_GEM_Digi_etaPartition, b_GEM_Digi_bx;
  float b_GEM_Digi_eta,    b_GEM_Digi_phi;
  TTree *t_GEM_rec;
  int   b_GEM_RecHit_region, b_GEM_RecHit_station, b_GEM_RecHit_ring, b_GEM_RecHit_chamber, b_GEM_RecHit_layer, b_GEM_RecHit_etaPartition, b_GEM_RecHit_bx;
  float b_GEM_RecHit_eta,    b_GEM_RecHit_phi; 

  /* DT */
  TTree *t_DT_digi;
  int   b_DT_Digi_sector, b_DT_Digi_station, b_DT_Digi_wheel, b_DT_Digi_chamber,  b_DT_Digi_superLayer, b_DT_Digi_layer, b_DT_Digi_wire;
  float b_DT_Digi_eta,    b_DT_Digi_phi;
  TTree *t_DT_rec;
  int   b_DT_RecHit_sector, b_DT_RecHit_station, b_DT_RecHit_wheel, b_DT_RecHit_chamber, b_DT_RecHit_superLayer, b_DT_RecHit_layer, b_DT_RecHit_wire;
  float b_DT_RecHit_eta,    b_DT_RecHit_phi;

  /* RPC */
  TTree *t_RPC_digi;
  int   b_RPC_Digi_region, b_RPC_Digi_sector, b_RPC_Digi_subSector, b_RPC_Digi_station, b_RPC_Digi_ring, b_RPC_Digi_chamber, b_RPC_Digi_layer, b_RPC_Digi_roll, b_RPC_Digi_isIRPC;
  float b_RPC_Digi_eta,    b_RPC_Digi_phi;
  TTree *t_RPC_rec;
  int   b_RPC_RecHit_region, b_RPC_RecHit_sector, b_RPC_RecHit_subSector, b_RPC_RecHit_station, b_RPC_RecHit_ring, b_RPC_RecHit_chamber, b_RPC_RecHit_layer, b_RPC_RecHit_roll, b_RPC_RecHit_isIRPC;
  float b_RPC_RecHit_eta,    b_RPC_RecHit_phi;
};

MuonDetHitAnalyser::MuonDetHitAnalyser(const edm::ParameterSet& iConfig)
{ 
  me0Digis_     = consumes<ME0DigiCollection>(iConfig.getParameter<edm::InputTag>("me0Digis"));
  me0RecHits_   = consumes<ME0RecHitCollection>(iConfig.getParameter<edm::InputTag>("me0RecHits"));

  csc2DRecHits_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("csc2DRecHits"));

  gemDigis_     = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigis"));
  gemRecHits_   = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));

  dtDigis_      = consumes<DTDigiCollection>(iConfig.getParameter<edm::InputTag>("dtDigis"));
  dtRecHits_    = consumes<DTRecHitCollection>(iConfig.getParameter<edm::InputTag>("dtRecHits"));

  rpcDigis_     = consumes<RPCDigiCollection>(iConfig.getParameter<edm::InputTag>("rpcDigis"));
  rpcRecHits_   = consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("rpcRecHits"));

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nME0Digis",     &b_nME0Digis,     "nME0Digis/I");
  t_event->Branch("nME0RecHits",   &b_nME0RecHits,   "nME0RecHits/I");
  t_event->Branch("nCSC2DRecHits", &b_nCSC2DRecHits, "nCSC2DRecHits/I");
  t_event->Branch("nGEMDigis",     &b_nGEMDigis,     "nGEMDigis/I");
  t_event->Branch("nGEMRecHits",   &b_nGEMRecHits,   "nGEMRecHits/I");
  t_event->Branch("nDTDigis",      &b_nDTDigis,      "nDTDigis/I");
  t_event->Branch("nDTRecHits",    &b_nDTRecHits,    "nDTRecHits/I");
  t_event->Branch("nRPCDigis",     &b_nRPCDigis,     "nRPCDigis/I");
  t_event->Branch("nRPCRecHits",   &b_nRPCRecHits,   "nRPCRecHits/I");

  /*ME0*/
  t_ME0_digi = fs->make<TTree>("ME0_Hit", "ME0_Hit");
  t_ME0_digi->Branch("Digi_eta",          &b_ME0_Digi_eta,          "Digi_eta/F");
  t_ME0_digi->Branch("Digi_phi",          &b_ME0_Digi_phi,          "Digi_phi/F");
  t_ME0_digi->Branch("Digi_bx",           &b_ME0_Digi_bx,           "Digi_bx/I");
  t_ME0_digi->Branch("Digi_region",       &b_ME0_Digi_region,       "Digi_region/I");
  t_ME0_digi->Branch("Digi_chamber",      &b_ME0_Digi_chamber,      "Digi_chamber/I");
  t_ME0_digi->Branch("Digi_layer",        &b_ME0_Digi_layer,        "Digi_layer/I");
  t_ME0_digi->Branch("Digi_etaPartition", &b_ME0_Digi_etaPartition, "Digi_etaPartition/I");

  t_ME0_rec = fs->make<TTree>("ME0_RecHit", "ME0_RecHit");
  t_ME0_rec->Branch("RecHit_eta",          &b_ME0_RecHit_eta,          "RecHit_eta/F");
  t_ME0_rec->Branch("RecHit_phi",          &b_ME0_RecHit_phi,          "RecHit_phi/F");
  t_ME0_rec->Branch("RecHit_region",       &b_ME0_RecHit_region,       "RecHit_region/I");
  t_ME0_rec->Branch("RecHit_chamber",      &b_ME0_RecHit_chamber,      "RecHit_chaber/I");
  t_ME0_rec->Branch("RecHit_layer",        &b_ME0_RecHit_layer,        "RecHit_layer/I");
  t_ME0_rec->Branch("RecHit_etaPartition", &b_ME0_RecHit_etaPartition, "RecHit_etaParition/I");

  /*CSC*/
  t_CSC_rec = fs->make<TTree>("CSC_2DRecHit", "CSC_2DRecHit");
  t_CSC_rec->Branch("RecHit_eta",     &b_CSC_2DRecHit_eta,     "RecHit_eta/F");
  t_CSC_rec->Branch("RecHit_phi",     &b_CSC_2DRecHit_phi,     "RecHit_phi/F");
  t_CSC_rec->Branch("RecHit_region",  &b_CSC_2DRecHit_region,  "RecHit_region/I");
  t_CSC_rec->Branch("RecHit_station", &b_CSC_2DRecHit_station, "RecHit_station/I");
  t_CSC_rec->Branch("RecHit_ring",    &b_CSC_2DRecHit_ring,    "RecHit_ring/I");
  t_CSC_rec->Branch("RecHit_chamber", &b_CSC_2DRecHit_chamber, "RecHit_chaber/I");
  t_CSC_rec->Branch("RecHit_layer",   &b_CSC_2DRecHit_layer,   "RecHit_layer/I");

  /*GEM*/
  t_GEM_digi = fs->make<TTree>("GEM_Hit", "GEM_Hit");
  t_GEM_digi->Branch("Digi_eta",          &b_GEM_Digi_eta,          "Digi_eta/F");
  t_GEM_digi->Branch("Digi_phi",          &b_GEM_Digi_phi,          "Digi_phi/F");
  t_GEM_digi->Branch("Digi_bx",           &b_GEM_Digi_bx,           "Digi_bx/I");
  t_GEM_digi->Branch("Digi_region",       &b_GEM_Digi_region,       "Digi_region/I");
  t_GEM_digi->Branch("Digi_station",      &b_GEM_Digi_station,      "Digi_station/I");
  t_GEM_digi->Branch("Digi_ring",         &b_GEM_Digi_ring,         "Digi_ring/I");
  t_GEM_digi->Branch("Digi_chamber",      &b_GEM_Digi_chamber,      "Digi_chamber/I");
  t_GEM_digi->Branch("Digi_layer",        &b_GEM_Digi_layer,        "Digi_layer/I");
  t_GEM_digi->Branch("Digi_etaPartition", &b_GEM_Digi_etaPartition, "Digi_etaPartition/I");

  t_GEM_rec = fs->make<TTree>("GEM_RecHit", "GEM_RecHit");
  t_GEM_rec->Branch("RecHit_eta",          &b_GEM_RecHit_eta,          "RecHit_eta/F");
  t_GEM_rec->Branch("RecHit_phi",          &b_GEM_RecHit_phi,          "RecHit_phi/F");
  t_GEM_rec->Branch("RecHit_bx",           &b_GEM_RecHit_bx,           "RecHit_bx/I");
  t_GEM_rec->Branch("RecHit_region",       &b_GEM_RecHit_region,       "RecHit_region/I");
  t_GEM_rec->Branch("RecHit_station",      &b_GEM_RecHit_station,      "RecHit_station/I");
  t_GEM_rec->Branch("RecHit_ring",         &b_GEM_RecHit_ring,         "RecHit_ring/I");
  t_GEM_rec->Branch("RecHit_chamber",      &b_GEM_RecHit_chamber,      "RecHit_chaber/I");
  t_GEM_rec->Branch("RecHit_layer",        &b_GEM_RecHit_layer,        "RecHit_layer/I");
  t_GEM_rec->Branch("RecHit_etaPartition", &b_GEM_RecHit_etaPartition, "RecHit_etaParition/I");

  /*DT*/
  t_DT_digi = fs->make<TTree>("DT_Hit", "DT_Hit");
  t_DT_digi->Branch("Digi_eta",        &b_DT_Digi_eta,        "Digi_eta/F");
  t_DT_digi->Branch("Digi_phi",        &b_DT_Digi_phi,        "Digi_phi/F");
  t_DT_digi->Branch("Digi_sector",     &b_DT_Digi_sector,     "Digi_sector/I");
  t_DT_digi->Branch("Digi_station",    &b_DT_Digi_station,    "Digi_station/I");
  t_DT_digi->Branch("Digi_wheel",      &b_DT_Digi_wheel,      "Digi_wheel/I");
  t_DT_digi->Branch("Digi_chamber",    &b_DT_Digi_chamber,    "Digi_chamber/I");
  t_DT_digi->Branch("Digi_superLayer", &b_DT_Digi_superLayer, "Digi_superLayer/I");
  t_DT_digi->Branch("Digi_layer",      &b_DT_Digi_layer,      "Digi_layer/I");
  t_DT_digi->Branch("Digi_wire",       &b_DT_Digi_wire,       "Digi_wire/I");

  t_DT_rec = fs->make<TTree>("DT_RecHit", "DT_RecHit");
  t_DT_rec->Branch("RecHit_eta",        &b_DT_RecHit_eta,        "RecHit_eta/F");
  t_DT_rec->Branch("RecHit_phi",        &b_DT_RecHit_phi,        "RecHit_phi/F");
  t_DT_rec->Branch("RecHit_sector",     &b_DT_RecHit_sector,     "RecHit_sector/I");
  t_DT_rec->Branch("RecHit_station",    &b_DT_RecHit_station,    "RecHit_station/I");
  t_DT_rec->Branch("RecHit_wheel",      &b_DT_RecHit_wheel,      "RecHit_wheel/I");
  t_DT_rec->Branch("RecHit_chamber",    &b_DT_RecHit_chamber,    "RecHit_chaber/I");
  t_DT_rec->Branch("RecHit_superLayer", &b_DT_RecHit_superLayer, "RecHit_superLayer/I");
  t_DT_rec->Branch("RecHit_layer",      &b_DT_RecHit_layer,      "RecHit_layer/I");
  t_DT_rec->Branch("RecHit_wire",       &b_DT_RecHit_wire,       "RecHit_wire/I");

  /*RPC*/
  t_RPC_digi = fs->make<TTree>("RPC_Hit", "RPC_Hit");
  t_RPC_digi->Branch("Digi_eta",       &b_RPC_Digi_eta,       "Digi_eta/F");
  t_RPC_digi->Branch("Digi_phi",       &b_RPC_Digi_phi,       "Digi_phi/F");
  t_RPC_digi->Branch("Digi_isIRPC",    &b_RPC_Digi_isIRPC,    "Digi_isIRPC/I");
  t_RPC_digi->Branch("Digi_region",    &b_RPC_Digi_region,    "Digi_region/I");
  t_RPC_digi->Branch("Digi_sector",    &b_RPC_Digi_sector,    "Digi_sector/I");
  t_RPC_digi->Branch("Digi_subSector", &b_RPC_Digi_subSector, "Digi_subSector/I");
  t_RPC_digi->Branch("Digi_station",   &b_RPC_Digi_station,   "Digi_station/I");
  t_RPC_digi->Branch("Digi_ring",      &b_RPC_Digi_ring,      "Digi_ring/I");
  t_RPC_digi->Branch("Digi_chamber",   &b_RPC_Digi_chamber,   "Digi_chamber/I");
  t_RPC_digi->Branch("Digi_layer",     &b_RPC_Digi_layer,     "Digi_layer/I");
  t_RPC_digi->Branch("Digi_roll",      &b_RPC_Digi_roll,      "Digi_roll/I");
  t_RPC_digi->Branch("Digi_isIRPC",    &b_RPC_Digi_isIRPC,    "Digi_isIRPC/I");
  t_RPC_digi->Branch("Digi_eta",       &b_RPC_Digi_eta,       "Digi_eta/F");
  t_RPC_digi->Branch("Digi_phi",       &b_RPC_Digi_phi,       "Digi_phi/F");

  t_RPC_rec = fs->make<TTree>("RPC_RecHit", "RPC_RecHit");
  t_RPC_rec->Branch("RecHit_eta",       &b_RPC_RecHit_eta,       "RecHit_eta/F");
  t_RPC_rec->Branch("RecHit_phi",       &b_RPC_RecHit_phi,       "RecHit_phi/F");
  t_RPC_rec->Branch("RecHit_isIRPC",    &b_RPC_RecHit_isIRPC,    "RecHit_isIRPC/I");
  t_RPC_rec->Branch("RecHit_region",    &b_RPC_RecHit_region,    "RecHit_region/I");
  t_RPC_rec->Branch("RecHit_sector",    &b_RPC_RecHit_sector,    "RecHit_sector/I");
  t_RPC_rec->Branch("RecHit_subSector", &b_RPC_RecHit_subSector, "RecHit_subSector/I");
  t_RPC_rec->Branch("RecHit_station",   &b_RPC_RecHit_station,   "RecHit_station/I");
  t_RPC_rec->Branch("RecHit_ring",      &b_RPC_RecHit_ring,      "RecHit_ring/I");
  t_RPC_rec->Branch("RecHit_chamber",   &b_RPC_RecHit_chamber,   "RecHit_chaber/I");
  t_RPC_rec->Branch("RecHit_layer",     &b_RPC_RecHit_layer,     "RecHit_layer/I");
  t_RPC_rec->Branch("RecHit_roll",      &b_RPC_RecHit_roll,      "RecHit_roll/I");
}

MuonDetHitAnalyser::~MuonDetHitAnalyser()
{
}

void
MuonDetHitAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<ME0Geometry> hME0Geom;
  iSetup.get<MuonGeometryRecord>().get(hME0Geom);
  const ME0Geometry* ME0Geometry_ = &*hME0Geom;

  /* ME0 Geometry */
  edm::Handle<ME0DigiCollection> me0Digis;
  iEvent.getByToken(me0Digis_, me0Digis);

  edm::Handle<ME0RecHitCollection> me0RecHits;
  iEvent.getByToken(me0RecHits_, me0RecHits);

  /* CSC Geometry */
  edm::ESHandle<CSCGeometry> hCSCGeom;
  iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  const CSCGeometry* CSCGeometry_ = &*hCSCGeom;

  edm::Handle<CSCRecHit2DCollection> csc2DRecHits;
  iEvent.getByToken(csc2DRecHits_, csc2DRecHits);

  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<GEMDigiCollection> gemDigis;
  iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  /* DT Geometry */
  edm::ESHandle<DTGeometry> hDTGeom;
  iSetup.get<MuonGeometryRecord>().get(hDTGeom);
  const DTGeometry* DTGeometry_ = &*hDTGeom;

  edm::Handle<DTDigiCollection> dtDigis;
  iEvent.getByToken(dtDigis_, dtDigis);

  edm::Handle<DTRecHitCollection> dtRecHits;
  iEvent.getByToken(dtRecHits_, dtRecHits);

  /* RPC Geometry */
  edm::ESHandle<RPCGeometry> hRPCGeom;
  iSetup.get<MuonGeometryRecord>().get(hRPCGeom);
  const RPCGeometry* RPCGeometry_ = &*hRPCGeom;

  edm::Handle<RPCDigiCollection> rpcDigis;
  iEvent.getByToken(rpcDigis_, rpcDigis);

  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHits_, rpcRecHits);

  initValue();

  hEvents_->Fill(1);

  /* ME0 */
  for (auto roll : ME0Geometry_->etaPartitions()) {
    ME0DetId rId = roll->id();
    // ME0 digi
    auto digisRange = me0Digis->get(rId);
    auto me0Digi = digisRange.first;
    for (auto hit = me0Digi; hit != digisRange.second; ++hit) {
      auto strip = hit->strip();
      auto digiLp = roll->centreOfStrip(strip);
      auto gp = roll->toGlobal(digiLp);
      b_ME0_Digi_eta          = gp.eta();
      b_ME0_Digi_phi          = gp.phi();
      b_ME0_Digi_bx           = hit->bx();
      b_ME0_Digi_region       = rId.region(); // Region id: 0 for Barrel Not in use, +/-1 For +/- Endcap 
      b_ME0_Digi_chamber      = rId.chamber();
      b_ME0_Digi_layer        = rId.layer(); 
      b_ME0_Digi_etaPartition = rId.roll(); 
      t_ME0_digi->Fill();
      b_nME0Digis++;
    }
    // ME0 RecHit
    auto recRange = me0RecHits->get(rId);
    auto me0Rec = recRange.first;
    for (auto rec = me0Rec; rec != recRange.second; ++rec) {
      auto recLd = rec->localPosition();
      auto gp = roll->toGlobal(recLd);
      b_ME0_RecHit_eta          = gp.eta();
      b_ME0_RecHit_phi          = gp.phi();
      b_ME0_RecHit_region       = rId.region(); 
      b_ME0_RecHit_chamber      = rId.chamber();
      b_ME0_RecHit_layer        = rId.layer(); 
      b_ME0_RecHit_etaPartition = rId.roll(); 
      b_nME0RecHits++;
      t_ME0_rec->Fill();
    }
  }
  for ( auto me0Hit : *me0RecHits ) {
    const auto detId = me0Hit.me0Id();
    auto region  = std::to_string(detId.region());
    auto ch      = std::to_string(detId.chamber());
    auto iEta    = std::to_string(detId.roll());
    auto station = std::to_string(detId.station());
    auto ring    = std::to_string(0);//std::to_string(detId.ring());
    auto layer   = std::to_string(detId.layer());
    const string rollName = "Region_"+region+"_station_"+station+"_ring_"+ring+"_CH_"+ch+"_layer_"+layer+"_iEta_"+iEta;
    const int idx = hME0Area_->GetXaxis()->FindBin(rollName.c_str());
    hME0Counts_->Fill(idx);
  }

  /* CSC */
  for (auto ly : CSCGeometry_->layers()) {
    CSCDetId lId = ly->id();
    // CSC rechit
    auto recRange = csc2DRecHits->get(lId);
    auto cscRec = recRange.first;
    for (auto rec = cscRec; rec != recRange.second; ++rec) {
      auto recLd = rec->localPosition();
      auto gp = ly->toGlobal(recLd);
      b_CSC_2DRecHit_eta     = gp.eta();
      b_CSC_2DRecHit_phi     = gp.phi();
      b_CSC_2DRecHit_region  = lId.endcap();
      b_CSC_2DRecHit_station = lId.station();
      b_CSC_2DRecHit_ring    = lId.ring();
      b_CSC_2DRecHit_chamber = lId.chamber();
      b_CSC_2DRecHit_layer   = lId.layer();
      b_nCSC2DRecHits++;
      t_CSC_rec->Fill();
    }
  }
  for ( auto cscHit : *csc2DRecHits ) {
    const auto detId = cscHit.cscDetId();
    auto endcap  = std::to_string(detId.endcap()); // +1 : forward , -1 : backward
    auto ch      = std::to_string(detId.chamber());
    auto ring    = std::to_string(detId.ring());
    auto station = std::to_string(detId.station());
    auto layer   = std::to_string(detId.layer());
    const string name = "Endcap_"+endcap+"_station_"+station+"_ring_"+ring+"_CH_"+ch+"_layer_"+layer;
    const int idx = hCSCArea_->GetXaxis()->FindBin(name.c_str());
    hCSCCounts_->Fill(idx);
  }


  /* GEM */
  for (auto roll : GEMGeometry_->etaPartitions()) {
    GEMDetId rId = roll->id();
    // GEM digi
    auto digisRange = gemDigis->get(rId);
    auto gemDigi = digisRange.first;
    for (auto hit = gemDigi; hit != digisRange.second; ++hit) {
      auto strip = hit->strip();
      auto digiLp = roll->centreOfStrip(strip);
      auto gp = roll->toGlobal(digiLp);
      b_GEM_Digi_eta          = gp.eta();
      b_GEM_Digi_phi          = gp.phi();
      b_GEM_Digi_bx           = hit->bx();
      b_GEM_Digi_region       = rId.region();
      b_GEM_Digi_station      = rId.station();
      b_GEM_Digi_ring         = rId.ring();
      b_GEM_Digi_chamber      = rId.chamber();
      b_GEM_Digi_layer        = rId.layer();
      b_GEM_Digi_etaPartition = rId.roll();
      t_GEM_digi->Fill();
      b_nGEMDigis++;
    }
    // GEM rechit
    auto recRange = gemRecHits->get(rId);
    auto gemRec = recRange.first;
    for (auto rec = gemRec; rec != recRange.second; ++rec) {
      auto recLd = rec->localPosition();
      auto gp = roll->toGlobal(recLd);
      b_GEM_RecHit_eta          = gp.eta();
      b_GEM_RecHit_phi          = gp.phi();
      b_GEM_RecHit_bx           = rec->BunchX();
      b_GEM_RecHit_region       = rId.region();
      b_GEM_RecHit_station      = rId.station();
      b_GEM_RecHit_ring         = rId.ring();
      b_GEM_RecHit_chamber      = rId.chamber();
      b_GEM_RecHit_layer        = rId.layer();
      b_GEM_RecHit_etaPartition = rId.roll();
      b_nGEMRecHits++;
      t_GEM_rec->Fill();
    }
  }
  for ( auto gemHit : *gemRecHits ) {
    const auto detId = gemHit.gemId();
    auto region  = std::to_string(detId.region());
    auto ch      = std::to_string(detId.chamber());
    auto iEta    = std::to_string(detId.roll());
    auto station = std::to_string(detId.station());
    auto ring    = std::to_string(detId.ring());
    auto layer   = std::to_string(detId.layer());
    const string rollName = "Region_"+region+"_station_"+station+"_ring_"+ring+"_CH_"+ch+"_layer_"+layer+"_iEta_"+iEta;

    const int idx = hGEMArea_->GetXaxis()->FindBin(rollName.c_str());
    hGEMCounts_->Fill(idx);
  }

  /* DT */
  for (auto ly : DTGeometry_->layers()) {
    DTLayerId lId = ly->id();
    //auto ch = ly->chamber();
    //DTChamberId cId = ch->id();
    // DT digi
    auto digisRange = dtDigis->get(lId);
    auto dtDigi = digisRange.first;
    for (auto hit = dtDigi; hit != digisRange.second; ++hit) {
      b_DT_Digi_sector     = lId.sector();
      b_DT_Digi_station    = lId.station();
      b_DT_Digi_wheel      = lId.wheel();
      //b_DT_Digi_chamber    = cId.chamber();
      b_DT_Digi_superLayer = lId.superLayer(); /// Superlayers are numbered 1 (phi), 2 (Z), 3 (phi)
      b_DT_Digi_layer      = lId.layer();
      b_DT_Digi_wire       = hit->wire();
      t_DT_digi->Fill();
      b_nDTDigis++;
    }
    // DT rechit
    auto recRange = dtRecHits->get(lId);
    auto dtRec = recRange.first;
    for (auto rec = dtRec; rec != recRange.second; ++rec) {
      auto recLd = rec->localPosition();
      auto gp = ly->toGlobal(recLd);
      auto wireId = rec->wireId();
      b_DT_RecHit_eta        = gp.eta();
      b_DT_RecHit_phi        = gp.phi();
      b_DT_RecHit_sector     = lId.sector();
      b_DT_RecHit_station    = lId.station();
      b_DT_RecHit_wheel      = lId.wheel();
      //b_DT_RecHit_chamber    = cId.chamber()
      b_DT_RecHit_superLayer = lId.superLayer();
      b_DT_RecHit_layer      = lId.layer();
      b_DT_RecHit_wire       = wireId.wire();
      b_nDTRecHits++;
      t_DT_rec->Fill();
    }
  }

  /* RPC */
  for (auto roll : RPCGeometry_->rolls()) {
    RPCDetId rId = roll->id();
    //auto ch = roll->chamber();
    //auto cId = ch->id();
    // RPC digi
    auto digisRange = rpcDigis->get(rId);
    auto rpcDigi = digisRange.first;
    for (auto hit = rpcDigi; hit != digisRange.second; ++hit) {
      auto strip = hit->strip();
      auto digiLp = roll->centreOfStrip(strip);
      auto gp = roll->toGlobal(digiLp);
      b_RPC_Digi_eta       = gp.eta();
      b_RPC_Digi_phi       = gp.phi();
      b_RPC_Digi_isIRPC    = (int)roll->isIRPC();
      b_RPC_Digi_region    = rId.region();    // region : 0 for barrel +/-1 for +/- endcap
      b_RPC_Digi_sector    = rId.sector();    // Sector id: the group of chambers at same phi (and increasing r) 
      b_RPC_Digi_subSector = rId.subsector(); // subSector id : some sectors are divided along the phi direction in subsectors (from 1 to 4 in Barrel, from 1 to 6 in Endcap) 
      b_RPC_Digi_station   = rId.station();   // Station id : For Barrel: the four groups of chambers at same r (distance from beam axis) and increasing phi / For Endcap: the three groups of chambers at same z (distance from interaction point), i.e. the disk
      b_RPC_Digi_ring      = rId.ring();      // Ring id: Wheel number in Barrel (from -2 to +2) Ring Number in Endcap (from 1 to 3) Ring has a different meaning in Barrel and Endcap! In Barrel it is wheel, in Endcap Ring has a different meaning in
      //b_RPC_Digi_chamber   = cId.chamber();
      b_RPC_Digi_layer     = rId.layer();     // Layer id: each station can have two layers of chambers: layer 1 is the inner chamber and layer 2 is the outer chamber (when present) // Only in Barrel: RB1 and RB
      b_RPC_Digi_roll      = rId.roll();           // Roll id  (also known as eta partition): each chamber is divided along the strip direction in - two or three parts (rolls) for Barrel and two, three or four parts for endcap
      t_RPC_digi->Fill();
      b_nRPCDigis++;
    }
    // RPC rechit
    auto recRange = rpcRecHits->get(rId);
    auto rpcRec = recRange.first;
    for (auto rec = rpcRec; rec != recRange.second; ++rec) {
      auto recLd = rec->localPosition();
      auto gp = roll->toGlobal(recLd);
      b_RPC_RecHit_eta       = gp.eta();
      b_RPC_RecHit_phi       = gp.phi();
      b_RPC_RecHit_isIRPC    = (int)roll->isIRPC();
      b_RPC_RecHit_region    = rId.region();
      b_RPC_RecHit_sector    = rId.sector();
      b_RPC_RecHit_subSector = rId.subsector();
      b_RPC_RecHit_station   = rId.station();
      b_RPC_RecHit_ring      = rId.ring();
      //b_RPC_RecHit_chamber   = cId.chamber();
      b_RPC_RecHit_layer     = rId.layer();
      b_RPC_RecHit_roll      = rId.roll();
      b_nRPCRecHits++;
      t_RPC_rec->Fill();
    }
  }
  for ( auto rpcHit : *rpcRecHits ) {
    const auto detId = rpcHit.rawId();
    const string rollName = RPCGeomServ(detId).name();
    const int idx = hRPCArea_->GetXaxis()->FindBin(rollName.c_str());
    hRPCCounts_->Fill(idx);
  }
  t_event->Fill();
}

void MuonDetHitAnalyser::beginJob(){}
void MuonDetHitAnalyser::endJob(){}

void MuonDetHitAnalyser::beginRun(const edm::Run& run, const edm::EventSetup& iSetup){//(Run const& run, EventSetup const&){

  hEvents_ = fs->make<TH1D>("hEvents", "hEvents;Types;Number of Events", 1, 1, 2);

  hME0Area_   = fs->make<TProfile>("hME0Area", "Roll area;Roll Name;Area [cm^{2}]", 5000, 1, 5001, 0, 100000);
  hME0Counts_ = fs->make<TH1D>("hME0Counts", "Counts;Roll Index;Number of RecHits", 5000, 1, 5001);

  hGEMArea_   = fs->make<TProfile>("hGEMArea", "Roll area;Roll Name;Area [cm^{2}]", 5000, 1, 5001, 0, 100000);
  hGEMCounts_ = fs->make<TH1D>("hGEMCounts", "Counts;Roll Index;Number of RecHits", 5000, 1, 5001);

  hCSCArea_   = fs->make<TProfile>("hCSCArea", "Roll area;Roll Name;Area [cm^{2}]", 5000, 1, 5001, 0, 100000);
  hCSCCounts_ = fs->make<TH1D>("hCSCCounts", "Counts;Roll Index;Number of RecHits", 5000, 1, 5001);

  hRPCArea_   = fs->make<TProfile>("hRPCArea", "Roll area;Roll Name;Area [cm^{2}]", 5000, 1, 5001, 0, 100000);
  hRPCCounts_ = fs->make<TH1D>("hRPCCounts", "Counts;Roll Index;Number of RecHits", 5000, 1, 5001);


  edm::ESHandle<ME0Geometry> me0Geom;
  iSetup.get<MuonGeometryRecord>().get(me0Geom);

  edm::ESHandle<GEMGeometry> gemGeom;
  iSetup.get<MuonGeometryRecord>().get(gemGeom);

  edm::ESHandle<CSCGeometry> cscGeom;
  iSetup.get<MuonGeometryRecord>().get(cscGeom);

  edm::ESHandle<RPCGeometry> rpcGeom;
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);

  /* ME0 */
  int iME0 = 0;
  for ( auto roll : me0Geom->etaPartitions() ) {
    ++iME0;
    const auto detId = roll->id();
    auto region  = std::to_string(detId.region());
    auto ch      = std::to_string(detId.chamber());
    auto iEta    = std::to_string(detId.roll());
    auto station = std::to_string(detId.station()); // for now, always 1
    auto ring    = std::to_string(0);//std::to_string(detId.ring());
    auto layer   = std::to_string(detId.layer());
    const string rollName = "Region_"+region+"_station_"+station+"_ring_"+ring+"_CH_"+ch+"_layer_"+layer+"_iEta_"+iEta;
    const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*>(&(roll->topology())));
    const float striplength(top_->stripLength());
    const float pitch(roll->pitch());
    auto nstrip = top_->nstrips();
//    std::cout << " ME0 : " << rollName.c_str() << " / area : " << striplength*pitch*nstrip << std::endl;    
    hME0Area_->GetXaxis()->SetBinLabel(iME0, rollName.c_str());
    hME0Area_->Fill(iME0, striplength*pitch*nstrip);
  }

  /* GEM */
  int iGEM = 0;
  for ( auto roll : gemGeom->etaPartitions() ) {
    ++iGEM;
    const auto detId = roll->id();
    auto region  = std::to_string(detId.region());
    auto ch      = std::to_string(detId.chamber());
    auto iEta    = std::to_string(detId.roll());
    auto station = std::to_string(detId.station());
    auto ring    = std::to_string(detId.ring());
    auto layer   = std::to_string(detId.layer());
    const string rollName = "Region_"+region+"_station_"+station+"_ring_"+ring+"_CH_"+ch+"_layer_"+layer+"_iEta_"+iEta;
    const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*>(&(roll->topology())));
    const float striplength(top_->stripLength());
    const float pitch(roll->pitch());
    auto nstrip = top_->nstrips();    
//    std::cout << " GEM : " << rollName.c_str() << " / area : " << striplength*pitch*nstrip << std::endl;    
    hGEMArea_->GetXaxis()->SetBinLabel(iGEM, rollName.c_str());
    hGEMArea_->Fill(iGEM, striplength*pitch*nstrip);
  }

  /* CSC */
  int iCSC = 0;
  for (auto ly : cscGeom->layers()) {
    ++iCSC;
    auto detId = ly->id();
    auto endcap  = std::to_string(detId.endcap()); // +1 : forward , -1 : backward
    auto ch      = std::to_string(detId.chamber());
    auto ring    = std::to_string(detId.ring());
    auto station = std::to_string(detId.station());
    auto layer   = std::to_string(detId.layer());
    const string name = "Endcap_"+endcap+"_station_"+station+"_ring_"+ring+"_CH_"+ch+"_layer_"+layer;
    const RadialStripTopology* top_(dynamic_cast<const RadialStripTopology*>(&(ly->topology())));
    const float striplength(top_->stripLength());
    const float pitch(ly->geometry()->stripPitch());
    const float nstrip = top_->nstrips();
    //std::cout << " CSC : " << name << " / area : " << striplength*pitch*nstrip << std::endl;
    hCSCArea_->GetXaxis()->SetBinLabel(iCSC, name.c_str());
    hCSCArea_->Fill(iCSC, striplength*pitch*nstrip);
  }

  /* RPC */
  int iRPC = 0;
  for ( const RPCRoll* roll : rpcGeom->rolls() ) {
    ++iRPC;
    const auto detId = roll->id();
    const string rollName = RPCGeomServ(detId).name();
    const double width = roll->surface().bounds().width();
    const double height = roll->surface().bounds().length();
    //std::cout << " RPC : " << rollName.c_str() << " / area : " << width*height << std::endl;
    hRPCArea_->GetXaxis()->SetBinLabel(iRPC, rollName.c_str());
    hRPCArea_->Fill(iRPC, width*height);
  }

}
void MuonDetHitAnalyser::endRun(Run const&, EventSetup const&){}

void MuonDetHitAnalyser::initValue() {
  b_nME0Digis = 0; b_nME0RecHits   = 0;
                   b_nCSC2DRecHits = 0;
  b_nGEMDigis = 0; b_nGEMRecHits   = 0;
  b_nDTDigis  = 0; b_nDTRecHits    = 0; 
  b_nRPCDigis = 0; b_nRPCRecHits   = 0;  

  /*ME0 digi*/
  b_ME0_Digi_region = -9; b_ME0_Digi_chamber = -9; b_ME0_Digi_layer = -9; b_ME0_Digi_etaPartition = -9; b_ME0_Digi_bx = -999;
  b_ME0_Digi_eta    = -9; b_ME0_Digi_phi     = -9;
  /*ME0 rechit*/
  b_ME0_RecHit_region = -9; b_ME0_RecHit_chamber = -9; b_ME0_RecHit_layer = -9; b_ME0_RecHit_etaPartition = -9;
  b_ME0_RecHit_eta    = -9; b_ME0_RecHit_phi = -9;

  /*CSC rechit*/
  b_CSC_2DRecHit_region = -9; b_CSC_2DRecHit_station = -9; b_CSC_2DRecHit_ring = -9; b_CSC_2DRecHit_chamber = -9; b_CSC_2DRecHit_layer = -9; 
  b_CSC_2DRecHit_eta    = -9; b_CSC_2DRecHit_phi     = -9;

  /*GEM digi*/
  b_GEM_Digi_region = -9; b_GEM_Digi_station = -9; b_GEM_Digi_ring = -9; b_GEM_Digi_chamber = -9; b_GEM_Digi_layer = -9; b_GEM_Digi_etaPartition = -9; b_GEM_Digi_bx = -999;
  b_GEM_Digi_eta    = -9; b_GEM_Digi_phi     = -9;
  /*GEM rechit*/
  b_GEM_RecHit_region = -9; b_GEM_RecHit_station = -9; b_GEM_RecHit_ring = -9; b_GEM_RecHit_chamber = -9; b_GEM_RecHit_layer = -9; b_GEM_RecHit_etaPartition = -9; b_GEM_RecHit_bx = -999;
  b_GEM_RecHit_eta    = -9; b_GEM_RecHit_phi     = -9;;

  /* DT digi*/
  b_DT_Digi_sector = -9; b_DT_Digi_station = -9; b_DT_Digi_wheel = -9; b_DT_Digi_chamber = -9; b_DT_Digi_superLayer = -9; b_DT_Digi_layer = -9; b_DT_Digi_wire = -9;
  b_DT_Digi_eta    = -9; b_DT_Digi_phi     = -9;
  /* DT rechit*/
  b_DT_RecHit_sector = -9; b_DT_RecHit_station = -9; b_DT_RecHit_wheel = -9; b_DT_RecHit_chamber = -9; b_DT_RecHit_superLayer = -9; b_DT_RecHit_layer = -9; b_DT_RecHit_wire = -9;
  b_DT_RecHit_eta    = -9; b_DT_RecHit_phi     = -9;

  /*RPC digi*/
  b_RPC_Digi_region = -9; b_RPC_Digi_sector = -9; b_RPC_Digi_subSector = -9; b_RPC_Digi_station = -9; b_RPC_Digi_ring = -9; b_RPC_Digi_chamber = -9; b_RPC_Digi_layer = -9; b_RPC_Digi_roll = -9;
  b_RPC_Digi_eta    = -9; b_RPC_Digi_phi    = -9;
  b_RPC_Digi_isIRPC = -1;
  /*RPC rechit*/
  b_RPC_RecHit_region = -9; b_RPC_RecHit_sector = -9; b_RPC_RecHit_subSector = -9; b_RPC_RecHit_station = -9; b_RPC_RecHit_ring = -9; b_RPC_RecHit_chamber = -9; b_RPC_RecHit_layer = -9; b_RPC_RecHit_roll = -9;
  b_RPC_RecHit_eta    = -9; b_RPC_RecHit_phi    = -9;
  b_RPC_RecHit_isIRPC = -1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonDetHitAnalyser);
