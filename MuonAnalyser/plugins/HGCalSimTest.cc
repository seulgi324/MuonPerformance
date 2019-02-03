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

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class HGCalSimTest : public edm::EDAnalyzer {
public:
  explicit HGCalSimTest(const edm::ParameterSet&);
  ~HGCalSimTest();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  void initValue();

  // ----------member data ---------------------------
  edm::EDGetTokenT<ME0DigiCollection>        me0Digis_;
  edm::EDGetTokenT<ME0SegmentCollection>     me0Segments_;
  edm::EDGetTokenT<ME0RecHitCollection>      me0RecHits_;

  edm::EDGetTokenT<CSCSegmentCollection>     cscSegments_;
  edm::EDGetTokenT<CSCRecHit2DCollection>    csc2DRecHits_;

  edm::EDGetTokenT<GEMDigiCollection>        gemDigis_;
  edm::EDGetTokenT<GEMSegmentCollection>     gemSegments_;
  edm::EDGetTokenT<GEMRecHitCollection>      gemRecHits_;

  edm::EDGetTokenT<DTDigiCollection>         dtDigis_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dt4DSegments_;
  edm::EDGetTokenT<DTRecHitCollection>       dtRecHits_;

  edm::EDGetTokenT<RPCDigiCollection>        rpcDigis_;
  edm::EDGetTokenT<RPCRecHitCollection>      rpcRecHits_;

  edm::Service<TFileService> fs;

  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nME0Digis, b_nME0Segments, b_nME0RecHits, b_nCSCSegments, b_nCSC2DRecHits, b_nGEMDigis, b_nGEMSegments, b_nGEMRecHits, b_nDTDigis, b_nDT4DSegments, b_nDTRecHits, b_nRPCDigis, b_nRPCRecHits;

  /* ME0 */
  TTree *t_ME0_digi;
  int b_ME0_Digi_chamber, b_ME0_Digi_layer, b_ME0_Digi_etaPartition;
  TTree *t_ME0_seg;
  int b_ME0_Seg_chamber, b_ME0_Seg_layer, b_ME0_Seg_nRecHits;
  float b_ME0_Seg_eta;
  TTree *t_ME0_rec;
  int b_ME0_RecHit_chamber, b_ME0_RecHit_layer, b_ME0_RecHit_etaPartition;
  float b_ME0_RecHit_eta; 

  /* CSC */
  TTree *t_CSC_seg;
  int b_CSC_Seg_chamber, b_CSC_Seg_layer, b_CSC_Seg_station, b_CSC_Seg_ring, b_CSC_Seg_nRecHits;
  float b_CSC_Seg_eta;
  TTree *t_CSC_rec;
  int b_CSC_2DRecHit_chamber, b_CSC_2DRecHit_layer, b_CSC_2DRecHit_station, b_CSC_2DRecHit_ring, b_CSC_2DRecHit_etaPartition;
  float b_CSC_2DRecHit_eta; 

  /* GEM */
  TTree *t_GEM_digi;
  int b_GEM_Digi_chamber, b_GEM_Digi_layer, b_GEM_Digi_station, b_GEM_Digi_ring, b_GEM_Digi_etaPartition;
  TTree *t_GEM_seg;
  int b_GEM_Seg_chamber, b_GEM_Seg_layer, b_GEM_Seg_station, b_GEM_Seg_ring, b_GEM_Seg_nRecHits;
  float b_GEM_Seg_eta;
  TTree *t_GEM_rec;
  int b_GEM_RecHit_chamber, b_GEM_RecHit_layer, b_GEM_RecHit_station, b_GEM_RecHit_ring, b_GEM_RecHit_etaPartition;
  float b_GEM_RecHit_eta; 

  /* DT */
  TTree *t_DT_digi;
  int b_DT_Digi_chamber, b_DT_Digi_layer, b_DT_Digi_superLayer, b_DT_Digi_wheel, b_DT_Digi_sector, b_DT_Digi_station;
  TTree *t_DT_seg;
  int b_DT_4DSeg_chamber, b_DT_4DSeg_wheel, b_DT_4DSeg_sector, b_DT_4DSeg_station, b_DT_4DSeg_nRecHits;
  float b_DT_4DSeg_eta;
  TTree *t_DT_rec;
  int b_DT_RecHit_chamber, b_DT_RecHit_layer, b_DT_RecHit_superLayer, b_DT_RecHit_wheel, b_DT_RecHit_sector, b_DT_RecHit_station;
  float b_DT_RecHit_eta;

  /* RPC */
  TTree *t_RPC_digi;
  int b_RPC_Digi_chamber, b_RPC_Digi_layer, b_RPC_Digi_station, b_RPC_Digi_ring, b_RPC_Digi_roll, b_RPC_Digi_sector, b_RPC_Digi_subSector;
  int b_RPC_Digi_isIRPC;
  TTree *t_RPC_rec;
  int b_RPC_RecHit_chamber, b_RPC_RecHit_layer, b_RPC_RecHit_station, b_RPC_RecHit_ring, b_RPC_RecHit_roll, b_RPC_RecHit_sector, b_RPC_RecHit_subSector;
  int b_RPC_RecHit_isIRPC;
  float b_RPC_RecHit_eta;

};

HGCalSimTest::HGCalSimTest(const edm::ParameterSet& iConfig)
{ 
  me0Digis_     = consumes<ME0DigiCollection>(iConfig.getParameter<edm::InputTag>("me0Digis"));
  me0Segments_  = consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0Segments"));
  me0RecHits_   = consumes<ME0RecHitCollection>(iConfig.getParameter<edm::InputTag>("me0RecHits"));

  cscSegments_  = consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"));
  csc2DRecHits_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("csc2DRecHits"));

  gemDigis_     = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigis"));
  gemSegments_  = consumes<GEMSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemSegments"));
  gemRecHits_   = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));

  dtDigis_      = consumes<DTDigiCollection>(iConfig.getParameter<edm::InputTag>("dtDigis"));
  dt4DSegments_ = consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dt4DSegments"));
  dtRecHits_    = consumes<DTRecHitCollection>(iConfig.getParameter<edm::InputTag>("dtRecHits"));

  rpcDigis_     = consumes<RPCDigiCollection>(iConfig.getParameter<edm::InputTag>("rpcDigis"));
  rpcRecHits_   = consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("rpcRecHits"));


  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nME0Digis",     &b_nME0Digis,     "nME0Digis/I");
  t_event->Branch("nME0Segments",  &b_nME0Segments,  "nME0Segments/I");
  t_event->Branch("nME0RecHits",   &b_nME0RecHits,   "nME0RecHits/I");
  t_event->Branch("nCSCSegments",  &b_nCSCSegments,  "nCSCSegments/I");
  t_event->Branch("nCSC2DRecHits", &b_nCSC2DRecHits, "nCSC2DRecHits/I");
  t_event->Branch("nGEMDigis",     &b_nGEMDigis,     "nGEMDigis/I");
  t_event->Branch("nGEMSegments",  &b_nGEMSegments,  "nGEMSegments/I");
  t_event->Branch("nGEMRecHits",   &b_nGEMRecHits,   "nGEMRecHits/I");
  t_event->Branch("nDTDigis",      &b_nDTDigis,      "nDTDigis/I");
  t_event->Branch("nDT4DSegments", &b_nDT4DSegments, "nDT4DSegments/I");
  t_event->Branch("nDTRecHits",    &b_nDTRecHits,    "nDTRecHits/I");
  t_event->Branch("nRPCDigis",     &b_nRPCDigis,     "nRPCDigis/I");
  t_event->Branch("nRPCRecHits",   &b_nRPCRecHits,   "nRPCRecHits/I");

  /*ME0*/
  t_ME0_digi = fs->make<TTree>("ME0_Hit", "ME0_Hit");
  t_ME0_digi->Branch("Digi_chamber",      &b_ME0_Digi_chamber,      "Digi_chamber/I");
  t_ME0_digi->Branch("Digi_layer",        &b_ME0_Digi_layer,        "Digi_layer/I");
  t_ME0_digi->Branch("Digi_etaPartition", &b_ME0_Digi_etaPartition, "Digi_etaPartition/I");

  t_ME0_seg = fs->make<TTree>("ME0_Segment", "ME0_Segment");
  t_ME0_seg->Branch("Seg_eta",      &b_ME0_Seg_eta,      "Seg_eta/F");
  t_ME0_seg->Branch("Seg_chamber",  &b_ME0_Seg_chamber,  "Seg_chamber/I");
  t_ME0_seg->Branch("Seg_layer",    &b_ME0_Seg_layer,    "Seg_layer/I");
  t_ME0_seg->Branch("Seg_nRecHits", &b_ME0_Seg_nRecHits, "Seg_nRecHits/I");

  t_ME0_rec = fs->make<TTree>("ME0_RecHit", "ME0_RecHit");
  t_ME0_rec->Branch("RecHit_chamber",      &b_ME0_RecHit_chamber,      "RecHit_chaber/I");
  t_ME0_rec->Branch("RecHit_layer",        &b_ME0_RecHit_layer,        "RecHit_layer/I");
  t_ME0_rec->Branch("RecHit_etaPartition", &b_ME0_RecHit_etaPartition, "RecHit_etaParition/I");
  t_ME0_rec->Branch("RecHit_eta",          &b_ME0_RecHit_eta,          "RecHit_eta/F");

  /*CSC*/
  t_CSC_seg = fs->make<TTree>("CSC_Segment", "CSC_Segment");
  t_CSC_seg->Branch("Seg_eta",      &b_CSC_Seg_eta,      "Seg_eta/F");
  t_CSC_seg->Branch("Seg_chamber",  &b_CSC_Seg_chamber,  "Seg_chamber/I");
  t_CSC_seg->Branch("Seg_layer",    &b_CSC_Seg_layer,    "Seg_layer/I");
  t_CSC_seg->Branch("Seg_nRecHits", &b_CSC_Seg_nRecHits, "Seg_nRecHits/I");
  t_CSC_seg->Branch("Seg_station",  &b_CSC_Seg_station,  "Seg_station/I");
  t_CSC_seg->Branch("Seg_ring",     &b_CSC_Seg_ring,     "Seg_ring/I");

  t_CSC_rec = fs->make<TTree>("CSC_2DRecHit", "CSC_2DRecHit");
  t_CSC_rec->Branch("RecHit_chamber",      &b_CSC_2DRecHit_chamber,      "RecHit_chaber/I");
  t_CSC_rec->Branch("RecHit_layer",        &b_CSC_2DRecHit_layer,        "RecHit_layer/I");
  t_CSC_rec->Branch("RecHit_eta",          &b_CSC_2DRecHit_eta,          "RecHit_eta/F");
  t_CSC_rec->Branch("RecHit_station",      &b_CSC_2DRecHit_station,      "RecHit_station/I");
  t_CSC_rec->Branch("RecHit_ring",         &b_CSC_2DRecHit_ring,         "RecHit_ring/I");

  /*GEM*/
  t_GEM_digi = fs->make<TTree>("GEM_Hit", "GEM_Hit");
  t_GEM_digi->Branch("Digi_chamber",      &b_GEM_Digi_chamber,      "Digi_chamber/I");
  t_GEM_digi->Branch("Digi_layer",        &b_GEM_Digi_layer,        "Digi_layer/I");
  t_GEM_digi->Branch("Digi_etaPartition", &b_GEM_Digi_etaPartition, "Digi_etaPartition/I");
  t_GEM_digi->Branch("Digi_station",      &b_GEM_Digi_station,      "Digi_station/I");
  t_GEM_digi->Branch("Digi_ring",         &b_GEM_Digi_ring,         "Digi_ring/I");

  t_GEM_seg = fs->make<TTree>("GEM_Segment", "GEM_Segment");
  t_GEM_seg->Branch("Seg_eta",      &b_GEM_Seg_eta,      "Seg_eta/F");
  t_GEM_seg->Branch("Seg_chamber",  &b_GEM_Seg_chamber,  "Seg_chamber/I");
  t_GEM_seg->Branch("Seg_layer",    &b_GEM_Seg_layer,    "Seg_layer/I");
  t_GEM_seg->Branch("Seg_nRecHits", &b_GEM_Seg_nRecHits, "Seg_nRecHits/I");
  t_GEM_seg->Branch("Seg_station",  &b_GEM_Seg_station,  "Seg_station/I");
  t_GEM_seg->Branch("Seg_ring",     &b_GEM_Seg_ring,     "Seg_ring/I");

  t_GEM_rec = fs->make<TTree>("GEM_RecHit", "GEM_RecHit");
  t_GEM_rec->Branch("RecHit_chamber",      &b_GEM_RecHit_chamber,      "RecHit_chaber/I");
  t_GEM_rec->Branch("RecHit_layer",        &b_GEM_RecHit_layer,        "RecHit_layer/I");
  t_GEM_rec->Branch("RecHit_etaPartition", &b_GEM_RecHit_etaPartition, "RecHit_etaParition/I");
  t_GEM_rec->Branch("RecHit_eta",          &b_GEM_RecHit_eta,          "RecHit_eta/F");
  t_GEM_rec->Branch("RecHit_station",      &b_GEM_RecHit_station,      "RecHit_station/I");
  t_GEM_rec->Branch("RecHit_ring",         &b_GEM_RecHit_ring,         "RecHit_ring/I");

  /*DT*/
  t_DT_digi = fs->make<TTree>("DT_Hit", "DT_Hit");
  t_DT_digi->Branch("Digi_chamber",      &b_DT_Digi_chamber,      "Digi_chamber/I");
  t_DT_digi->Branch("Digi_layer",        &b_DT_Digi_layer,        "Digi_layer/I");
  t_DT_digi->Branch("Digi_superLayer",   &b_DT_Digi_superLayer,   "Digi_superLayer/I");
  t_DT_digi->Branch("Digi_wheel",        &b_DT_Digi_wheel,        "Digi_wheel/I");
  t_DT_digi->Branch("Digi_sector",       &b_DT_Digi_sector,       "Digi_sector/I");
  t_DT_digi->Branch("Digi_station",      &b_DT_Digi_station,      "Digi_station/I");

  t_DT_seg = fs->make<TTree>("DT_4DSegment", "DT_4DSegment");
  t_DT_seg->Branch("Seg_eta",        &b_DT_4DSeg_eta,        "Seg_eta/F");
  t_DT_seg->Branch("Seg_chamber",    &b_DT_4DSeg_chamber,    "Seg_chamber/I");
  t_DT_seg->Branch("Seg_nRecHits",   &b_DT_4DSeg_nRecHits,   "Seg_nRecHits/I");
  t_DT_seg->Branch("Seg_wheel",      &b_DT_4DSeg_wheel,      "Seg_wheel/I");
  t_DT_seg->Branch("Seg_sector",     &b_DT_4DSeg_sector,     "Seg_sector/I");
  t_DT_seg->Branch("Seg_station",    &b_DT_4DSeg_station,    "Seg_station/I");

  t_DT_rec = fs->make<TTree>("DT_RecHit", "DT_RecHit");
  t_DT_rec->Branch("RecHit_chamber",      &b_DT_RecHit_chamber,      "RecHit_chaber/I");
  t_DT_rec->Branch("RecHit_layer",        &b_DT_RecHit_layer,        "RecHit_layer/I");
  t_DT_rec->Branch("RecHit_superLayer",   &b_DT_RecHit_superLayer,   "RecHit_superLayer/I");
  t_DT_rec->Branch("RecHit_eta",          &b_DT_RecHit_eta,          "RecHit_eta/F");
  t_DT_rec->Branch("RecHit_wheel",        &b_DT_RecHit_wheel,        "RecHit_wheel/I");
  t_DT_rec->Branch("RecHit_sector",       &b_DT_RecHit_sector,       "RecHit_sector/I");
  t_DT_rec->Branch("RecHit_station",      &b_DT_RecHit_station,      "RecHit_station/I");

  /*RPC*/
  t_RPC_digi = fs->make<TTree>("RPC_Hit", "RPC_Hit");
  t_RPC_digi->Branch("Digi_chamber",      &b_RPC_Digi_chamber,      "Digi_chamber/I");
  t_RPC_digi->Branch("Digi_layer",        &b_RPC_Digi_layer,        "Digi_layer/I");
  t_RPC_digi->Branch("Digi_station",      &b_RPC_Digi_station,      "Digi_station/I");
  t_RPC_digi->Branch("Digi_ring",         &b_RPC_Digi_ring,         "Digi_ring/I");
  t_RPC_digi->Branch("Digi_roll",         &b_RPC_Digi_roll,         "Digi_roll/I");
  t_RPC_digi->Branch("Digi_sector",       &b_RPC_Digi_sector,       "Digi_sector/I");
  t_RPC_digi->Branch("Digi_subSector",    &b_RPC_Digi_subSector,    "Digi_subSector/I");
  t_RPC_digi->Branch("Digi_isIRPC",       &b_RPC_Digi_isIRPC,       "Digi_isIRPC/I");

  t_RPC_rec = fs->make<TTree>("RPC_RecHit", "RPC_RecHit");
  t_RPC_rec->Branch("RecHit_chamber",      &b_RPC_RecHit_chamber,      "RecHit_chaber/I");
  t_RPC_rec->Branch("RecHit_layer",        &b_RPC_RecHit_layer,        "RecHit_layer/I");
  t_RPC_rec->Branch("RecHit_eta",          &b_RPC_RecHit_eta,          "RecHit_eta/F");
  t_RPC_rec->Branch("RecHit_station",      &b_RPC_RecHit_station,      "RecHit_station/I");
  t_RPC_rec->Branch("RecHit_ring",         &b_RPC_RecHit_ring,         "RecHit_ring/I");
  t_RPC_rec->Branch("RecHit_roll",         &b_RPC_RecHit_roll,         "RecHit_roll/I");
  t_RPC_rec->Branch("RecHit_sector",       &b_RPC_RecHit_sector,       "RecHit_sector/I");
  t_RPC_rec->Branch("RecHit_subSector",    &b_RPC_RecHit_subSector,    "RecHit_subSector/I");
  t_RPC_rec->Branch("RecHit_isIRPC",       &b_RPC_RecHit_isIRPC,       "RecHit_isIRPC/I");

}

HGCalSimTest::~HGCalSimTest()
{
}

void
HGCalSimTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<ME0Geometry> hME0Geom;
  iSetup.get<MuonGeometryRecord>().get(hME0Geom);
  const ME0Geometry* ME0Geometry_ = &*hME0Geom;

  /* ME0 Geometry */
  edm::Handle<ME0DigiCollection> me0Digis;
  iEvent.getByToken(me0Digis_, me0Digis);

  edm::Handle<ME0SegmentCollection> me0Segments;
  iEvent.getByToken(me0Segments_, me0Segments);

  edm::Handle<ME0RecHitCollection> me0RecHits;
  iEvent.getByToken(me0RecHits_, me0RecHits);

  /* CSC Geometry */
  edm::ESHandle<CSCGeometry> hCSCGeom;
  iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  const CSCGeometry* CSCGeometry_ = &*hCSCGeom;

  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegments_, cscSegments);

  edm::Handle<CSCRecHit2DCollection> csc2DRecHits;
  iEvent.getByToken(csc2DRecHits_, csc2DRecHits);

  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<GEMDigiCollection> gemDigis;
  iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMSegmentCollection> gemSegments;
  iEvent.getByToken(gemSegments_, gemSegments);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  /* DT Geometry */
  edm::ESHandle<DTGeometry> hDTGeom;
  iSetup.get<MuonGeometryRecord>().get(hDTGeom);
  const DTGeometry* DTGeometry_ = &*hDTGeom;

  edm::Handle<DTDigiCollection> dtDigis;
  iEvent.getByToken(dtDigis_, dtDigis);

  edm::Handle<DTRecSegment4DCollection> dt4DSegments;
  iEvent.getByToken(dt4DSegments_, dt4DSegments);

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

  for (auto ch : ME0Geometry_->chambers()) {
    /* ME0 seg */
    ME0DetId cId = ch->id();
    auto segsRange = me0Segments->get(cId);
    b_ME0_Seg_chamber = cId.chamber();
    b_ME0_Seg_layer   = cId.layer();
    auto me0Seg = segsRange.first; 
    for (auto seg = me0Seg; seg != segsRange.second; ++seg) {
      auto segLd = me0Seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_ME0_Seg_eta = gp.eta();
      b_ME0_Seg_nRecHits = me0Seg->nRecHits();
      std::cout << b_ME0_Seg_chamber << " ==> seg x : " << segLd.x() << " seg y : " << segLd.y() << " seg eta : " << gp.eta() << " nRecHits : " << b_ME0_Seg_nRecHits << std::endl;
      b_nME0Segments++;
      t_ME0_seg->Fill();
    }
    for (auto ly : ch->layers()) {
      for (auto roll : ly->etaPartitions()) {
        ME0DetId rId = roll->id();
        int roll_ = rId.roll();
        int chamber = rId.chamber();
        int layer = rId.layer();
        /* ME0 digi */
        auto digisRange = me0Digis->get(rId);
        auto me0Digi = digisRange.first;
        b_ME0_Digi_chamber = chamber;
        b_ME0_Digi_layer = layer;
        for (auto hit = me0Digi; hit != digisRange.second; ++hit) {
          b_ME0_Digi_etaPartition = roll_;
          t_ME0_digi->Fill();
          b_nME0Digis++;
        }
        /* ME0 RecHit */
        auto recRange = me0RecHits->get(rId);
        auto me0Rec = recRange.first;
        b_ME0_RecHit_chamber = chamber;
        b_ME0_RecHit_layer = layer;
        for (auto rec = me0Rec; rec != recRange.second; ++rec) {
          auto recLd = me0Rec->localPosition();
          auto gp = roll->toGlobal(recLd);
          b_ME0_RecHit_eta = gp.eta();
          b_ME0_RecHit_etaPartition = roll_;
          std::cout << b_ME0_RecHit_chamber << " ==> seg x : " << recLd.x() << " seg y : " << recLd.y() << " seg eta : " << gp.eta() << " nME0RecHits : " << b_nME0RecHits << std::endl;
          b_nME0RecHits++;
          t_ME0_rec->Fill();
        }
      }
    }
  }

  for (auto ch : CSCGeometry_->chambers()) {
    /* CSC seg */
    CSCDetId cId = ch->id();
    b_CSC_Seg_station      = cId.station();
    b_CSC_2DRecHit_station = cId.station();
    b_CSC_Seg_ring         = cId.ring();
    b_CSC_2DRecHit_ring    = cId.ring();
    b_CSC_Seg_chamber      = cId.chamber();
    b_CSC_Seg_layer        = cId.layer();
    auto segsRange = cscSegments->get(cId);
    auto cscSeg = segsRange.first;
    for (auto seg = cscSeg; seg != segsRange.second; ++seg) { 
      auto segLd = cscSeg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_CSC_Seg_eta = gp.eta();
      b_CSC_Seg_nRecHits = cscSeg->nRecHits();
      b_nCSCSegments++;
      t_CSC_seg->Fill();
    }
    for (auto ly : ch->layers()) {
      CSCDetId lId = ly->id();
      /* CSC rechit */
      b_CSC_2DRecHit_chamber = lId.chamber();
      b_CSC_2DRecHit_layer   = lId.layer();
      auto recRange = csc2DRecHits->get(lId);
      auto cscRec = recRange.first;
      for (auto rec = cscRec; rec != recRange.second; ++rec) {
        auto recLd = cscRec->localPosition();
        auto gp = ly->toGlobal(recLd);
        b_CSC_2DRecHit_eta = gp.eta();
        b_nCSC2DRecHits++;
        t_CSC_rec->Fill();
      }
    }
  }

  for (auto ch : GEMGeometry_->chambers()) {
    /* GEM seg */
    GEMDetId cId = ch->id();
    b_GEM_Digi_station   = cId.station();
    b_GEM_Seg_station    = cId.station();
    b_GEM_RecHit_station = cId.station();
    b_GEM_Digi_ring      = cId.ring();
    b_GEM_Seg_ring       = cId.ring();
    b_GEM_RecHit_ring    = cId.ring();
    b_GEM_Seg_chamber    = cId.chamber();
    b_GEM_Seg_layer      = cId.layer();
    auto segsRange = gemSegments->get(cId);
    auto gemSeg = segsRange.first;
    for (auto seg = gemSeg; seg != segsRange.second; ++seg) {
      auto segLd = gemSeg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_GEM_Seg_eta = gp.eta();
      b_GEM_Seg_nRecHits = gemSeg->nRecHits();
      b_nGEMSegments++;
      t_GEM_seg->Fill();
    }
    for (auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();
      int roll_ = rId.roll();
      b_GEM_Digi_chamber   = rId.chamber();
      b_GEM_RecHit_chamber = rId.chamber();
      b_GEM_Digi_layer     = rId.layer();
      b_GEM_RecHit_layer   = rId.layer();
      /* GEM digi */
      auto digisRange = gemDigis->get(rId);
      auto gemDigi = digisRange.first;
      for (auto hit = gemDigi; hit != digisRange.second; ++hit) {
        b_GEM_Digi_etaPartition = roll_;
        t_GEM_digi->Fill();
        b_nGEMDigis++;
      }
      /* GEM rechit */
      auto recRange = gemRecHits->get(rId);
      auto gemRec = recRange.first;
      for (auto rec = gemRec; rec != recRange.second; ++rec) {
        auto recLd = gemRec->localPosition();
        auto gp = roll->toGlobal(recLd);
        b_GEM_RecHit_eta = gp.eta();
        b_GEM_RecHit_etaPartition = roll_;
        b_nGEMRecHits++;
        t_GEM_rec->Fill();
      }
    }
  }

  for (auto ch : DTGeometry_->chambers()) {
    /* DT seg */
    DTChamberId cId = ch->id();
    b_DT_4DSeg_station       = cId.station();
    b_DT_4DSeg_sector        = cId.sector();
    b_DT_4DSeg_wheel         = cId.wheel();
    //b_DT_4DSeg_chamber       = ch->id();
    auto segsRange = dt4DSegments->get(cId);
    auto dt4DSeg = segsRange.first;
    for (auto seg = dt4DSeg; seg != segsRange.second; ++seg) {
      auto segLd = dt4DSeg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_DT_4DSeg_eta = gp.eta();
      b_DT_4DSeg_nRecHits = dt4DSeg->recHits().size();
      b_nDT4DSegments++;
      t_DT_seg->Fill();
    }
  }
  for (auto ly : DTGeometry_->layers()) {
    DTLayerId lId = ly->id();
    b_DT_Digi_station      = lId.station();
    b_DT_RecHit_station    = lId.station();
    b_DT_Digi_sector       = lId.sector();
    b_DT_RecHit_sector     = lId.sector();
    b_DT_Digi_wheel        = lId.wheel();
    b_DT_RecHit_wheel      = lId.wheel();
    b_DT_Digi_layer        = lId.layer();
    b_DT_RecHit_layer      = lId.layer();
    b_DT_Digi_superLayer   = lId.superLayer();
    b_DT_RecHit_superLayer = lId.superLayer();
    /* DT digi */
    auto digisRange = dtDigis->get(lId);
    auto dtDigi = digisRange.first;
    for (auto hit = dtDigi; hit != digisRange.second; ++hit) { 
      t_DT_digi->Fill();
      b_nDTDigis++;
    }
    /* DT rechit */
    auto recRange = dtRecHits->get(lId);
    auto dtRec = recRange.first;
    for (auto rec = dtRec; rec != recRange.second; ++rec) {
      auto recLd = dtRec->localPosition();
      auto gp = ly->toGlobal(recLd);
      b_DT_RecHit_eta = gp.eta();
      b_nDTRecHits++;
      t_DT_rec->Fill();
    }
  }

  for (auto ch : RPCGeometry_->chambers()) {
    RPCDetId cId = ch->id();
    b_RPC_Digi_station     = cId.station();
    b_RPC_RecHit_station   = cId.station();
    b_RPC_Digi_ring        = cId.ring();
    b_RPC_RecHit_ring      = cId.ring();
    //b_RPC_Digi_chamber     = cId.chamber();
    //b_RPC_RecHit_chamber   = cId.chamber();
    b_RPC_Digi_sector      = cId.sector();
    b_RPC_RecHit_sector    = cId.sector();
    b_RPC_Digi_subSector   = cId.subsector();
    b_RPC_RecHit_subSector = cId.subsector();

    for (auto roll : ch->rolls()) {
      b_RPC_Digi_isIRPC   = (int)roll->isIRPC();
      b_RPC_RecHit_isIRPC = (int)roll->isIRPC();
      RPCDetId rId = roll->id();
      int roll_ = rId.roll();
      b_RPC_Digi_layer   = rId.layer();
      b_RPC_RecHit_layer = rId.layer();
      /* RPC digi */
      auto digisRange = rpcDigis->get(rId);
      auto rpcDigi = digisRange.first;
      for (auto hit = rpcDigi; hit != digisRange.second; ++hit) {
        b_RPC_Digi_roll = roll_;
        t_RPC_digi->Fill();
        b_nRPCDigis++;
      }
      /* RPC rechit */
      auto recRange = rpcRecHits->get(rId);
      auto rpcRec = recRange.first;
      for (auto rec = rpcRec; rec != recRange.second; ++rec) {
        auto recLd = rpcRec->localPosition();
        auto gp = roll->toGlobal(recLd);
        b_RPC_RecHit_eta = gp.eta();
        b_RPC_RecHit_roll = roll_;
        b_nRPCRecHits++;
        t_RPC_rec->Fill();
      }
    }
  }
  t_event->Fill();
}

void HGCalSimTest::beginJob(){}
void HGCalSimTest::endJob(){}

void HGCalSimTest::beginRun(Run const& run, EventSetup const&){
}
void HGCalSimTest::endRun(Run const&, EventSetup const&){}

void HGCalSimTest::initValue() {
  b_nME0Digis = 0; b_nME0Segments  = 0;  b_nME0RecHits = 0;
                   b_nCSCSegments  = 0;  b_nCSC2DRecHits = 0;
  b_nGEMDigis = 0; b_nGEMSegments  = 0;  b_nGEMRecHits = 0;
  b_nDTDigis  = 0; b_nDT4DSegments = 0; b_nDTRecHits = 0; 
  b_nRPCDigis = 0;                      b_nRPCRecHits = 0;  

  /*ME0 digi*/
  b_ME0_Digi_chamber = -9; b_ME0_Digi_layer = -9; b_ME0_Digi_etaPartition = -9;
  /*ME0 seg*/
  b_ME0_Seg_chamber = -9; b_ME0_Seg_layer = -9; b_ME0_Seg_nRecHits = -9;
  b_ME0_Seg_eta = -9;
  /*ME0 rechit*/
  b_ME0_RecHit_etaPartition = -9;
  b_ME0_RecHit_eta = -9;

  /*CSC seg*/
  b_CSC_Seg_chamber = -9; b_CSC_Seg_layer = -9; b_CSC_Seg_station = -9; b_CSC_Seg_ring = -9; b_CSC_Seg_nRecHits = -9;
  b_CSC_Seg_eta = -9;
  /*CSC rechit*/
  b_CSC_2DRecHit_chamber = -9; b_CSC_2DRecHit_layer = -9; b_CSC_2DRecHit_station = -9; b_CSC_2DRecHit_ring = -9; b_CSC_2DRecHit_etaPartition = -9;
  b_CSC_2DRecHit_eta = -9;

  /*GEM digi*/
  b_GEM_Digi_chamber = -9; b_GEM_Digi_layer = -9; b_GEM_Digi_station = -9; b_GEM_Digi_ring = -9; b_GEM_Digi_etaPartition = -9;
  /*GEM seg*/
  b_GEM_Seg_chamber = -9; b_GEM_Seg_layer = -9; b_GEM_Seg_station = -9; b_GEM_Seg_ring = -9; b_GEM_Seg_nRecHits = -9;
  b_GEM_Seg_eta = -9;
  /*GEM rechit*/
  b_GEM_RecHit_chamber = -9; b_GEM_RecHit_layer = -9; b_GEM_RecHit_station = -9; b_GEM_RecHit_ring = -9; b_GEM_RecHit_etaPartition = -9;
  b_GEM_RecHit_eta = -9;

  /* DT digi*/
  b_DT_Digi_chamber = -9; b_DT_Digi_layer = -9; b_DT_Digi_superLayer = -9; b_DT_Digi_wheel = -9; b_DT_Digi_sector = -9; b_DT_Digi_station = -9;
  /* DT seg*/
  b_DT_4DSeg_chamber = -9; b_DT_4DSeg_wheel = -9; b_DT_4DSeg_sector = -9; b_DT_4DSeg_station = -9; b_DT_4DSeg_nRecHits = -9;
  b_DT_4DSeg_eta = -9;
  /* DT rechit*/
  b_DT_RecHit_chamber = -9; b_DT_RecHit_layer = -9; b_DT_RecHit_superLayer = -9; b_DT_RecHit_wheel = -9; b_DT_RecHit_sector = -9; b_DT_RecHit_station = -9;
  b_DT_RecHit_eta = -9;

  /*RPC digi*/
  b_RPC_Digi_chamber = -9; b_RPC_Digi_layer = -9; b_RPC_Digi_station = -9; b_RPC_Digi_ring = -9; b_RPC_Digi_roll = -9; b_RPC_Digi_sector = -9; b_RPC_Digi_subSector = -9;
  b_RPC_Digi_isIRPC = -1;
  /*RPC rechit*/
  b_RPC_RecHit_chamber = -9; b_RPC_RecHit_layer = -9; b_RPC_RecHit_station = -9; b_RPC_RecHit_ring = -9; b_RPC_RecHit_roll = -9; b_RPC_RecHit_sector = -9; b_RPC_RecHit_subSector = -9;
  b_RPC_RecHit_eta = -9;
  b_RPC_RecHit_isIRPC = -1;

}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalSimTest);
