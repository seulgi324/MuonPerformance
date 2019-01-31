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

#include "DataFormats/GEMDigi/interface/ME0DigiCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

//#include "DataFormats/CSCRecHit/interface/CSCDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
//#include "Geometry/CSCGeometry/interface/CSCEtaPartition.h"
//#include "Geometry/CSCGeometry/interface/CSCEtaPartitionSpecs.h"

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
  edm::EDGetTokenT<ME0DigiCollection>    me0Digis_;
  edm::EDGetTokenT<ME0SegmentCollection> me0Segments_;
  edm::EDGetTokenT<ME0RecHitCollection>  me0RecHits_;

  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::EDGetTokenT<CSCRecHit2DCollection>  cscRecHits_;

  edm::EDGetTokenT<GEMDigiCollection>    gemDigis_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegments_;
  edm::EDGetTokenT<GEMRecHitCollection>  gemRecHits_;

  edm::Service<TFileService> fs;


  TTree *t_event;
  int b_run, b_lumi, b_latency;
  int b_nME0Digis, b_nME0Segments, b_nME0RecHits, b_nCSCSegments, b_nCSCRecHits, b_nGEMDigis, b_nGEMSegments, b_nGEMRecHits;

  /* me0 */
  TTree *t_me0_digi;
  int b_me0_Digi_firstStrip, b_me0_Digi_nStrips, b_me0_Digi_chamber, b_me0_Digi_layer, b_me0_Digi_etaPartition;

  TTree *t_me0_seg;
  int b_me0_Seg_chamber, b_me0_Seg_layer, b_me0_Seg_nRecHits;
  float b_me0_Seg_eta;

  TTree *t_me0_rec;
  int b_me0_RecHit_chamber, b_me0_RecHit_layer, b_me0_RecHit_etaPartition;
  float b_me0_RecHit_eta; 

  /* me11 */
  TTree *t_me11_digi;
  int b_me11_Digi_firstStrip, b_me11_Digi_nStrips, b_me11_Digi_chamber, b_me11_Digi_layer, b_me11_Digi_etaPartition;

  TTree *t_me11_seg;
  int b_me11_Seg_chamber, b_me11_Seg_layer, b_me11_Seg_nRecHits;
  float b_me11_Seg_eta;

  TTree *t_me11_rec;
  int b_me11_RecHit_chamber, b_me11_RecHit_layer, b_me11_RecHit_etaPartition;
  float b_me11_RecHit_eta; 

  /* ge11 */
  TTree *t_ge11_digi;
  int b_ge11_Digi_firstStrip, b_ge11_Digi_nStrips, b_ge11_Digi_chamber, b_ge11_Digi_layer, b_ge11_Digi_etaPartition;

  TTree *t_ge11_seg;
  int b_ge11_Seg_chamber, b_ge11_Seg_layer, b_ge11_Seg_nRecHits;
  float b_ge11_Seg_eta;

  TTree *t_ge11_rec;
  int b_ge11_RecHit_chamber, b_ge11_RecHit_layer, b_ge11_RecHit_etaPartition;
  float b_ge11_RecHit_eta; 
};

HGCalSimTest::HGCalSimTest(const edm::ParameterSet& iConfig)
{ 
  me0Digis_    = consumes<ME0DigiCollection>(iConfig.getParameter<edm::InputTag>("me0Digis"));
  me0Segments_ = consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0Segments"));
  me0RecHits_  = consumes<ME0RecHitCollection>(iConfig.getParameter<edm::InputTag>("me0RecHits"));

  cscSegments_ = consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"));
  cscRecHits_  = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("csc2DRecHits"));

  gemDigis_    = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigis"));
  gemSegments_ = consumes<GEMSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemSegments"));
  gemRecHits_  = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));


  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nME0Digis",    &b_nME0Digis,    "nME0Digis/I");
  t_event->Branch("nME0Segments", &b_nME0Segments, "nME0Segments/I");
  t_event->Branch("nME0RecHits",  &b_nME0RecHits,  "nME0RecHits/I");
  t_event->Branch("nCSCSegments", &b_nCSCSegments, "nCSCSegments/I");
  t_event->Branch("nCSCRecHits",  &b_nCSCRecHits,  "nCSCRecHits/I");
  t_event->Branch("nGEMDigis",    &b_nGEMDigis,    "nGEMDigis/I");
  t_event->Branch("nGEMSegments", &b_nGEMSegments, "nGEMSegments/I");
  t_event->Branch("nGEMRecHits",  &b_nGEMRecHits,  "nGEMRecHits/I");

  /*me0*/
  t_me0_digi = fs->make<TTree>("ME0_Hit", "ME0_Hit");
  t_me0_digi->Branch("Digi_firstStrip",   &b_me0_Digi_firstStrip,   "Digi_firstStrip/I");
  t_me0_digi->Branch("Digi_nStrips",      &b_me0_Digi_nStrips,      "Digi_nStrips/I");
  t_me0_digi->Branch("Digi_chamber",      &b_me0_Digi_chamber,      "Digi_chamber/I");
  t_me0_digi->Branch("Digi_layer",        &b_me0_Digi_layer,        "Digi_layer/I");
  t_me0_digi->Branch("Digi_etaPartition", &b_me0_Digi_etaPartition, "Digi_etaPartition/I");

  t_me0_seg = fs->make<TTree>("ME0_Segment", "ME0_Segment");
  t_me0_seg->Branch("Seg_eta",      &b_me0_Seg_eta,      "Seg_eta/F");
  t_me0_seg->Branch("Seg_chamber",  &b_me0_Seg_chamber,  "Seg_chamber/I");
  t_me0_seg->Branch("Seg_layer",    &b_me0_Seg_layer,    "Seg_layer/I");
  t_me0_seg->Branch("Seg_nRecHits", &b_me0_Seg_nRecHits, "Seg_nRecHits/I");

  t_me0_rec = fs->make<TTree>("ME0_RecHit", "ME0_RecHit");
  t_me0_rec->Branch("RecHit_chamber",      &b_me0_RecHit_chamber,      "RecHit_chaber/I");
  t_me0_rec->Branch("RecHit_layer",        &b_me0_RecHit_layer,        "RecHit_layer/I");
  t_me0_rec->Branch("RecHit_etaPartition", &b_me0_RecHit_etaPartition, "RecHit_etaParition/I");
  t_me0_rec->Branch("RecHit_eta",          &b_me0_RecHit_eta,          "RecHit_eta/F");

  /*me11*/
  t_me11_digi = fs->make<TTree>("ME11_Hit", "ME11_Hit");
  t_me11_digi->Branch("Digi_firstStrip",   &b_me11_Digi_firstStrip,   "Digi_firstStrip/I");
  t_me11_digi->Branch("Digi_nStrips",      &b_me11_Digi_nStrips,      "Digi_nStrips/I");
  t_me11_digi->Branch("Digi_chamber",      &b_me11_Digi_chamber,      "Digi_chamber/I");
  t_me11_digi->Branch("Digi_layer",        &b_me11_Digi_layer,        "Digi_layer/I");
  t_me11_digi->Branch("Digi_etaPartition", &b_me11_Digi_etaPartition, "Digi_etaPartition/I");

  t_me11_seg = fs->make<TTree>("ME11_Segment", "ME11_Segment");
  t_me11_seg->Branch("Seg_eta",      &b_me11_Seg_eta,      "Seg_eta/F");
  t_me11_seg->Branch("Seg_chamber",  &b_me11_Seg_chamber,  "Seg_chamber/I");
  t_me11_seg->Branch("Seg_layer",    &b_me11_Seg_layer,    "Seg_layer/I");
  t_me11_seg->Branch("Seg_nRecHits", &b_me11_Seg_nRecHits, "Seg_nRecHits/I");

  t_me11_rec = fs->make<TTree>("ME11_RecHit", "ME11_RecHit");
  t_me11_rec->Branch("RecHit_chamber",      &b_me11_RecHit_chamber,      "RecHit_chaber/I");
  t_me11_rec->Branch("RecHit_layer",        &b_me11_RecHit_layer,        "RecHit_layer/I");
  t_me11_rec->Branch("RecHit_etaPartition", &b_me11_RecHit_etaPartition, "RecHit_etaParition/I");
  t_me11_rec->Branch("RecHit_eta",          &b_me11_RecHit_eta,          "RecHit_eta/F");

  /*ge11*/
  t_ge11_digi = fs->make<TTree>("GE11_Hit", "GE11_Hit");
  t_ge11_digi->Branch("Digi_firstStrip",   &b_ge11_Digi_firstStrip,   "Digi_firstStrip/I");
  t_ge11_digi->Branch("Digi_nStrips",      &b_ge11_Digi_nStrips,      "Digi_nStrips/I");
  t_ge11_digi->Branch("Digi_chamber",      &b_ge11_Digi_chamber,      "Digi_chamber/I");
  t_ge11_digi->Branch("Digi_layer",        &b_ge11_Digi_layer,        "Digi_layer/I");
  t_ge11_digi->Branch("Digi_etaPartition", &b_ge11_Digi_etaPartition, "Digi_etaPartition/I");

  t_ge11_seg = fs->make<TTree>("GE11_Segment", "GE11_Segment");
  t_ge11_seg->Branch("Seg_eta",      &b_ge11_Seg_eta,      "Seg_eta/F");
  t_ge11_seg->Branch("Seg_chamber",  &b_ge11_Seg_chamber,  "Seg_chamber/I");
  t_ge11_seg->Branch("Seg_layer",    &b_ge11_Seg_layer,    "Seg_layer/I");
  t_ge11_seg->Branch("Seg_nRecHits", &b_ge11_Seg_nRecHits, "Seg_nRecHits/I");

  t_ge11_rec = fs->make<TTree>("GE11_RecHit", "GE11_RecHit");
  t_ge11_rec->Branch("RecHit_chamber",      &b_ge11_RecHit_chamber,      "RecHit_chaber/I");
  t_ge11_rec->Branch("RecHit_layer",        &b_ge11_RecHit_layer,        "RecHit_layer/I");
  t_ge11_rec->Branch("RecHit_etaPartition", &b_ge11_RecHit_etaPartition, "RecHit_etaParition/I");
  t_ge11_rec->Branch("RecHit_eta",          &b_ge11_RecHit_eta,          "RecHit_eta/F");


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

  edm::Handle<CSCRecHit2DCollection> cscRecHits;
  iEvent.getByToken(cscRecHits_, cscRecHits);

  /* Gem Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<GEMDigiCollection> gemDigis;
  iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMSegmentCollection> gemSegments;
  iEvent.getByToken(gemSegments_, gemSegments);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  initValue();

  for (auto ch : ME0Geometry_->chambers()) {
    /* me0 seg */
    ME0DetId cId = ch->id();
    auto segsRange = me0Segments->get(cId);
    b_me0_Seg_chamber = cId.chamber();
    auto me0Seg = segsRange.first; 
    for (auto seg = me0Seg; seg != segsRange.second; ++seg) {
      auto segLd = me0Seg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_me0_Seg_eta = gp.eta();
      b_me0_Seg_nRecHits = me0Seg->nRecHits();
      std::cout << b_me0_Seg_chamber << " ==> seg x : " << segLd.x() << " seg y : " << segLd.y() << " seg eta : " << gp.eta() << " nRecHits : " << b_me0_Seg_nRecHits << std::endl;
      b_nME0Segments++;
      t_me0_seg->Fill();
    }
    for (auto ly : ch->layers()) {
      for (auto roll : ly->etaPartitions()) {
        ME0DetId rId = roll->id();
        int roll_ = rId.roll();
        int chamber = rId.chamber();
        int layer = rId.layer();
        /* me0 digi */
        auto digisRange = me0Digis->get(rId);
        auto me0Digi = digisRange.first;
        b_me0_Digi_chamber = chamber;
        b_me0_Digi_layer = layer;
        for (auto hit = me0Digi; hit != digisRange.second; ++hit) {
          int strip = hit->strip();
          b_me0_Digi_firstStrip = strip;
          b_me0_Digi_etaPartition = roll_;
          t_me0_digi->Fill();
          b_nME0Digis++;
        }
        /* me0 RecHit */
        auto recRange = me0RecHits->get(rId);
        auto me0Rec = recRange.first;
        b_me0_RecHit_chamber = chamber;
        b_me0_RecHit_layer = layer;
        for (auto rec = me0Rec; rec != recRange.second; ++rec) {
          auto recLd = me0Rec->localPosition();
          auto gp = roll->toGlobal(recLd);
          b_me0_RecHit_eta = gp.eta();
          b_me0_RecHit_etaPartition = roll_;
          std::cout << b_me0_RecHit_chamber << " ==> seg x : " << recLd.x() << " seg y : " << recLd.y() << " seg eta : " << gp.eta() << " nME0RecHits : " << b_nME0RecHits << std::endl;
          b_nME0RecHits++;
          t_me0_rec->Fill();
        }
      }
    }
  }

  for (auto ch : CSCGeometry_->chambers()) {
    /* me11 seg */
    CSCDetId cId = ch->id();
    if (cId.ring() != 1 || cId.station() != 1) continue; // choose me11
    auto segsRange = cscSegments->get(cId);
    b_me11_Seg_chamber = cId.chamber();
    auto cscSeg = segsRange.first;
    for (auto seg = cscSeg; seg != segsRange.second; ++seg) { 
      auto segLd = cscSeg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_me11_Seg_eta = gp.eta();
      b_me11_Seg_nRecHits = cscSeg->nRecHits();
      b_nCSCSegments++;
      t_me11_seg->Fill();
    }
    for (auto ly : ch->layers()) {
      CSCDetId lId = ly->id();
      int chamber = lId.chamber();
      int layer = lId.layer();
      /* me11 digi */
      /* me11 rechit */
      auto recRange = cscRecHits->get(lId);
      auto cscRec = recRange.first;
      b_me11_RecHit_chamber = chamber;
      b_me11_RecHit_layer = layer;
      for (auto rec = cscRec; rec != recRange.second; ++rec) {
        auto recLd = cscRec->localPosition();
        auto gp = ly->toGlobal(recLd);
        b_me11_RecHit_eta = gp.eta();
        b_nCSCRecHits++;
        t_me11_rec->Fill();
      }
    }
  }

  for (auto ch : GEMGeometry_->chambers()) {
    /* ge11 seg */
    GEMDetId cId = ch->id();
    if (cId.station() != 1 || (cId.ring() != 1 && cId.ring() != 2)) continue; // choose ge11
    auto segsRange = gemSegments->get(cId);
    b_ge11_Seg_chamber = cId.chamber();
    auto gemSeg = segsRange.first;
    for (auto seg = gemSeg; seg != segsRange.second; ++seg) {
      auto segLd = gemSeg->localPosition();
      auto gp = ch->toGlobal(segLd);
      b_ge11_Seg_eta = gp.eta();
      b_ge11_Seg_nRecHits = gemSeg->nRecHits();
      b_nGEMSegments++;
      t_ge11_seg->Fill();
    }
    for (auto roll : ch->etaPartitions()) {
      GEMDetId rId = roll->id();
      int roll_ = rId.roll();
      int chamber = rId.chamber();
      int layer = rId.layer();
      /* ge11 digi */
      auto digisRange = gemDigis->get(rId);
      auto gemDigi = digisRange.first;
      b_ge11_Digi_chamber = chamber;
      b_ge11_Digi_layer = layer;
      for (auto hit = gemDigi; hit != digisRange.second; ++hit) {
        int strip = hit->strip();
        b_ge11_Digi_firstStrip = strip;
        b_ge11_Digi_etaPartition = roll_;
        t_ge11_digi->Fill();
        b_nGEMDigis++;
      }
      /* ge11 rechit */
      auto recRange = gemRecHits->get(rId);
      auto gemRec = recRange.first;
      b_ge11_RecHit_chamber = chamber;
      b_ge11_RecHit_layer = layer;
      for (auto rec = gemRec; rec != recRange.second; ++rec) {
        auto recLd = gemRec->localPosition();
        auto gp = roll->toGlobal(recLd);
        b_ge11_RecHit_eta = gp.eta();
        b_ge11_RecHit_etaPartition = roll_;
        b_nGEMRecHits++;
        t_ge11_rec->Fill();
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
  b_nME0Digis = 0; b_nME0Segments = 0; b_nME0RecHits = 0;
  b_nCSCSegments = 0; b_nCSCRecHits = 0;
  b_nGEMDigis = 0; b_nGEMSegments = 0; b_nGEMRecHits = 0;
 
  /* me0 digi*/
  b_me0_Digi_firstStrip = -1; b_me0_Digi_nStrips = -1; b_me0_Digi_chamber = -1; b_me0_Digi_layer = -1; b_me0_Digi_etaPartition = -1;
  /* me0 seg*/
  b_me0_Seg_chamber = -1, b_me0_Seg_layer = -1, b_me0_Seg_nRecHits = -1;
  b_me0_Seg_eta = -9;
  /*me0 rechit*/
  b_me0_RecHit_etaPartition = -1;
  b_me0_RecHit_eta = -9;

  /* me11 digi*/
  b_me11_Digi_firstStrip = -1; b_me11_Digi_nStrips = -1; b_me11_Digi_chamber = -1; b_me11_Digi_layer = -1; b_me11_Digi_etaPartition = -1;
  /* me11 seg*/
  b_me11_Seg_chamber = -1, b_me11_Seg_layer = -1, b_me11_Seg_nRecHits = -1;
  b_me11_Seg_eta = -9;
  /*me11 rechit*/
  b_me11_RecHit_etaPartition = -1;
  b_me11_RecHit_eta = -9;

  /* ge11 digi*/
  b_ge11_Digi_firstStrip = -1; b_ge11_Digi_nStrips = -1; b_ge11_Digi_chamber = -1; b_ge11_Digi_layer = -1; b_ge11_Digi_etaPartition = -1;
  /* ge11 seg*/
  b_ge11_Seg_chamber = -1, b_ge11_Seg_layer = -1, b_ge11_Seg_nRecHits = -1;
  b_ge11_Seg_eta = -9;
  /*ge11 rechit*/
  b_ge11_RecHit_etaPartition = -1;
  b_ge11_RecHit_eta = -9;

}
//define this as a plug-in
DEFINE_FWK_MODULE(HGCalSimTest);
