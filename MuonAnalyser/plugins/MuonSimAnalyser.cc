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

//#include "DQMServices/Core/interface/DQMStore.h"
//#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
//#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// GEM
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
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

#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
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

class MuonSimAnalyser : public edm::EDAnalyzer {
public:
  explicit MuonSimAnalyser(const edm::ParameterSet&);
  ~MuonSimAnalyser();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  void initValue();

  // ----------member data ---------------------------
  edm::ParameterSet cfg_;
  edm::EDGetToken simHitsToken_;
  edm::EDGetToken simTracksToken_;
  edm::EDGetToken simVerticesToken_;
  
  edm::Service<TFileService> fs;

  TTree *t_event;
  int b_nGEMSimHits;

  /* GEM */
  TTree *t_GEM_simhit;
  int   b_GEM_SimHit_region, b_GEM_SimHit_station, b_GEM_SimHit_ring, b_GEM_SimHit_chamber, b_GEM_SimHit_layer, b_GEM_SimHit_etaPartition, b_GEM_SimHit_pdgId;
  float b_GEM_SimHit_pt,     b_GEM_SimHit_eta,     b_GEM_SimHit_phi; 
};

MuonSimAnalyser::MuonSimAnalyser(const edm::ParameterSet& iConfig)
{ 

  std::string simInputLabel_ = iConfig.getUntrackedParameter<std::string>("simInputLabel");

  simHitsToken_ = consumes<edm::PSimHitContainer>(edm::InputTag(simInputLabel_,"MuonGEMHits"));
  simTracksToken_ = consumes< edm::SimTrackContainer >(iConfig.getParameter<edm::InputTag>("simTrackCollection"));
  simVerticesToken_ = consumes< edm::SimVertexContainer >(iConfig.getParameter<edm::InputTag>("simVertexCollection"));
  cfg_ = iConfig;

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nGEMSimHits",   &b_nGEMSimHits,   "nGEMSimHits/I");

  /*GEM*/
  t_GEM_simhit = fs->make<TTree>("GEM_SimHit", "GEM_SimHit");
  t_GEM_simhit->Branch("SimHit_pdgId",        &b_GEM_SimHit_pdgId,        "SimHit_pdgId/I");
  t_GEM_simhit->Branch("SimHit_pt",           &b_GEM_SimHit_pt,           "SimHit_pt/F");
  t_GEM_simhit->Branch("SimHit_eta",          &b_GEM_SimHit_eta,          "SimHit_eta/F");
  t_GEM_simhit->Branch("SimHit_phi",          &b_GEM_SimHit_phi,          "SimHit_phi/F");
  t_GEM_simhit->Branch("SimHit_region",       &b_GEM_SimHit_region,       "SimHit_region/I");
  t_GEM_simhit->Branch("SimHit_station",      &b_GEM_SimHit_station,      "SimHit_station/I");
  t_GEM_simhit->Branch("SimHit_ring",         &b_GEM_SimHit_ring,         "SimHit_ring/I");
  t_GEM_simhit->Branch("SimHit_chamber",      &b_GEM_SimHit_chamber,      "SimHit_chaber/I");
  t_GEM_simhit->Branch("SimHit_layer",        &b_GEM_SimHit_layer,        "SimHit_layer/I");
  t_GEM_simhit->Branch("SimHit_etaPartition", &b_GEM_SimHit_etaPartition, "SimHit_etaParition/I");

}

MuonSimAnalyser::~MuonSimAnalyser()
{
}

void
MuonSimAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<edm::PSimHitContainer> simHits;
  edm::Handle<edm::SimTrackContainer> simTracks;
  edm::Handle<edm::SimVertexContainer> simVertices;
  iEvent.getByToken(simHitsToken_, simHits);
  iEvent.getByToken(simTracksToken_, simTracks);
  iEvent.getByToken(simVerticesToken_, simVertices);
  //if ( !simHits.isValid() || !simTracks.isValid() || !simVertices.isValid()) return;
  if (!simHits.isValid()) return;

  initValue();

  const edm::PSimHitContainer & sim_hits = *simHits.product();

  /* GEM */
  for (auto& sh : sim_hits){
    GEMDetId det_id = sh.detUnitId();
    auto vec = sh.momentumAtEntry();
    b_GEM_SimHit_pdgId        = sh.particleType();
    b_GEM_SimHit_pt           = vec.perp();
    b_GEM_SimHit_eta          = vec.eta();
    b_GEM_SimHit_phi          = sh.phiAtEntry();

    b_GEM_SimHit_region       = det_id.region();
    b_GEM_SimHit_station      = det_id.station();
    b_GEM_SimHit_ring         = det_id.ring();
    b_GEM_SimHit_chamber      = det_id.chamber();
    b_GEM_SimHit_layer        = det_id.layer();
    b_GEM_SimHit_etaPartition = det_id.roll();

    ++b_nGEMSimHits;
    t_GEM_simhit->Fill();
  }
  t_event->Fill();
}

void MuonSimAnalyser::beginJob(){}
void MuonSimAnalyser::endJob(){}

void MuonSimAnalyser::beginRun(const edm::Run& run, const edm::EventSetup& iSetup){//(Run const& run, EventSetup const&){
}
void MuonSimAnalyser::endRun(Run const&, EventSetup const&){}

void MuonSimAnalyser::initValue() {
  b_nGEMSimHits = -1;
  /*GEM */
  b_GEM_SimHit_region = -9; b_GEM_SimHit_station = -9; b_GEM_SimHit_ring = -9; b_GEM_SimHit_chamber = -9; b_GEM_SimHit_layer = -9; b_GEM_SimHit_etaPartition = -9; b_GEM_SimHit_pdgId = -99;
  b_GEM_SimHit_pt     = -9; b_GEM_SimHit_eta     = -9; b_GEM_SimHit_phi  = -9;

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonSimAnalyser);
