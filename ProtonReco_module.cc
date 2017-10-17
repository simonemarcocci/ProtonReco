#ifndef ProtonReco_Module
#define ProtonReco_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
 #include "lardataobj/RawData/TriggerData.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH3.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#define DEBUG 1
//========================================================================

namespace recoanalysis{

class ProtonReco : public art::EDAnalyzer {
public:

    explicit ProtonReco(fhicl::ParameterSet const& pset);
    virtual ~ProtonReco();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

    void processMC(const art::Event& evt, bool &isFiducial);
    void processData(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
    double truthLength( const simb::MCParticle *MCparticle );
    bool insideFV(double vertex[4]);
    void doEfficiencies();
    void clear_vectors();
    void AnaTree_HitsPurity (std::vector< art::Ptr<recob::Hit> > const& , Int_t&, Float_t& , double&);

private:

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fMCTrackModuleLabel;
    std::string fTrackModuleLabel;
    std::string fVertexModuleLabel;
    std::string fPFParticleModuleLabel;
    std::string fDigitModuleLabel;
    std::string fTrackingAlgorithmLabel;
    int         fNeutrinoPDGcode;
    int		fLeptonPDGcode;
    double      fMaxNeutrinoE;
    double      fMaxLeptonP;
    bool    	fisNeutrinoInt;
	
    //output tree and variables
    TTree* ftracks;
    TFile* ffile;

    //neutrino or general info on the event 
    bool fFV;
    std::map<Int_t,Int_t> fparticle_count;
    double fnu_E;
    double fnu_x;
    double fnu_y;
    double fnu_z;
    int fn_particles; //which is the size of the vectors below
    int fn_fake_reco; //number of reco tracks which do not belong to any MC truth track
    int fcount_pr; //number of protons in the event
    int fev;
    int frun;
    int fsubrun;

    //info on the MC particle
    std::vector<float> flength;
    std::vector<float> fstartT;
    std::vector<float> fstart_x;
    std::vector<float> fstart_y;
    std::vector<float> fstart_z;
    std::vector<float> fend_x;
    std::vector<float> fend_y;
    std::vector<float> fend_z;
    std::vector<int> fn_steps;
    std::vector<int> fpdg;
    std::vector<int> fmother_id;
   std::vector<int > fg4_id; 
   std::vector<float> fp0;//initial momentum
   std::vector<float> fp0x;
   std::vector<float> fp0y;
   std::vector<float> fp0z;
   std::vector<float> fkinE;
   std::vector<float> fcostheta_muon;
   std::vector<bool> fis_leading;

   //info coming from the tracking algorithm - when there is mc truth
   std::vector<bool> fis_tracked;
   std::vector<bool> fis_mismatched; //it says if the MC truth assignment in different planes is different (possible hint for wrong or problematic backtracking)
   std::vector<float> fcostheta_muon_reco; //MCtruth info
   std::vector<int> fpdg_reco;
   std::vector<float> flength_reco;
   std::vector<float> freco_momentum_mcs; //(GeV) MCS
   std::vector<float> freco_momentum_mcs_llhd; //(GeV) MCS LLHD
   std::vector<float> freco_momentum_range; //MeV
   std::vector<float> fpurity;
   std::vector<float> fcompleteness;
   std::vector<double> fnhits;
   std::vector<float> freco_kinE;
   
   //info coming from the tracking algorithm - when there is NO mc truth
   std::vector<bool> ffake_is_tracked;
   std::vector<bool> ffake_is_mismatched;
   std::vector<float> ffake_costheta_muon_reco; //MCtruth info
   std::vector<int> ffake_pdg_reco;
   std::vector<float> ffake_length_reco;
   std::vector<float> ffake_reco_momentum_mcs; //(GeV) MCS
   std::vector<float> ffake_reco_momentum_mcs_llhd; //(GeV) MCS LLHD
   std::vector<float> ffake_reco_momentum_range; //MeV
   std::vector<float> ffake_purity;
   std::vector<float> ffake_completeness;
   std::vector<double> ffake_nhits;
   std::vector<float> ffake_kinE;

    int    MC_isCC;
    int    MC_incoming_PDG;
    double MC_incoming_P[4];
    double MC_vertex[4];
    double MC_lepton_startMomentum[4];

    int    MC_leading_protonID;
    int    MC_leading_PionPlusID;
    int    MC_leading_PionMinusID;
    int    MC_leptonID;
    int    MC_kaonID;
    int    MC_michelID;

    double MC_leptonP;
    double MC_leading_PionPlusP;
    double MC_leading_ProtonP;
    double MC_leading_PionMinusP;
    double MC_kaonP;
    double MC_michelP;


    TH3D* h_vertex_distribution;
    TH1D* h_track_length;
    TH1D* h_direction_y;
    TH1D* h_direction_theta;
    TH1D* h_direction_phi;
    TH1D* h_n_pfparticles;
    TH1D* h_n_vertexes;
    TH1D* h_n_tracks;
    TH1D* h_n_hits_track;
    TH1D* h_n_hits;
    TH1D* h_cos_recotruth;
    TH1D* h_dis_start_start;

    TH1D *h_Ev_den;
    TH1D *h_Ev_num;
    TH1D *h_Pmu_den;
    TH1D *h_Pmu_num;
    TH1D *h_theta_den;
    TH1D *h_theta_num;
    TH1D *h_Pproton_den;
    TH1D *h_Pproton_num;
    TH1D *h_Ppion_plus_den; 
    TH1D *h_Ppion_plus_num; 
    TH1D *h_Ppion_minus_den; 
    TH1D *h_Ppion_minus_num; 

    TH1D *h_Efrac_lepton;     
    TH1D *h_Ecomplet_lepton;     
    TH1D *h_Efrac_proton;     
    TH1D *h_Ecomplet_proton;     
    TH1D *h_Efrac_pion_plus;     
    TH1D *h_Ecomplet_pion_plus;     
    TH1D *h_Efrac_pion_minus;     
    TH1D *h_Ecomplet_pion_minus;     
    TH1D *h_trackRes_lepton;
    TH1D *h_trackRes_proton;
    TH1D *h_trackRes_pion_plus;
    TH1D *h_trackRes_pion_minus;

    TH1D *h_muon_length;
    TH1D *h_proton_length;
    TH1D *h_pionp_length;
    TH1D *h_pionm_length;
    TH1D *h_muonwtrk_length;
    TH1D *h_protonwtrk_length;
    TH1D *h_pionpwtrk_length;
    TH1D *h_pionmwtrk_length;

    TEfficiency* h_Eff_Ev = 0;
    TEfficiency* h_Eff_Pmu = 0;
    TEfficiency* h_Eff_theta = 0;
    TEfficiency* h_Eff_Pproton = 0;
    TEfficiency* h_Eff_Ppion_plus = 0;
    TEfficiency* h_Eff_Ppion_minus = 0;

    TEfficiency* h_Eff_Lmuon = 0;
    TEfficiency* h_Eff_Lproton = 0;
    TEfficiency* h_Eff_Lpion_plus = 0;
    TEfficiency* h_Eff_Lpion_minus = 0;


    //nucleon decay histograms
    TH1D *h_Pkaon_den;
    TH1D *h_Pkaon_num; 
    TH1D *h_Pmichel_e_den;
    TH1D *h_Pmichel_e_num;
    TH1D *h_Efrac_kaon;
    TH1D *h_Ecomplet_kaon;
    TH1D *h_trackRes_kaon; 
    TH1D *h_Efrac_michel;
    TH1D *h_Ecomplet_michel;
    TH1D *h_trackRes_michel;
    TH1D *h_kaon_length;
    TH1D *h_michel_length;
    TH1D *h_kaonwtrk_length;
    TH1D *h_michelwtrk_length;
    TEfficiency* h_Eff_Pkaon =0;
    TEfficiency* h_Eff_Pmichel =0;
    TEfficiency* h_Eff_Lkaon = 0;
    TEfficiency* h_Eff_Lmichel =0;


    float fFidVolCutX;
    float fFidVolCutY;
    float fFidVolCutZ;

    float fFidVolXmin;
    float fFidVolXmax;
    float fFidVolYmin;
    float fFidVolYmax;
    float fFidVolZmin;
    float fFidVolZmax;

    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;
    art::ServiceHandle<geo::Geometry> geom;

 }; // class ProtonReco


//========================================================================
ProtonReco::ProtonReco(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
ProtonReco::~ProtonReco(){
  //destructor
}
//========================================================================
void ProtonReco::reconfigure(fhicl::ParameterSet const& p){

    fMCTruthModuleLabel  = p.get<std::string>("MCTruthModuleLabel");
    fMCTrackModuleLabel  = p.get<std::string>("MCTrackModuleLabel");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fCalorimetryModuleLabel    = p.get<std::string>("CalorimetryModuleLabel");
    fTrackingAlgorithmLabel    = p.get<std::string>("TrackingAlgorithmLabel");
    fVertexModuleLabel    = p.get<std::string>("VertexModuleLabel");
    fPFParticleModuleLabel    = p.get<std::string>("PFParticleModuleLabel");
    fDigitModuleLabel    = p.get<std::string>("DigitModuleLabel");
    fNeutrinoPDGcode     = p.get<int>("NeutrinoPDGcode");
    fFidVolCutX          = p.get<float>("FidVolCutX");
    fFidVolCutY          = p.get<float>("FidVolCutY");
    fFidVolCutZ          = p.get<float>("FidVolCutZ");
    
}

void ProtonReco::AnaTree_HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe){
	trackid = -1;
	purity = -1;
	
	art::ServiceHandle<cheat::BackTracker> bt;
	std::map<int,double> trkide;
	 
	for(size_t h = 0; h < hits.size(); ++h){
	     art::Ptr<recob::Hit> hit = hits[h];
	     std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(hit);
#if DEBUG == 1
	     std::cout << "hit time " << hit->PeakTime() << " peak amplitude " << hit->PeakAmplitude() <<std::endl; 
	     std::cout << "integral " << hit->Integral() << " multiplicity " << hit->Multiplicity() <<std::endl; 
	     //std::cout << "ide size " << ides.size() << std::endl;
	     std::cout << "eveID size  " << eveIDs.size() << std::endl; 
#endif
	     //std::vector<sim::IDE> ides;
	     //bt->HitToSimIDEs(hits[h],ides);
	     for(size_t e = 0; e < eveIDs.size(); ++e){
	       trkide[eveIDs[e].trackID] += eveIDs[e].energy;
	     }
	}
	 
	maxe = -1;
	double tote = 0;
	for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	    tote += ii->second;
	    if ((ii->second)>maxe){
	       maxe = ii->second;
	       trackid = ii->first;
	    }
	}
	 
#if DEBUG == 1
	std::cout << "the total energy of this reco track is: " << tote << std::endl;
	std::cout << "track id " << trackid << std::endl;
#endif	 
	if (tote>0){
	purity = maxe/tote;
	}
} //the MC track releasing the most energy in the hits is associated to the reco track

void ProtonReco::clear_vectors(){
    //info on the MC particle
    flength.clear();
    fstartT.clear();
    fstart_x.clear();
    fstart_y.clear();
    fstart_z.clear();
    fend_x.clear();
    fend_y.clear();
    fend_z.clear();
  fn_steps.clear();
  fpdg.clear();
  fmother_id.clear();
  fg4_id.clear(); 
   fp0.clear();//initial momentum
   fp0x.clear();
   fp0y.clear();
   fp0z.clear();
   fkinE.clear();
   fcostheta_muon.clear();
  fis_leading.clear();

   //info coming from the tracking algorithm - when there is mc truth
  fis_tracked.clear();
  fis_mismatched.clear();
   fcostheta_muon_reco.clear(); //MCtruth info
 fpdg.clear();
   flength.clear();
   freco_momentum_mcs.clear(); //(GeV) MCS
   freco_momentum_mcs_llhd.clear(); //(GeV) MCS LLHD
   freco_momentum_range.clear(); //MeV
   fpurity.clear();
   fcompleteness.clear();
   fnhits.clear();
   fkinE.clear();
   
   //info coming from the tracking algorithm - when there is NO mc truth
  ffake_is_tracked.clear();
  ffake_is_mismatched.clear();
   ffake_costheta_muon_reco.clear(); //MCtruth info
 ffake_pdg_reco.clear();
   ffake_length_reco.clear();
   ffake_reco_momentum_mcs.clear(); //(GeV) MCS
   ffake_reco_momentum_mcs_llhd.clear(); //(GeV) MCS LLHD
   ffake_reco_momentum_range.clear(); //MeV
   ffake_purity.clear();
   ffake_completeness.clear();
   ffake_nhits.clear();
   ffake_kinE.clear();
}
//========================================================================
void ProtonReco::beginJob(){
  
// Get geometry.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  std::cout<<"Fiducial volume:"<<"\n"
	   <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
	   <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
	   <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";

  art::ServiceHandle<art::TFileService> tfs;

  h_vertex_distribution = tfs->make<TH3D>("reco_vertex","reco_vertex",100, fFidVolXmin, fFidVolXmax, 100, fFidVolYmin, fFidVolYmax, 500, fFidVolZmin, fFidVolZmax);
  h_track_length = tfs->make<TH1D>("track_length","track_length", 10000, 0,1000);
  h_direction_y = tfs->make<TH1D>("zenit_direction","zenit_direction",314,0,3.14);
  h_direction_theta = tfs->make<TH1D>("theta","theta",314,0,3.14);
  h_direction_phi = tfs->make<TH1D>("phi","phi",628,0,6.28);
  h_n_pfparticles = tfs->make<TH1D>("n_pfparticles","n_pfparticles",31,-0.5,30.5);
  h_n_vertexes = tfs->make<TH1D>("n_vertexes","n_vertexes",31,-0.5,30.5);
  h_n_tracks = tfs->make<TH1D>("n_tracks","n_tracks",31,-0.5,30.5);
  h_n_hits_track = tfs->make<TH1D>("n_hits_track","n_hits_track",31,-0.5,30.5);
  h_n_hits = tfs->make<TH1D>("n_hits","n_hits",31,-0.5,30.5);
    
  h_cos_recotruth = tfs->make<TH1D>("cos_recotruth","cos_recotruth",100,-1,1);
  h_dis_start_start = tfs->make<TH1D>("dis_start_start","dis_start_start",1000,0,10);

  /*
  h_Ev_den = tfs->make<TH1D>("h_Ev_den","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
  h_Ev_num = tfs->make<TH1D>("h_Ev_num","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
  h_Pmu_den = tfs->make<TH1D>("h_Pmu_den","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,E_bins);
  h_Pmu_num = tfs->make<TH1D>("h_Pmu_num","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,E_bins);
  h_theta_den = tfs->make<TH1D>("h_theta_den","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
  h_theta_num = tfs->make<TH1D>("h_theta_num","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
  h_Pproton_den = tfs->make<TH1D>("h_Pproton_den","Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Pproton_num = tfs->make<TH1D>("h_Pproton_num","Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Ppion_plus_den = tfs->make<TH1D>("h_Ppion_plus_den", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ppion_plus_num = tfs->make<TH1D>("h_Ppion_plus_num", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ppion_minus_den = tfs->make<TH1D>("h_Ppion_minus_den", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ppion_minus_num = tfs->make<TH1D>("h_Ppion_minus_num", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ev_den->Sumw2();
  h_Ev_num->Sumw2();
  h_Pmu_den->Sumw2();
  h_Pmu_num->Sumw2();
  h_theta_den->Sumw2();
  h_theta_num->Sumw2();
  h_Pproton_den->Sumw2();
  h_Pproton_num->Sumw2();
  h_Ppion_plus_den->Sumw2();
  h_Ppion_plus_num->Sumw2();
  h_Ppion_minus_den->Sumw2();
  h_Ppion_minus_num->Sumw2();

  h_Efrac_lepton = tfs->make<TH1D>("h_Efrac_lepton","Efrac Lepton; Track Purity;",60,0,1.2);
  h_Ecomplet_lepton = tfs->make<TH1D>("h_Ecomplet_lepton","Ecomplet Lepton; Track Completeness;",60,0,1.2);
  h_Efrac_proton = tfs->make<TH1D>("h_Efrac_proton","Efrac Proton; Track Purity;",60,0,1.2);
  h_Ecomplet_proton = tfs->make<TH1D>("h_Ecomplet_proton","Ecomplet Proton; Track Completeness;",60,0,1.2);
  h_Efrac_pion_plus = tfs->make<TH1D>("h_Efrac_pion_plus","Efrac Pion +; Track Purity;",60,0,1.2);
  h_Ecomplet_pion_plus = tfs->make<TH1D>("h_Ecomplet_pion_plus","Ecomplet Pion +; Track Completeness;",60,0,1.2);
  h_Efrac_pion_minus = tfs->make<TH1D>("h_Efrac_pion_minus","Efrac Pion -; Track Purity;",60,0,1.2);
  h_Ecomplet_pion_minus = tfs->make<TH1D>("h_Ecomplet_pion_minus","Ecomplet Pion -; Track Completeness;",60,0,1.2);
  h_trackRes_lepton = tfs->make<TH1D>("h_trackRes_lepton", "Muon Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_proton = tfs->make<TH1D>("h_trackRes_proton", "Proton Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_pion_plus = tfs->make<TH1D>("h_trackRes_pion_plus", "Pion + Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_pion_minus = tfs->make<TH1D>("h_trackRes_pion_minus", "Pion - Residual; Truth length - Reco length (cm);",200,-100,100);
  h_Efrac_lepton->Sumw2();
  h_Ecomplet_lepton->Sumw2();
  h_Efrac_proton->Sumw2();
  h_Ecomplet_proton->Sumw2();
  h_Efrac_pion_plus->Sumw2();
  h_Ecomplet_pion_plus->Sumw2();
  h_Efrac_pion_minus->Sumw2();
  h_Ecomplet_pion_minus->Sumw2();
  h_trackRes_lepton->Sumw2();
  h_trackRes_proton->Sumw2();
  h_trackRes_pion_plus->Sumw2();
  h_trackRes_pion_minus->Sumw2();

  h_muon_length = tfs->make<TH1D>("h_muon_length","Muon Length;  Muon Truth Length (cm)",40,0,100);
  h_proton_length = tfs->make<TH1D>("h_proton_length","Proton Length; Proton Truth Length (cm)",40,0,100);
  h_pionp_length = tfs->make<TH1D>("h_pionp_length","Pion + Length; Pion^{+} Truth Length (cm)",40,0,100);
  h_pionm_length = tfs->make<TH1D>("h_pionm_length","Pion - Length; Pion^{-} Truth Length (cm)",40,0,100);

  h_muonwtrk_length = tfs->make<TH1D>("h_muonwtrk_length","Muon Length; Muon Truth Length (cm)",40,0,100);
  h_protonwtrk_length = tfs->make<TH1D>("h_protonwtrk_length","Proton Length; Proton Truth Length (cm)",40,0,100);
  h_pionpwtrk_length = tfs->make<TH1D>("h_pionpwtrk_length","Pion + Length; Pion^{+} Truth Length (cm)",40,0,100);
  h_pionmwtrk_length = tfs->make<TH1D>("h_pionmwtrk_length","Pion - Length; Pion^{-} Truth Length (cm)",40,0,100);

  h_muon_length->Sumw2();
  h_muonwtrk_length->Sumw2();
  h_proton_length->Sumw2();
  h_protonwtrk_length->Sumw2();
  h_pionp_length->Sumw2();
  h_pionpwtrk_length->Sumw2();
  h_pionm_length->Sumw2();
  h_pionmwtrk_length->Sumw2();

  h_Pkaon_den = tfs->make<TH1D>("h_Pkaon_den","Kaon; Kaon Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Pkaon_num = tfs->make<TH1D>("h_Pkaon_num","Kaon; Kaon Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Pmichel_e_den = tfs->make<TH1D>("h_Pmichel_e_den","Michel Electron; Michele e Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Pmichel_e_num = tfs->make<TH1D>("h_Pmichel_e_num","Michel Electron; Michele e Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Pkaon_den->Sumw2();
  h_Pkaon_num->Sumw2();
  h_Pmichel_e_den->Sumw2(); 
  h_Pmichel_e_num->Sumw2(); 
  h_Efrac_kaon = tfs->make<TH1D>("h_Efrac_kaon","Efrac Kaon; Track Purity;",60,0,1.2);
  h_Ecomplet_kaon = tfs->make<TH1D>("h_Ecomplet_kaon","Ecomplet Kaon; Track Completeness;",60,0,1.2);
  h_trackRes_kaon = tfs->make<TH1D>("h_trackRes_kaon","Kaon Residual; Truth length - Reco length (cm);",200,-100,100);
  h_Efrac_michel = tfs->make<TH1D>("h_Efrac_michel","Efrac Michel; Track Energy fraction;",60,0,1.2);
  h_Ecomplet_michel = tfs->make<TH1D>("h_Ecomplet_michel","Ecomplet Michel; Track Completeness;",60,0,1.2);
  h_trackRes_michel = tfs->make<TH1D>("h_trackRes_michel","Michel Residual; Truth length - Reco length (cm);",200,-100,100);
  h_kaon_length = tfs->make<TH1D>("h_kaon_length","Kaon Length; Kaon Truth Length (cm)",40,0,100);
  h_kaonwtrk_length = tfs->make<TH1D>("h_kaonwtrk_length","Kaon Length; Kaon Truth Length (cm)",40,0,100);
  h_michel_length = tfs->make<TH1D>("h_michel_length","Michel Length; Michel e Truth Length (cm)",40,0,100);
  h_michelwtrk_length = tfs->make<TH1D>("h_michelwtrk_length","Michel Length; Michel e Truth Length (cm)",40,0,100);

  h_Efrac_kaon->Sumw2();
  h_Ecomplet_kaon->Sumw2();
  h_trackRes_kaon->Sumw2();
  h_Efrac_michel->Sumw2();
  h_Ecomplet_michel->Sumw2();
  h_trackRes_michel->Sumw2();
  h_kaon_length->Sumw2();
  h_kaonwtrk_length->Sumw2();
  h_michel_length->Sumw2();
  h_michelwtrk_length->Sumw2();
*/

 ftracks = tfs->make<TTree>("tracks", "tracks"); //ntuple
 ftracks->Branch("ev",&fev);
 ftracks->Branch("run",&frun);
 ftracks->Branch("subrun",&fsubrun);
 ftracks->Branch("nu_E",&fnu_E);
 ftracks->Branch("nu_x",&fnu_x);
 ftracks->Branch("nu_y",&fnu_y);
 ftracks->Branch("nu_z",&fnu_z);
 ftracks->Branch("count_pr",&fcount_pr);
 ftracks->Branch("particle_count",&fparticle_count);
 ftracks->Branch("FV",&fFV);
 ftracks->Branch("n_particles",&fn_particles);
 ftracks->Branch("n_fake_reco",&fn_fake_reco);


    //info on the MC particle
    ftracks->Branch("length",&flength);
    ftracks->Branch("startT",&fstartT);
    ftracks->Branch("start_x",&fstart_x);
    ftracks->Branch("start_y",&fstart_y);
    ftracks->Branch("start_z",&fstart_z);
    ftracks->Branch("end_x",&fend_x);
    ftracks->Branch("end_y",&fend_y);
    ftracks->Branch("end_z",&fend_z);
    ftracks->Branch("steps",&fn_steps);
    ftracks->Branch("pdg",&fpdg);
    ftracks->Branch("mother_id",&fmother_id);
   ftracks->Branch("g4_id",&fg4_id); //first index is t
   ftracks->Branch("p0",&fp0);//initial momentum
   ftracks->Branch("p0x",&fp0x);
   ftracks->Branch("p0y",&fp0y);
   ftracks->Branch("p0z",&fp0z);
   ftracks->Branch("kinE",&fkinE);
   ftracks->Branch("costheta_muon",&fcostheta_muon_reco);
   ftracks->Branch("is_leading",&fis_leading);

   //info coming from the tracking algorithm - when there is mc truth
   ftracks->Branch("is_tracked",&fis_tracked);
   ftracks->Branch("is_mismatched",&fis_mismatched);
   ftracks->Branch("length_reco",&flength);
   ftracks->Branch("reco_momentum_mcs",&freco_momentum_mcs); //(GeV) MCS
   ftracks->Branch("reco_momentum_mcs_llhd",&freco_momentum_mcs_llhd); //(GeV) MCS LLHD
   ftracks->Branch("reco_momentum_range",&freco_momentum_range); //MeV
   ftracks->Branch("purity",&fpurity);
   ftracks->Branch("completeness",&fcompleteness);
   ftracks->Branch("nhits",&fnhits);
   ftracks->Branch("reco_kinE",&freco_kinE);
   
   //info coming from the tracking algorithm - when there is NO mc truth
   ftracks->Branch("fake_is_tracked",&ffake_is_tracked);
   ftracks->Branch("fake_is_mismatched",&ffake_is_mismatched);
   ftracks->Branch("fake_length_reco",&ffake_length_reco);
   ftracks->Branch("fake_reco_momentum_mcs",&ffake_reco_momentum_mcs); //(GeV) MCS
   ftracks->Branch("fake_reco_momentum_mcs_llhd",&ffake_reco_momentum_mcs_llhd); //(GeV) MCS LLHD
   ftracks->Branch("fake_reco_momentum_range",&ffake_reco_momentum_range); //MeV
   ftracks->Branch("fake_purit",&ffake_purity);
   ftracks->Branch("fake_completeness",&ffake_completeness);
   ftracks->Branch("fake_nhits",&ffake_nhits);
   ftracks->Branch("fake_kinE",&ffake_kinE);
  
   std::cout << "TTree BOOKED." << std::endl;
}
//========================================================================
void ProtonReco::endJob(){
  //doEfficiencies();
  std::cout << "END!" << std::endl;

}
//========================================================================
void ProtonReco::beginRun(const art::Run& /*run*/){
  mf::LogInfo("ProtonReco")<<"begin run..."<<std::endl;
}
//========================================================================
void ProtonReco::analyze( const art::Event& event ){
    bool isFiducial = false;
cout << "ANALYZE" << endl;    
    if (event.isRealData()) 
	    processData(event, isFiducial) ;
    else
	    processMC(event, isFiducial);
}
//========================================================================
void ProtonReco::processMC( const art::Event& event, bool &isFiducial){
//	std::cout << "processing MC...." << endl;

    art::ServiceHandle<sim::LArG4Parameters> larParameters;
 //   std::cout << "electrons to GeV " <<  1./larParameters->GeVToElectrons() <<std::endl;

    art::ServiceHandle<cheat::BackTracker> bt;

    //!save neutrino's interaction info
    //neutrino frames
    art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> MCeventlist;
    if (event.getByLabel(fMCTruthModuleLabel, MCtruthHandle))
    art::fill_ptr_vector(MCeventlist, MCtruthHandle);

    //wire truth information
    art::Handle<std::vector<sim::SimChannel>> simChannelHandle;
    std::vector<art::Ptr<sim::SimChannel> > channellist;
    if (event.getByLabel(fMCTruthModuleLabel ,simChannelHandle)) 
    art::fill_ptr_vector(channellist, simChannelHandle);

    //mc track info
    art::Handle< std::vector<sim::MCTrack> > mctrackh;
    event.getByLabel(fMCTrackModuleLabel, mctrackh);
    //map of MC particles
    std::map< int, const simb::MCParticle* > particleMap;


 //loop su frame di neutrino
 for ( auto const& neutrino_event : (MCeventlist) ) {
    
      if( !neutrino_event->NeutrinoSet() ) //if there is no nu info, there was an error
		 std::cout << "NO NEUTRINO INFO!! ERROR!!!!" << std::endl;
      
      simb::MCNeutrino nu = neutrino_event->GetNeutrino();
      if( nu.CCNC() != 0 ) continue; //we want only CC events!!
      simb::MCParticle neutrino = nu.Nu(); //get the neutrino
      
      if (nu.Nu().PdgCode()!=fNeutrinoPDGcode) continue; //we want only nu mu

      //const TLorentzVector& nu_momentum = neutrino.Momentum(0); //store initial momentum
      //std::cout << nu_momentum.X() <<std::endl;
      //const TLorentzVector& vertex =neutrino.Position(0);  //store initial position -> si puo' salvare tutti gli steps actually
      //fnu_x = vertex.X();
      //fnu_y = vertex.Y();
      //fnu_z = vertex.Z();

     // * info on vertexes
      art::Handle< std::vector<recob::Vertex> > vertexListHandle;
      std::vector<art::Ptr<recob::Vertex> >  vertexlist;
      if ( event.getByLabel(fVertexModuleLabel , vertexListHandle) )
	      art::fill_ptr_vector(vertexlist, vertexListHandle);

     
      //loop on vertexes
      for ( auto const & vxt : vertexlist) {
	double xyz[3];
	vxt->XYZ(xyz);
	h_vertex_distribution->Fill(xyz[0], xyz[1], xyz[2]);
      }

      frun = event.run();
      fsubrun = event.subRun();
      fev = event.event();

      //GEANT particle info 
      const sim::ParticleList& plist = bt->ParticleList();
      
      double muon_p = 0;
      double muon_px = 0;
      double muon_py = 0;
      double muon_pz = 0;
      double muon_startx = 0;
      double muon_starty = 0;
      double muon_startz = 0;

      int n_mc_particles = 0;
      for ( auto const& geant_particle : plist ) { //look for the muon first
	if (geant_particle.second->PdgCode()!=13) continue;
		if (muon_p!=0) {
			cout << "ERROR!!!! Muon should be 1!!!!" << endl;
		}
	muon_p = geant_particle.second->Momentum().Rho();
	muon_px = geant_particle.second->Momentum().X();
	muon_py = geant_particle.second->Momentum().Y();
	muon_pz = geant_particle.second->Momentum().Z();
	muon_startx = geant_particle.second->Position().X();
	muon_starty = geant_particle.second->Position().Y();
	muon_startz = geant_particle.second->Position().Z();
#if DEBUG == 1
	cout << " p " << muon_p << endl;
	cout << " modulus " << sqrt(muon_px*muon_px + muon_py*muon_py + muon_pz*muon_pz) << endl;
#endif
      }

      long total = 0 ;
      for ( auto const& geant_particle : plist ) {
	  
	  //save only GENIE output and not the secondaries
	  if (geant_particle.second->Mother() > 0) continue;
	  if (geant_particle.second->StatusCode()!=1 ) continue;
	
	  n_mc_particles++;
#if DEBUG == 1
	      std::cout << ">>>>> Particle # " << geant_particle.first <<std::endl;
	      std::cout << "Track ID " << geant_particle.second->TrackId() <<std::endl;
	      std::cout << "PDG " << geant_particle.second->PdgCode() << " status " << geant_particle.second->StatusCode() <<std::endl;
	      std::cout << "Mother " << geant_particle.second->Mother() << " " << geant_particle.second->Trajectory().size() <<std::endl;
	      std::cout << "N tr points " << geant_particle.second->NumberTrajectoryPoints() << " trajectory points " << geant_particle.second->Trajectory().size() <<std::endl;
	      std::cout << "Start X " << geant_particle.second->Position().X() << " Start Y " << geant_particle.second->Position().Y() << " Start Z " << geant_particle.second->Position().Z() <<std::endl;
	      std::cout << "End X " << geant_particle.second->EndPosition().X() << " End Y " << geant_particle.second->EndPosition().Y() << " End Z " << geant_particle.second->EndPosition().Z() <<std::endl;
	      std::cout << "Start PX " << geant_particle.second->Momentum().X() << " Start PY " << geant_particle.second->Momentum().Y() << " Start PZ " << geant_particle.second->Momentum().Z() <<std::endl;
	      std::cout << "End PX " << geant_particle.second->EndMomentum().X() << " End PY " << geant_particle.second->EndMomentum().Y() << " End PZ " << geant_particle.second->EndMomentum().Z() <<std::endl;
#endif
	      total+=geant_particle.second->NumberTrajectoryPoints();
	//count particle types
	if (fparticle_count.find( geant_particle.second->PdgCode() ) == fparticle_count.end()) //not found
	     fparticle_count[ geant_particle.second->PdgCode() ] = 1;
	else
	     fparticle_count[geant_particle.second->PdgCode()] = fparticle_count[geant_particle.second->PdgCode()] +1;
	
	flength.push_back( geant_particle.second->Trajectory().TotalLength() );
	fstartT.push_back( geant_particle.second->T() );
	fstart_x.push_back ( geant_particle.second->Position().X() );
	fstart_y.push_back ( geant_particle.second->Position().Y() );
	fstart_z.push_back ( geant_particle.second->Position().Z() );
	fend_x.push_back ( geant_particle.second->EndPosition().X() );
	fend_y.push_back ( geant_particle.second->EndPosition().Y() );
	fend_z.push_back ( geant_particle.second->EndPosition().Z() );
	fn_steps.push_back ( geant_particle.second->Trajectory().size() );
	fmother_id.push_back ( geant_particle.second->Mother() );
	fpdg.push_back( geant_particle.second->PdgCode() );	
	fg4_id.push_back(geant_particle.second->TrackId());
	fp0.push_back( geant_particle.second->Momentum().Rho());
	fp0x.push_back( geant_particle.second->Momentum().X() );
	fp0y.push_back( geant_particle.second->Momentum().Y() );
	fp0z.push_back( geant_particle.second->Momentum().Z() );
	
	//loop on other tracked_particles and decide if this is the leading based on initial kinetic energy. This must be done before filling the kinE vector for the current particle!
	bool is_leading = true;
	float current_kinE = geant_particle.second->E() - geant_particle.second->Mass() ;
	for (unsigned jj=0; jj< fkinE.size(); jj++) { //loop on previous particles
	if ( fpdg[jj] == geant_particle.second->PdgCode() ) {
	    if ( fkinE[jj] > current_kinE )
	       is_leading = false;
	    else if ( fkinE[jj] < current_kinE )
	       fis_leading[jj] = false;
		}
	}
	fis_leading.push_back( is_leading );

	//now write current kin E
	fkinE.push_back( geant_particle.second->E() - geant_particle.second->Mass() );

	if ( geant_particle.second->Momentum().Rho() != 0) {
	     fcostheta_muon.push_back( muon_px*geant_particle.second->Momentum().X()/geant_particle.second->Momentum().Rho()/muon_p
			     	+ muon_py*geant_particle.second->Momentum().Y()/geant_particle.second->Momentum().Rho()/muon_p
				+ muon_pz*geant_particle.second->Momentum().Z()/geant_particle.second->Momentum().Rho()/muon_p );
	} else {
	     fcostheta_muon.push_back(-2);
	}
	
      
	//add dummy entries for the reco variables. They will be possibly updated later, if a matching reco track is found
	fis_tracked.push_back(0);
	fis_mismatched.push_back(0);
	flength.push_back(-1);
	freco_momentum_mcs.push_back(-1);
	freco_momentum_mcs_llhd.push_back(-1);
	freco_momentum_range.push_back(-1);
	fpurity.push_back(-1);
	fcompleteness.push_back(-1);
	fnhits.push_back(-1);
	freco_kinE.push_back(-1);
      
      } //for all geant_particle
      fn_particles = flength.size();

      std::cout << "TOTAL NUMBER OF MCparticles " << total << std::endl;
      //go with matching with reco track!
    art::Handle<std::vector<recob::Track>> TrackHandle;
    std::vector<art::Ptr<recob::Track> > reco_tracks;
    if (event.getByLabel(fTrackingAlgorithmLabel ,TrackHandle)) 
    art::fill_ptr_vector( reco_tracks, TrackHandle);

    art::FindManyP<recob::Hit>  fmht(TrackHandle, event, fTrackingAlgorithmLabel);
  
    art::FindManyP<anab::Calorimetry> fmcal( TrackHandle, event, fCalorimetryModuleLabel); //handle for calorimetry


    for (unsigned n_track = 0; n_track < reco_tracks.size(); n_track++ ) { //sort of loop on tracks
   
    std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(n_track);
    std::vector<art::Ptr<recob::Hit> > hits_planes[3]; //vector of hits for the 3 planes
	
    Int_t trkid_truth[3];
    Float_t purity_truth[3];
    //int pdg_truth[3];

    	for (auto & hit : allHits) {
   
	   if ( hit->WireID().Plane <  3){
		  hits_planes[hit->WireID().Plane].push_back(hit);
		} else {
			std::cout << "ERROR! Weird plane number " << hit->WireID().Plane << std::endl;
	   }
	}
	double maxe = 0;
	
	for (size_t ipl = 0; ipl < 3; ++ipl){ //loop on all planes, evaluate purity on single planes
	maxe = 0;
	
	AnaTree_HitsPurity(hits_planes[ipl], trkid_truth[ipl], purity_truth[ipl], maxe);
	   
		if ( trkid_truth[ipl] > 0 ){ //if there is matching
		double tote = 0;
		std::vector<sim::IDE> vide( bt->TrackIDToSimIDE( trkid_truth[ipl] ) );
		for (const sim::IDE& ide: vide) {
			tote += ide.energy;
	//	std::cout << "ide energy " << ide.energy << " plane " <<  ipl << std::endl;
		}
		
		//pdg_truth[ipl] = particle->PdgCode();
		purity_truth[ipl] = maxe/(tote/3); //average on the single plane
	//	std::cout << "purity plane " << ipl << " " <<  purity_truth[ipl] << std::endl;
		}
	}
	
	maxe = 0; //evaluate global purity and completeness
	Int_t trk_g4id = 0;
	Float_t trk_purity = 0;
	double trk_completeness = 0;
	AnaTree_HitsPurity(allHits, trk_g4id, trk_purity, maxe);
   
	float totenergy = 0;
	for (auto & hit : allHits) {
	std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(hit);
		for(size_t e = 0; e < eveIDs.size(); ++e){
	    		if (eveIDs[e].trackID==trk_g4id) totenergy += eveIDs[e].energy;
		}
	}

	bool mismatch = false;
	for (size_t ipl = 0; ipl < 3; ++ipl) {
		if (trkid_truth[ipl]!=trk_g4id) mismatch = true;
	}

	if (totenergy) trk_completeness = maxe/totenergy;
	
	//fill info for this truck in the vectors
	
	Int_t index_MC_particle = -1;
	//look for the MC truth particle
	for (Int_t i=0; i<fn_particles; i++) {
	if (fg4_id[i] == trk_g4id && index_MC_particle != -1 ) std::cout << "more than 1 matching particle!! ERROR!!!" << std::endl;
	if (fg4_id[i] == trk_g4id && index_MC_particle == -1) index_MC_particle = i;
	}

	if (index_MC_particle==-1) { //the MC particle is not found (!)
	ffake_is_tracked.push_back(true);
	ffake_is_mismatched.push_back(mismatch);
        ffake_length_reco.push_back( reco_tracks[n_track]->Trajectory().Length() );

	trkf::TrackMomentumCalculator trkm; //track momentum calculator
	trkm.SetMinLength(0); //minimum track length for momentum calculation
	ffake_reco_momentum_mcs.push_back( trkm.GetMomentumMultiScatterChi2( reco_tracks[n_track] ));
	ffake_reco_momentum_mcs_llhd.push_back( trkm.GetMomentumMultiScatterLLHD( reco_tracks[n_track] ));
	ffake_reco_momentum_range.push_back( trkm.GetTrackMomentum( reco_tracks[n_track]->Trajectory().Length(), 13 )); //if you don't know the particle, assume it's a muon (??? probably stupid, but...)
	ffake_purity.push_back( trk_purity );
	ffake_completeness.push_back( trk_completeness );

	//calorimetry info
	if (fmcal.isValid()){ 
	std::vector< art::Ptr<anab::Calorimetry> > calos = fmcal.at(n_track);
	ffake_kinE.push_back( calos[geo::kCollection]->KineticEnergy() );
	ffake_nhits.push_back( calos[geo::kCollection]->dEdx().size() );
	}

	} else { //particle is found!

  	fis_tracked[index_MC_particle] = 1;
	flength[index_MC_particle] = reco_tracks[n_track]->Trajectory().Length();
	trkf::TrackMomentumCalculator trkm; //track momentum calculator
	trkm.SetMinLength(0); //minimum track length for momentum calculation
	freco_momentum_mcs[index_MC_particle] = trkm.GetMomentumMultiScatterChi2( reco_tracks[n_track] );
	freco_momentum_mcs_llhd[index_MC_particle] = trkm.GetMomentumMultiScatterLLHD( reco_tracks[n_track] ) ;
	freco_momentum_range[index_MC_particle] = trkm.GetTrackMomentum( reco_tracks[n_track]->Trajectory().Length(), fpdg[index_MC_particle] ); //use info on pdg
	fpurity[index_MC_particle] = trk_purity;
	fcompleteness[index_MC_particle] = trk_completeness;
	fis_mismatched[index_MC_particle] = mismatch;
	//calorimetry
	if (fmcal.isValid()){
	std::vector< art::Ptr<anab::Calorimetry> > calos = fmcal.at(n_track);
	fnhits[index_MC_particle] = calos[geo::kCollection]->dEdx().size();
	freco_kinE[index_MC_particle] = calos[geo::kCollection]->KineticEnergy();
	}
      
      } //save reconstruction info
    
    //loop on tracks

	//relative to the beginning of the track!
	h_track_length->Fill( reco_tracks[n_track]->Trajectory().Length() );
	h_direction_y->Fill(reco_tracks[n_track]->Trajectory().ZenithAngle() );
	h_direction_theta->Fill(reco_tracks[n_track]->Trajectory().Theta() );
	h_direction_phi->Fill(reco_tracks[n_track]->Trajectory().Phi() );
    }
	
    //now loop on MC particles and try to characterize non-reco tracks (position to 

      if ( unsigned(plist.size()) != flength.size() ) {
	      //std::cout << "plist.size() " << plist.size()  << " flength.size()" << " " << flength.size() << std::endl;
	      //std::cout << "ERROR!!!! Mismatch in number of MC particles and recorded MC particles" << std::endl;
	      //exit(-1);
      }

      for ( unsigned j=0; j<flength.size(); j++) {
	if ( fis_tracked[j] ) continue; //if is tracked, it's not interesting :)
	
	if (fpdg[j] != 2212) {
		std::cout << "Something different than a proton!" << std::endl;
		continue;
	}
	
	double dis1 = sqrt( pow( fstart_x[j]-muon_startx,2) + pow( fstart_y[j]-muon_starty,2) + pow( fstart_z[j]-muon_startz,2)); 
	h_dis_start_start->Fill(dis1);	
	double costheta = 1./muon_p * 1./fp0[j] * ( muon_px*fp0x[j] + muon_py*fp0y[j] + muon_pz*fp0z[j]);
	h_cos_recotruth ->Fill(costheta);

      }

      ftracks->Fill();
      clear_vectors();
 }//for all MC events     
      
 }

//========================================================================
void ProtonReco::processData( const art::Event& event, bool &isFiducial){
    art::ServiceHandle<sim::LArG4Parameters> larParameters;

    // Declare object-ID-to-PFParticleID maps so we can assign hasPFParticle and PFParticleID to the tracks, showers, vertices.
//    std::map<Short_t, Short_t> trackIDtoPFParticleIDMap, vertexIDtoPFParticleIDMap, showerIDtoPFParticleIDMap;

// * trigger information
      art::Handle< std::vector<raw::Trigger> > triggerListHandle;
      std::vector<art::Ptr<raw::Trigger> > triggerlist;
      if (event.getByLabel(fDigitModuleLabel, triggerListHandle))
	art::fill_ptr_vector(triggerlist, triggerListHandle);
/*
      if (triggerlist.size()){
	      std:: cout << "trigger number " << triggerlist[0]->TriggerNumber() << std::endl;
	      std::cout << "trigger time " << triggerlist[0]->TriggerTime() << std::endl;
	      std::cout << "trigger beam gate time " << triggerlist[0]->BeamGateTime() << std::endl;
	      std::cout << "trigger trigger bits " << triggerlist[0]->TriggerBits() << std::endl;
	      //print out trigger type
	      if ( triggerlist[0]-> Triggered(trigger::kTriggerBNB))
		      std::cout << "BNB trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kPMTTriggerBeam))
		      std::cout << "PMT beam trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kPMTTriggerCosmic))
		      std::cout << "PMT cosmic trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kTriggerPC))
		      std::cout << "PC trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kTriggerEXT))
		      std::cout << "EXT trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kTriggerNuMI))
		      std::cout << "NuMI trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kFakeBeam))
		      std::cout << "Fake Beam trigger!" << std::endl;
	      if ( triggerlist[0]-> Triggered(trigger::kFakeGate))
		      std::cout << "Fake Gate trigger!" << std::endl;
      }
*/

//software trigger
  //SW Trigger
      art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
      event.getByLabel("swtrigger", softwareTriggerHandle);
      if (softwareTriggerHandle.isValid()){ 
         if (softwareTriggerHandle->getNumberOfAlgorithms() == 1) {
             std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
             std::cout << "SW trigger name: " << algoNames[0] << " passed " << softwareTriggerHandle->passedAlgo(algoNames[0]) <<std::endl;
           // _is_swtriggered = (softwareTriggerHandle->passedAlgo(algoNames[0]) ? 1 : 0);
	   }
	 } else {
      cout << "software trigger is INVALID" << endl;
	}


    // PFParticle
    // * PFParticles
    lar_pandora::PFParticleVector pfparticlelist;
    lar_pandora::PFParticlesToClusters pfParticleToClusterMap;
    lar_pandora::LArPandoraHelper::CollectPFParticles(event, fPFParticleModuleLabel, pfparticlelist, pfParticleToClusterMap);

    lar_pandora::PFParticleVector pfparticle_finalstatelist;
    lar_pandora::LArPandoraHelper::SelectFinalStatePFParticles (pfparticlelist, pfparticle_finalstatelist);

    if ( pfparticle_finalstatelist.size() != pfparticlelist.size() )
	    cout << "WARNING!! different size: total = " << pfparticlelist.size() << " final state: " << pfparticle_finalstatelist.size() << endl;

    /*for ( auto const & pfparticle : pfparticle_finalstatelist ) { //loop on pf particles
    }*/

    h_n_pfparticles->Fill(pfparticle_finalstatelist.size());

     // * info on vertexes
      art::Handle< std::vector<recob::Vertex> > vertexListHandle;
      std::vector<art::Ptr<recob::Vertex> >  vertexlist;
      if ( event.getByLabel(fVertexModuleLabel , vertexListHandle) )
	      art::fill_ptr_vector(vertexlist, vertexListHandle);
	
      h_n_vertexes->Fill(vertexlist.size());
     
      //loop on vertexes
      for ( auto const & vxt : vertexlist) {
	double xyz[3];
	vxt->XYZ(xyz);
	h_vertex_distribution->Fill(xyz[0], xyz[1], xyz[2]);
      }

      frun = event.run();
      fsubrun = event.subRun();
      fev = event.event();

    art::Handle<std::vector<recob::Track>> TrackHandle;
    std::vector<art::Ptr<recob::Track> > reco_tracks;
    if (event.getByLabel(fTrackingAlgorithmLabel ,TrackHandle)) 
    art::fill_ptr_vector( reco_tracks, TrackHandle);

    art::FindManyP<recob::Hit>  fmht(TrackHandle, event, fTrackingAlgorithmLabel);
  
    art::FindManyP<anab::Calorimetry> fmcal( TrackHandle, event, fCalorimetryModuleLabel); //handle for calorimetry

    h_n_tracks->Fill(reco_tracks.size());
	
    long nhits = 0;

    for (unsigned n_track = 0; n_track < reco_tracks.size(); n_track++ ) { //sort of loop on tracks
   
    std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(n_track);
    std::vector<art::Ptr<recob::Hit> > hits_planes[3]; //vector of hits for the 3 planes
	
    h_n_hits_track->Fill(allHits.size());
    nhits+= allHits.size();
    
    //loop on tracks

	//relative to the beginning of the track!
	h_track_length->Fill( reco_tracks[n_track]->Trajectory().Length() );
	h_direction_y->Fill(reco_tracks[n_track]->Trajectory().ZenithAngle() );
	h_direction_theta->Fill(reco_tracks[n_track]->Trajectory().Theta() );
	h_direction_phi->Fill(reco_tracks[n_track]->Trajectory().Phi() );
    }

    h_n_hits->Fill(nhits);

 }

/*
    //!save FS particles

    double tmp_leadingPionPlusE = 0.0;
    double tmp_leadingPionMinusE = 0.0;
    double tmp_leadingProtonE  = 0.0; 
  
    simb::MCParticle *MClepton = NULL; 
    simb::MCParticle *MCproton = NULL;
    simb::MCParticle *MCpion_plus = NULL;
    simb::MCParticle *MCpion_minus = NULL;
    simb::MCParticle *MCkaon = NULL;
    simb::MCParticle *MCmichel = NULL;
 
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index

    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
       particle = ipar->second;
       if( particle->PdgCode() == fLeptonPDGcode && particle->Mother() == 0 ){  //primary lepton
         const TLorentzVector& lepton_momentum =particle->Momentum(0);
         lepton_momentum.GetXYZT(MC_lepton_startMomentum);
         MC_leptonID = particle->TrackId();
         MC_leptonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
         MClepton = particle;
       }
       if( particle->Mother() == 0 ){   //save primary particle i.e. from the neutrino interaction
         //save leading pion and proton
         if( particle->PdgCode() == 2212 ){
           if(particle->Momentum().E() > tmp_leadingProtonE){
             tmp_leadingProtonE = particle->Momentum().E();
             MC_leading_protonID = particle->TrackId();          
             MC_leading_ProtonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             MCproton = particle;
           } 
         }
         else if( particle->PdgCode() == 211 ){
           if(particle->Momentum().E() > tmp_leadingPionPlusE){
             tmp_leadingPionPlusE = particle->Momentum().E();
             MC_leading_PionPlusID = particle->TrackId();          
             MC_leading_PionPlusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             MCpion_plus = particle;
           } 
         }
         else if( particle->PdgCode() == -211 ){
           if(particle->Momentum().E() > tmp_leadingPionMinusE){
             tmp_leadingPionMinusE = particle->Momentum().E();
             MC_leading_PionMinusID = particle->TrackId();          
             MC_leading_PionMinusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             MCpion_minus = particle;
           } 
         }
         i++; //paticle index
       }
       
       //=======================================================================================
       //add Nucleon decay stuff and particle cannon 
       //=======================================================================================
       if(!fisNeutrinoInt ){
         if( particle->Mother() == 0 ){
           const TLorentzVector& positionStart = particle->Position(0); 
           positionStart.GetXYZT(MC_vertex); 
         }
         if( particle->PdgCode() == 321 ){   //save primary Kaon
           MC_kaonID = particle->TrackId();    
           MC_kaonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
           MCkaon = particle;
         }   
         else if( particle->PdgCode() == fLeptonPDGcode ){  // Particle cannon muon
           MC_leptonID = particle->TrackId();
           MC_leptonP =  sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2)); 
           MClepton = particle;
         }
         else if( particle->PdgCode() == 2212 ){
           if(particle->Momentum().E() > tmp_leadingProtonE){
             tmp_leadingProtonE = particle->Momentum().E();
             MC_leading_protonID = particle->TrackId();          
             MC_leading_ProtonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             MCproton = particle;
           } 
         }
         else if( particle->PdgCode() == 211 ){
           if(particle->Momentum().E() > tmp_leadingPionPlusE){
             tmp_leadingPionPlusE = particle->Momentum().E();
             MC_leading_PionPlusID = particle->TrackId();          
             MC_leading_PionPlusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             MCpion_plus = particle;
           } 
         }
         else if( particle->PdgCode() == -211 ){
           if(particle->Momentum().E() > tmp_leadingPionMinusE){
             tmp_leadingPionMinusE = particle->Momentum().E();
             MC_leading_PionMinusID = particle->TrackId();          
             MC_leading_PionMinusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             MCpion_minus = particle;
           } 
         }
         else if( particle->Process() =="Decay" && particle->PdgCode() == -11){  // michel electron from muon decay
           MC_michelID = particle->TrackId();
           MC_michelP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
           MCmichel = particle;
         }

       }
    } 
    //===================================================================
    //Saving denominator histograms  
    //===================================================================
    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
    double Pv  = sqrt(pow(MC_incoming_P[0],2)+pow(MC_incoming_P[1],2)+pow(MC_incoming_P[2],2));
    double theta_mu = acos((MC_incoming_P[0]*MC_lepton_startMomentum[0] + MC_incoming_P[1]*MC_lepton_startMomentum[1] +MC_incoming_P[2]*MC_lepton_startMomentum[2])/(Pv*MC_leptonP) );
    theta_mu *= (180.0/3.14159);
    double truth_lengthLepton = truthLength(MClepton); 
    double proton_length = truthLength(MCproton);
    double pion_plus_length = truthLength(MCpion_plus);
    double pion_minus_length = truthLength(MCpion_minus);
    double kaonLength = truthLength(MCkaon);
    double michelLength = truthLength(MCmichel);

    //save CC events within the fiducial volume with the favorite neutrino flavor 
    if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
       if( MClepton ){
         h_Ev_den->Fill(MC_incoming_P[3]);
         h_Pmu_den->Fill(MC_leptonP);
         h_theta_den->Fill(theta_mu);
	 h_muon_length->Fill(truth_lengthLepton);
       }
       if( MCproton ){
         h_Pproton_den->Fill(MC_leading_ProtonP);
         h_proton_length->Fill(proton_length);
       }
       if( MCpion_plus ){
         h_Ppion_plus_den->Fill( MC_leading_PionPlusP);
         h_pionp_length->Fill(pion_plus_length);
       }
       if( MCpion_minus ){
         h_Ppion_minus_den->Fill( MC_leading_PionMinusP);
         h_pionm_length->Fill(pion_minus_length);
       }
    } 
  
    //save events for Nucleon decay and particle cannon
    if(!fisNeutrinoInt ){
      if( MClepton ){
         h_Pmu_den->Fill(MC_leptonP);
	 h_muon_length->Fill(truth_lengthLepton);
       }
       if( MCkaon ){
         h_Pkaon_den->Fill(MC_kaonP);
         h_kaon_length->Fill(kaonLength);
       }
       if( MCproton ){
         h_Pproton_den->Fill(MC_leading_ProtonP);
         h_proton_length->Fill(proton_length);
       }
       if( MCpion_plus ){
         h_Ppion_plus_den->Fill( MC_leading_PionPlusP);
         h_pionp_length->Fill(pion_plus_length);
       }
       if( MCpion_minus ){
         h_Ppion_minus_den->Fill( MC_leading_PionMinusP);
         h_pionm_length->Fill(pion_minus_length);
       }
       if( MCmichel ){
         h_Pmichel_e_den->Fill(MC_michelP);
	 h_michel_length->Fill(michelLength);
       }
    }
 
    //========================================================================
    // Reco stuff, once we have selected a MC Particle let's find out if there is a track associated 
    //========================================================================
    art::Handle<std::vector<recob::Track>> trackListHandle;
    if(! event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track> > tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);
    int n_recoTrack = tracklist.size();
   
    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    if( n_recoTrack == 0 ){
      LOG_DEBUG("ProtonReco")<<"There are no reco tracks... bye";
      return; 
    }
    LOG_DEBUG("ProtonReco")<<"Found this many reco tracks "<<n_recoTrack;
   
 
    double Efrac_lepton =0.0;
    double Ecomplet_lepton =0.0;
    double Efrac_proton =0.0;
    double Ecomplet_proton =0.0;
    double Efrac_pionplus =0.0;
    double Ecomplet_pionplus =0.0;
    double Efrac_pionminus =0.0;
    double Ecomplet_pionminus =0.0;
    double Efrac_kaon =0.0;
    double Ecomplet_kaon =0.0;
    double Efrac_michel =0.0;
    double Ecomplet_michel =0.0;
    double trackLength_lepton =0.0;
    double trackLength_proton =0.0;
    double trackLength_pion_plus =0.0;
    double trackLength_pion_minus =0.0;
    double trackLength_kaon =0.0;
    double trackLength_michel =0.0;
    const simb::MCParticle *MClepton_reco = NULL; 
    const simb::MCParticle *MCproton_reco = NULL;
    const simb::MCParticle *MCpion_plus_reco = NULL;
    const simb::MCParticle *MCpion_minus_reco = NULL;
    const simb::MCParticle *MCkaon_reco = NULL;
    const simb::MCParticle *MCmichel_reco = NULL;

    std::vector<art::Ptr<recob::Hit>> tmp_all_trackHits = track_hits.at(0);  
    std::vector<art::Ptr<recob::Hit>> all_hits;
    art::Handle<std::vector<recob::Hit>> hithandle;
    if(event.get(tmp_all_trackHits[0].id(), hithandle))  art::fill_ptr_vector(all_hits, hithandle);
 
    for(int i=0; i<n_recoTrack; i++) {
       art::Ptr<recob::Track> track = tracklist[i];
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       double tmpEfrac = 0;
       double tmpEcomplet =0;
       const simb::MCParticle *particle;
       truthMatcher( all_hits,  all_trackHits, particle, tmpEfrac, tmpEcomplet );
       if (!particle) continue;
       //std::cout<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac<<std::endl;
       if(  (particle->PdgCode() == fLeptonPDGcode) && (particle->TrackId() == MC_leptonID) ){
         //save the best track ... based on completeness if there is more than one track 
         //if( tmpEfrac > Efrac_lepton ){ ///this was base on purity 
         if( tmpEcomplet > Ecomplet_lepton ){
           Ecomplet_lepton = tmpEcomplet;
           Efrac_lepton = tmpEfrac;
           trackLength_lepton = track->Length(); 
           MClepton_reco = particle;
         }
       }
       else if( (particle->PdgCode() == 2212) && (particle->TrackId() == MC_leading_protonID) ){
         //save the best track ... based on completeness if there is more than one track 
         if( tmpEcomplet > Ecomplet_proton ){
           Ecomplet_proton = tmpEcomplet;
           Efrac_proton = tmpEfrac;
           trackLength_proton = track->Length();
           MCproton_reco = particle;
         }
       }
       else if( (particle->PdgCode() == 211) && (particle->TrackId() == MC_leading_PionPlusID) ){
         //save the best track ... based on completeness if there is more than one track 
         if( tmpEcomplet > Ecomplet_pionplus ){
           Ecomplet_pionplus = tmpEcomplet;
           Efrac_pionplus = tmpEfrac;
           trackLength_pion_plus = track->Length();
           MCpion_plus_reco = particle;
         }
       }
       else if( (particle->PdgCode() == -211) && (particle->TrackId() == MC_leading_PionMinusID) ){
         //save the best track ... based on completeness if there is more than one track 
         if( tmpEcomplet > Ecomplet_pionminus ){
           Ecomplet_pionminus = tmpEcomplet;
           Efrac_pionminus = tmpEfrac;
           trackLength_pion_minus = track->Length();
           MCpion_minus_reco = particle;
         }
       }
       //kaon from nucleon decay
       else if( (particle->PdgCode() == 321) && (particle->TrackId() == MC_kaonID) ){
         //save the best track ... based on completeness if there is more than one track 
         if( tmpEcomplet > Ecomplet_kaon ){
           Ecomplet_kaon = tmpEcomplet;
           Efrac_kaon = tmpEfrac;
           trackLength_kaon = track->Length();
           MCkaon_reco = particle;
         }
       }
       //michel from nucleon decay
       else if( (particle->PdgCode() == -11) && (particle->TrackId() == MC_michelID) ){
         //save the best track ... based on completeness if there is more than one track 
         if( tmpEcomplet > Ecomplet_michel ){
           Ecomplet_michel = tmpEcomplet;
           Efrac_michel = tmpEfrac;
           trackLength_michel = track->Length();
           MCmichel_reco = particle;
         }
       } 

    }

    double Reco_LengthRes =  truth_lengthLepton-trackLength_lepton;
    double Reco_LengthResProton = proton_length-trackLength_proton; 
    double Reco_LengthResPionPlus = pion_plus_length-trackLength_pion_plus; 
    double Reco_LengthResPionMinus = pion_minus_length-trackLength_pion_minus;

    if( MClepton_reco && MClepton  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){ 
        h_Pmu_num->Fill(MC_leptonP);
        h_Ev_num->Fill(MC_incoming_P[3]);
        h_theta_num->Fill(theta_mu);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_Ecomplet_lepton->Fill(Ecomplet_lepton);
        h_trackRes_lepton->Fill(Reco_LengthRes);  
        h_muonwtrk_length->Fill(truth_lengthLepton);
      }
    }
    if( MCproton_reco && MCproton ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
        h_Pproton_num->Fill(MC_leading_ProtonP);     
        h_Efrac_proton->Fill(Efrac_proton);
        h_Ecomplet_proton->Fill(Ecomplet_proton);
        h_trackRes_proton->Fill(Reco_LengthResProton);       
	h_protonwtrk_length->Fill(proton_length);
      }
    }
    if( MCpion_plus_reco && MCpion_plus ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
        h_Ppion_plus_num->Fill(MC_leading_PionPlusP);     
        h_Efrac_pion_plus->Fill(Efrac_pionplus);
        h_Ecomplet_pion_plus->Fill(Ecomplet_pionplus);
        h_trackRes_pion_plus->Fill(Reco_LengthResPionPlus);
	h_pionpwtrk_length->Fill(pion_plus_length);
      }
    }
    if( MCpion_minus_reco && MCpion_minus  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ) {
        h_Ppion_minus_num->Fill(MC_leading_PionMinusP);     
        h_Efrac_pion_minus->Fill(Efrac_pionminus);
        h_Ecomplet_pion_minus->Fill(Ecomplet_pionminus);
        h_trackRes_pion_minus->Fill(Reco_LengthResPionMinus);
	h_pionmwtrk_length->Fill(pion_minus_length);
      }
    }
    //Non neutrino events 
    //=========================================================
    if(!fisNeutrinoInt ){
      if( MClepton_reco && MClepton  ){
        h_Pmu_num->Fill(MC_leptonP);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_Ecomplet_lepton->Fill(Ecomplet_lepton);
        h_trackRes_lepton->Fill(Reco_LengthRes);
	h_muonwtrk_length->Fill(truth_lengthLepton);
      }
      if( MCkaon_reco && MCkaon ){
        h_Pkaon_num->Fill(MC_kaonP);
        h_Efrac_kaon->Fill(Efrac_kaon);
        h_Ecomplet_kaon->Fill(Ecomplet_kaon);
        h_trackRes_kaon->Fill(kaonLength-trackLength_kaon);
	h_kaonwtrk_length->Fill(kaonLength);
      }
      if( MCproton_reco && MCproton ){
        h_Pproton_num->Fill(MC_leading_ProtonP);     
        h_Efrac_proton->Fill(Efrac_proton);
        h_Ecomplet_proton->Fill(Ecomplet_proton);
        h_trackRes_proton->Fill(Reco_LengthResProton);       
	h_protonwtrk_length->Fill(proton_length);
      }
      if( MCpion_plus_reco && MCpion_plus ){
        h_Ppion_plus_num->Fill(MC_leading_PionPlusP);     
        h_Efrac_pion_plus->Fill(Efrac_pionplus);
        h_Ecomplet_pion_plus->Fill(Ecomplet_pionplus);
        h_trackRes_pion_plus->Fill(Reco_LengthResPionPlus);
	h_pionpwtrk_length->Fill(pion_plus_length);
      }
      if( MCpion_minus_reco && MCpion_minus  ){
        h_Ppion_minus_num->Fill(MC_leading_PionMinusP);     
        h_Efrac_pion_minus->Fill(Efrac_pionminus);
        h_Ecomplet_pion_minus->Fill(Ecomplet_pionminus);
        h_trackRes_pion_minus->Fill(Reco_LengthResPionMinus);
	h_pionmwtrk_length->Fill(pion_minus_length);
      }
      if( MCmichel_reco && MCmichel ){
        h_Pmichel_e_num->Fill(MC_michelP);
        h_Efrac_michel->Fill(Efrac_michel);
        h_Ecomplet_michel->Fill(Ecomplet_michel);
        h_trackRes_michel->Fill(michelLength-trackLength_michel);
	h_michelwtrk_length->Fill(michelLength);
      }

    }
*/    
//}
//========================================================================
void ProtonReco::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

    //std::cout<<"truthMatcher..."<<std::endl;
    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < track_hits.size(); ++j){
       art::Ptr<recob::Hit> hit = track_hits[j];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t k = 0; k < TrackIDs.size(); k++){
          trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
       }            
    }
    double E_em =0.0;
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 
    //!since we are looking for muons/pions/protons this should be enough 
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         partial_E = ii->second;
         max_E = ii->second;
         TrackID = ii->first;
         if( TrackID < 0 ) E_em += ii->second;
       }
    } 
    MCparticle = bt->TrackIDToParticle(TrackID);

    //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of 
    //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle
    //we don't want to track gammas or any other EM activity 
    if( TrackID < 0 ) return;

    //Efrac = (partial_E+E_em)/total_E;
    Efrac = (partial_E)/total_E;

    //completeness
    double totenergy =0;
    for(size_t k = 0; k < all_hits.size(); ++k){
       art::Ptr<recob::Hit> hit = all_hits[k];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t l = 0; l < TrackIDs.size(); ++l){
          if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;
       }
    } 
    Ecomplet = partial_E/totenergy;
}
//========================================================================
double ProtonReco::truthLength( const simb::MCParticle *MCparticle ){
   //calculate the truth length considering only the part that is inside the TPC
   //Base on a peace of code from dune/TrackingAna/TrackingEfficiency_module.cc

   if( !MCparticle ) return -999.0;
   const int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
   double *TPCLengthHits = new double[numberTrajectoryPoints];
   int FirstHit=0, LastHit=0;
   double TPCLength = 0.0;
   bool BeenInVolume = false;

   for(int MCHit=0; MCHit < numberTrajectoryPoints; ++MCHit) {
      const TLorentzVector& tmpPosition= MCparticle->Position(MCHit);
      double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      if (MCHit!=0) TPCLengthHits[MCHit] = sqrt( pow( (MCparticle->Vx(MCHit-1)-MCparticle->Vx(MCHit)),2)+ pow( (MCparticle->Vy(MCHit-1)-MCparticle->Vy(MCHit)),2)+ pow( (MCparticle->Vz(MCHit-1)-MCparticle->Vz(MCHit)),2));
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if(tpcid.isValid) {
        // -- Check if hit is within drift window...
        geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
        geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
        double XPlanePosition      = tpc.PlaneLocation(0)[0];
        double DriftTimeCorrection = fabs( tmpPosition[0] - XPlanePosition ) / XDriftVelocity;
        double TimeAtPlane         = MCparticle->T() + DriftTimeCorrection; 
        if( TimeAtPlane < detprop->TriggerOffset() || TimeAtPlane > detprop->TriggerOffset() + WindowSize ) continue;
        LastHit = MCHit;
        if( !BeenInVolume ) {
	  BeenInVolume = true;
          FirstHit = MCHit;
	}
      }		
   }
   for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) TPCLength += TPCLengthHits[Hit];
   return TPCLength;
}
//========================================================================
bool ProtonReco::insideFV( double vertex[4]){ 

     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     if (x>fFidVolXmin && x<fFidVolXmax&&
	 y>fFidVolYmin && y<fFidVolYmax&&
	 z>fFidVolZmin && z<fFidVolZmax)
       return true;
     else
       return false;
}
//========================================================================
void ProtonReco::doEfficiencies(){

   art::ServiceHandle<art::TFileService> tfs;

   if(TEfficiency::CheckConsistency(*h_Ev_num,*h_Ev_den)){
     h_Eff_Ev = tfs->make<TEfficiency>(*h_Ev_num,*h_Ev_den);
     TGraphAsymmErrors *grEff_Ev = h_Eff_Ev->CreateGraph();
     grEff_Ev->Write("grEff_Ev");
     h_Eff_Ev->Write("h_Eff_Ev");
   }
   if(TEfficiency::CheckConsistency(*h_Pmu_num,*h_Pmu_den)){ 
     h_Eff_Pmu = tfs->make<TEfficiency>(*h_Pmu_num,*h_Pmu_den);
     TGraphAsymmErrors *grEff_Pmu = h_Eff_Pmu->CreateGraph();
     grEff_Pmu->Write("grEff_Pmu");
     h_Eff_Pmu->Write("h_Eff_Pmu");
   }
   if(TEfficiency::CheckConsistency(*h_theta_num,*h_theta_den)){
     h_Eff_theta = tfs->make<TEfficiency>(*h_theta_num,*h_theta_den);
     TGraphAsymmErrors *grEff_theta = h_Eff_theta->CreateGraph();
     grEff_theta->Write("grEff_theta");
     h_Eff_theta->Write("h_Eff_theta");
   }
   if(TEfficiency::CheckConsistency(*h_Pproton_num,*h_Pproton_den)){
     h_Eff_Pproton = tfs->make<TEfficiency>(*h_Pproton_num,*h_Pproton_den);
     TGraphAsymmErrors *grEff_Pproton = h_Eff_Pproton->CreateGraph();
     grEff_Pproton->Write("grEff_Pproton");
     h_Eff_Pproton->Write("h_Eff_Pproton");
   }
   if(TEfficiency::CheckConsistency(*h_Ppion_plus_num,*h_Ppion_plus_den)){
     h_Eff_Ppion_plus = tfs->make<TEfficiency>(*h_Ppion_plus_num,*h_Ppion_plus_den);
     TGraphAsymmErrors *grEff_Ppion_plus = h_Eff_Ppion_plus->CreateGraph();
     grEff_Ppion_plus->Write("grEff_Ppion_plus");
     h_Eff_Ppion_plus->Write("h_Eff_Ppion_plus");
   }
   if(TEfficiency::CheckConsistency(*h_Ppion_minus_num,*h_Ppion_minus_den)){
     h_Eff_Ppion_minus = tfs->make<TEfficiency>(*h_Ppion_minus_num,*h_Ppion_minus_den);
     TGraphAsymmErrors *grEff_Ppion_minus = h_Eff_Ppion_minus->CreateGraph();
     grEff_Ppion_minus->Write("grEff_Ppion_minus");
     h_Eff_Ppion_minus->Write("h_Eff_Ppion_minus");
   }
   if(TEfficiency::CheckConsistency(*h_muonwtrk_length,*h_muon_length)){
     h_Eff_Lmuon = tfs->make<TEfficiency>(*h_muonwtrk_length,*h_muon_length);
     TGraphAsymmErrors *grEff_Lmuon = h_Eff_Lmuon->CreateGraph();
     grEff_Lmuon->Write("grEff_Lmuon");
     h_Eff_Lmuon->Write("h_Eff_Lmuon");
   }
   if(TEfficiency::CheckConsistency(*h_protonwtrk_length,*h_proton_length)){
     h_Eff_Lproton = tfs->make<TEfficiency>(*h_protonwtrk_length,*h_proton_length);
     TGraphAsymmErrors *grEff_Lproton = h_Eff_Lproton->CreateGraph();
     grEff_Lproton->Write("grEff_Lproton");
     h_Eff_Lproton->Write("h_Eff_Lproton");
   }
   if(TEfficiency::CheckConsistency(*h_pionpwtrk_length,*h_pionp_length)){
     h_Eff_Lpion_plus = tfs->make<TEfficiency>(*h_pionpwtrk_length,*h_pionp_length);
     TGraphAsymmErrors *grEff_Lpion_plus = h_Eff_Lpion_plus->CreateGraph();
     grEff_Lpion_plus->Write("grEff_Lpion_plus");
     h_Eff_Lpion_plus->Write("h_Eff_Lpion_plus");
   }
   if(TEfficiency::CheckConsistency(*h_pionpwtrk_length,*h_pionp_length)){
     h_Eff_Lpion_minus = tfs->make<TEfficiency>(*h_pionmwtrk_length,*h_pionm_length);
     TGraphAsymmErrors *grEff_Lpion_minus = h_Eff_Lpion_minus->CreateGraph();
     grEff_Lpion_minus->Write("grEff_Lpion_minus");
     h_Eff_Lpion_minus->Write("h_Eff_Lpion_minus");
   }
   if(TEfficiency::CheckConsistency(*h_Pkaon_num,*h_Pkaon_den)){
       h_Eff_Pkaon = tfs->make<TEfficiency>(*h_Pkaon_num,*h_Pkaon_den);
       TGraphAsymmErrors *grEff_Pkaon = h_Eff_Pkaon->CreateGraph();
       grEff_Pkaon->Write("grEff_Pkaon");
       h_Eff_Pkaon->Write("h_Eff_Pkaon");
   }
   if(TEfficiency::CheckConsistency(*h_kaonwtrk_length,*h_kaon_length)){
       h_Eff_Lkaon = tfs->make<TEfficiency>(*h_kaonwtrk_length,*h_kaon_length);
       TGraphAsymmErrors *grEff_Lkaon = h_Eff_Lkaon->CreateGraph();
       grEff_Lkaon->Write("grEff_Lkaon");
       h_Eff_Lkaon->Write("h_Eff_Lkaon");
   }
   if(TEfficiency::CheckConsistency(*h_Pmichel_e_num,*h_Pmichel_e_den)){
       h_Eff_Pmichel = tfs->make<TEfficiency>(*h_Pmichel_e_num,*h_Pmichel_e_den);
       TGraphAsymmErrors *grEff_Pmichel = h_Eff_Pmichel->CreateGraph();
       grEff_Pmichel->Write("grEff_Pmichel");
       h_Eff_Pmichel->Write("h_Eff_Pmichel");
   }
   if(TEfficiency::CheckConsistency(*h_michelwtrk_length,*h_michel_length)){
       h_Eff_Lmichel = tfs->make<TEfficiency>(*h_michelwtrk_length,*h_michel_length);
       TGraphAsymmErrors *grEff_Lmichel = h_Eff_Lmichel->CreateGraph();
       grEff_Lmichel->Write("grEff_Lmichel");
       h_Eff_Lmichel->Write("h_Eff_Lmichel");
   }

}
//========================================================================
DEFINE_ART_MODULE(ProtonReco)

} 

#endif // ProtonReco_Module
