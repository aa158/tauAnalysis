// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminantFunctions.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "helpers.h"
#include "candidateAuxFunctions.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class dataMCcomp : public edm::EDAnalyzer {
	public:
		explicit dataMCcomp(const edm::ParameterSet&);
		~dataMCcomp();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	        void setWeightInfo(const edm::Event& event);
	        double getEventWeight(const edm::Event& event);

		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;
		edm::EDGetTokenT<pat::MuonCollection> muonToken_;
		edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;

		edm::Service<TFileService> fs;

		std::string tauID_;
		
		bool isMC_;

		TTree* tree;
		Float_t tauPt_;
		Float_t tauEta_;
		Int_t oldDMF_;
		Int_t newDMF_;
		Int_t tauIndex_;
		Int_t passDiscr_;
		Int_t nvtx_;
		Int_t goodReco_;
		Int_t tauDecayMode_;
		float tauMass_;
		float muMass_;
		int tauChargePt_;
		int tauChargeDirect_;
		int tauChargeSum_;
		int muChargeDirect_;
		double maxDR_;
		int eventCount_;
		double weight_;
		double initSumWeights_;
		double fidSumWeights_;
		unsigned int nProcessed_;
		unsigned int nPass_;
		double crossSection_;
		
		virtual void endJob() override;
		
};

dataMCcomp::dataMCcomp(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	genEventInfoToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator")))
{
	tauID_    = iConfig.getParameter<std::string>("tauID");
	isMC_ = iConfig.getParameter<bool>("isMC");

	tree = fs->make<TTree>("Ntuple", "Ntuple");
	tree->Branch("tauPt", &tauPt_,"tauPt_/F");
	tree->Branch("tauEta", &tauEta_,"tauEta_/F");
	tree->Branch("tauIndex", &tauIndex_,"tauIndex_/I");
	tree->Branch("passDiscr", &passDiscr_,"passDiscr_/I");
	tree->Branch("oldDMF", &oldDMF_,"oldDMF_/I");
	tree->Branch("newDMF",&newDMF_,"newDMF_/I");
	tree->Branch("nvtx",&nvtx_,"nvtx_/I");
	tree->Branch("goodReco",&goodReco_,"goodReco_/I");
	tree->Branch("tauDecayMode",&tauDecayMode_,"tauDecayMode_/I");
	tree->Branch("tauMass",&tauMass_,"tauMass_/I");
	tree->Branch("muMass",&muMass_,"muMass_/I");
	tree->Branch("tauChargePt",&tauChargePt_,"tauChargePt_/I");
	tree->Branch("tauChargeDirect",&tauChargeDirect_,"tauChargeDirect_/I");
	tree->Branch("tauChargeSum",&tauChargeSum_,"tauChargeSum_/I");
	tree->Branch("muChargeDirect",&muChargeDirect_,"muChargeDirect_/I");
	tree->Branch("eventCount",&eventCount_,"eventCount_/I");
	tree->Branch("weight",&weight_,"weight_/I");
	maxDR_ = 0.3;
	initSumWeights_ = 0;
	fidSumWeights_ = 0;
	nPass_ = 0;
	crossSection_ = iConfig.getUntrackedParameter<double>("xSec", -1);

}

dataMCcomp::~dataMCcomp()
{
}

	void
dataMCcomp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//if (isMC_) {
	//std::cout << "we analyzing" << std::endl;
	nProcessed_++;
	if (isMC_) {initSumWeights_ += getEventWeight(iEvent);}
		eventCount_=1;
		edm::Handle<reco::VertexCollection> vertices;
		iEvent.getByToken(vtxToken_, vertices);
		if (vertices->empty()) return;
		const reco::Vertex &PV = vertices->front();
		nvtx_=vertices->size();
		edm::Handle<pat::TauCollection> taus;
		edm::Handle<pat::MuonCollection> mus;
		iEvent.getByToken(tauToken_, taus);
		iEvent.getByToken(muonToken_, mus);
		muChargeDirect_=-99;

		const std::vector<pat::Muon> * muons = mus.product();
		unsigned int nbMuon =  muons->size();
		std::auto_ptr<pat::MuonCollection> musID(new pat::MuonCollection);
		int musSize = 0;
		for(unsigned i = 0 ; i < nbMuon; i++){
			pat::Muon muon(muons->at(i));
			float muIso = (muon.chargedHadronIso()+std::max(muon.photonIso()+muon.neutralHadronIso()-(0.5*(muon.puChargedHadronIso())),0.0))/(muon.pt());
			////std::cout << "muon.chargedHadronIso(): " << muon.chargedHadronIso() << "; ";
			////std::cout << "second term: " << std::max(muon.photonIso()+muon.neutralHadronIso()-(0.5*(muon.puChargedHadronIso())),0.0) << std::endl;
			////std::cout << muIso << std::endl;
		 //////std::cout<<"muIso: "<<muIso<<std::endl;
			//float muIso = (muon.pfIsolationR03().sumChargedHadronPt + max(
		 //muon.pfIsolationR03().sumNeutralHadronPt +
		 //muon.pfIsolationR03().sumPhotonEt - 
		 //0.5 * muon.pfIsolationR03().sumPUPt, 0.0)) / muon.pt();
      			bool goodGlob = muon.isGlobalMuon() && 
                      			muon.globalTrack()->normalizedChi2() < 3 && 
                      			muon.combinedQuality().chi2LocalPosition < 12 && 
                      			muon.combinedQuality().trkKink < 20; 
      			bool isMedium = muon::isLooseMuon(muon) && 
                      			muon.innerTrack()->validFraction() > 0.8 && 
                      			muon::segmentCompatibility(muon) > (goodGlob ? 0.303 : 0.451); 
			int muId = isMedium;
			/* 
			if (muon.isLooseMuon()&&(((muon.isGlobalMuon()&&muon.globalTrack()->normalizedChi2()<3&&muon.combinedQuality().chi2LocalPosition<12&&muon.combinedQuality().trkKink<20)&&(muon.innerTrack()->validFraction()>=0.8&&muon.segmentCompatibility()>=0.303))||(!(muon.isGlobalMuon()&&muon.globalTrack()->normalizedChi2()<3&&muon.combinedQuality().chi2LocalPosition<12&&muon.combinedQuality().trkKink<20)&&(muon.innerTrack()->validFraction()>=0.8&&muon.segmentCompatibility()>=0.451))))
			{
				muId=1;
			}
			*/
			muon.addUserFloat("dBRelIso",muIso);
			////std::cout << "does it has? " << muon.hasUserFloat("dBRelIso") << " and it's " << muon.userFloat("dBRelIso") << std::endl;
		
			muon.addUserInt("mediumID",muId);
			musSize++;
			musID->push_back(muon);
		}

		////std::cout << "cand mus: " << musSize << std::endl;

		tauPt_=-999;
		tauEta_=-999;
		tauIndex_=-1;
		passDiscr_=0;
		tauDecayMode_=-1;
		oldDMF_=0;
		newDMF_=0;
		tauMass_=0;
		tauChargePt_=-99;
		tauChargeDirect_=-99;
		tauChargeSum_=-99;
		goodReco_=0;
		muMass_=0;
		int itau = -1;
		for (const pat::Tau &tau : *taus) {
			tauPt_=tau.pt();
			tauEta_=tau.eta();
			tauMass_=tau.mass(); //eV
			//////std::cout << "about to access charge for reco tau" << std::endl;
			tauChargePt_ = highPtCharge(tau);
			//////std::cout << "stored Pt charge" << std::endl;
			tauChargeDirect_=tau.charge();
			tauChargeSum_ = sumCharge(tau); // sumCharge(tau);
			//////std::cout << "stored sum charge" << std::endl;
			itau++;
			oldDMF_=tau.tauID("decayModeFinding"); // this is the old DMF; strictly tighter than new DMF
			newDMF_=tau.tauID("decayModeFindingNewDMs");
			bool pass_Loose = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>.5;
			tauDecayMode_=tau.decayMode();
			bool singleMu = 0;
				//	closevtx = goodVertex(tau, PV);
			for(pat::Muon muon : *musID){
				////std::cout << "muon pt: " << muon.pt() << "; ";
				////std::cout << "muon eta: " << muon.eta() << "; ";
				////std::cout << "delta R: " << reco::deltaR(tau,muon) << "; ";
				////std::cout << "medID: " << muon.userInt("mediumID") << "; ";
				////std::cout << "iso: " << muon.userFloat("dBRelIso") << std::endl;
				if (muon.pt()>20 && muon.eta() < 2.1 && reco::deltaR(tau,muon) > .5 && muon.userInt("mediumID")>.5 && muon.userFloat("dBRelIso")<.1 && !muon.hasUserInt("Used")){
					////std::cout << "good muon" << std::endl;
					singleMu = 1;
					muMass_ = (tau.p4() + muon.p4()).mass();
					muChargeDirect_ = muon.charge();
					muon.addUserInt("Used",1);
					break;
				}
			}
			if (tauID_ == "decayModeFinding") {
				// //////std::cout << "Old DMF\n";
				if (tauPt_ > 20 && abs(tauEta_)<2.3 && pass_Loose && oldDMF_>.5 && singleMu) {
					//tauPt_=tau.pt();
					goodReco_=1;
					passDiscr_=tau.tauID(tauID_);
					nPass_++;
					if (isMC_) {setWeightInfo(iEvent);}
					tree->Fill();
					////////std::cout << "new size of used_taus is " << used_taus.size() << std::endl;
				}
			}
			else {
				if (tauPt_ > 20 && abs(tauEta_)<2.3 && tau.tauID(tauID_)>.5 && newDMF_>.5 && singleMu) {
					//tauPt_=tau.pt();
					goodReco_=1;
					passDiscr_=tau.tauID(tauID_);
					nPass_++;
					if (isMC_) {setWeightInfo(iEvent);}
					tree->Fill();
					////////std::cout << "new size of used_taus is " << used_taus.size() << std::endl;
				} // end if tau passes criteria
			}
			//std::cout << "filling tree" << std::endl;
			eventCount_=0;
		} // end tau for loop
	//} // end if it's a mc
	//else {
	//} // end if data
}

double dataMCcomp::getEventWeight(const edm::Event& event) {
    edm::Handle<GenEventInfoProduct> genEventInfo;
    event.getByToken(genEventInfoToken_, genEventInfo);
    return genEventInfo->weights()[0];
}

void
dataMCcomp::setWeightInfo(const edm::Event& event) {
    weight_ = getEventWeight(event);
    fidSumWeights_ += weight_;
}

void 
dataMCcomp::endJob() 
{    
    TTree* metaData = fs->make<TTree>("MetaData", "MetaData");
    metaData->Branch("nProcessedEvents", &nProcessed_);
    metaData->Branch("nPass", &nPass_);
    metaData->Branch("inputXSection", &crossSection_);
    metaData->Branch("initSumWeights", &initSumWeights_);
    metaData->Branch("fidSumWeights", &fidSumWeights_);
    float fidXSec = crossSection_ * fidSumWeights_/initSumWeights_;
    metaData->Branch("fidXSection", &fidXSec);
    metaData->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(dataMCcomp);


