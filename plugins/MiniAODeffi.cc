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

class MiniAODeffi : public edm::EDAnalyzer {
	public:
		explicit MiniAODeffi(const edm::ParameterSet&);
		~MiniAODeffi();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;

		std::string tauID_;
		
		bool isMC_;

		TTree* tree;
		Float_t tauPt_;
		Float_t genPt_;
		Float_t ratioPt_;
		Float_t tauEta_;
		Float_t genEta_;
		Int_t oldDMF_;
		Int_t newDMF_;
		Int_t tauIndex_;
		Int_t passDiscr_;
		Int_t genMatchedTau_;
		Int_t nvtx_;
		Int_t goodReco_;
		Int_t tauDecayMode_;
		Int_t genDecayMode_;
		float tauMass_;
		float genMass_;
		int genCharge_;
		int tauChargePt_;
		int tauChargeDirect_;
		int tauChargeSum_;
		double maxDR_;
};

MiniAODeffi::MiniAODeffi(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
{

	tauID_    = iConfig.getParameter<std::string>("tauID");
	isMC_ = iConfig.getParameter<bool>("isMC");

	edm::Service<TFileService> fs;

	tree = fs->make<TTree>("Ntuple", "Ntuple");
	tree->Branch("tauPt", &tauPt_,"tauPt_/F");
	tree->Branch("genPt", &genPt_,"genPt_/F");
	tree->Branch("ratioPt",&ratioPt_,"ratioPt_/F");
	tree->Branch("genEta", &genEta_,"genEta_/F");
	tree->Branch("tauEta", &tauEta_,"tauEta_/F");
	tree->Branch("tauIndex", &tauIndex_,"tauIndex_/I");
	tree->Branch("passDiscr", &passDiscr_,"passDiscr_/I");
	tree->Branch("oldDMF", &oldDMF_,"oldDMF_/I");
	tree->Branch("newDMF",&newDMF_,"newDMF_/I");
	tree->Branch("genMatchedTau", &genMatchedTau_,"genMatchedTau_/I");
	tree->Branch("nvtx",&nvtx_,"nvtx_/I");
	tree->Branch("goodReco",&goodReco_,"goodReco_/I");
	tree->Branch("tauDecayMode",&tauDecayMode_,"tauDecayMode_/I");
	tree->Branch("genDecayMode",&genDecayMode_,"genDecayMode_/I");
	tree->Branch("tauMass",&tauMass_,"tauMass_/I");
	tree->Branch("genMass",&genMass_,"genMass_/I");
	tree->Branch("genCharge",&genCharge_,"genCharge_/I");
	tree->Branch("tauChargePt",&tauChargePt_,"tauChargePt_/I");
	tree->Branch("tauChargeDirect",&tauChargeDirect_,"tauChargeDirect_/I");
	tree->Branch("tauChargeSum",&tauChargeSum_,"tauChargeSum_/I");
	maxDR_ = 0.3;


}

MiniAODeffi::~MiniAODeffi()
{
}

	void
MiniAODeffi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if (isMC_) {
		edm::Handle<reco::VertexCollection> vertices;
		iEvent.getByToken(vtxToken_, vertices);
		const reco::Vertex &PV = vertices->front();
		nvtx_=vertices->size();
		edm::Handle<pat::TauCollection> taus;
		iEvent.getByToken(tauToken_, taus);
		std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollectionMiniAOD(iEvent);
		genMatchedTau_=0;
		tauPt_=-999;
		genPt_=-999;
		ratioPt_=0;
		tauEta_=-999;
		genEta_=-999;
		tauIndex_=-1;
		passDiscr_=0;
		goodReco_=0;
		tauDecayMode_=-1;
		genDecayMode_=-1;
		oldDMF_=0;
		newDMF_=0;
		tauMass_=0;
		genMass_=0;
		genCharge_=-99;
		tauChargePt_=-99;
		tauChargeDirect_=-99;
		tauChargeSum_=-99;
		int tau_position=-1;
		for (size_t i = 0; i < GenObjects.size(); ++i) {
			ratioPt_=0;
			genPt_ = getVisMomentum(GenObjects[i]).Pt();
			genEta_=GenObjects[i]->eta();
			genMass_=getVisMass(GenObjects[i]); //eV
			//std::cout << "getting gencharge" << std::endl;
			genCharge_=GenObjects[i]->charge();
			//std::cout << "got gencharge" << std::endl;
			tauIndex_=-1;
			passDiscr_=0;
			goodReco_=0;
			tau_position++;
			genMatchedTau_=1;
			oldDMF_=0;
			newDMF_=0;
			genDecayMode_=-1;
			tauDecayMode_=-1;
			// add vertex requirement
			// also plot the vertex variable
			// also do old dm for comparison
			// closevtx = goodVertex(GenObjects[i], PV);
			if (genPt_ > 20 && abs(genEta_)<2.3) {
				genMatchedTau_=1;
				tauIndex_=tau_position;
				int itau = -1;
				genDecayMode_=genDecayMode(GenObjects[i]);
				std::vector<int> used_taus;
				for (const pat::Tau &tau : *taus) {
					tauPt_=tau.pt();
					ratioPt_=tauPt_/genPt_;
					tauEta_=tau.eta();
					tauMass_=tau.mass(); //eV
					//std::cout << "about to access charge for reco tau" << std::endl;
					tauChargePt_ = highPtCharge(tau);
					//std::cout << "stored Pt charge" << std::endl;
					tauChargeDirect_=tau.charge();
					tauChargeSum_ = sumCharge(tau); // sumCharge(tau);
					//std::cout << "stored sum charge" << std::endl;
					itau++;
					oldDMF_=tau.tauID("decayModeFinding"); // this is the old DMF; strictly tighter than new DMF
					newDMF_=tau.tauID("decayModeFindingNewDMs");
					double deltaR = reco::deltaR(tau, *GenObjects[i]);
					bool pass_Loose = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>.5;
					tauDecayMode_=tau.decayMode();
					if (std::find(used_taus.begin(), used_taus.end(), itau) != used_taus.end()) {
						////std::cout << "found a used tau, number " << itau << std::endl;
						continue;
					}
				//	closevtx = goodVertex(tau, PV);
					if (tauID_ == "decayModeFinding") {
						// //std::cout << "Old DMF\n";
						if (tauPt_ > 20 && abs(tauEta_)<2.3 && pass_Loose && deltaR<maxDR_ && oldDMF_>.5) {
							//tauPt_=tau.pt();
							goodReco_=1;
							passDiscr_=tau.tauID(tauID_);
							used_taus.push_back(itau);
							////std::cout << "new size of used_taus is " << used_taus.size() << std::endl;
							break;
						}
					}
					else if (tauID_ == "kOneProng0PiZero") {
						if (genDecayMode_!=0) {
							genMatchedTau_=0;
							break;
						}
						if (tauPt_ > 20 && abs(tauEta_)<2.3 && pass_Loose && deltaR<maxDR_ && tauDecayMode_==0 && newDMF_>.5) {
                                                        goodReco_=1;
                                                        passDiscr_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
							used_taus.push_back(itau);
                                                        break;
						}
					}
					else if (tauID_ == "kOneProng1PiZero") {
						if (genDecayMode_!=1) {
							genMatchedTau_=0;
							break;
						}
                                                if (tauPt_ > 20 && abs(tauEta_)<2.3 && pass_Loose && deltaR<maxDR_ && tauDecayMode_==1 && newDMF_>.5) {
                                                        goodReco_=1;
                                                        passDiscr_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
							used_taus.push_back(itau);
                                                        break;
						}
					}
					else if (tauID_ == "kOneProng2PiZero") {
						if (genDecayMode_!=2) {
							genMatchedTau_=0;
							break;
						}
                                                if (tauPt_ > 20 && abs(tauEta_)<2.3 && pass_Loose && deltaR<maxDR_ && tauDecayMode_==2 && newDMF_>.5) {
                                                        goodReco_=1;
                                                        passDiscr_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
							used_taus.push_back(itau);
                                                        break;
						}
					}
					else {
						if (tauPt_ > 20 && abs(tauEta_)<2.3 && tau.tauID(tauID_)>.5 && deltaR<maxDR_ && newDMF_>.5) {
							//tauPt_=tau.pt();
							goodReco_=1;
							passDiscr_=tau.tauID(tauID_);
							used_taus.push_back(itau);
							////std::cout << "new size of used_taus is " << used_taus.size() << std::endl;
							break;
						} // end if tau passes criteria
					}
				} // end tau for loop
				tree->Fill();
			} //end if gen tau matches critera
		} //end gen tau for loop
	} // end if it's a mc
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODeffi);


