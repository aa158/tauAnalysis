// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "helpers.h"

reco::PFJetRef getJetRef(const reco::PFTau& tau) {
	if (tau.jetRef().isNonnull())
		return tau.jetRef();
	else if (tau.pfTauTagInfoRef()->pfjetRef().isNonnull())
		return tau.pfTauTagInfoRef()->pfjetRef();
	else throw cms::Exception("cant find jet ref");
}

reco::PFJetRef getJetRef(const pat::Tau tau) {
        if (tau.pfJetRef().isNonnull())
                return tau.pfJetRef();
        else throw cms::Exception("cant find jet ref");
}

bool genMatchingMiniAOD(const pat::Tau tau, std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
	bool tau_match=false;
	for (size_t i = 0; i < GenPart.size(); ++i) {
		if (abs(GenPart[i]->pdgId())==15){
			double deltaR = reco::deltaR(tau, *GenPart[i]);
			if (deltaR < maxDR) {
				tau_match=true;
			}
		}
	}
	return tau_match;
}

std::vector<const reco::GenParticle*> getGenParticleCollectionMiniAOD(const edm::Event& evt) {
	std::vector<const reco::GenParticle*> output;
	edm::Handle< std::vector<reco::GenParticle> > handle;
	evt.getByLabel("prunedGenParticles", handle);
	// Loop over objects in current collection
	for (size_t j = 0; j < handle->size(); ++j) {
		const reco::GenParticle& object = handle->at(j);
		if(abs(object.pdgId()) == 15) output.push_back(&object);
	}
	return output;
}

std::vector<const reco::GenParticle*> getGenEleCollectionMiniAOD(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("prunedGenParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                if(abs(object.pdgId()) == 11) output.push_back(&object);
        }
        return output;
}

std::vector<const reco::GenParticle*> getGenMuCollectionMiniAOD(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("prunedGenParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                if(abs(object.pdgId()) == 13) output.push_back(&object);
        }
        return output;
}

// Get collection of generator particles with status 2
std::vector<const reco::GenParticle*> getGenParticleCollection(const edm::Event& evt) {
	std::vector<const reco::GenParticle*> output;
	edm::Handle< std::vector<reco::GenParticle> > handle;
	evt.getByLabel("genParticles", handle);
	// Loop over objects in current collection
	for (size_t j = 0; j < handle->size(); ++j) {
		const reco::GenParticle& object = handle->at(j);
		//if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
		if(abs(object.pdgId()) == 15) output.push_back(&object);
	}
	return output;
}

std::vector<const reco::GenParticle*> getGenEleCollection(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("genParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                //if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
                if(abs(object.pdgId()) == 11) output.push_back(&object);
        }
        return output;
}

std::vector<const reco::GenParticle*> getGenMuCollection(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("genParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                //if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
                if(abs(object.pdgId()) == 13) output.push_back(&object);
        }
        return output;
}


// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
//const reco::GenParticle* findBestGenMatch1(const reco::PFTau TagTauObj,
/*
const pat::Tau* findBestTauMatch(const reco::GenParticle* GenPart, edm::Handle<pat::TauCollection> taus, double maxDR) {
        const pat::Tau *output = NULL;
        double bestDeltaR = maxDR;
	for (const pat::Tau &tau : *taus) {
                double deltaR = reco::deltaR(tau, *GenPart);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                output = tau;
                                bestDeltaR = deltaR;
                        }
                }
        }
        return output;
}
*/

const reco::GenParticle* findBestGenMatch(const reco::PFTau& tauObj,
		std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
	const reco::GenParticle* output = NULL;
	double bestDeltaR = maxDR;
	for (size_t i = 0; i < GenPart.size(); ++i) {
		double deltaR = reco::deltaR(tauObj, *GenPart[i]);
		if (deltaR < maxDR) {
			if (deltaR < bestDeltaR) {
				output = GenPart[i];
				bestDeltaR = deltaR;
			}
		}
	}
	return output;
}

const reco::GenParticle* findBestGenMatch(const pat::Tau& tauObj,
                std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
        const reco::GenParticle* output = NULL;
        double bestDeltaR = maxDR;
        for (size_t i = 0; i < GenPart.size(); ++i) {
                double deltaR = reco::deltaR(tauObj, *GenPart[i]);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                output = GenPart[i];
                                bestDeltaR = deltaR;
                        }
                }
        }
        return output;
}

int findBestGenMatchIndex(const pat::Tau& tauObj,
                std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
        double bestDeltaR = maxDR;
	int index = -1;
        for (size_t i = 0; i < GenPart.size(); ++i) {
                double deltaR = reco::deltaR(tauObj, *GenPart[i]);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                bestDeltaR = deltaR;
                        }
                }
        }
	return index;
}

const pat::Jet* findBestJetMatch(const pat::Tau& tauObj,
                std::vector<const pat::Jet*>& jet_denom_vec, double maxDR) {
        const pat::Jet* output = NULL;
        double bestDeltaR = maxDR;
        for (size_t i = 0; i < jet_denom_vec.size(); ++i) {
                double deltaR = reco::deltaR(tauObj, *jet_denom_vec[i]);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                output = jet_denom_vec[i];
                                bestDeltaR = deltaR;
                        }
                }
        }
        return output;
}

bool isLooseJet(const reco::PFJet jet){
	bool loose = true;
	if (jet.neutralHadronEnergyFraction() >= 0.99) loose = false;
	if (jet.neutralEmEnergyFraction() >= 0.99) loose = false;
	if (jet.numberOfDaughters() <= 1) loose = false; //getPFConstitutents broken in miniAOD
	if (std::abs(jet.eta()) < 2.4) {
		if (jet.chargedHadronEnergyFraction() == 0) loose = false;
		if (jet.chargedHadronMultiplicity() == 0) loose = false;
		if (jet.chargedEmEnergyFraction() >= 0.99) loose = false;
	}
	return loose;
}
bool isMediumJet(const reco::PFJet jet){
	bool medium = true;
	if (jet.neutralHadronEnergyFraction() >= 0.95) medium = false;
	if (jet.neutralEmEnergyFraction() >= 0.95) medium = false;
	if (jet.numberOfDaughters() <= 1) medium = false; //getPFConstitutents broken in miniAOD
	if (std::abs(jet.eta()) < 2.4) {
		if (jet.chargedHadronEnergyFraction() == 0) medium = false;
		if (jet.chargedHadronMultiplicity() == 0) medium = false;
		if (jet.chargedEmEnergyFraction() >= 0.99) medium = false;
	}
	return medium;
}

bool isTightJet(const reco::PFJet jet){
	bool tight = true;
	if (jet.neutralHadronEnergyFraction() >= 0.90) tight = false;
	if (jet.neutralEmEnergyFraction() >= 0.90) tight = false;
	if (jet.numberOfDaughters() <= 1) tight = false; //getPFConstitutents broken in miniAOD
	if (std::abs(jet.eta()) < 2.4) {
		if (jet.chargedHadronEnergyFraction() == 0) tight = false;
		if (jet.chargedHadronMultiplicity() == 0) tight = false;
		if (jet.chargedEmEnergyFraction() >= 0.99) tight = false;
	}
	return tight;
}
bool isLooseJet(const pat::Jet jet){
        bool loose = true;
        if (jet.neutralHadronEnergyFraction() >= 0.99) loose = false;
        if (jet.neutralEmEnergyFraction() >= 0.99) loose = false;
        if (jet.numberOfDaughters() <= 1) loose = false; //getPFConstitutents broken in miniAOD
        if (std::abs(jet.eta()) < 2.4) {
                if (jet.chargedHadronEnergyFraction() == 0) loose = false;
                if (jet.chargedHadronMultiplicity() == 0) loose = false;
                if (jet.chargedEmEnergyFraction() >= 0.99) loose = false;
        }
        return loose;
}
bool isMediumJet(const pat::Jet jet){
        bool medium = true;
        if (jet.neutralHadronEnergyFraction() >= 0.95) medium = false;
        if (jet.neutralEmEnergyFraction() >= 0.95) medium = false;
        if (jet.numberOfDaughters() <= 1) medium = false; //getPFConstitutents broken in miniAOD
        if (std::abs(jet.eta()) < 2.4) {
                if (jet.chargedHadronEnergyFraction() == 0) medium = false;
                if (jet.chargedHadronMultiplicity() == 0) medium = false;
                if (jet.chargedEmEnergyFraction() >= 0.99) medium = false;
        }
        return medium;
}

bool isTightJet(const pat::Jet jet){
        bool tight = true;
        if (jet.neutralHadronEnergyFraction() >= 0.90) tight = false;
        if (jet.neutralEmEnergyFraction() >= 0.90) tight = false;
        if (jet.numberOfDaughters() <= 1) tight = false; //getPFConstitutents broken in miniAOD
        if (std::abs(jet.eta()) < 2.4) {
                if (jet.chargedHadronEnergyFraction() == 0) tight = false;
                if (jet.chargedHadronMultiplicity() == 0) tight = false;
                if (jet.chargedEmEnergyFraction() >= 0.99) tight = false;
        }
        return tight;
}

bool goodVertex(const pat::Tau tau, const reco::Vertex PV){
	bool close = false;
	if (abs(tau.vertex().z() - PV.z()) < .2) {
		float dxy = sqrt(pow(tau.vertex().x() - PV.x(),2) +pow(tau.vertex().y() - PV.y(),2));
		if (dxy < .045) {
			close = true;
		}
	}
	return close;
}

bool goodVertex(const reco::GenParticle* gen, const reco::Vertex PV){
	bool close = false;
	if (abs(gen->vertex().z() - PV.z()) < .2) {
		float dxy = sqrt(pow(gen->vertex().x() - PV.x(),2) + pow(gen->vertex().y() - PV.y(),2));
		if (dxy < .045) {
			close = true;
		}
	}
	return close;
}

const reco::GenParticle* getGenTau(const pat::Tau& patTau)
{
  std::vector<reco::GenParticleRef> associatedGenParticles = patTau.genParticleRefs();
  for ( std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin();
	it != associatedGenParticles.end(); ++it ) {
    if ( it->isAvailable() ) {
      const reco::GenParticleRef& genParticle = (*it);
      double deltaR = reco::deltaR(patTau, *genParticle);
      
      if (abs(genParticle->pdgId()) == 15 &&  deltaR<.3 && genParticle->pt()>20 && abs(genParticle->eta())<2.3) return genParticle.get();
    }
  }
  return 0;
}

void countDecayProducts(const reco::GenParticle* genParticle,
			  int& numElectrons, int& numElecNeutrinos, int& numMuons, int& numMuNeutrinos, 
			  int& numChargedHadrons, int& numPi0s, int& numOtherNeutralHadrons, int& numPhotons, int& chargeSum)
  {
    //std::cout << " genParticle: pdgId = " << genParticle->pdgId() << std::endl;

    int absPdgId = TMath::Abs(genParticle->pdgId());
    int status   = genParticle->status();
    int charge   = genParticle->charge();

    if      ( absPdgId == 111 ) ++numPi0s;
    else if ( status   ==   1 ) {
      chargeSum = chargeSum + charge;
      if      ( absPdgId == 11 ) ++numElectrons;
      else if ( absPdgId == 12 ) ++numElecNeutrinos;
      else if ( absPdgId == 13 ) ++numMuons;
      else if ( absPdgId == 14 ) ++numMuNeutrinos;
      else if ( absPdgId == 15 ) { 
	edm::LogError ("countDecayProducts")
	  << "Found tau lepton with status code 1 !!";
	return; 
      }
      else if ( absPdgId == 16 ) return; // no need to count tau neutrinos
      else if ( absPdgId == 22 ) ++numPhotons;
      else if ( charge   !=  0 ) ++numChargedHadrons;
      else                       ++numOtherNeutralHadrons;
    } else {
      unsigned numDaughters = genParticle->numberOfDaughters();
      for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
	const reco::GenParticle* daughter = genParticle->daughterRef(iDaughter).get();
	
	countDecayProducts(daughter, 
			   numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
			   numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons, chargeSum);
      }
    }
  }

std::string getGenTauDecayMode(const reco::GenParticle* genParticle) 
  {
//--- determine generator level tau decay mode
//
//    NOTE: 
//        (1) function implements logic defined in PhysicsTools/JetMCUtils/src/JetMCTag::genTauDecayMode
//            for different type of argument 
//        (2) this implementation should be more robust to handle cases of tau --> tau + gamma radiation
//
  
    //std::cout << "<MCEmbeddingValidationAnalyzer::getGenTauDecayMode>:" << std::endl;

    int numElectrons           = 0;
    int numElecNeutrinos       = 0;
    int numMuons               = 0;
    int numMuNeutrinos         = 0; 
    int numChargedHadrons      = 0;
    int numPi0s                = 0; 
    int numOtherNeutralHadrons = 0;
    int numPhotons             = 0;
    int chargeSum	       = 0;
    
    countDecayProducts(genParticle,
		       numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
		       numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons, chargeSum);
    
    if      ( numElectrons == 1 && numElecNeutrinos == 1 ) return std::string("tauDecaysElectron");
    else if ( numMuons     == 1 && numMuNeutrinos   == 1 ) return std::string("tauDecayMuon");
    
    switch ( numChargedHadrons ) {
    case 1 : 
      if ( numOtherNeutralHadrons != 0 ) return std::string("tauDecayOther");
      switch ( numPi0s ) {
      case 0:
	return std::string("oneProng0Pi0");
      case 1:
	return std::string("oneProng1Pi0");
      case 2:
	return std::string("oneProng2Pi0");
      case 3:
	return std::string("oneProng3Pi0");
      case 4:
	return std::string("oneProng4Pi0");
      default:
	return std::string("tauDecayOther");
      }
    case 2 :
	if ( numOtherNeutralHadrons != 0 ) return std::string("tauDecayOther");
        switch ( numPi0s) {
	case 0:
	    return std::string("twoProng0Pi0");
	case 1:
	    return std::string("twoProng1Pi0");
	case 2:
	    return std::string("twoProng2Pi0");
	case 3:
	    return std::string("twoProng3Pi0");
	case 4:
	    return std::string("twoProng4Pi0");
	default:
	    return std::string("tauDecayOther");
	}
    case 3 : 
      if ( numOtherNeutralHadrons != 0 ) return std::string("tauDecayOther");
      switch ( numPi0s ) {
      case 0:
	return std::string("threeProng0Pi0");
      case 1:
	return std::string("threeProng1Pi0");
      case 2:
	return std::string("threeProng2Pi0");
      case 3:
	return std::string("threeProng3Pi0");
      case 4:
	return std::string("threeProng4Pi0");
      default:
	return std::string("tauDecayOther");
      }
    default:
      return std::string("tauDecayOther");
    }
  }

const int genDecayMode(const reco::GenParticle* genTau)
{
    std::string decayModeString = getGenTauDecayMode(genTau);
    if (decayModeString == "oneProng0Pi0") {return 0;}
    else if (decayModeString == "oneProng1Pi0") {return 1;}
    else if (decayModeString == "oneProng2Pi0") {return 2;}
    else if (decayModeString == "oneProng3Pi0") {return 3;}
    else if (decayModeString == "oneProng4Pi0") {return 4;}
    else if (decayModeString == "twoProng0Pi0") {return 5;}
    else if (decayModeString == "twoProng1Pi0") {return 6;}
    else if (decayModeString == "twoProng2Pi0") {return 7;}
    else if (decayModeString == "twoProng3Pi0") {return 8;}
    else if (decayModeString == "twoProng4Pi0") {return 9;}
    else if (decayModeString == "threeProng0Pi0") {return 10;}
    else if (decayModeString == "threeProng1Pi0") {return 11;}
    else if (decayModeString == "threeProng2Pi0") {return 12;}
    else if (decayModeString == "threeProng3Pi0") {return 13;}
    else if (decayModeString == "threeProng4Pi0") {return 14;}
    else if (decayModeString == "tauDecaysElectron") {return 15;}
    else if (decayModeString == "tauDecayMuon") {return 16;}
    else if (decayModeString == "tauDecayOther") {return 17;}
    else {return -1;}
}

const int highPtCharge(const pat::Tau &recoTau)
{
	pat::Tau tau = recoTau;
	//float highPt = -999;
	int charge = -99;
	charge = tau.leadChargedHadrCand().get()->charge();
	//std::cout << "charge is " << charge << std::endl;
/*	tau.embedSignalTauChargedHadronCandidates();
	//tau.embedSignalTauCandidates();
	std::cout << "num cands is " << tau.signalTauChargedHadronCandidates().size() << std::endl;
	for (size_t itrk = 0; itrk < tau.signalTauChargedHadronCandidates().size(); ++itrk) {
	//	std::cout << "cand pt is " << tau.signalTauChargedHadronCandidates()[itrk].pt() << std::endl;
		if (tau.signalTauChargedHadronCandidates()[itrk].pt() > highPt)
		{
			charge = tau.signalTauChargedHadronCandidates()[itrk].charge();
			highPt = tau.signalTauChargedHadronCandidates()[itrk].pt();
		}
	}
*/	return charge;
}

const int sumCharge(const pat::Tau &recoTau)
{
        pat::Tau tau = recoTau;
	int chargeSum = 0;
	reco::CandidatePtrVector signalTracks = recoTau.signalCands();
	for ( unsigned iTrack = 0; iTrack < signalTracks.size(); ++iTrack ) {
		chargeSum = chargeSum + signalTracks[iTrack]->charge();
	}

	/*
	if (recoTau.decayMode() == 5) {
		for (unsigned iTrack = 0; iTrack < signalTracks.size(); ++iTrack ) {
			//std::cout << "pgID of sigcand is " << TMath::Abs(signalTracks[iTrack]->pdgId()) << "with charge " << signalTracks[iTrack]->charge() << std::endl;
		}
		std::cout << "total charge is " << chargeSum << std::endl;
	}*/

    return chargeSum;
}
