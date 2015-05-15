
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
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TauTrigMatch.h"
#include "helpers.h"

reco::PFJetRef getJetRef(const reco::PFTau& tau) {
	if (tau.jetRef().isNonnull())
		return tau.jetRef();
	else if (tau.pfTauTagInfoRef()->pfjetRef().isNonnull())
		return tau.pfTauTagInfoRef()->pfjetRef();
	else throw cms::Exception("cant find jet ref");
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
		if(object.status() == 2) output.push_back(&object);
	}
	return output;
}


// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
//const reco::GenParticle* findBestGenMatch1(const reco::PFTau TagTauObj,
const reco::GenParticle* findBestGenMatch(const reco::PFTau& TagTauObj,
		std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
	const reco::GenParticle* output = NULL;
	double bestDeltaR = -1;
	for (size_t i = 0; i < GenPart.size(); ++i) {
		double deltaR = reco::deltaR(TagTauObj, *GenPart[i]);
		if (deltaR < maxDR) {
			if (!output || deltaR < bestDeltaR) {
				output = GenPart[i];
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

