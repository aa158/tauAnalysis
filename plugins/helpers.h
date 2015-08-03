/*
 * =====================================================================================
 *
 *       Filename:  Helpers.h
 *
 *    Description:  Common UCT functions.
 *
 *         Author:  M. Cepeda, S. Dasu, E. Friis
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#ifndef HELPERS_W9QK6HND
#define HELPERS_W9QK6HND
//MINIAOD
std::vector<const reco::GenParticle*> getGenParticleCollectionMiniAOD(const edm::Event& evt);
std::vector<const reco::GenParticle*> getGenEleCollectionMiniAOD(const edm::Event& evt);
std::vector<const reco::GenParticle*> getGenMuCollectionMiniAOD(const edm::Event& evt);
bool genMatchingMiniAOD(const pat::Tau tau, std::vector<const reco::GenParticle*>& GenPart, double maxDR);
//AODSIM
reco::PFJetRef getJetRef(const reco::PFTau& tau);
std::vector<const reco::GenParticle*> getGenParticleCollection(const edm::Event& evt);	
const reco::GenParticle* findBestGenMatch(const reco::PFTau& TagTauObj,std::vector<const reco::GenParticle*>& GenPart, double maxDR);
bool isLooseJet(const reco::PFJet jet);
bool isMediumJet(const reco::PFJet jet);
bool isTightJet(const reco::PFJet jet);
bool isLooseJet(const pat::Jet jet);
bool isMediumJet(const pat::Jet jet);
bool isTightJet(const pat::Jet jet);
const reco::GenParticle* findBestGenMatch(const pat::Tau& tauObj,std::vector<const reco::GenParticle*>& GenPart, double maxDR);
const pat::Jet* findBestJetMatch(const pat::Tau& tauObj,std::vector<const pat::Jet*>& jet_denom_vec, double maxDR);
int findBestGenMatchIndex(const pat::Tau& tauObj,std::vector<const reco::GenParticle*>& GenPart, double maxDR);
bool goodVertex(const pat::Tau tau, const reco::Vertex PV);
bool goodVertex(const reco::GenParticle* gen, const reco::Vertex PV);
const pat::Tau* findBestTauMatch(const reco::GenParticle GenPart, edm::Handle<pat::TauCollection> taus, double maxDR);
const reco::GenParticle* getGenTau(const pat::Tau& patTau);
void countDecayProducts(const reco::GenParticle* genParticle,
                          int& numElectrons, int& numElecNeutrinos, int& numMuons, int& numMuNeutrinos,
                          int& numChargedHadrons, int& numPi0s, int& numOtherNeutralHadrons, int& numPhotons);
std::string getGenTauDecayMode(const reco::GenParticle* genParticle);
const int genDecayMode(const reco::GenParticle* genTau);
const int highPtCharge(const pat::Tau &tau);
const int sumCharge(const pat::Tau &recoTau);
#endif

