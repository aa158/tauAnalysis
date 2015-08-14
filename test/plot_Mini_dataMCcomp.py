'''
Usage:python plot.py RootFile.root labelselection,[optional]

Script to make some quick efficiency plots to test ntuple generation.


Author: L. Dodd, UW Madison

'''
import sys
from subprocess import Popen
from sys import argv, exit, stdout, stderr

import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

######## File #########
if len(argv) < 4:
   print 'Usage:python plot.py WJets.root DYJets.root singleMuon.root labelselection,[optional]'
   exit()

WJetsFile = argv[1]
DYJetsFile = argv[2]
singleMuonFile = argv[3]
ntuple_WJets = ROOT.TFile(WJetsFile)
ntuple_DYJets = ROOT.TFile(DYJetsFile)
ntuple_singleMuon = ROOT.TFile(singleMuonFile)

######## LABEL & SAVE WHERE #########

saveWhere='~/myAnalysis/CMSSW_7_4_0/src/RecoTauTag/tauAnalysis/outputs/'



#####################################
#Get Effi NTUPLE                 #
#####################################

byLooseCmbIso3_WJets = ntuple_WJets.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byMedCmbIso3_WJets = ntuple_WJets.Get("byMediumCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byTightCmbIso3_WJets = ntuple_WJets.Get("byTightCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
OldDMF_WJets = ntuple_WJets.Get("decayModeFinding/Ntuple")
NewDMF_WJets = ntuple_WJets.Get("decayModeFindingNewDMs/Ntuple")
WJets_meta = ntuple_WJets.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/MetaData")

byLooseCmbIso3_DYJets = ntuple_DYJets.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byMedCmbIso3_DYJets = ntuple_DYJets.Get("byMediumCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byTightCmbIso3_DYJets = ntuple_DYJets.Get("byTightCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
OldDMF_DYJets = ntuple_DYJets.Get("decayModeFinding/Ntuple")
NewDMF_DYJets = ntuple_DYJets.Get("decayModeFindingNewDMs/Ntuple")
DYJets_meta = ntuple_DYJets.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/MetaData")

byLooseCmbIso3_singleMuon = ntuple_singleMuon.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byMedCmbIso3_singleMuon = ntuple_singleMuon.Get("byMediumCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byTightCmbIso3_singleMuon = ntuple_singleMuon.Get("byTightCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
OldDMF_singleMuon = ntuple_singleMuon.Get("decayModeFinding/Ntuple")
NewDMF_singleMuon = ntuple_singleMuon.Get("decayModeFindingNewDMs/Ntuple")
singleMuon_meta = ntuple_singleMuon.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/MetaData")

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot_mass(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def produce_mass(ntuple, variable, PtCut,selection,binning, filename,color,xaxis,title,norm):
    hist = make_plot_mass(ntuple, variable, selection,binning,xaxis, title)

    rebin = 1
    fcolor = color # ROOT.kGreen+1
    lcolor = color
    sf = norm
    fillStyle = 1
    markerStyle = 21

    hist.Rebin( rebin )
    hist.SetFillColor( fcolor )
    #hist.SetLineColor( lcolor )
    hist.SetMarkerStyle(markerStyle)
    hist.SetMarkerColor(fcolor)
    #hist.SetLineWidth( 2 )
    hist.Scale( sf )
    max_hist = hist.GetMaximum()
    #hist.GetYaxis().SetRangeUser(0,1.2*max_hist)

    return hist

def make_stacked_comp(ntuple1,legend1,ntuple2, legend2, data, legend3, variable, PtCut,
                        selection1,selection2, selection3,
			cross1, cross2, cross3,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    hs = ROOT.THStack("hs","Stacked 1D Histograms")

    ntuple1.Draw("eventCount>>hist1", "", "goff")
    events1 = ROOT.gDirectory.Get("hist1").GetMean() * ROOT.gDirectory.Get("hist1").GetEntries()
    #print events1
    ntuple2.Draw("eventCount>>hist2", "", "goff")
    events2 = ROOT.gDirectory.Get("hist2").GetMean() * ROOT.gDirectory.Get("hist2").GetEntries()
    #print events2
    data.Draw("eventCount>>hist3", "", "goff")
    events3 = ROOT.gDirectory.Get("hist3").GetMean() * ROOT.gDirectory.Get("hist3").GetEntries()
    #print events3

    norm1 = 61526.7 * 40 / 24151270
    norm2 = 6025 * 40 / 28825132
    norm3 = 1
    
    l1 = produce_mass(ntuple1,variable, PtCut,selection1,binning, filename,ROOT.kMagenta-3,xaxis,title,norm1)
    l2 = produce_mass(ntuple2,variable, PtCut,selection2,binning, filename,ROOT.kBlue-9, xaxis, title,norm2)

    l3 = produce_mass(data,variable, PtCut,selection3,binning, filename,ROOT.kRed+3, xaxis, title,norm3)

    hs.Add(l1)
    hs.Add(l2)

    if (hs.GetMaximum() > l3.GetMaximum()):
        frame.SetMaximum(1.2*hs.GetMaximum())
    else :
        frame.SetMaximum(1.2*l3.GetMaximum())

    #frame.SetMaximum(1)
    #frame.SetMinimum(0)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    hs.Draw('same')
    l3.Draw('e1psame')
    #canvas.SetLogy()
    legend = ROOT.TLegend(0.1,0.75,0.28,0.88, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "f")
    legend.AddEntry(l2,legend2, "f")
    legend.AddEntry(l3,legend3, "pe")
    legend.Draw('sames')
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_plot_mass_sum(tree, variable1, variable2, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "(%s+%s)/1000000000>>htemp(%s)" % (variable1, variable2, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def produce_mass_sum(ntuple, variable1, variable2, PtCut,selection,binning, filename,color,xaxis,title,norm):
    hist = make_plot_mass_sum(ntuple, variable1, variable2, selection,binning,xaxis, title)

    rebin = 1
    fcolor = color # ROOT.kGreen+1
    lcolor = color
    sf = norm
    fillStyle = 1
    markerStyle = 21

    hist.Rebin( rebin )
    hist.SetFillColor( fcolor )
    #hist.SetLineColor( lcolor )
    hist.SetMarkerStyle(markerStyle)
    hist.SetMarkerColor(fcolor)
    #hist.SetLineWidth( 2 )
    hist.Scale( sf )
    max_hist = hist.GetMaximum()
    #hist.GetYaxis().SetRangeUser(0,1.2*max_hist)

    return hist

def make_stacked_comp_sum(ntuple1,legend1,ntuple2, legend2, data, legend3, variable1, variable2, PtCut,
                        selection1,selection2, selection3,
                        cross1, cross2, cross3,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    hs = ROOT.THStack("hs","Stacked 1D Histograms")

    ntuple1.Draw("eventCount>>shist1", "", "goff")
    events1 = ROOT.gDirectory.Get("shist1").GetMean() * ROOT.gDirectory.Get("shist1").GetEntries()
    #print events1
    ntuple2.Draw("eventCount>>shist2", "", "goff")
    events2 = ROOT.gDirectory.Get("shist2").GetMean() * ROOT.gDirectory.Get("shist2").GetEntries()
    #print events2
    data.Draw("eventCount>>shist3", "", "goff")
    events3 = ROOT.gDirectory.Get("shist3").GetMean() * ROOT.gDirectory.Get("shist3").GetEntries()
    #print events3

    norm1 = events1 / cross1 * lumi
    norm2 = events2 / cross2 * lumi
    norm3 = events3 / cross3 * lumi

    l1 = produce_mass_sum(ntuple1,variable1, variable2, PtCut,selection1,binning, filename,ROOT.kMagenta-3,xaxis,title,norm1)
    l2 = produce_mass_sum(ntuple2,variable1, variable2, PtCut,selection2,binning, filename,ROOT.kBlue-9, xaxis, title,norm2)

    l3 = produce_mass_sum(data,variable1, variable2, PtCut,selection3,binning, filename,ROOT.kRed+3, xaxis, title,norm3)

    hs.Add(l1)
    hs.Add(l2)

    if (hs.GetMaximum() > l3.GetMaximum()):
        frame.SetMaximum(1.2*hs.GetMaximum())
    else :
        frame.SetMaximum(1.2*l3.GetMaximum())

    #frame.SetMaximum(1)
    #frame.SetMinimum(0)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    hs.Draw('same')
    l3.Draw('e1psame')
    #canvas.SetLogy()
    legend = ROOT.TLegend(0.1,0.75,0.28,0.88, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "f")
    legend.AddEntry(l2,legend2, "f")
    legend.AddEntry(l3,legend3, "pe")
    legend.Draw('sames')
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

make_stacked_comp(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'tauMass', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
			.1125,.03367,1,
                        [36,900000000,1200000000],
			'tau_mass_comp',
                        'Tau Mass for Single Muon', 'Tau Mass (GeV)', 'Events')

make_stacked_comp_sum(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'tauMass','muMass', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
                        .1125,.03367,1,
                        [36,1.9,2.2],
                        'tauMu_mass_comp',
                        'Tau+Muon Mass for Single Muon', 'Tau+Muon Mass (GeV)', 'Events')

make_stacked_comp(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'tauPt', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
                        .1125,.03367,1,
                        [50,0,450],
                        'tau_pT_comp',
                        'Tau pT for Single Muon', 'Tau pT (GeV)', 'Events')

make_stacked_comp(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'tauEta', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
                        .1125,.03367,1,
                        [30,-2.4,2.4],
                        'tau_eta_comp',
                        'Tau Eta for Single Muon', 'Tau Eta', 'Events')

make_stacked_comp(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'nvtx', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
                        .1125,.03367,1,
                        [36,0,35],
                        'tau_nvtx_comp',
                        'Tau Nvtx for Single Muon', 'Nvtx', 'Events')

make_stacked_comp(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'tauChargeDirect', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
                        .1125,.03367,1,
                        [10,-3,3],
                        'tau_charge_comp',
                        'Tau Charge for Single Muon', 'Tau Charge (tau.charge())', 'Events')

make_stacked_comp(byLooseCmbIso3_WJets,"W+Jets",byLooseCmbIso3_DYJets, "DY Jets", byLooseCmbIso3_singleMuon, "data",'tauDecayMode', 20,
                        "goodReco==1","goodReco==1", "goodReco==1",
                        .1125,.03367,1,
                        [100,-10,10],
                        'tau_decayMode_comp',
                        'Tau Decay Mode for Single Muon', 'Tau Decay Mode', 'Events')
