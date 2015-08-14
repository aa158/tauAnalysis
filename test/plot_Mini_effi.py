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

import time

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root labelselection,[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/myAnalysis/CMSSW_7_4_0/src/RecoTauTag/tauAnalysis/outputs/'+argvselection,[2]+'_'
else:
   saveWhere='~/myAnalysis/CMSSW_7_4_0/src/RecoTauTag/tauAnalysis/outputs/'



#####################################
#Get Effi NTUPLE                 #
#####################################

byLooseCmbIso3 = ntuple_file.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byMedCmbIso3 = ntuple_file.Get("byMediumCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byTightCmbIso3 = ntuple_file.Get("byTightCombinedIsolationDeltaBetaCorr3Hits/Ntuple")

ntrlIsoPtSum = ntuple_file.Get("neutralIsoPtSum/Ntuple")
puCorrPtSum = ntuple_file.Get("puCorrPtSum/Ntuple")
MuLoose3 = ntuple_file.Get("againstMuonLoose3/Ntuple")
MuTight3 = ntuple_file.Get("againstMuonTight3/Ntuple")
EleVLooseMVA5 = ntuple_file.Get("againstElectronVLooseMVA5/Ntuple")
EleLooseMVA5 = ntuple_file.Get("againstElectronLooseMVA5/Ntuple")
EleMediumMVA5 = ntuple_file.Get("againstElectronMediumMVA5/Ntuple")

OldDMF = ntuple_file.Get("decayModeFinding/Ntuple")
NewDMF = ntuple_file.Get("decayModeFindingNewDMs/Ntuple")

kOneProng0PiZero = ntuple_file.Get("kOneProng0PiZero/Ntuple")
kOneProng1PiZero = ntuple_file.Get("kOneProng1PiZero/Ntuple")
kOneProng2PiZero = ntuple_file.Get("kOneProng2PiZero/Ntuple")

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def make_efficiency(denom, num):
    ''' Make an efficiency graph with the style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.kBlack)
    return eff

def make_num(ntuple, variable,PtCut,selection,binning):
    num = make_plot(
        ntuple, variable,
        selection,#"goodReco==1",
        binning
    )
    return num
def make_denom(ntuple, variable,PtCut,selection,binning):
    denom = make_plot(
        ntuple, variable,
        selection,#"genMatchedTau==1",
        binning
    )
    return denom

def produce_efficiency(ntuple, variable, PtCut,selection_num,selection_denom,binning, filename,color):
    denom = make_denom(ntuple, variable,PtCut,selection_denom,binning)
    num = make_num(ntuple,variable,PtCut,selection_num,binning)
    l1 = make_efficiency(denom,num)
    l1.SetMarkerColor(color)
    return l1

def compare_efficiencies(ntuple1,legend1,ntuple2,legend2, variable, PtCut,
			selection_num1,selection_denom1,selection_num2,selection_denom2,
			binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,selection_num1,selection_denom1,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,selection_num2,selection_denom2,binning, filename,ROOT.kBlue-9)
    frame.SetMaximum(l1.GetMaximum())
    frame.SetMinimum(0)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    legend = ROOT.TLegend(0.5, 0.1, 0.89, 0.2, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def compare_3efficiencies(ntuple1,legend1,ntuple2, legend2,ntuple3, legend3, variable, PtCut, 
			selection_num1,selection_denom1,selection_num2,selection_denom2,selection_num3,selection_denom3,
			binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,selection_num1,selection_denom1,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,selection_num2,selection_denom2,binning, filename,ROOT.kBlue-9)
    l3 = produce_efficiency(ntuple3,variable, PtCut,selection_num3,selection_denom3,binning, filename,ROOT.kRed+3)
    frame.SetMaximum(1)
    frame.SetMinimum(0)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    l3.Draw('pesame')
    legend = ROOT.TLegend(0.1,0.1,0.48,0.3, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.AddEntry(l3,legend3, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def compare_4efficiencies(ntuple1,legend1,ntuple2, legend2,ntuple3, legend3,ntuple4,legend4, variable, PtCut,
                        selection_num1,selection_denom1,selection_num2,selection_denom2,selection_num3,selection_denom3,selection_num4,selection_denom4,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,selection_num1,selection_denom1,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,selection_num2,selection_denom2,binning, filename,ROOT.kBlue-9)
    l3 = produce_efficiency(ntuple3,variable, PtCut,selection_num3,selection_denom3,binning, filename,ROOT.kRed+3)
    l4 = produce_efficiency(ntuple4,variable, PtCut,selection_num4,selection_denom4,binning, filename,ROOT.kRed-3)

    frame.SetMaximum(l1.GetMaximum())
    frame.SetMinimum(0)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    l3.Draw('pesame')
    l4.Draw('pesame')
    legend = ROOT.TLegend(0.5, 0.1, 0.89, 0.4, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.AddEntry(l3,legend3, "pe")
    legend.AddEntry(l4,legend4, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def produce_ratio(ntuple,x_var,y_var,selection,binning,filename,color,xaxis,yaxis,title):
    draw_string = "%s:%s>>Graph(%s)" % (y_var, x_var, ", ".join(str(x) for x in binning))
    ntuple.Draw(draw_string, selection, "prof goff")
    ratio = ROOT.gDirectory.Get("Graph").Clone()
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.5)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.GetXaxis().SetTitle(xaxis)
    ratio.GetYaxis().SetTitle(yaxis)
    ratio.SetTitle(title)
    return ratio

def plot_ratio(ntuple,legend, x_var, y_var,
			selection,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_ratio(ntuple,x_var,y_var,selection,binning,filename,ROOT.kMagenta-3,xaxis,yaxis,title)
    frame.SetMaximum(1.05)
    frame.SetMinimum(.95)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    legend = ROOT.TLegend(0.5, 0.1, 0.89, 0.4, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    #legend.AddEntry(l1,legend, "pe")
    #legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_plot_mass(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s/1000000>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def produce_mass(ntuple, variable, PtCut,selection,binning, filename,color,xaxis,title):
    hist = make_plot_mass(ntuple, variable, selection,binning,xaxis, title)

    rebin = 1
    fcolor = 0 # ROOT.kGreen+1
    lcolor = color
    sf = 1
    fillStyle = 1

    hist.Rebin( rebin )
    hist.SetFillColor( fcolor )
    hist.SetLineColor( lcolor )
    hist.SetLineWidth( 2 )
    hist.Scale( sf )
    max_hist = hist.GetMaximum()
    #hist.GetYaxis().SetRangeUser(0,1.2*max_hist)

    hist.SetMarkerColor(color)
    return hist

def compare_4masses(ntuple1,legend1,ntuple2, legend2,ntuple3, legend3,ntuple4,legend4, variable, PtCut,
			selection1,selection2,selection3,selection4,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_mass(ntuple1,variable, PtCut,selection1,binning, filename,ROOT.kMagenta-3,xaxis,title)
    l2 = produce_mass(ntuple2,variable, PtCut,selection2,binning, filename,ROOT.kBlue-9, xaxis, title)
    l3 = produce_mass(ntuple3,variable, PtCut,selection3,binning, filename,ROOT.kRed+3, xaxis, title)
    l4 = produce_mass(ntuple4,variable, PtCut,selection4,binning, filename,ROOT.kRed-3, xaxis, title)
    #frame.SetMaximum(1)
    #frame.SetMinimum(0)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('hist')
    l2.Draw('histsame')
    l3.Draw('histsame')
    l4.Draw('histsame')
    legend = ROOT.TLegend(0.7,0.75,0.88,0.88, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "f")
    legend.AddEntry(l2,legend2, "f")
    legend.AddEntry(l3,legend3, "f")
    legend.AddEntry(l4,legend4, "f")
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.Draw('sames')
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def compare_2masses(ntuple1,legend1,ntuple2, legend2, variable, PtCut,
                        selection1,selection2,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_mass(ntuple1,variable, PtCut,selection1,binning, filename,ROOT.kMagenta-3,xaxis,title)
    l2 = produce_mass(ntuple2,variable, PtCut,selection2,binning, filename,ROOT.kBlue-9, xaxis, title)

    if (l1.GetMaximum() > l2.GetMaximum()):
        frame.SetMaximum(1.2*l1.GetMaximum())
    else :
        frame.SetMaximum(1.2*l2.GetMaximum())

    #frame.SetMaximum(1)
    #frame.SetMinimum(0)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('ehist')
    l2.Draw('ehistsame')
    legend = ROOT.TLegend(0.1,0.75,0.28,0.88, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "f")
    legend.AddEntry(l2,legend2, "f")
    legend.Draw('sames')
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_sumPlot_mass(tree, variable1, variable2, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "(%s+%s)/1000000>>htemp(%s)" % (variable1,variable2, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def sum_mass(ntuple,variable1,variable2, PtCut,selection,binning, filename,color, xaxis, title):
    hist = make_sumPlot_mass(ntuple, variable1, variable2, selection,binning,xaxis, title)

    rebin = 1
    fcolor = 0 # ROOT.kGreen+1
    lcolor = color
    sf = 1
    fillStyle = 1

    hist.Rebin( rebin )
    hist.SetFillColor( fcolor )
    hist.SetLineColor( lcolor )
    hist.SetLineWidth( 2 )
    hist.Scale( sf )
    max_hist = hist.GetMaximum()
    #hist.GetYaxis().SetRangeUser(0,1.2*max_hist)

    hist.SetMarkerColor(color)
    return hist

def reconstruct_genMass(ntuple,legend1, legend2, legend3, variable1, variable2, PtCut,
                        selection,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_mass(ntuple,variable1, PtCut,selection,binning, filename,ROOT.kMagenta-3,xaxis,title)
    l2 = produce_mass(ntuple,variable2, PtCut,selection,binning, filename,ROOT.kBlue-9, xaxis, title)
    l3 = sum_mass(ntuple,variable1,variable2, PtCut,selection,binning, filename,ROOT.kRed+3, xaxis, title)

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
    l1.Draw('ehist')
    l2.Draw('ehistsame')
    l3.Draw('ehistsame')
    legend = ROOT.TLegend(0.1,0.75,0.28,0.9, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "f")
    legend.AddEntry(l2,legend2, "f")
    legend.AddEntry(l3,legend3, "f")
    legend.Draw('sames')
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def make_charge_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s/genCharge>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def produce_charge(ntuple, variable, PtCut,selection,binning, filename,color,xaxis,title):
    hist = make_charge_plot(ntuple, variable, selection,binning,xaxis, title)

    rebin = 1
    fcolor = 0 #ROOT.kGreen+1
    lcolor = color
    sf = 1
    fillStyle = 1

    hist.Rebin( rebin )
    hist.SetFillColor( fcolor )
    hist.SetLineColor( lcolor )
    hist.SetLineWidth( 2 )
    hist.Scale( sf )
    max_hist = hist.GetMaximum()
    #hist.GetYaxis().SetRangeUser(0,1.2*max_hist)
    hist.SetOption("E")

    hist.SetMarkerColor(color)
    return hist

def compare_2charge(ntuple1,legend1,ntuple2, legend2, variable, PtCut,
                        selection1,selection2,
                        binning, filename,
                        title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_charge(ntuple1,variable, PtCut,selection1,binning, filename,ROOT.kMagenta-3,xaxis,title)
    l2 = produce_charge(ntuple2,variable, PtCut,selection2,binning, filename,ROOT.kBlue-9, xaxis, title)

    if (l1.GetMaximum() > l2.GetMaximum()):
	frame.SetMaximum(1.2*l1.GetMaximum())
    else :
	frame.SetMaximum(1.2*l2.GetMaximum())

    #frame.SetMaximum(1)
    #frame.SetMinimum(0)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('ehist')
    l2.Draw('ehistsame')
    legend = ROOT.TLegend(0.7,0.75,0.88,0.88, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "f")
    legend.AddEntry(l2,legend2, "f")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','genPt', 20,
		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
		    [30, 0, 450],#variable, ptcut, binning
                    'iso_effi_pT',#filename
                    "Tau Efficiency",#title
                    "gen Tau p_{T}",#xaxis
                    "efficiency" #yaxis
)

#compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','genPt',20,
#		    "goodReco==1","genMatchedTau==1",[30,0,450],"goodReco==1","genMatchedTau==1",
#		    'PtSum_effi_pT',
#		    "Tau Efficiency",
#		    "gen Tau p_{T} (GeV)",
#		    "efficiency"
#)

#compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','genPt',20,
#		    "goodReco==1","genMatchedTau==1",[30,0,450],"goodReco==1","genMatchedTau==1",
#                    'Mu_effi_pT',
#                    "Tau Efficiency",
#                    "gen Tau p_{T} (GeV)",
#                    "efficiency"
#)

#compare_3efficiencies(EleVLooseMVA5,'againstElectronVLooseMVA5',EleLooseMVA5,'againstElectronLooseMVA5',EleMediumMVA5,'againstElectronMediumMVA5','genPt',20,
#		    "goodReco==1","genMatchedTau==1",[30,0,450],"goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#                    'Ele_effi_pT',
#                    "Tau Efficiency",
#                    "gen Tau p_{T} (GeV)",
#                    "efficiency"
#)

compare_efficiencies(byLooseCmbIso3,'New DMF (loose)',OldDMF,'Old DMF (tight)','genPt',20,
		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
		    [30,0,450],
                    'DMF_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_4efficiencies(byLooseCmbIso3,'kOneProng1PiZero (new DMF)',OldDMF,'kOneProng1PiZero (old DMF)',byLooseCmbIso3,'kThreeProng0PiZero (new DMF)',OldDMF,'kThreeProng0PiZero (old DMF)','genPt',20,
		    "goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
		    [30,0,450],
		    'decayMode_DMF_effi_pT',
		    "Tau Efficiency",
		    "gen Tau p_{T} (GeV)",
		    "efficiency"
)

compare_3efficiencies(byLooseCmbIso3,'kOneProng0PiZero',byLooseCmbIso3,'kOneProng1PiZero',byLooseCmbIso3,'kThreeProng0PiZero','genPt',20,
                    "goodReco==1 && genDecayMode==0","genMatchedTau==1 && genDecayMode==0","goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
                    [30,0,450],
                    'decayMode_looseNewDMF_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_3efficiencies(OldDMF,'kOneProng0PiZero',OldDMF,'kOneProng1PiZero',OldDMF,'kThreeProng0PiZero','genPt',20,
                    "goodReco==1 && genDecayMode==0","genMatchedTau==1 && genDecayMode==0","goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
                    [30,0,450],
                    'decayMode_looseOldDMF_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_3efficiencies(NewDMF,'kOneProng0PiZero',NewDMF,'kOneProng1PiZero',NewDMF,'kThreeProng0PiZero','genPt',20,
                    "goodReco==1 && genDecayMode==0","genMatchedTau==1 && genDecayMode==0","goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
                    [30,0,450],
                    'decayMode_nocutNewDMF_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_efficiencies(byLooseCmbIso3,'kOneProng2PiZero',byLooseCmbIso3,'kThreeProng1PiZero','genPt',20,
                    "goodReco==1 && genDecayMode==2","genMatchedTau==1 && genDecayMode==2","goodReco==1 && genDecayMode==11","genMatchedTau==1 && genDecayMode==11",
                    [30,0,450],
                    'decayMode_badModesLoose_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_efficiencies(OldDMF,'kOneProng2PiZero',OldDMF,'kThreeProng1PiZero','genPt',20,
                    "goodReco==1 && genDecayMode==2","genMatchedTau==1 && genDecayMode==2","goodReco==1 && genDecayMode==11","genMatchedTau==1 && genDecayMode==11",
                    [30,0,450],
                    'decayMode_badModesOld_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_efficiencies(NewDMF,'kOneProng2PiZero',NewDMF,'kThreeProng1PiZero','genPt',20,
                    "goodReco==1 && genDecayMode==2","genMatchedTau==1 && genDecayMode==2","goodReco==1 && genDecayMode==11","genMatchedTau==1 && genDecayMode==11",
                    [30,0,450],
                    'decayMode_badModesNew_effi_pT',
                    "Tau Efficiency",
                    "gen Tau p_{T} (GeV)",
                    "efficiency"
)

compare_efficiencies(byLooseCmbIso3,'Loose Iso Cut (New DMF)',NewDMF,'No Iso Cut (New DMF)','genPt',20,
		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
		    [30,0,450],
		    'isono_effi_pT',
		    "Tau Efficiency",
		    "gen Tau p_{T} (GeV)",
		    "efficiency"
)

plot_ratio(byLooseCmbIso3,'Loose Iso Cut (New DMF)', 'genPt', 'ratioPt',
                        "goodReco==1",
                        [30,0,450], 'ratio_pT',
                        "Pt ratio", "generated P_{T} [GeV]", "<P_{T}^{rec}/P_{T}^{gen}>")

#compare_4masses(byLooseCmbIso3,'kOneProng1PiZero (new DMF)',OldDMF,'kOneProng1PiZero (old DMF)',byLooseCmbIso3,'kThreeProng0PiZero (new DMF)',OldDMF,'kThreeProng0PiZero (old DMF)',
#                    'tauMass',20,
#                    "goodReco==1","genMatchedTau==1 && genDecayMode==1","goodReco==1","genMatchedTau==1 && genDecayMode==1","goodReco==1","genMatchedTau==1 && genDecayMode==10","goodReco==1","genMatchedTau==1 && genDecayMode==10",
#                    [30,0,450],
#                    'decayMode_DMF_effi_pT',
#                    "Tau Efficiency",
#                    "gen Tau p_{T} (GeV)",
#                    "efficiency"
#)

## eta plots
#compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','genEta', 20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#		    [20,-2.4,2.4],#variable, ptcut, binning
#                    'iso_effi_eta',#filename
#                    "Tau Efficiency",#title
#                    "gen Tau Eta",#xaxis
#                    "efficiency" #yaxis             
#)

#compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','genEta',20,
#	            "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#		    [20,-2.4,2.4],
#                    'PtSum_effi_eta',
#                    "Tau Efficiency",
#                    "gen Tau Eta",
#                    "efficiency"
#)
#
#compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','genEta',20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#		    [20,-2.4,2.4],
#                    'Mu_effi_eta',
#                    "Tau Efficiency",
#                    "gen Tau Eta",
#                    "efficiency"
#)
#
#compare_3efficiencies(EleVLooseMVA5,'againstElectronVLooseMVA5',EleLooseMVA5,'againstElectronLooseMVA5',EleMediumMVA5,'againstElectronMediumMVA5','genEta',20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#		    [20,-2.4,2.4],
#                    'Ele_effi_eta',
#                    "Tau Efficiency",
#                    "gen Tau Eta",
#                    "efficiency"
#)
#compare_4efficiencies(byLooseCmbIso3,'kOneProng0PiZero',byLooseCmbIso3,'kOneProng1PiZero',byLooseCmbIso3,'kTwoProng0PiZero',byLooseCmbIso3,'kThreeProng0PiZero',
#		    'genEta',20,
#		    "goodReco==1 && genDecayMode==0","genMatchedTau==1 && genDecayMode==0","goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==5","genMatchedTau==1 && genDecayMode==5","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
#		    [30,-2.4,2.4],
#                    'decayMode_effi_eta',
#                    "Tau Efficiency",
#                    "gen Tau Eta",
#                    "efficiency"
#)

#compare_efficiencies(byLooseCmbIso3,'Loose Iso Cut (New DMF)',NewDMF,'No Iso Cut (New DMF)','genEta',20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#	            [30,-2.4,2.4],
#                    'isono_effi_eta',
#                    "Tau Efficiency",
#                    "gen Tau Eta",
#                    "efficiency"
#)

## nvtx plots
compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','nvtx', 20,
		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
		    [20,0,35],#variable, ptcut, binning
                    'iso_effi_nvtx',#filename
                    "Tau Efficiency",#title
                    "N_{vtx}",#xaxis
                    "efficiency" #yaxis             
)

#compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','nvtx',20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#                    'PtSum_effi_nvtx',
#                    "Tau Efficiency",
#                    "N_{vtx}",
#                    "efficiency"
#)

#compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','nvtx',20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#		    [20,0,35],
#                    'Mu_effi_nvtx',
#                    "Tau Efficiency",
#                    "N_{vtx}",
#                    "efficiency"
#)

#compare_3efficiencies(EleVLooseMVA5,'againstElectronVLooseMVA5',EleLooseMVA5,'againstElectronLooseMVA5',EleMediumMVA5,'againstElectronMediumMVA5','nvtx',20,
#		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
#		    [20,0,35],
#                    'Ele_effi_nvtx',
#                    "Tau Efficiency",
#                    "N_{vtx}",
#                    "efficiency"
#)

compare_4efficiencies(byLooseCmbIso3,'kOneProng0PiZero',byLooseCmbIso3,'kOneProng1PiZero',byLooseCmbIso3,'kTwoProng0PiZero',byLooseCmbIso3,'kThreeProng0PiZero',
		    'nvtx',20,
		    "goodReco==1 && genDecayMode==0","genMatchedTau==1 && genDecayMode==0","goodReco==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1","goodReco==1 && genDecayMode==5","genMatchedTau==1 && genDecayMode==5","goodReco==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
		    [30,0,35],
                    'decayMode_effi_nvtx',
                    "Tau Efficiency",
                    "N_{vtx}",
                    "efficiency"
)

compare_efficiencies(byLooseCmbIso3,'Loose Iso Cut (New DMF)',NewDMF,'No Iso Cut (New DMF)','nvtx',20,
		    "goodReco==1","genMatchedTau==1","goodReco==1","genMatchedTau==1",
		    [30,0,35],
                    'isono_effi_nvtx',
                    "Tau Efficiency",
                    "N_{vtx}",
                    "efficiency"
)

## mass plots
# mass for one prong 0 pi0, one prong 1 pi0, two prong n pi0, three prong n pi0

compare_4masses(byLooseCmbIso3,'kOneProng0PiZero',byLooseCmbIso3,'kOneProng1PiZero',byLooseCmbIso3,'kTwoProng0PiZero',byLooseCmbIso3,'kThreeProng0PiZero','tauMass', 20,
			"goodReco==1 && tauDecayMode==0","goodReco==1 && tauDecayMode==1", "goodReco==1 && tauDecayMode==5", "goodReco==1 && tauDecayMode == 10", 
                        [101,970,1070],
			'newDMF_mass',
                        "Tau Mass: New DMF",
			"Mass (MeV)",
			"Events"
)

compare_2masses(byLooseCmbIso3,'Reco Tau (New DMF)',OldDMF,'Reco Tau (Old DMF)',
                        'tauMass', 20,
                        "goodReco==1 && tauDecayMode==0","goodReco==1 && tauDecayMode==0",
                        [101,970,1070],
                        '1P0pi_tauMass',
                        "Reco Tau Mass: kOneProng0PiZero",
                        "Reco Tau Mass (MeV)",
                        "Events"
)

#compare_2masses(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (tight)',
#                        'tauMass', 20,
#                        "goodReco==1 && tauDecayMode==1","goodReco==1 && tauDecayMode==1",
#                        [101,970,1070],
#                        '1P1pi_tauMass',
#                        "Reco Tau Mass: kOneProng1PiZero",
#                        "Reco Tau Mass (MeV)",
#                        "Events"
#)

compare_2masses(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)','tauMass', 20,
                        "goodReco==1 && tauDecayMode==5","goodReco==1 && tauDecayMode==5",
                        [101,970,1070],
                        '2P0pi_tauMass',
                        "Reco Tau Mass: kTwoProng0PiZero",
                        "Reco Tau Mass (MeV)",
                        "Events"
)

compare_2masses(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauMass', 20,
                        "goodReco==1 && tauDecayMode==10","goodReco==1 && tauDecayMode==10",
                        [101,970,1070],
                        '3P0pi_tauMass',
                        "Reco Tau Mass: kThreeProng0PiZero",
                        "Reco Tau Mass (MeV)",
                        "Events"
)

compare_2masses(byLooseCmbIso3,'Reco Tau (New DMF)',OldDMF,'Reco Tau (Old DMF)',
                        'genMass', 20,
                        "genMatchedTau==1 && genDecayMode==0","genMatchedTau==1 && genDecayMode==0",
                        [101,970,1070],
                        '1P0pi_genMass',
                        "Gen Tau Mass: kOneProng0PiZero",
                        "Gen Mass (MeV)",
                        "Events"
)

compare_2masses(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (tight)',
                        'genMass', 20,
                        "genMatchedTau==1 && genDecayMode==1","genMatchedTau==1 && genDecayMode==1",
                        [101,970,1070],
                        '1P1pi_genMass',
                        "Gen Tau Mass: kOneProng1PiZero",
                        "Gen Mass (MeV)",
                        "Events"
)

compare_2masses(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'genMass', 20,
                        "genMatchedTau==1 && genDecayMode==5","genMatchedTau==1 && genDecayMode==5",
                        [101,970,1070],
                        '2P0pi_genMass',
                        "Gen Mass: kTwoProng0PiZero",
                        "Gen Mass (MeV)",
                        "Events"
)

compare_2masses(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'genMass', 20,
                        "genMatchedTau==1 && genDecayMode==10","genMatchedTau==1 && genDecayMode==10",
                        [101,970,1070],
                        '3P0pi_genMass',
                        "Gen Tau Mass: kThreeProng0PiZero",
                        "Gen Mass (MeV)",
                        "Events"
)

compare_2masses(byLooseCmbIso3,'Gen Tau Mass',OldDMF,'Gen Tau Mass',
                        'genMass', 20,
                        "genMatchedTau==1","genMatchedTau==1",
                        [101,970,1070],
                        'all_genMass',
                        "Gen Tau Mass: kThreeProng0PiZero",
                        "Gen Mass (MeV)",
                        "Events"
)

#reconstruct_genMass(byLooseCmbIso3,'Visible Mass', 'Neutrino Mass', 'Total Mass', 'genMass', 'genNuMass', 20,
#                       	"genMatchedTau==1",
#                        [101,970,1070],
#			'genMass_sum',
#                        "Gen Tau Mass (loose iso cut)", "Mass (MeV)", "Events")

## Charge plots

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargePt', 20,
                        "goodReco==1 && tauDecayMode==5","goodReco==1 && tauDecayMode==5",
                        [50,-6,6],
                        '2P0pi_chargePt',
                        "Tau Charge (High Pt Track): kTwoProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargeSum', 20,
                        "goodReco==1 && tauDecayMode==5","goodReco==1 && tauDecayMode==5",
                        [50,-6,6],
                        '2P0pi_chargeSum',
                        "Tau Charge (Sum of Tracks): kTwoProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargeDirect', 20,
                        "goodReco==1 && tauDecayMode==5","goodReco==1 && tauDecayMode==5",
                        [50,-6,6],
                        '2P0pi_chargeDirect',
                        "Tau Charge (tau.charge()): kTwoProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargePt', 20,
                        "goodReco==1 && tauDecayMode==0","goodReco==1 && tauDecayMode==0",
                        [50,-6,6],
                        '1P0pi_chargePt',
                        "Tau Charge (High Pt Track): kOneProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargeSum', 20,
                        "goodReco==1 && tauDecayMode==0","goodReco==1 && tauDecayMode==0",
                        [50,-6,6],
                        '1P0pi_chargeSum',
                        "Tau Charge (Sum of Tracks): kOneProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargeDirect', 20,
                        "goodReco==1 && tauDecayMode==0","goodReco==1 && tauDecayMode==0",
                        [50,-6,6],
                        '1P0pi_chargeDirect',
                        "Tau Charge (tau.charge()): kOneProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

#compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
#                        'tauChargePt', 20,
#                        "goodReco==1 && tauDecayMode==1","goodReco==1 && tauDecayMode==1",
#                        [50,-6,6],
#                        '1P1pi_chargePt',
#                        "Tau Charge (High Pt Track): kOneProng1PiZero",
#                        "Tau Charge / Gen Charge",
#                        "Events"
#)

#compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
#                        'tauChargeSum', 20,
#                        "goodReco==1 && tauDecayMode==1","goodReco==1 && tauDecayMode==1",
#                        [50,-6,6],
#                        '1P1pi_chargeSum',
#                        "Tau Charge (Sum of Tracks): kOneProng1PiZero",
#                        "Tau Charge / Gen Charge",
#                        "Events"
#)

#compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
#                        'tauChargeDirect', 20,
#                        "goodReco==1 && tauDecayMode==1","goodReco==1 && tauDecayMode==1",
#                        [50,-6,6],
#                        '1P1pi_chargeDirect',
#                        "Tau Charge (tau.charge()): kOneProng1PiZero",
#                        "Tau Charge / Gen Charge",
#                        "Events"
#)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargePt', 20,
                        "goodReco==1 && tauDecayMode==10","goodReco==1 && tauDecayMode==10",
                        [50,-6,6],
                        '3P0pi_chargePt',
                        "Tau Charge (High Pt Track): kThreeProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

time.sleep(1)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargeSum', 20,
                        "goodReco==1 && tauDecayMode==10","goodReco==1 && tauDecayMode==10",
                        [50,-6,6],
                        '3P0pi_chargeSum',
                        "Tau Charge (Sum of Tracks): kThreeProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'New DMF (Loose)',OldDMF,'Old DMF (Tight)',
                        'tauChargeDirect', 20,
                        "goodReco==1 && tauDecayMode==10","goodReco==1 && tauDecayMode==10",
                        [50,-6,6],
                        '3P0pi_chargeDirect',
                        "Tau Charge (tau.charge()): kThreeProng0PiZero",
                        "Tau Charge / Gen Charge",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'highest pt track',byLooseCmbIso3,'sum of tracks',
                        'tauDecayMode', 20,
                        "abs(tauChargePt)>1","abs(tauChargePt)>1",
                        [100,-10,10],
                        'decayMode_forCharge2',
                        "Tau decay mode",
                        "Tau decay mode",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'OneProng0PiZero',byLooseCmbIso3,'all',
                        'tauPt', 20,
                        "goodReco==1  && (tauChargeDirect/genCharge) == -1","goodReco==1 && (tauChargeDirect/genCharge)==-1",
                        [100,-450,450],
                        'bad_charges_byPt',
                        "Tau decay mode",
                        "Tau decay mode",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'OneProng0PiZero',byLooseCmbIso3,'all',
                        'tauEta', 20,
                        "goodReco==1  && (tauChargeDirect/genCharge) == -1","goodReco==1 && (tauChargeDirect/genCharge)==-1",
                        [100,-2.4,2.4],
                        'bad_charges_byEta',
                        "Tau decay mode",
                        "Tau decay mode",
                        "Events"
)

compare_2charge(byLooseCmbIso3,'OneProng0PiZero',byLooseCmbIso3,'all',
                        'nvtx', 20,
                        "goodReco==1  && (tauChargeDirect/genCharge) == -1","goodReco==1 && (tauChargeDirect/genCharge)==-1",
                        [100,-35,35],
                        'bad_charges_byNvtx',
                        "Tau decay mode",
                        "Tau decay mode",
                        "Events"
)
