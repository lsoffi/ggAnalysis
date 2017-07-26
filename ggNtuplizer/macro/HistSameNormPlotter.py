#! /usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-c", "--colors", dest="colors", type=int, default=[2,4,8,1,6,9,3,7,5], nargs='*', help="Histgram Colors")
parser.add_argument("-f", "--filelist", dest="filelist", default=["root://cmsxrootd.fnal.gov//store/user/drberry/HATS/ggntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root","root://cmsxrootd.fnal.gov//store/user/drberry/HATS/ggntuples/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","root://cmsxrootd.fnal.gov//store/user/drberry/HATS/ggntuples/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root"], nargs='*', help="List of input ROOT files")
parser.add_argument("-o", "--options", dest="options", default="norm", help="Cuts placed on TTree")
parser.add_argument("-v", "--variable", dest="variable", default=["elePt"], nargs='*', help="Variable to plot")
parser.add_argument("-t", "--ttree", dest="ttree", default="ggNtuplizer/EventTree", help="TTree Name")
parser.add_argument("-x", "--cuts", dest="cuts", default="", help="Cuts placed on TTree")
parser.add_argument("--outputtype", dest="outputtype", default=".png", help="Outputfile type: .png")
parser.add_argument("--logx", dest="logx", action="store_true", default=False, help="Set TCanvas to LogX")
parser.add_argument("--logy", dest="logy", action="store_true", default=False, help="Set TCanvas to LogY")
parser.add_argument("--xsize", dest="xsize", type=int, default=800, help="X width of TCanvas")
parser.add_argument("--ysize", dest="ysize", type=int, default=800, help="Y width of TCanvas")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(000000)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.UseCurrentStyle()

can = ROOT.TCanvas("Plots","Plots",args.xsize,args.ysize)
if args.logx: can.SetLogx()
if args.logy: can.SetLogy()
print args.logy
for var in args.variable:
    index = 0
    for filename in args.filelist:
        tfile = ROOT.TFile.Open(filename)
        ttree = tfile.Get(args.ttree)
        ttree.SetLineColor(args.colors[index])
        ttree.SetLineWidth(2)
        if index==0:
            ttree.Draw(var,args.cuts,args.options)
        else:
            ttree.Draw(var,args.cuts,args.options+" same")
        index += 1
    can.SaveAs(var+args.outputtype)
