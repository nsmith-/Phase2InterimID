import ROOT
import array
import glob
import math
import re
import sys
from bdtCommon import BarrelIDConfig, EndcapIDConfig, makeReader


ROOT.gStyle.SetOptDate(0)
ROOT.gStyle.SetHistLineWidth(2)
def fixLegend(leg, opt):
    for p in leg.GetListOfPrimitives():
        p.SetOption(opt)


def makeROC(filesig, filebkg, cuts, truecut, bkgcut):
    x, y, exl, exh, eyl, eyh = map(lambda _: array.array('d', [0.]*len(cuts)), xrange(6))

    file = ROOT.TFile.Open(filesig)
    tree = file.Get("ntupler/photons")
    hdummy = ROOT.TH1D("dummy", "", 2, 0, 2)
    for i, cut in enumerate(cuts):
        hdummy.Reset()
        tree.Draw("(%s) >>dummy" % cut, truecut, "goff")
        nPass, nAll = hdummy.GetBinContent(2), hdummy.Integral()
        y[i] = nPass / nAll
        eyl[i] = y[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False)
        eyh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) - y[i]

    file = ROOT.TFile.Open(filebkg)
    tree = file.Get("ntupler/photons")
    nevents = tree.GetEntries()
    hdummy = ROOT.TH1D("dummy", "", 2, 0, 2)
    for i, cut in enumerate(cuts):
        hdummy.Reset()
        tree.Draw("(%s) >>dummy" % cut, bkgcut, "goff")
        nPass, nAll = hdummy.GetBinContent(2), hdummy.Integral()
        x[i] = nPass / nevents
        exl[i] = x[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) * nAll / nevents
        exh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) * nAll / nevents - x[i]

    gROC = ROOT.TGraphAsymmErrors(len(cuts), x, y, exl, exh, eyl, eyh)
    gROC.GetXaxis().SetTitle("Background Photons / Event")
    gROC.GetYaxis().SetTitle("Signal Photon Efficiency")
    #gROC.SetMarkerSize(0.6)
    #gROC.SetLineStyle(ROOT.kDashed)
    return gROC



infiles = [
    "/eos/cms/store/user/ncsmith/932phoID/RelValH125GGgluonfusion_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/171017_222812/0000/output_4.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/171017_222826/0000/output_1.root",
]
sigFile = infiles[0]
bkgFile = infiles[1]

makeReader(BarrelIDConfig)
makeReader(EndcapIDConfig)

mmin, mmax, mn = -0.15, 0.35, 12
mva = lambda cutfcn: [cutfcn+' > %f' % (mmin+mmax*x/float(mn-1)) for x in range(mn)]

rocBar = makeROC(sigFile, bkgFile, mva(BarrelIDConfig.evalFcn), BarrelIDConfig.trueCut, BarrelIDConfig.bkgCut)
rocBar.SetNameTitle("barrel", "Barrel")
rocBar.SetMarkerStyle(ROOT.kFullSquare)
rocBar.SetMarkerColor(ROOT.kBlue-3)
rocBar.SetLineColor(ROOT.kBlue-3)

rocEnd = makeROC(sigFile, bkgFile, mva(EndcapIDConfig.evalFcn), EndcapIDConfig.trueCut, EndcapIDConfig.bkgCut)
rocEnd.SetNameTitle("endcap", "Endcap")
rocEnd.SetMarkerStyle(ROOT.kFullSquare)
rocEnd.SetMarkerColor(ROOT.kBlue-3)
rocEnd.SetLineColor(ROOT.kBlue-3)

cRes = ROOT.TCanvas()
#cRes.SetLogx(True)
mg = ROOT.TMultiGraph()
mg.Add(rocBar, "ple")
mg.Add(rocEnd, "ple")
mg.Draw("alpe")
mg.GetXaxis().SetTitle("Background Photons / Event")
mg.GetYaxis().SetTitle("Signal Photon Efficiency")
l = cRes.BuildLegend(0.5, 0.2, 0.8, 0.5)
fixLegend(l, "ple")
# l.SetHeader("")
#tdr.formatHisto(mg)
mg.GetXaxis().SetRangeUser(0., 0.1)
mg.GetYaxis().SetRangeUser(0., 1.2)
mg.GetXaxis().SetTitleOffset(1.4)
#tdr.drawCMS()
#tdr.drawEnPu('200')

