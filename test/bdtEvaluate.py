import ROOT
import array
import glob
import math
import re
import sys
import bdtCommon


ROOT.gStyle.SetOptDate(0)
ROOT.gStyle.SetHistLineWidth(2)
def fixLegend(leg, opt):
    for p in leg.GetListOfPrimitives():
        p.SetOption(opt)


def makeROC(trees, idConfig):
    npts = 100
    x, y, exl, exh, eyl, eyh = map(lambda _: array.array('d', [0.]*npts), xrange(6))

    ROOT.gROOT.cd()
    hdummy = ROOT.TH1D("dummy", "", npts, -1., 1.)
    for tree in trees:
        ROOT.gROOT.cd()
        tree.Draw("(%s) >>+dummy" % idConfig.name, idConfig.trueCut, "goff")

    for i in range(npts):
        nPass, nAll = hdummy.Integral(i+1, -1), hdummy.Integral()
        y[i] = nPass / nAll
        eyl[i] = y[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False)
        eyh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) - y[i]

    hdummy.Reset()
    nevents = 0.
    for tree in trees:
        ROOT.gROOT.cd()
        tree.Draw("(%s) >>+dummy" % idConfig.name, idConfig.bkgCut, "goff")
        nevents += tree.GetEntries()

    for i in range(npts):
        nPass, nAll = hdummy.Integral(i+1, -1), hdummy.Integral()
        x[i] = nPass / nevents
        exl[i] = x[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) * nAll / nevents
        exh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) * nAll / nevents - x[i]

    gROC = ROOT.TGraphAsymmErrors(len(x), x, y, exl, exh, eyl, eyh)
    gROC.GetXaxis().SetTitle("Background Photons / Event")
    gROC.GetYaxis().SetTitle("Signal Photon Efficiency")
    gROC.SetMarkerSize(0.6)
    gROC.SetLineWidth(2)
    #gROC.SetLineStyle(ROOT.kDashed)
    return gROC


idConfig = filter(lambda x: x.name==sys.argv[-1], bdtCommon.idconfigs)
if len(idConfig)==1:
    idConfig = idConfig[0]
else:
    raise Exception("Check idConfigs in bdtCommon.py")

filenames = [
    "/data/ncsmith/932phoID_round2/GJets.root",
    "/data/ncsmith/932phoID_round2/DiPhotonSherpa.root",
    "/data/ncsmith/932phoID_round2/QCD_1.root",
    "/data/ncsmith/932phoID_round2/DY2J.root",
]
files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
for i, tree in enumerate(trees):
    tree.AddFriend("tree", idConfig.friendTreeFile(filenames[i]))


fout = ROOT.TFile.Open("plots_%s.root" % idConfig.name, "recreate")

roc = makeROC(trees[0:1], idConfig)
roc.SetNameTitle(idConfig.name, idConfig.name)
fout.cd()
roc.Write()

gencats = {
    "promptUnconvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && !(abs(gen_conversionZ[gedReco_iGen])<300 && gen_conversionRho[gedReco_iGen]<120)",
    "promptConvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && (abs(gen_conversionZ[gedReco_iGen])<300 && gen_conversionRho[gedReco_iGen]<120)",
    "nonpromptPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && !gen_isPromptFinalState[gedReco_iGen]",
    "ele": "gedReco_iGen>=0 && abs(gen_id[gedReco_iGen]) == 11",
    "fake": "gedReco_iGen<0",
}
for name, cut in gencats.iteritems():
    cut = idConfig.trainingCut + "&&" + cut
    if type(idConfig) is bdtCommon.EndcapIDConfig:
        cut = cut.replace("gedReco", "localReco")
    hmva = ROOT.TH1F("mva_"+name, "%s;BDT output;Normalized" % name, 100, -1, 1)
    for tree in trees:
        tree.Draw("({0}) >>+mva_{1}".format(idConfig.name, name), cut, "goff")
    hmva.Write()
