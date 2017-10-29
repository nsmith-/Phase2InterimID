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


def makeROC(trees, idConfig, trueCut, bkgCut):
    npts = 100
    x, y, exl, exh, eyl, eyh = map(lambda _: array.array('d', [0.]*npts), xrange(6))

    hdummy = ROOT.TH1D("dummy", "", npts, -1., 1.)
    for tree in trees:
        tree.Draw("(%s) >>+dummy" % idConfig.name, trueCut, "goff")

    for i in range(npts):
        nPass, nAll = hdummy.Integral(i+1, -1), hdummy.Integral()
        y[i] = nPass / nAll
        eyl[i] = y[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False)
        eyh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) - y[i]

    hdummy.Reset()
    nevents = 0.
    for tree in trees:
        tree.Draw("(%s) >>+dummy" % idConfig.name, bkgCut, "goff")
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


def makeEff(trees, histDef, ncut, dcut):
    htmpnum = histDef.Clone("tmpnum")
    htmpdenom = histDef.Clone("tmpdenom")
    for tree in trees:
        tree.Draw(histDef.GetTitle()+">>tmpnum", "(%s)&&(%s)" % (dcut, ncut), "goff")
        tree.Draw(histDef.GetTitle()+">>tmpdenom", "%s" % dcut, "goff")
    hEff = ROOT.TEfficiency(htmpnum, htmpdenom)
    return hEff


idConfig = filter(lambda x: x.name==sys.argv[-1], bdtCommon.idconfigs)
if len(idConfig)==1:
    idConfig = idConfig[0]
else:
    raise Exception("Check idConfigs in bdtCommon.py")

filenames = bdtCommon.allInputFiles
files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
for i, tree in enumerate(trees):
    tree.AddFriend("tree", idConfig.friendTreeFile(filenames[i]))


fout = ROOT.TFile.Open("plots_%s.root" % idConfig.name, "recreate")


# ROC curves
rocGJets = makeROC(trees[0:1], idConfig, idConfig.trueCut, idConfig.bkgCut)
rocGJets.SetNameTitle("gjetsROC", idConfig.name)
rocGJets.Write()
noe = idConfig.bkgCut + "&& (gedReco_iGen<0||abs(gen_id[gedReco_iGen])!=11)"
juste = idConfig.bkgCut + "&& gedReco_iGen>=0 && abs(gen_id[gedReco_iGen])==11 && gen_isPromptFinalState[gedReco_iGen]"
if isinstance(idConfig, bdtCommon.EndcapIDConfig):
    noe = noe.replace("gedReco", "localReco")
    juste = juste.replace("gedReco", "localReco")
rocNoE = makeROC(trees, idConfig, idConfig.trueCut, noe)
rocNoE.SetNameTitle("noeROC", idConfig.name)
rocNoE.Write()
rocJustE = makeROC(trees, idConfig, idConfig.trueCut, juste)
rocJustE.SetNameTitle("justeROC", idConfig.name)
rocJustE.Write()


# MVA & input values for gen categories
gencats = {
    "promptUnconvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) != 11 && !(abs(gen_conversionZ[gedReco_iGen])<300 && gen_conversionRho[gedReco_iGen]<120)",
    "promptConvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) != 11 && (abs(gen_conversionZ[gedReco_iGen])<300 && gen_conversionRho[gedReco_iGen]<120)",
    "promptPhoEleMom": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) == 11",
    "nonpromptPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && !gen_isPromptFinalState[gedReco_iGen]",
    "ele": "gedReco_iGen>=0 && abs(gen_id[gedReco_iGen]) == 11",
    "fake": "gedReco_iGen<0",
}
for name, cut in gencats.iteritems():
    cut = idConfig.trainingCut + "&&" + cut
    if isinstance(idConfig, bdtCommon.EndcapIDConfig):
        cut = cut.replace("gedReco", "localReco")
    hmva = ROOT.TH1F("mva_"+name, "%s;BDT output;Normalized" % name, 100, -1, 1)
    for tree in trees:
        tree.Draw("({0}) >>+mva_{1}".format(idConfig.name, name), cut, "goff")
    hmva.Write()


# Gen to x efficiencies
effpt_def = ROOT.TH1D("effpt_def", "gen_pt;Gen. p_{T}^{#gamma} [GeV];Efficiency", 20, 10, 150)
effeta_def = ROOT.TH1D("effeta_def", "abs(gen_eta);Gen. |#eta^{#gamma}|;Efficiency", 32, 0, 3.2)
denom = "gen_id == 22 && gen_isPromptFinalState && abs(gen_parentId)!=11"
if isinstance(idConfig, bdtCommon.EndcapIDConfig):
    denom += " && abs(gen_eta)>1.5"
else:
    denom += " && abs(gen_eta)<1.4"

efficiencies = {
    "reco": "gen_gedRecoDeltaR < 0.1",
    "reco0p3": "gen_gedRecoDeltaR < 0.3",
    "idmva0": "gen_gedRecoDeltaR < 0.1 && %s[gen_iGedReco] > 0." % idConfig.name,
    "idmva0p1": "gen_gedRecoDeltaR < 0.1 && %s[gen_iGedReco] > 0.1" % idConfig.name,
}
for name, cut in efficiencies.iteritems():
    if isinstance(idConfig, bdtCommon.EndcapIDConfig):
        cut = cut.replace("gedReco", "localReco")
        cut = cut.replace("GedReco", "LocalReco")
    eff = makeEff(trees[0:2], effpt_def, cut, denom)
    eff.SetNameTitle("effPt_%s" % name, name)
    eff.Write()
    eff2 = makeEff(trees[0:2], effeta_def, cut, denom)
    eff2.SetNameTitle("effEta_%s" % name, name)
    eff2.Write()


