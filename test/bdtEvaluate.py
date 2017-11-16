import ROOT
import array
import glob
import math
import re
import sys
import bdtCommon

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

def makeROC(trees, idConfig, trueCut, bkgCut, bkgPerEvent=True):
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
        if bkgPerEvent:
            x[i] = nPass / nevents
            exl[i] = x[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) * nAll / nevents
            exh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) * nAll / nevents - x[i]
        else:
            x[i] = nPass / nAll
            exl[i] = x[i] - ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False)
            exh[i] = ROOT.TEfficiency.ClopperPearson(nAll, nPass, 0.68, False) - x[i]

    gROC = ROOT.TGraphAsymmErrors(len(x), x, y, exl, exh, eyl, eyh)
    gROC.GetXaxis().SetTitle("Background / Event" if bkgPerEvent else "Background Efficiency")
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


def effectiveSigma(hres, debug=None):
    if hres.GetEffectiveEntries() == 0. or hres.Integral() == 0.:
        return (-1., 0.)
    quantiles = array.array('d', [0., 0.68])
    qvals = array.array('d', [0., 0.])
    npts = hres.GetNbinsX()
    def evalpt(i):
        quantiles[0] = (1.-0.68)*i/float(npts-1)
        quantiles[1] = quantiles[0] + 0.68
        hres.GetQuantiles(2, qvals, quantiles)
    best = (0, float('inf'))
    for i in range(npts):
        evalpt(i)
        if qvals[1]-qvals[0] < best[1]:
            best = (i, qvals[1]-qvals[0])
    evalpt(best[0])
    xlo, xhi = qvals
    effectiveSigma = (xhi-xlo)/(xhi+xlo)
    errEstimate = effectiveSigma / math.sqrt(2*hres.GetEffectiveEntries())

    if debug:
        cold = ROOT.gPad.Pad()
        ROOT.gStyle.SetOptStat(111111)
        c = ROOT.TCanvas("tmp")
        hres.SetStats(True)
        hres.Draw()
        lsigma = ROOT.TLine()
        lsigma.SetLineWidth(3)
        lsigma.SetLineColorAlpha(ROOT.kBlue, 0.7)
        lsigma.DrawLine(xlo, hres.GetMinimum(), xlo, hres.GetMaximum())
        lsigma.DrawLine(xhi, hres.GetMinimum(), xhi, hres.GetMaximum())
        lsigma.SetLineColorAlpha(ROOT.kRed, 0.7)
        lsigma.DrawLine((xhi+xlo)/2., hres.GetMinimum(), (xhi+xlo)/2., hres.GetMaximum())
        c.Print(debug)
        ROOT.gStyle.SetOptStat(0)
        if cold:
            ROOT.SetOwnership(cold, False)
            cold.cd()
        c = None

    return (effectiveSigma, errEstimate)


def resolutionPlot(trees, name, eoetrue, cut, ptbinning, debugName):
    nbins = len(ptbinning)-1
    x, y, exl, exh, eyl, eyh = map(lambda _: array.array('d', [0.]*nbins), xrange(6))
    resbins = array.array('d', [0.5+i/199. for i in range(200)])
    hres = ROOT.TH2D("hres2d", "dummy;Pt;E/Etrue;Counts", len(ptbinning)-1, ptbinning, len(resbins)-1, resbins)
    for tree in trees[:2]:
        tree.Draw(eoetrue+">>+hres2d", cut, "goff")
    for i in range(nbins):
        x[i] = (ptbinning[i+1]+ptbinning[i])/2.
        exl[i] = (ptbinning[i+1]-ptbinning[i])/2.
        exh[i] = exl[i]
        effSig, effSigErr = effectiveSigma(hres.ProjectionY("tmptmp", i+1, i+1), debug="resolutions/%s_%s_ptbin%d.pdf" % (debugName, name, i))
        y[i] = effSig
        eyl[i] = effSigErr
        eyh[i] = effSigErr
    gSigma = ROOT.TGraphAsymmErrors(len(x), x, y, exl, exh, eyl, eyh)
    gSigma.SetNameTitle(name, name)
    gSigma.GetXaxis().SetTitle("Generated photon energy (GeV)")
    gSigma.GetYaxis().SetTitle("Energy resolution, #sigma_{eff}(E)/E")
    gSigma.SetMarkerSize(0.8)
    return gSigma


idConfig = filter(lambda x: x.name==sys.argv[-1], bdtCommon.idconfigs)
if len(idConfig)==1:
    idConfig = idConfig[0]
else:
    raise Exception("Check idConfigs in bdtCommon.py")

filenames = bdtCommon.allInputFiles
if "noPU" in idConfig.name:
    filenames = idConfig.inputFiles
files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
for i, tree in enumerate(trees):
    tree.AddFriend("tree", idConfig.friendTreeFile(filenames[i]))


fout = ROOT.TFile.Open("plots_%s.root" % idConfig.name, "recreate")


# ROC curves
rocGJets = makeROC(trees[0:1], idConfig, idConfig.trueCut, idConfig.bkgCut)
rocGJets.SetNameTitle("gjetsROC", idConfig.name)
rocGJets.Write()
rocGJetsRej = makeROC(trees[0:1], idConfig, idConfig.trueCut, idConfig.bkgCut, False)
rocGJetsRej.SetNameTitle("gjetsROCrej", idConfig.name)
rocGJetsRej.Write()
hiptcut = "&&localReco_pt>50." if isinstance(idConfig, bdtCommon.EndcapIDConfig) else "&&gedReco_pt>50."
rocGJetsHighPt = makeROC(trees[0:1], idConfig, idConfig.trueCut+hiptcut, idConfig.bkgCut+hiptcut)
rocGJetsHighPt.SetNameTitle("gjetsHiPtROC", idConfig.name)
rocGJetsHighPt.Write()
rocGJetsHighPtRej = makeROC(trees[0:1], idConfig, idConfig.trueCut+hiptcut, idConfig.bkgCut+hiptcut, False)
rocGJetsHighPtRej.SetNameTitle("gjetsHiPtROCrej", idConfig.name)
rocGJetsHighPtRej.Write()
rocGJetsNoPre = makeROC(trees[0:1], idConfig, idConfig.trueCut.replace(idConfig.preselection+" && ", ""), idConfig.bkgCut.replace(idConfig.preselection+" && ", ""))
rocGJetsNoPre.SetNameTitle("gjetsROC_nopre", idConfig.name)
rocGJetsNoPre.Write()
rocGamma = makeROC(trees[0:2], idConfig, idConfig.trueCut, idConfig.bkgCut)
rocGamma.SetNameTitle("gammaROC", idConfig.name)
rocGamma.Write()
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
if True:
    gencats = {
        "promptUnconvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) != 11 && !(abs(gen_conversionZ[gedReco_iGen])<290 && gen_conversionRho[gedReco_iGen]<120)",
        "promptConvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) != 11 && (abs(gen_conversionZ[gedReco_iGen])<290 && gen_conversionRho[gedReco_iGen]<120)",
        "promptPhoEleMom": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) == 11",
        "nonpromptPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && !gen_isPromptFinalState[gedReco_iGen]",
        "phoPiZero": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_parentId[gedReco_iGen]==111",
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


# Define working points
def findWPs(roc):
    wp85, wp95 = None, None
    for i in range(roc.GetN()):
        if roc.GetY()[i] <= 0.85 and not wp85:
            wp85 = ROOT.TParameter("double")("wp85", -1. + 2.*i/float(roc.GetN()))
        if roc.GetY()[i] <= 0.95 and not wp95:
            wp95 = ROOT.TParameter("double")("wp95", -1. + 2.*i/float(roc.GetN()))
    if not wp85:
        wp85 = ROOT.TParameter("double")("wp85", -1.)
    if not wp95:
        wp95 = ROOT.TParameter("double")("wp95", -1.)
    return (wp85, wp95)

wp85, wp95 = findWPs(rocGJets)
wp85.Write()
wp95.Write()


# Gen to x efficiencies
ptbinning = array.array('d', [10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 250., 300., 350., 400., 450., 500., 750., 1000.])
effpt_def = ROOT.TH1D("effpt_def", "gen_pt;Generated p_{T}^{#gamma} [GeV];Efficiency", len(ptbinning)-1, ptbinning)
effeta_def = ROOT.TH1D("effeta_def", "abs(gen_eta);Generated |#eta^{#gamma}|;Efficiency", 32, 0, 3.2)
denom_pt = "gen_id == 22 && gen_isPromptFinalState && abs(gen_parentId)!=11"
denom_eta = denom_pt
pre = "1."
if isinstance(idConfig, bdtCommon.EndcapIDConfig):
    denom_pt += " && abs(gen_eta)>1.6 && abs(gen_eta)<2.9"
    denom_eta += " && gen_pt>30."
    pre = "abs(localReco_depthCompatibility[gen_iLocalReco]) < 20. && localReco_seedEnergyFH[gen_iLocalReco]/localReco_seedEnergyEE[gen_iLocalReco] < 20."
else:
    denom_pt += " && abs(gen_eta)<1.4"
    denom_eta += " && gen_pt>30."
efficiencies = {
    "reco": "gen_gedRecoDeltaR < 0.1",
    "recoPre": "gen_gedRecoDeltaR < 0.1 && " + pre,
    "idmva85": "gen_gedRecoDeltaR < 0.1 && %s && %s[gen_iGedReco] > %f" % (pre, idConfig.name, wp85.GetVal()),
    "idmva95": "gen_gedRecoDeltaR < 0.1 && %s && %s[gen_iGedReco] > %f" % (pre, idConfig.name, wp95.GetVal()),
}
for name, cut in efficiencies.iteritems():
    if isinstance(idConfig, bdtCommon.EndcapIDConfig):
        cut = cut.replace("gedReco", "localReco")
        cut = cut.replace("GedReco", "LocalReco")
    eff = makeEff(trees[0:2], effpt_def, cut, denom_pt)
    eff.SetNameTitle("effPt_%s" % name, name)
    eff.Write()
    eff2 = makeEff(trees[0:2], effeta_def, cut, denom_eta)
    eff2.SetNameTitle("effEta_%s" % name, name)
    eff2.Write()
    # eff3 = makeEff(trees[0:2], effpt_def, cut, denom_pt+" && !(abs(gen_conversionZ)<290 && gen_conversionRho<120)")
    # eff3.SetNameTitle("effPt_%s_unconverted" % name, name)
    # eff3.Write()


# Bkg efficiencies
reco = "localReco" if isinstance(idConfig, bdtCommon.EndcapIDConfig) else "gedReco"
bkgeffpt_def = ROOT.TH1D("bkgeffpt_def", "%s_pt;p_{T}^{#gamma} [GeV];Efficiency" % reco, len(ptbinning)-1, ptbinning)
bkgeffeta_def = ROOT.TH1D("bkgeffeta_def", "abs(%s_eta);|#eta^{#gamma}|;Efficiency" % reco, 32, 0, 3.2)
denom_pt = "!(%s)" % idConfig.trueDef
denom_eta = "!(%s)" % idConfig.trueDef
if isinstance(idConfig, bdtCommon.EndcapIDConfig):
    denom_pt += " && abs(localReco_scEta)>1.6 && abs(localReco_scEta)<2.9"
    denom_eta += " && localReco_pt>30."
else:
    denom_pt += " && abs(gedReco_scEta)<1.4"
    denom_eta += " && gedReco_pt>30."
efficiencies = {
    "idmva85": "%s && %s > %f" % (pre, idConfig.name, wp85.GetVal()),
    "idmva95": "%s && %s > %f" % (pre, idConfig.name, wp95.GetVal()),
}
for name, cut in efficiencies.iteritems():
    if isinstance(idConfig, bdtCommon.EndcapIDConfig):
        cut = cut.replace("gedReco", "localReco")
        cut = cut.replace("GedReco", "LocalReco")
    eff = makeEff(trees[0:3], bkgeffpt_def, cut, denom_pt)
    eff.SetNameTitle("bkgEffPt_%s" % name, name)
    eff.Write()
    eff2 = makeEff(trees[0:3], bkgeffeta_def, cut, denom_eta)
    eff2.SetNameTitle("bkgEffEta_%s" % name, name)
    eff2.Write()


# Resolution
ebinning = array.array('d', [10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 240., 280., 320., 360., 400., 450., 500., 550., 600., 700., 800., 900., 1000., 1200., 1400., 1600., 1800., 2000., 2500., 3000.])
if isinstance(idConfig, bdtCommon.EndcapIDConfig):
    idcut = idConfig.trueCut + "&& %s > %f" % (idConfig.name, wp95.GetVal())
    etabins = {
        "eta1p7to2p0": "&& abs(localReco_eta)>1.7 && abs(localReco_eta)<2.0",
        "eta2p5to2p8": "&& abs(localReco_eta)>2.5 && abs(localReco_eta)<2.8",
        "eta1p6to2p8": "&& abs(localReco_eta)>1.6 && abs(localReco_eta)<2.8",
    }
    for ebn, ebc in etabins.iteritems():
        cut = idcut + ebc
        for energy in ["scRawEnergy", "seedEnergy", "seedOrigEnergy"]:
            if energy == "seedOrigEnergy" and "noPU" in idConfig.name:
                continue
            resAllE = resolutionPlot(trees[:2], "res_allPhotons_%s_%s" % (energy, ebn), "localReco_%s/(gen_pt[localReco_iGen]*cosh(gen_eta[localReco_iGen])):gen_pt[localReco_iGen]*cosh(gen_eta[localReco_iGen])" % energy, cut, ebinning, idConfig.name)
            resAllE.Write()
            resConvE = resolutionPlot(trees[:2], "res_unconverted_%s_%s" % (energy, ebn), "localReco_%s/(gen_pt[localReco_iGen]*cosh(gen_eta[localReco_iGen])):gen_pt[localReco_iGen]*cosh(gen_eta[localReco_iGen])" % energy, cut+" && !(abs(gen_conversionZ[localReco_iGen])<290 && gen_conversionRho[localReco_iGen]<120)", ebinning, idConfig.name)
            resConvE.Write()
else:
    cut = idConfig.trueCut + " && %s > %f" % (idConfig.name, wp95.GetVal())
    resAll = resolutionPlot(trees[:2], "res_allPhotons",   "gedReco_energy_nmax15/(gen_pt[gedReco_iGen]*cosh(gen_eta[gedReco_iGen])):gen_pt[gedReco_iGen]*cosh(gen_eta[gedReco_iGen])", cut, ebinning, idConfig.name)
    resConv = resolutionPlot(trees[:2], "res_unconverted", "gedReco_energy_nmax15/(gen_pt[gedReco_iGen]*cosh(gen_eta[gedReco_iGen])):gen_pt[gedReco_iGen]*cosh(gen_eta[gedReco_iGen])", cut+" && !(abs(gen_conversionZ[gedReco_iGen])<290 && gen_conversionRho[gedReco_iGen]<120)", ebinning, idConfig.name)
    resAll.Write()
    resConv.Write()


# Reco to x efficiencies
if False:
    recoeffpt_def  = ROOT.TH1D("recoeffpt_def", "gedReco_pt;p_{T}^{#gamma} [GeV];Efficiency", len(ptbinning)-1, ptbinning)
    recoeffeta_def = ROOT.TH1D("recoeffeta_def", "abs(gedReco_eta);|#eta^{#gamma}|;Efficiency", 30, 0, 3.)
    denom = "gedReco_pt > 10. && abs(gedReco_eta)<1.4"
    if isinstance(idConfig, bdtCommon.EndcapIDConfig):
        denom = "localReco_pt > 10. && abs(localReco_eta)>1.5 && abs(localReco_eta)<2.9"
        recoeffpt_def.SetTitle("localReco_pt")
        recoeffeta_def.SetTitle("abs(localReco_eta)")

    denom += " && " + idConfig.trueDef

    efficiencies = {
        "idmva85": "%s > %f" % (idConfig.name, wp85.GetVal()),
        "idmva95": "%s > %f" % (idConfig.name, wp95.GetVal()),
    }
    for name, cut in efficiencies.iteritems():
        eff = makeEff(trees[0:2], recoeffpt_def, cut, denom)
        eff.SetNameTitle("recoEffPt_%s" % name, name)
        eff.Write()
        eff2 = makeEff(trees[0:2], recoeffeta_def, cut, denom)
        eff2.SetNameTitle("recoEffEta_%s" % name, name)
        eff2.Write()

