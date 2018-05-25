import ROOT
import array
import glob
import math
import re
import sys
import bdtCommon

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True


def makeEff(trees, histDef, ncut, dcut):
    dcute = dcut.replace("gedReco", "localReco")
    ncutb = ncut[0]
    ncute = ncut[1]
    htmpnum = histDef.Clone("tmpnum")
    htmpnum.Sumw2()
    htmpdenom = histDef.Clone("tmpdenom")
    eta = "seedEta"
    for tree in trees:
        tree.Draw("gedReco_energy_nmax15/cosh(gedReco_{0})      :abs(gedReco_{0})   >>+tmpnum".format(eta), "(%s)&&(%s)" % (dcut, ncutb), "goff")
        tree.Draw("localReco_seedOrigEnergy/cosh(localReco_{0}) :abs(localReco_{0}) >>+tmpnum".format(eta), "(%s)&&(%s)" % (dcute, ncute), "goff")
        tree.Draw("gedReco_energy_nmax15/cosh(gedReco_{0})      :abs(gedReco_{0})   >>+tmpdenom".format(eta), "%s" % dcut, "goff")
        tree.Draw("localReco_seedOrigEnergy/cosh(localReco_{0}) :abs(localReco_{0}) >>+tmpdenom".format(eta), "%s" % dcute, "goff")
    heff = histDef.Clone("tmpeff")
    heff.Divide(htmpnum, htmpdenom, 1, 1, "b")
    return heff


idBarrel = filter(lambda x: x.name == "barrelV4", bdtCommon.idconfigs)[0]
idEndcap = filter(lambda x: x.name == "endcapV4", bdtCommon.idconfigs)[0]

filenames = idBarrel.inputFiles
files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
for i, tree in enumerate(trees):
    tree.AddFriend("tree", idBarrel.friendTreeFile(filenames[i]))
    tree.AddFriend("tree", idEndcap.friendTreeFile(filenames[i]))

categories = {
    "promptPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen]",
    "promptUnconvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) != 11 && !(abs(gen_conversionZ[gedReco_iGen])<290 && gen_conversionRho[gedReco_iGen]<120)",
    "promptConvPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen] && abs(gen_parentId[gedReco_iGen]) != 11 && (abs(gen_conversionZ[gedReco_iGen])<290 && gen_conversionRho[gedReco_iGen]<120)",
    "nonpromptPho": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && !gen_isPromptFinalState[gedReco_iGen]",
    "phoPiZero": "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_parentId[gedReco_iGen]==111",
    "ele": "gedReco_iGen>=0 && abs(gen_id[gedReco_iGen]) == 11",
    "fake": "gedReco_iGen<0",
}


ptbins = array.array('d', [0, 10, 20, 30, 40, 60, 90, 120, 160, 240, 300, 500, 1000])
etabins = array.array('d', [0, 0.8, 1.4, 1.48, 1.6, 2.0, 2.5, 2.8, 2.9, 3.0])
fout = ROOT.TFile.Open("delphes.root", "recreate")
ptEtaDef = ROOT.TH2D("ptEtaDef", "Pt-Eta efficiency;Photon |#eta|;Photon p_{T} [GeV];Efficiency", len(etabins)-1, etabins, len(ptbins)-1, ptbins)

effs = []


genCut = "gen_id == 22 && gen_isPromptFinalState"
recoCut = "(gen_gedRecoDeltaR < 0.1 && abs(gedReco_energy_nmax15[gen_iGedReco]/(gen_pt*cosh(gen_eta))-1) < 0.5) || (gen_localRecoDeltaR < 0.1 && abs(localReco_seedEnergy[gen_iLocalReco]/(gen_pt*cosh(gen_eta))-1) < 0.5)"
htmpnum = ptEtaDef.Clone("tmpnumreco")
htmpnum.Sumw2()
htmpdenom = ptEtaDef.Clone("tmpdenomreco")
for tree in trees:
    tree.Draw("gen_pt:abs(gen_eta)>>+tmpnumreco", "(%s)&&(%s)" % (genCut, recoCut), "goff")
    tree.Draw("gen_pt:abs(gen_eta)>>+tmpdenomreco", "%s" % genCut, "goff")
geneff = ptEtaDef.Clone("eff_reconstruction")
geneff.Divide(htmpnum, htmpdenom, 1, 1, "b")
effs.append(geneff)


looseIdBarrel = "(%s>=%f)" % (idBarrel.name, 0.)
tightIdBarrel = "(%s>=%f)" % (idBarrel.name,  0.56)
looseIdEndcap = "(%s>=%f)" % (idEndcap.name, 0.2)
tightIdEndcap = "(%s>=%f)" % (idEndcap.name, 0.68)

for cat, cut in categories.iteritems():
    eff = makeEff(trees, ptEtaDef, (looseIdBarrel, looseIdEndcap), cut)
    eff.SetName("eff_%s_looseID" % cat)
    effs.append(eff)

    eff = makeEff(trees, ptEtaDef, (tightIdBarrel, tightIdEndcap), cut)
    eff.SetName("eff_%s_tightID" % cat)
    effs.append(eff)


for eff in effs:
    eff.SetDirectory(fout)
    eff.Write()
