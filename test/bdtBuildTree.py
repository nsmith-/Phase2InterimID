import ROOT
import os
import array
import warnings
from bdtCommon import BarrelIDConfig, EndcapIDConfig

# https://root-forum.cern.ch/t/creating-converter-when-using-ttreeformula/13845/2
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )

idConfig = BarrelIDConfig
# idConfig = EndcapIDConfig

filenames = [
    "/eos/cms/store/user/ncsmith/932phoID/RelValH125GGgluonfusion_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/171017_222812/0000/output_1.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValH125GGgluonfusion_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/171017_222812/0000/output_2.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValH125GGgluonfusion_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/171017_222812/0000/output_3.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValH125GGgluonfusion_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/171017_222812/0000/output_4.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/171017_222826/0000/output_1.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/171017_222826/0000/output_2.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/171017_222826/0000/output_3.root",
    "/eos/cms/store/user/ncsmith/932phoID/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/171017_222826/0000/output_4.root",
]


files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
if not all(var.checkTree(tree) for tree in trees for var in idConfig.varmap):
    exit()


ROOT.gROOT.cd()
hreweight_tmp = idConfig.hreweight_def.Clone("hreweight_tmp")
for tree in trees:
    ROOT.gROOT.cd()
    tree.Draw("%s:%s>>+hreweight_tmp" % (idConfig.reweightvar2.formula, idConfig.reweightvar1.formula), idConfig.trueCut, "goff")
ROOT.gROOT.cd()
hreweight = hreweight_tmp.Clone("hreweight")
hreweight_tmp.Reset()
for tree in trees:
    ROOT.gROOT.cd()
    tree.Draw("%s:%s>>hreweight_tmp" % (idConfig.reweightvar2.formula, idConfig.reweightvar1.formula), idConfig.trainingCut, "goff")
hreweight.Divide(hreweight_tmp)

fout = ROOT.TFile("%sin.root" % idConfig.name, "recreate")
hreweight.Write()
tout = ROOT.TTree("tree", "flat like a pancake")
for var in idConfig.varmap:
    if var is idConfig.trainvar:
        continue
    var.makeBranch(tout)
weight = array.array('f', [0.])
tout.Branch("weight", weight, "weight/f")

for tree in trees:
    fm = ROOT.TTreeFormulaManager()
    for var in idConfig.varmap:
        var.setTree(tree, fm)
    fm.Sync()

    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        for i in range(fm.GetNdata()):
            for var in idConfig.varmap:
                var.fill(i)
            if not idConfig.trainvar.val(): 
                continue
            if idConfig.truevar.val():
                weight[0] = 1.
            else:
                weight[0] = hreweight.GetBinContent(hreweight.FindBin(idConfig.reweightvar1.val(), idConfig.reweightvar2.val()))
                if weight[0] == 0.:
                    weight[0] = 1.
            tout.Fill()
            if tout.GetEntries() % 1000 == 0:
                print "Done with % 7d entries" % tout.GetEntries()


tout.Write()
fout.Close()

