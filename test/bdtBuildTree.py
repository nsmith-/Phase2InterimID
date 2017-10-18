import ROOT
import array
import os
import warnings
# https://root-forum.cern.ch/t/creating-converter-when-using-ttreeformula/13845/2
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )


class Variable:
    def __init__(self, name, formula):
        self.name = name
        self.formula = formula
        self._mem = array.array('f', [0.])

    def makeBranch(self, tree):
        tree.Branch(self.name, self._mem, self.name+"/f")

    def checkTree(self, tree):
        fclass = ROOT.TTreeFormula(self.name, self.formula, tree)
        if fclass.GetNdim() == 1:
            return True
        print "Failed to make variable %s for %r" % (self.name, tree)
        return False

    def setTree(self, tree, fm):
        self._fclass = ROOT.TTreeFormula(self.name, self.formula, tree)
        fm.Add(self._fclass)

    def fill(self, instance):
        self._mem[0] = self._fclass.EvalInstance(instance)

    def val(self):
        return self._mem[0]


trainingCut = "localReco_pt>25 && abs(localReco_eta)>1.5"
trueCut = trainingCut + " && localReco_iGen>=0 && gen_parentId[localReco_iGen] == 25"
varmap = [
    Variable("pt", "localReco_pt"),
    Variable("abseta", "abs(localReco_eta)"),
    Variable("isTrained", trainingCut),
    Variable("isTrue", trueCut),
    Variable("sigmaUU", "localReco_sigmaUU"),
    Variable("sigmaVV", "localReco_sigmaVV"),
    Variable("nLayers", "localReco_nLayers"),
    Variable("firstLayer", "localReco_firstLayer"),
    Variable("FHoverE", "localReco_energyFH/localReco_energyEE"),
    Variable("depthCompatibility", "localReco_depthCompatibility"),
    Variable("isoRing0", "localReco_caloIsoRing0"),
    Variable("isoRing1", "localReco_caloIsoRing1"),
    Variable("isoRing2", "localReco_caloIsoRing2"),
]
reweightvar1, reweightvar2, trainvar, truevar = varmap[0:4]
hreweight_tmp = ROOT.TH2D("hreweight_tmp", "Bkg to signal reweight;Photon p_{T} (GeV);Photon |#eta|", 15, 25, 100, 2, 1.5, 3.)

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
outputFilename = "tmvain.root"


files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
if not all(var.checkTree(tree) for tree in trees for var in varmap):
    exit()

for tree in trees:
    ROOT.gROOT.cd()
    tree.Draw("%s:%s>>+hreweight_tmp" % (reweightvar2.formula, reweightvar1.formula), trueCut, "goff")
ROOT.gROOT.cd()
hreweight = hreweight_tmp.Clone("hreweight")
hreweight_tmp.Reset()
for tree in trees:
    ROOT.gROOT.cd()
    tree.Draw("%s:%s>>hreweight_tmp" % (reweightvar2.formula, reweightvar1.formula), trainingCut, "goff")
hreweight.Divide(hreweight_tmp)

fout = ROOT.TFile(outputFilename, "recreate")
hreweight.Write()
tout = ROOT.TTree("tree", "flat like a pancake")
for var in varmap:
    if var is trainvar:
        continue
    var.makeBranch(tout)
weight = array.array('f', [0.])
tout.Branch("weight", weight, "weight/f")

for tree in trees:
    fm = ROOT.TTreeFormulaManager()
    for var in varmap:
        var.setTree(tree, fm)
    fm.Sync()

    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        for i in range(fm.GetNdata()):
            for var in varmap:
                var.fill(i)
            if not trainvar.val(): 
                continue
            if truevar.val():
                weight[0] = 1.
            else:
                weight[0] = hreweight.GetBinContent(hreweight.FindBin(reweightvar1.val(), reweightvar2.val()))
                if weight[0] == 0.:
                    weight[0] = 1.
            tout.Fill()
            if tout.GetEntries() % 1000 == 0:
                print "Done with % 7d entries" % tout.GetEntries()


tout.Write()
fout.Close()

