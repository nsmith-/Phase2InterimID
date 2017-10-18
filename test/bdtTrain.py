import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import array
import os

vars = [
    "sigmaUU",
    "sigmaVV",
    "nLayers",
    "firstLayer",
    "FHoverE",
    "depthCompatibility",
    "isoRing0",
    "isoRing1",
    "isoRing2",
]
truecut = "isTrue"
bkgcut = "!isTrue"
sigCut = ROOT.TCut(truecut)
bkgCut = ROOT.TCut(bkgcut)

fIn = ROOT.TFile.Open("tmvain.root")
tIn = fIn.Get("tree")
fOut = ROOT.TFile("tmvaout.root", "recreate")
# Transformations=I;D;P;G,D:
w = ROOT.TMVA.Factory("phase2photons", fOut, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification")
d = ROOT.TMVA.DataLoader()
for var in vars:
    d.AddVariable(var)
d.SetInputTrees(tIn, sigCut, bkgCut)
d.SetWeightExpression("weight")

#w.BookMethod(d, ROOT.TMVA.Types.kFisher, "Fisher", "H:!V")
#w.BookMethod(d, ROOT.TMVA.Types.kCuts, "Cuts", "FitMethod=SA")
w.BookMethod(d, ROOT.TMVA.Types.kBDT, "BDT", "NTrees=800:MaxDepth=3:nCuts=100")
w.TrainAllMethods()
w.TestAllMethods()
w.EvaluateAllMethods()

fOut.Close()
ROOT.TMVA.mvas("default", "tmvaout.root", ROOT.TMVA.kCompareType)
ROOT.gPad.Print("bdt_overtrain.pdf")

#g = ROOT.TMVA.TMVAGui("tmvaout.root")
