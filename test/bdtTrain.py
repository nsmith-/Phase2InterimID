import ROOT
import sys
from bdtCommon import BarrelIDConfig, EndcapIDConfig

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

idConfig = EndcapIDConfig
if sys.argv[-1] == 'b':
    idConfig = BarrelIDConfig

truecut = "isTrue"
bkgcut = "!isTrue"
sigCut = ROOT.TCut(truecut)
bkgCut = ROOT.TCut(bkgcut)

fIn = ROOT.TFile.Open("%sin.root" % idConfig.name)
tIn = fIn.Get("tree")
fOut = ROOT.TFile("%sout.root" % idConfig.name, "recreate")
# Transformations=I;D;P;G,D:
w = ROOT.TMVA.Factory(idConfig.name, fOut, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification")
d = ROOT.TMVA.DataLoader()
for var in idConfig.varmap:
    if var.train:
        d.AddVariable(var.name)
d.SetInputTrees(tIn, sigCut, bkgCut)
d.SetWeightExpression("weight")

#w.BookMethod(d, ROOT.TMVA.Types.kFisher, "Fisher", "H:!V")
#w.BookMethod(d, ROOT.TMVA.Types.kCuts, "Cuts", "FitMethod=SA")
w.BookMethod(d, ROOT.TMVA.Types.kBDT, "BDT", "NTrees=800:MaxDepth=3:nCuts=100")
w.TrainAllMethods()
w.TestAllMethods()
w.EvaluateAllMethods()

fOut.Close()
ROOT.TMVA.mvas("default", "%sout.root" % idConfig.name, ROOT.TMVA.kCompareType)
ROOT.gPad.Print("bdt_overtrain.pdf")

#g = ROOT.TMVA.TMVAGui("tmvaout.root")
