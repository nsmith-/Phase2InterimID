import ROOT
import sys
import bdtCommon

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

idConfig = filter(lambda x: x.name==sys.argv[-1], bdtCommon.idconfigs)
if len(idConfig)==1:
    idConfig = idConfig[0]
else:
    raise Exception("Check idconfigs in bdtCommon.py")

truecut = "isTrue"
bkgcut = "!isTrue"
sigCut = ROOT.TCut(truecut)
bkgCut = ROOT.TCut(bkgcut)

fIn = ROOT.TFile.Open("%sin.root" % idConfig.name)
tIn = fIn.Get("tree")
fOut = ROOT.TFile("%sout.root" % idConfig.name, "recreate")
hists = [k.ReadObj() for k in fIn.GetListOfKeys() if k.GetClassName == "TH2D"]
for h in hists:
    h.Write()
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
w.BookMethod(d, ROOT.TMVA.Types.kBDT, "BDT", idConfig.bdtSettings)
w.TrainAllMethods()
w.TestAllMethods()
w.EvaluateAllMethods()

fOut.Close()
ROOT.TMVA.mvas("default", "%sout.root" % idConfig.name, ROOT.TMVA.kCompareType)
ROOT.gPad.Print("bdt_overtrain_%s.pdf" % idConfig.name)

#g = ROOT.TMVA.TMVAGui("tmvaout.root")
