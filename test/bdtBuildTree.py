import ROOT
import array
import warnings
import sys
import bdtCommon

# https://root-forum.cern.ch/t/creating-converter-when-using-ttreeformula/13845/2
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )

idConfig = filter(lambda x: x.name==sys.argv[-1], bdtCommon.idconfigs)
if len(idConfig)==1:
    idConfig = idConfig[0]
else:
    raise Exception("Check idconfigs in bdtCommon.py")

filenames = idConfig.inputFiles

print "Checking files"
files = map(ROOT.TFile.Open, filenames)
trees = [f.Get("ntupler/photons") for f in files]
if not all(var.checkTree(tree) for tree in trees for var in idConfig.varmap):
    exit()


print "Making reweight numerator hist"
ROOT.gROOT.cd()
hreweight_tmp = idConfig.hreweight_def.Clone("hreweight_tmp")
for tree in trees:
    ROOT.gROOT.cd()
    tree.Draw("%s:%s>>+hreweight_tmp" % (idConfig.reweightvar2.formula, idConfig.reweightvar1.formula), idConfig.trueCut, "goff")
ROOT.gROOT.cd()
nEntries = hreweight_tmp.GetEntries()
hreweight_num = hreweight_tmp.Clone("hreweight_num")
hreweight_num.Scale(1./hreweight_num.Integral())
print "Total signal events:", nEntries
reweight_hists = []
for i, tree in enumerate(trees):
    print "Making reweight denominator for input tree", filenames[i]
    hreweight = hreweight_num.Clone("hreweight_fin%d" % i)
    hreweight_tmp.Reset()
    ROOT.gROOT.cd()
    tree.Draw("%s:%s>>hreweight_tmp" % (idConfig.reweightvar2.formula, idConfig.reweightvar1.formula), idConfig.bkgCut, "goff")
    print "File %s has %d bkg events" % (filenames[i], hreweight_tmp.GetEntries())
    nEntries += hreweight_tmp.GetEntries()
    # normalize so bkg from bkg-free files doesn't dominate
    hreweight_tmp.Scale(1./hreweight_tmp.Integral())
    hreweight.Divide(hreweight_tmp)
    if 'DY' in filenames[i]:
        hreweight.Scale(0.1)
    reweight_hists.append(hreweight)

print "Setting up output tree"
fout = ROOT.TFile("%sin.root" % idConfig.name, "recreate")
tout = ROOT.TTree("tree", "flat like a pancake")
for var in idConfig.varmap:
    if var is idConfig.trainvar:
        continue
    var.makeBranch(tout)
weight = array.array('f', [0.])
tout.Branch("weight", weight, "weight/f")

print "Filling, expect %d entries" % nEntries
for i, tree in enumerate(trees):
    hreweight = reweight_hists[i]
    fout.cd()
    hreweight.Write()

    fm = ROOT.TTreeFormulaManager()
    for var in idConfig.varmap:
        var.setTree(tree, fm)
    fm.Sync()

    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        for i in range(fm.GetNdata()):
            idConfig.trainvar.fill(i)
            if not idConfig.trainvar.val(): 
                continue
            for var in idConfig.varmap:
                var.fill(i)
            if idConfig.truevar.val():
                weight[0] = 1.
            else:
                weight[0] = hreweight.GetBinContent(hreweight.FindBin(idConfig.reweightvar1.val(), idConfig.reweightvar2.val()))
                if weight[0] == 0.:
                    weight[0] = 1.
            tout.Fill()
            if int(tout.GetEntries()*100./nEntries) > int((tout.GetEntries()-1)*100./nEntries):
                print "% 2d percent done" % int((tout.GetEntries()-1)*100./nEntries)


tout.Write()
#fout.Close()

