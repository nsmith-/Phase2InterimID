import ROOT
import array
import warnings
import sys
import bdtCommon

# https://root-forum.cern.ch/t/creating-converter-when-using-ttreeformula/13845/2
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )

idConfig = filter(lambda x: x.name==sys.argv[1], bdtCommon.idconfigs)
if len(idConfig)==1:
    idConfig = idConfig[0]
else:
    raise Exception("Check idconfigs in bdtCommon.py")

idConfig.makeReader()

filenames = bdtCommon.allInputFiles
if "_run2" in idConfig.name:
    filenames = idConfig.inputFiles
if "noPU" in idConfig.name:
    filenames = idConfig.inputFiles
if len(sys.argv) > 2:
    filenames = sys.argv[2:]

for filename in filenames:
    f = ROOT.TFile.Open(filename)
    tree = f.Get("ntupler/photons")
    if not all(var.checkTree(tree) for var in idConfig.varmap):
        print "File %s didn't have all the correct variables for this idconfig" % filename
        continue

    print "Setting up friend tree", idConfig.friendTreeFile(filename)
    fout = ROOT.TFile(idConfig.friendTreeFile(filename), "recreate")
    tout = ROOT.TTree("tree", "not flat")
    mvaval = ROOT.std.vector("float")()
    tout.Branch(idConfig.name, mvaval)
    fm = ROOT.TTreeFormulaManager()
    for var in idConfig.varmap:
        var.setTree(tree, fm)
    fm.Sync()

    nEntries = tree.GetEntries()
    print "Filling, expect %d entries" % nEntries
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        mvaval.clear()
        for i in range(fm.GetNdata()):
            for var in idConfig.varmap:
                var.fill(i)
            mvaval.push_back(idConfig.evalReader())
        tout.Fill()
        if int(tout.GetEntries()*100./nEntries) > int((tout.GetEntries()-1)*100./nEntries):
            print "% 2d percent done" % int((tout.GetEntries()-1)*100./nEntries)

    tout.Write()
    fout.Close()

