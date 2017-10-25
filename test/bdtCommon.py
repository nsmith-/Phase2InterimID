import ROOT
import array

class Variable:
    def __init__(self, name, formula, train=True):
        self.name = name
        self.formula = formula
        self.train = train
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


class IDConfig:
    def __init__(self, name, morevars, bdtSettings=""):
        self.name = name
        self.varmap = self._def_varmap + morevars
        self.bdtSettings = bdtSettings

    def makeReader(self):
        if not hasattr(self, "_reader"):
            self._reader = ROOT.TMVA.Reader("!Color:Silent")
            #self._reader.SetName(self.name+"Reader")
            # https://root.cern.ch/doc/v606/Reader_8h_source.html#l00126 (wtf?)
            iReader = ROOT.gROOT.GetListOfSpecials().GetSize()
            ROOT.gROOT.GetListOfSpecials().Add(self._reader)
            params = [v for v in self.varmap if v.train]
            for i, var in enumerate(params):
                self._reader.AddVariable(var.name, var._mem)
            self._reader.BookMVA(self.name, "default/weights/%s_BDT.weights.xml" % self.name)
            fcnArgs = ", ".join(("const double var%d" % i) for i in range(len(params)))
            fcnBodyVars = " ".join("in[{0}]=var{0};".format(i) for i in range(len(params)))
            evalFcn = 'double eval%s(%s){auto ptr = (TMVA::Reader*) gROOT->GetListOfSpecials()->At(%d); std::vector<double> in(%d); %s return ptr->EvaluateMVA(in, "%s", 0.);}' % (self.name, fcnArgs, iReader, len(params), fcnBodyVars, self.name)
            ROOT.gInterpreter.Declare(evalFcn)
            self.evalFcn = 'eval%s(%s)' % (self.name, ", ".join(p.formula for p in params))
            self.mvaVariable = Variable(self.name, self.evalFcn)
        return self._reader

    def friendTreeFile(self, filename):
        return "%s_%s.root" % (filename.replace(".root", ""), self.name)


class BarrelIDConfig(IDConfig):
    trainingCut = "gedReco_pt>25 && abs(gedReco_eta)<1.4"
    trueDef = "gedReco_iGen>=0 && gen_id[gedReco_iGen] == 22 && gen_isPromptFinalState[gedReco_iGen]"
    trueCut = trainingCut + " && (%s)" % trueDef
    bkgCut = trainingCut + " && !(%s)" % trueDef
    _def_varmap = [
        Variable("pt", "gedReco_pt", train=False),
        Variable("abseta", "abs(gedReco_eta)", train=False),
        Variable("isTrained", trainingCut, train=False),
        Variable("isTrue", trueCut, train=False),
    ]
    reweightvar1, reweightvar2, trainvar, truevar = _def_varmap[0:4]
    hreweight_def = ROOT.TH2D("hreweight_def_barrel", "Bkg to signal reweight;Photon p_{T} (GeV);Photon |#eta|", 25, 25, 150, 7, 0., 1.4)
    

class EndcapIDConfig(IDConfig):
    trainingCut = "localReco_pt>25 && abs(localReco_eta)>1.5 && localReco_nLayers > 0 && abs(localReco_depthCompatibility) < 10. && localReco_seedEnergyFH/localReco_seedEnergyEE < 10."
    trueDef = "localReco_iGen>=0 && gen_id[localReco_iGen] == 22 && gen_isPromptFinalState[localReco_iGen]"
    trueCut = trainingCut + " && (%s)" % trueDef
    bkgCut = trainingCut + " && !(%s)" % trueDef
    _def_varmap = [
        Variable("pt", "localReco_pt", train=False),
        Variable("abseta", "abs(localReco_eta)", train=False),
        Variable("isTrained", trainingCut, train=False),
        Variable("isTrue", trueCut, train=False),
    ]
    reweightvar1, reweightvar2, trainvar, truevar = _def_varmap[0:4]
    hreweight_def = ROOT.TH2D("hreweight_def_endcap", "Bkg to signal reweight;Photon p_{T} (GeV);Photon |#eta|", 25, 25, 150, 5, 1.5, 3.)


idconfigs = [
    BarrelIDConfig("barrel1", [
        Variable("full5x5_sigmaIetaIeta", "gedReco_full5x5_sigmaIetaIeta"),
        Variable("full5x5_sigmaIetaIphi", "gedReco_full5x5_sigmaIetaIphi"),
        Variable("full5x5_sigmaIphiIphi", "gedReco_full5x5_sigmaIphiIphi"),
        Variable("hadTowOverEm", "gedReco_hadTowOverEm"),
        Variable("chargedHadronIso", "gedReco_chargedHadronIso"),
        Variable("neutralHadronIso", "gedReco_neutralHadronIso"),
        Variable("photonIso", "gedReco_photonIso"),
        Variable("hasPixelSeed", "gedReco_hasPixelSeed"),
    ], "NTrees=800:MaxDepth=3:nCuts=100"),
    EndcapIDConfig("endcap1", [
        Variable("sigmaUU", "localReco_sigmaUU"),
        Variable("sigmaVV", "localReco_sigmaVV"),
        Variable("e4oEtot", "localReco_e4oEtot"),
        Variable("layerEfrac10", "localReco_layerEfrac10"),
        Variable("layerEfrac90", "localReco_layerEfrac90"),
        Variable("FHoverE", "localReco_seedEnergyFH/localReco_seedEnergyEE"),
        Variable("depthCompatibility", "localReco_depthCompatibility"),
        Variable("isoRing0", "localReco_caloIsoRing0"),
        Variable("isoRing1", "localReco_caloIsoRing1"),
        Variable("isoRing2", "localReco_caloIsoRing2"),
        Variable("isoRing3", "localReco_caloIsoRing3"),
        Variable("isoRing4", "localReco_caloIsoRing4"),
    ], "NTrees=800:MaxDepth=3:nCuts=100"),
    BarrelIDConfig("barrel2", [
        Variable("full5x5_sigmaIetaIeta", "gedReco_full5x5_sigmaIetaIeta"),
        Variable("full5x5_sigmaIetaIphi", "gedReco_full5x5_sigmaIetaIphi"),
        Variable("full5x5_sigmaIphiIphi", "gedReco_full5x5_sigmaIphiIphi"),
        Variable("etaWidth", "gedReco_etaWidth"),
        Variable("phiWidth", "gedReco_phiWidth"),
        Variable("full5x5_r9", "gedReco_full5x5_r9"),
        Variable("full5x5_s4", "gedReco_full5x5_s4"),
        Variable("hadTowOverEm", "gedReco_hadTowOverEm"),
        Variable("chargedHadronIso", "gedReco_chargedHadronIso"),
        Variable("neutralHadronIso", "gedReco_neutralHadronIso"),
        Variable("photonIso", "gedReco_photonIso"),
        Variable("hasPixelSeed", "gedReco_hasPixelSeed"),
        Variable("scRawEnergy", "gedReco_scRawEnergy"),
        Variable("trkSumPt", "gedReco_trkSumPtSolidConeDR04"),
        Variable("eVeto", "gedReco_conversionSafeElectronVeto"),
        Variable("rho", "rho"),
    ], "NTrees=2000:MaxDepth=3:nCuts=100"),
    EndcapIDConfig("endcap2", [
        Variable("sigmaUU", "localReco_sigmaUU"),
        Variable("sigmaVV", "localReco_sigmaVV"),
        Variable("e4oEtot", "localReco_e4oEtot"),
        Variable("layerEfrac10", "localReco_layerEfrac10"),
        Variable("layerEfrac90", "localReco_layerEfrac90"),
        Variable("FHoverE", "localReco_seedEnergyFH/localReco_seedEnergyEE"),
        Variable("depthCompatibility", "localReco_depthCompatibility"),
        Variable("isoRing0", "localReco_caloIsoRing0"),
        Variable("isoRing1", "localReco_caloIsoRing1"),
        Variable("isoRing2", "localReco_caloIsoRing2"),
        Variable("isoRing3", "localReco_caloIsoRing3"),
        Variable("isoRing4", "localReco_caloIsoRing4"),
        Variable("scEnergy", "localReco_scEnergy"),
        Variable("matchedTrackChi2", "localReco_matchedGsfChi2"),
        Variable("matchedTrackLostHits", "localReco_matchedGsfLostHits"),
        Variable("rho", "rho"),
    ], "NTrees=2000:MaxDepth=3:nCuts=100"),
]
