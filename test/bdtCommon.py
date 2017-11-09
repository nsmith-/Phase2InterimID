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
    def __init__(self, name, morevars, bdtSettings, inputFiles):
        self.name = name
        self.varmap = self._def_varmap + morevars
        self.bdtSettings = bdtSettings
        self.inputFiles = inputFiles

    def makeReader(self):
        if not hasattr(self, "_reader"):
            self._reader = ROOT.TMVA.Reader("!Color:Silent")
            params = [v for v in self.varmap if v.train]
            for i, var in enumerate(params):
                self._reader.AddVariable(var.name, var._mem)
            self._reader.BookMVA(self.name, "default/weights/%s_BDT.weights.xml" % self.name)
        return self._reader

    def evalReader(self):
        return self._reader.EvaluateMVA(self.name)

    def friendTreeFile(self, filename):
        return "%s_%s.root" % (filename.replace(".root", ""), self.name)

    def readerCppCode(self):
        params = [v for v in self.varmap if v.train]
        maxlen = max(len(p.name) for p in params)
        for i, v in enumerate(params):
            print '  idReader_->AddVariable("%s",%s &idVars_[%d]);' % (v.name, " "*(maxlen-len(v.name)), i)
        print
        for i, v in enumerate(params):
            print '  barrelVars_[%d] = 0.; // %s' % (i, v.formula)


class BarrelIDConfig(IDConfig):
    preselection = "1."
    trainingCut = preselection + " && gedReco_pt>25 && abs(gedReco_eta)<1.4 && (gedReco_iGen<0||abs(gen_parentId[gedReco_iGen])!=11)"
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
    hreweight_def = ROOT.TH2D("hreweight_def_barrel", "Bkg to signal reweight;Photon p_{T} (GeV);Photon |#eta|", 36, 20, 200, 14, 0., 1.4)
    

class EndcapIDConfig(IDConfig):
    preselection = "abs(localReco_depthCompatibility) < 20. && localReco_seedEnergyFH/localReco_seedEnergyEE < 20."
    trainingCut = preselection + " && localReco_pt>25 && abs(localReco_eta)>1.5 && (localReco_iGen<0||abs(gen_parentId[localReco_iGen])!=11)"
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
    hreweight_def = ROOT.TH2D("hreweight_def_endcap", "Bkg to signal reweight;Photon p_{T} (GeV);Photon |#eta|", 36, 20, 200, 15, 1.5, 3.)


class EndcapIDConfigRun2(EndcapIDConfig):
    preselection = "1."
    trainingCut = preselection + " && gedReco_pt>25 && abs(gedReco_eta)>1.5 && (gedReco_iGen<0||abs(gen_parentId[gedReco_iGen])!=11)"
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
    hreweight_def = ROOT.TH2D("hreweight_def_endcap_run2", "Bkg to signal reweight;Photon p_{T} (GeV);Photon |#eta|", 25, 25, 150, 5, 1.5, 2.5)


allInputFiles = [
    "/data/ncsmith/932phoID_round4/GJets.root",
    "/data/ncsmith/932phoID_round4/DiPhotonSherpa.root",
    "/data/ncsmith/932phoID_round4/QCD.root",
    "/data/ncsmith/932phoID_round4/DY2J.root",
]
allInputFilesPU0 = [
    "/data/ncsmith/932phoID_round3/GJetspu0.root",
    "/data/ncsmith/932phoID_round3/DiPhotonSherpapu0.root",
    "/data/ncsmith/932phoID_round3/QCDpu0.root",
    "/data/ncsmith/932phoID_round3/DY2Jpu0.root",
]

idconfigs = [
    EndcapIDConfig("endcapV4",
        [
            Variable("sigmaUU", "localReco_sigmaUU"),
            Variable("sigmaVV", "localReco_sigmaVV"),
            Variable("e4oEtot", "localReco_e4oEtot"),
            Variable("layerEfrac10", "localReco_layerEfrac10"),
            Variable("layerEfrac90", "localReco_layerEfrac90"),
            Variable("FHoverE", "localReco_seedEnergyFH/localReco_seedEnergyEE"),
            Variable("measuredDepth", "localReco_measuredDepth"),
            Variable("depthCompatibility", "localReco_depthCompatibility"),
            Variable("isoRing0", "localReco_caloIsoRing0"),
            Variable("isoRing1", "localReco_caloIsoRing1"),
            Variable("isoRing2", "localReco_caloIsoRing2"),
            Variable("isoRing3", "localReco_caloIsoRing3"),
            Variable("isoRing4", "localReco_caloIsoRing4"),
            Variable("scEnergy", "localReco_scRawEnergy"),
            Variable("matchedTrackChi2", "localReco_matchedGsfChi2"),
            Variable("matchedTrackHits", "localReco_matchedGsfHits"),
            Variable("matchedTrackLostHits", "localReco_matchedGsfLostHits"),
            Variable("rho", "rho"),
        ],
        "BoostType=Grad:Shrinkage=0.4:UseBaggedBoost=True:BaggedSampleFraction=0.6:NTrees=2000:MaxDepth=2:nCuts=200",
        allInputFiles
    ),
    BarrelIDConfig("barrelV4",
        [
            Variable("full5x5_sigmaIetaIeta", "gedReco_full5x5_sigmaIetaIeta"),
            Variable("full5x5_sigmaIetaIphi", "gedReco_full5x5_sigmaIetaIphi"),
            Variable("full5x5_sigmaIphiIphi", "gedReco_full5x5_sigmaIphiIphi"),
            Variable("etaWidth", "gedReco_etaWidth"),
            Variable("phiWidth", "gedReco_phiWidth"),
            Variable("full5x5_r9", "gedReco_full5x5_r9"),
            Variable("full5x5_s4", "gedReco_full5x5_s4"),
            Variable("hadronicOverEm", "gedReco_hadronicOverEm"),
            Variable("chargedHadronIso", "gedReco_chargedHadronIso"),
            Variable("neutralHadronIso", "gedReco_neutralHadronIso"),
            Variable("photonIso", "gedReco_photonIso"),
            Variable("hasPixelSeed", "gedReco_hasPixelSeed"),
            Variable("scRawEnergy", "gedReco_scRawEnergy"),
            Variable("scEta", "gedReco_scEta"),
            Variable("trkSumPt", "gedReco_trkSumPtSolidConeDR04"),
            Variable("eVeto", "gedReco_conversionSafeElectronVeto"),
            Variable("rho", "rho"),
        ],
        "BoostType=Grad:Shrinkage=0.4:UseBaggedBoost=True:BaggedSampleFraction=0.6:NTrees=2000:MaxDepth=2:nCuts=200",
        allInputFiles
    ),
    EndcapIDConfig("endcapV3noPU",
        [
            Variable("sigmaUU", "localReco_sigmaUU"),
            Variable("sigmaVV", "localReco_sigmaVV"),
            Variable("e4oEtot", "localReco_e4oEtot"),
            Variable("layerEfrac10", "localReco_layerEfrac10"),
            Variable("layerEfrac90", "localReco_layerEfrac90"),
            Variable("FHoverE", "localReco_seedEnergyFH/localReco_seedEnergyEE"),
            Variable("measuredDepth", "localReco_measuredDepth"),
            Variable("depthCompatibility", "localReco_depthCompatibility"),
            Variable("isoRing0", "localReco_caloIsoRing0"),
            Variable("isoRing1", "localReco_caloIsoRing1"),
            Variable("isoRing2", "localReco_caloIsoRing2"),
            Variable("isoRing3", "localReco_caloIsoRing3"),
            Variable("isoRing4", "localReco_caloIsoRing4"),
            Variable("scEnergy", "localReco_scRawEnergy"),
            Variable("matchedTrackChi2", "localReco_matchedGsfChi2"),
            Variable("matchedTrackHits", "localReco_matchedGsfHits"),
            Variable("matchedTrackLostHits", "localReco_matchedGsfLostHits"),
        ],
        "BoostType=Grad:Shrinkage=0.4:UseBaggedBoost=True:BaggedSampleFraction=0.6:NTrees=3000:MaxDepth=2:nCuts=200",
        allInputFilesPU0
    ),
    BarrelIDConfig("barrelV3noPU",
        [
            Variable("full5x5_sigmaIetaIeta", "gedReco_full5x5_sigmaIetaIeta"),
            Variable("full5x5_sigmaIetaIphi", "gedReco_full5x5_sigmaIetaIphi"),
            Variable("full5x5_sigmaIphiIphi", "gedReco_full5x5_sigmaIphiIphi"),
            Variable("etaWidth", "gedReco_etaWidth"),
            Variable("phiWidth", "gedReco_phiWidth"),
            Variable("full5x5_r9", "gedReco_full5x5_r9"),
            Variable("full5x5_s4", "gedReco_full5x5_s4"),
            Variable("hadronicOverEm", "gedReco_hadronicOverEm"),
            Variable("chargedHadronIso", "gedReco_chargedHadronIso"),
            Variable("neutralHadronIso", "gedReco_neutralHadronIso"),
            Variable("photonIso", "gedReco_photonIso"),
            Variable("hasPixelSeed", "gedReco_hasPixelSeed"),
            Variable("scRawEnergy", "gedReco_scRawEnergy"),
            Variable("scEta", "gedReco_scEta"),
            Variable("trkSumPt", "gedReco_trkSumPtSolidConeDR04"),
            Variable("eVeto", "gedReco_conversionSafeElectronVeto"),
        ],
        "BoostType=Grad:Shrinkage=0.4:UseBaggedBoost=True:BaggedSampleFraction=0.6:NTrees=3000:MaxDepth=2:nCuts=200",
        allInputFilesPU0
    ),
    BarrelIDConfig("barrelV5",
        [
            Variable("full5x5_sigmaIetaIeta", "gedReco_full5x5_sigmaIetaIeta"),
            Variable("full5x5_sigmaIetaIphi", "gedReco_full5x5_sigmaIetaIphi"),
            Variable("full5x5_sigmaIphiIphi", "gedReco_full5x5_sigmaIphiIphi"),
            Variable("etaWidth", "gedReco_etaWidth"),
            Variable("phiWidth", "gedReco_phiWidth"),
            Variable("full5x5_r9", "gedReco_full5x5_r9"),
            Variable("full5x5_s4", "gedReco_full5x5_s4"),
            Variable("hadronicOverEm", "gedReco_hadronicOverEm"),
            Variable("chargedHadronIso", "gedReco_chargedHadronIso"),
            Variable("neutralHadronIso", "gedReco_neutralHadronIso"),
            Variable("photonIso", "gedReco_photonIso"),
            Variable("scRawEnergy", "gedReco_scRawEnergy"),
            Variable("scAbsEta", "abs(gedReco_scEta)"),
            Variable("trkSumPt", "gedReco_trkSumPtSolidConeDR04"),
            Variable("rho", "rho"),
        ],
        "BoostType=Grad:Shrinkage=0.4:UseBaggedBoost=True:BaggedSampleFraction=0.6:NTrees=2000:MaxDepth=2:nCuts=200",
        allInputFiles[:3]
    ),
]
