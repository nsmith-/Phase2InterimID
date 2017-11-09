from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_phase2photons_round4'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'testPhase2PhotonTuples.py'

config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
#config.General.instance = 'preprod'

config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/%s/932phoID_round4' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

#config.Site.storageSite = 'T2_US_Wisconsin'
config.Site.storageSite = 'T2_CH_CERN'

pds_relval = [
    '/RelValH125GGgluonfusion_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/GEN-SIM-RECO',
    '/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-RECO',
    '/RelValSingleElectronPt15Eta1p7_2p7/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-RECO',
    '/RelValSingleGammaPt25Eta1p7_2p7/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-RECO',
    '/RelValSinglePiPt25Eta1p7_2p7/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-RECO',
    '/RelValTTbar_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-RECO',
    '/RelValZEE_14/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-RECO',
]

pds_photon = [
    '/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
    '/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
    '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
]

pds_photon_noPU = [
    '/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-noPU_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
    '/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/PhaseIITDRFall17DR-noPU_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRFall17DR-noPU_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
    '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17DR-noPU_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO',
]

pds_run2 = [
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM',
]

pds = pds_photon
# run2 = 50, photon = 30, photonpu0 = 10
config.Data.unitsPerJob = 30
#config.JobType.pyCfgParams = ["phase2=0"]
config.JobType.numCores = 4
# config.Data.totalUnits = 20*50

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    for i, pd in enumerate(pds):
        (_, primaryDS, conditions, dataTier) = pd.split('/')
        config.General.requestName = 'p2phoIDr3_%d_%s' % (i, primaryDS)
        config.Data.outputDatasetTag = conditions
        config.Data.inputDataset = pd
        if dataTier == 'MINIAODSIM':
            config.Data.useParent = True
        submit(config)
