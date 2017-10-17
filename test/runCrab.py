from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_phase2photons'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'testPhase2PhotonTuples.py'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/%s/932phoID' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
# config.Data.totalUnits = 20

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

pds = pds_relval

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
        config.General.requestName = 'p2phoID_%d_%s' % (i, primaryDS)
        config.Data.outputDatasetTag = conditions
        config.Data.inputDataset = pd
        if dataTier == 'MINIAODSIM':
            config.Data.useParent = True
        submit(config)
