j=Job(name='KShort')
myApp = GaudiExec()
myApp.directory = "./MooreDev_v53r4"
j.application = myApp
j.application.options = ['option.py']
j.application.platform = 'x86_64_v3-centos7-gcc11-opt+g'
bkPath = '/MC/Upgrade/Beam7000GeV-Upgrade-MagDown-Nu7.6-25ns-Pythia8/Sim10-Up03-OldP8Tuning/Reco-Up02/11104142/XDST'
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()
j.inputdata = data[0:]
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=25)
j.outputfiles = [LocalFile('output.root')]
j.postprocessors = TextMerger(files=['stdout'])
j.submit()
