from GaudiConf import IOHelper
from Configurables import DaVinci
#from StandardParticles import StdNoPIDsDownPions as Pions
from StandardParticles import StdLoosePions as Pions
from StandardParticles import StdLooseKaons as Kaons
from StandardParticles import StdLooseProtons as Protons

#from StandardParticles import StdLooseANNDownPions as Pions

from CommonParticles import StdLooseKs as Kss
Ks = Kss.StdLooseKsDD

from CommonParticles import StdLooseLambda as Lbd
Lb = Lbd.StdLooseLambdaDD


from PhysConf.Selections import SelectionSequence, CombineSelection, SimpleSelection
import GaudiConfUtils.ConfigurableGenerators as ConfigurableGenerators
from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *


sel_ks = CombineSelection(
		'Sel_Ks',
	#	ConfigurableGenerators.CombineParticles,
		[Pions],
		DecayDescriptor='KS0 -> pi+ pi-',
		DaughtersCuts = Ks.DaughtersCuts,
		CombinationCut = Ks.CombinationCut,
		MotherCut = Ks.MotherCut)

sel_lb = CombineSelection(
		'Sel_Lb',
	#	ConfigurableGenerators.CombineParticles,
		[Pions,Protons],
		DecayDescriptor="[Lambda0 -> p+ pi-]cc",
		DaughtersCuts = Lb.DaughtersCuts,
		CombinationCut = Lb.CombinationCut,
		MotherCut = Lb.MotherCut)

seq_ks = SelectionSequence('Ks_seq',TopSelection=sel_ks)
seq_lb = SelectionSequence('Lb_seq',TopSelection=sel_lb)

##--------------------------------------------
tool_list = [
		"TupleToolKinematic",
		"TupleToolPid",
		"TupleToolGeometry",
		"TupleToolPrimaries",
		"TupleToolTrackInfo",
		"TupleToolEventInfo",
		"TupleToolRecoStats",
		"TupleToolAngles",
#		"TupleToolMCTruth"
		]

dtt_ks = DecayTreeTuple('Ks2pipi')
print("out: ",seq_ks.outputLocations())
dtt_ks.Inputs =  seq_ks.outputLocations() 
dtt_ks.Decay = 'KS0 -> pi+ pi-'
dtt_ks.Branches =  {
	"h1" :	"KS0-> ^pi+ pi-",
	"h2" :	"KS0-> pi+ ^pi-"
}
dtt_ks.ToolList = tool_list
##--------------------------------------------

dtt_lb = DecayTreeTuple('Lambda2ppi')
print("out: ",seq_lb.outputLocations())
dtt_lb.Inputs =  seq_lb.outputLocations() 
dtt_lb.Decay = "[Lambda0 -> ^p+ ^pi-]CC"
dtt_lb.Branches =  {
	"h1" :	"[Lambda0 -> ^p+ pi-]CC",
	"h2" :	"[Lambda0 -> p+ ^pi-]CC"
}
##--------------------------------------------

DaVinci().InputType = 'DST'
DaVinci().TupleFile = 'V0_test_tupel.root'
DaVinci().PrintFreq = 1000
DaVinci().DataType = '2018'
DaVinci().HistogramFile = 'hist.root'
DaVinci().Simulation = False
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().EvtMax = 10000
DaVinci().DDDBtag = 'dddb-20171030-3'
DaVinci().CondDBtag='cond-20180202'
DaVinci().VerboseMessages = True
#DaVinci().RootInTES="/Event/Hlt2"
DaVinci().UserAlgorithms = [seq_ks.sequence(),seq_lb.sequence(),dtt_ks, dtt_lb]
#DaVinci().UserAlgorithms = [Pions,Ks,dt]

IOHelper().inputFiles([
# '00079432_00000058_1.full.dst'

# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision18/FULL.DST/00079432/0000/00079432_00000058_1.full.dst'
'/afs/cern.ch/user/s/skalavan/00079432_00000058_1.full.dst'
	
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000005_2.xdst'
	], clear = True)
