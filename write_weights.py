import ROOT

from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *

events = Events("/pnfs/desy.de/cms/tier2/store/user/tiroy/TT01j1l_HT800/TT01j1l_HT800/230824_060122/0000/GEN_LO_01j_102X_355.root")

handle = Handle('LHEEventProduct')

events.to(0)

events.getByLabel( "externalLHEProducer", handle )
lhe_product = handle.product()
output_file=ROOT.TFile.Open("WC_weights.root","RECREATE")
output_file.cd()
num=0
for weight in list(lhe_product.weights()):
    #print weight.id, weight.wgt
   # if 't' in weight.id:
	num+=1
	print weight.id, weight.wgt
	hist=ROOT.TH1D(weight.id,weight.id, 10000,0,1)
	hist.SetBinContent(1, weight.wgt)
#	hist.Sumw2()
	hist.Write()
print num
output_file.Close()
 
