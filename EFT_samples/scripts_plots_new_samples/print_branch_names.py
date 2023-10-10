import ROOT

file = ROOT.TFile.Open("root://eos.grid.vbc.ac.at//store/user/schoef/TT01j1lCA_HT500_v2/TT01j1lCA_HT500_v2/230918_163019/0000/GEN_LO_0j_102X_1.root")
tree = file.Get("Events")

for branch in tree.GetListOfBranches():
    print(branch.GetName())

file.Close()
