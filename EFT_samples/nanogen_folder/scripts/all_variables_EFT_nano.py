import ROOT
import os
import glob
from array import array
import subprocess
from DataFormats.FWLite import Events, Handle
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Config as edm
from XRootD import client
from XRootD.client.flags import DirListFlags, StatInfoFlags, OpenFlags, QueryCode



ROOT.gROOT.SetBatch(True)

totalEvents = 0

h_leptonPt = ROOT.TH1F("h_leptonPt", "Lepton pT; pT (GeV);Events", 1000, 0, 1000)
h_topPt = ROOT.TH1F("h_topPt", "Top Quark pT; pT (GeV);Events", 1000, 0, 3000)
h_antitopPt = ROOT.TH1F("h_antitopPt", "Anti-Top Quark p_{T}; p_{T} [GeV];Events", 1000, 0, 3000)
h_leptoneta = ROOT.TH1F("h_leptoneta", "eta; #eta;Events", 100, -5, 5)
h_leptonphi = ROOT.TH1F("h_leptonphi", "Azimuthal Angle; #phi;Events", 100, -ROOT.TMath.Pi(), ROOT.TMath.Pi())
h_invariantMass = ROOT.TH1F("h_invariantMass", "Invariant Mass; M (GeV);Events", 100, 0, 7000)
h_partonMultiplicity = ROOT.TH1F("h_partonMultiplicity", "Jet Multiplicity; N_{jets};Events", 20, 0, 100)
h_MET = ROOT.TH1F("hMET", "MET;MET (GeV);Events", 100, 0, 200)
h_bquark_pt = ROOT.TH1F("hbquarkPt", "b-quark pT;pT (GeV);Events", 150, 0, 1000)
h_bquark_eta = ROOT.TH1F("hbquarkEta", "b-quark #eta;#eta;Events", 100, -5, 5)
h_angle_top_antitop = ROOT.TH1F("h_angle", "Angle between top and antitop;Angle (radians);Events", 50, 0, ROOT.TMath.Pi())

h_decayChannel = ROOT.TH1F("h_decayChannel", "Top Decay Channels; Channel; Events", 2, 0, 2)
h_decayChannel.GetXaxis().SetBinLabel(1, "t -> W+b")
h_decayChannel.GetXaxis().SetBinLabel(2, "Other")


h_topMultiplicity = ROOT.TH1F("h_topMultiplicity", "Top Multiplicity; N_{top};Events", 5, 0, 5)
# h_antitopMultiplicity = ROOT.TH1F("h_antitopMultiplicity", "Anti-Top Multiplicity; N_{antitop};Events", 5, 0, 5)

h_missingParticles = ROOT.TH1F("h_missingParticles", "Missing Particles; Particle Type; Events", 4, 0, 4)
h_missingParticles.GetXaxis().SetBinLabel(1, "No Top")
h_missingParticles.GetXaxis().SetBinLabel(2, "No Anti-Top")
h_missingParticles.GetXaxis().SetBinLabel(3, "No Top & No Anti-Top")



bin_edges = [-16.5, -14.5, -12.5, -10.5, 10.5, 12.5, 14.5, 16.5]
h_leptonFlavor = ROOT.TH1F("h_leptonFlavor", "Lepton Flavor; PDG ID;Events", len(bin_edges)-1, array('d', bin_edges))

h_leptonFlavor.GetXaxis().SetBinLabel(1, "muon-")
h_leptonFlavor.GetXaxis().SetBinLabel(2, "electron-")
h_leptonFlavor.GetXaxis().SetBinLabel(4, "electron+")
h_leptonFlavor.GetXaxis().SetBinLabel(5, "muon+")


h_nonTopMotherJets = ROOT.TH1F("h_nonTopMotherJets", "Jets without Top as Mother; Count;Events", 10, 0, 50)
h_jetMultiplicity = ROOT.TH1F("h_jetMultiplicity", "Number of Jets per Event", 10, 0, 50)

h_topMother = ROOT.TH1F("h_topMother", "Mother of Top Quarks; Mother; Events", 3, 0, 3)
h_topMother.GetXaxis().SetBinLabel(1, "qq")
h_topMother.GetXaxis().SetBinLabel(2, "gg")
h_topMother.GetXaxis().SetBinLabel(3, "Other")

h_motherPdgId = ROOT.TH1F("h_motherPdgId", "PDG ID of Top's Mother;PDG ID;Counts", 23, -6, 22)

def analyze(filename):
    print("Processing file:", filename)
    
    file = ROOT.TFile.Open(filename)
    tree = file.Get("Events")
    
    global totalEvents
    totalEvents += tree.GetEntries()
    print("Number of events in file:", tree.GetEntries())
    
    relevant_pdgIds = {12, 14, 16, 24, 1, 2, 3, 4, 5, 6, 21, 11, 13, 15}
    
    tops = []
    
    for entry in tree:
        nGenPart = entry.nGenPart
        nGenJet = entry.nGenJet
        
        top_count = 0
        antitop_count = 0
        partons = []
        
        met_vector = ROOT.TLorentzVector()
        
        for i in range(nGenPart):
            pdgId = entry.GenPart_pdgId[i]
            pt = entry.GenPart_pt[i]
            eta = entry.GenPart_eta[i]
            phi = entry.GenPart_phi[i]
            mass = entry.GenPart_mass[i]
            mother_idx = entry.GenPart_genPartIdxMother[i]
            
            if abs(pdgId) in relevant_pdgIds:
            
                # Tops
                if abs(pdgId) == 6: 
                    daughters = [entry.GenPart_pdgId[j] for j in range(nGenPart) if entry.GenPart_genPartIdxMother[j] == i]

                    if 6 in [abs(d) for d in daughters]:
                        continue
                    
                    # mothers of the top and antitop
                    mother1_pdgId = entry.GenPart_pdgId[mother_idx] if mother_idx >= 0 else None
                    # Mother 2 is more complicated in NanoAOD. We assume one mother here for simplicity.
                    if mother1_pdgId:
                        h_motherPdgId.Fill(mother1_pdgId)
                        
                    w_quark_daughter = 24 in daughters
                    b_quark_daughter = 5 in daughters
                    
                    if not w_quark_daughter or not b_quark_daughter:
                        continue
                    
                    h_decayChannel.Fill(0)  # t -> W+b
                    
                    top_4vec = ROOT.TLorentzVector()
                    top_4vec.SetPtEtaPhiM(pt, eta, phi, mass)
                    tops.append(top_4vec)
                    
                    if pdgId == 6:
                        top_count += 1
                        h_topPt.Fill(pt)
                    elif pdgId == -6:
                        antitop_count += 1
                        h_antitopPt.Fill(pt)
                    
                    # leptons from W decay
                    leptons_from_W = [entry.GenPart_pdgId[j] for j in range(nGenPart) if entry.GenPart_genPartIdxMother[j] == i and abs(entry.GenPart_pdgId[j]) in [11, 13]]
                    if leptons_from_W:
                        lepton_pdgId = leptons_from_W[0]
                        h_leptonPt.Fill(entry.GenPart_pt[i])
                        h_leptoneta.Fill(eta)
                        h_leptonphi.Fill(phi)
                        h_leptonFlavor.Fill(lepton_pdgId)
                            
                    
                    if mother1_pdgId == 21:  # gg
                        h_topMother.Fill(1)
                    elif mother1_pdgId in [1,2,3,4,5]:  # qq
                        h_topMother.Fill(0)
                    else:
                        h_topMother.Fill(2)
            
                
                # Partons
                if abs(pdgId) in [1, 2, 3, 4, 5, 6, 21]:
                    partons.append(i)  # Store index of the parton for later reference

                # b-quarks
                if abs(pdgId) == 5:
                    b_vector = ROOT.TLorentzVector()
                    b_vector.SetPtEtaPhiM(pt, eta, phi, mass)
                    
                    if pdgId == 5:
                        h_bquark_pt.Fill(b_vector.Pt())
                        h_bquark_eta.Fill(b_vector.Eta())
            
                if abs(pdgId) in [12, 14, 16]:
                    neutrino = ROOT.TLorentzVector()
                    neutrino.SetPtEtaPhiM(pt, eta, phi, mass)
                    met_vector += neutrino
                
        h_MET.Fill(met_vector.Pt())
                
        h_partonMultiplicity.Fill(len(partons))
        
        # Fill histograms related to tops
        if top_count == 0:
            h_missingParticles.Fill(0)  # Filling "no top" bin

        if antitop_count == 0:
            h_missingParticles.Fill(1)  # Filling "no antitop" bin

        if top_count == 0 and antitop_count == 0:
            h_missingParticles.Fill(2)  # Filling "no top and antitop" bin

        # If we have both top and antitop, we can calculate the invariant mass and angle
        top_idx = next((idx for idx, pdg in enumerate(entry.GenPart_pdgId) if pdg == 6), None)
        antitop_idx = next((idx for idx, pdg in enumerate(entry.GenPart_pdgId) if pdg == -6), None)

        if top_idx is not None and antitop_idx is not None:
            top_4vec = ROOT.TLorentzVector()
            top_4vec.SetPtEtaPhiM(entry.GenPart_pt[top_idx], entry.GenPart_eta[top_idx], entry.GenPart_phi[top_idx], entry.GenPart_mass[top_idx])
            
            antitop_4vec = ROOT.TLorentzVector()
            antitop_4vec.SetPtEtaPhiM(entry.GenPart_pt[antitop_idx], entry.GenPart_eta[antitop_idx], entry.GenPart_phi[antitop_idx], entry.GenPart_mass[antitop_idx])

            ttbar = top_4vec + antitop_4vec
            h_invariantMass.Fill(ttbar.M())
            h_angle_top_antitop.Fill(top_4vec.Angle(antitop_4vec.Vect()))

        h_topMultiplicity.Fill(top_count)
        
        non_top_mother_jet_count_j = []
        jet_count = 0
        
        def deltaR(eta1, phi1, eta2, phi2):
            deta = eta1 - eta2
            dphi = abs(phi1 - phi2)
            if dphi > ROOT.TMath.Pi():
                dphi = 2*ROOT.TMath.Pi() - dphi
            return (deta*deta + dphi*dphi)**0.5
        
        GenJet_eta = entry.GenJet_eta
        GenJet_phi = entry.GenJet_phi
        non_top_mother_jet_count_j = []
        jet_count = 0
        
        for i in range(len(GenJet_eta)):
            jet_count += 1

            jet_eta = GenJet_eta[i]
            jet_phi = GenJet_phi[i]
            
            is_from_top = False
            for top in tops:  
                if deltaR(jet_eta, jet_phi, top.Eta(), top.Phi()) < 0.4: 
                    is_from_top = True
                    break

            if not is_from_top:
                non_top_mother_jet_count_j.append(i)
        
        h_jetMultiplicity.Fill(jet_count)
        h_nonTopMotherJets.Fill(len(non_top_mother_jet_count_j))
    # deltaR calculates the deltaR distance in eta - phi space. 
    # If this distance is less than a threshold (like 0.4, a typical jet size), we can say the jet is possibly from the top quark.         
        
    file.Close()

# url = "davs://dcache-cms-webdav-wan.desy.de:2880/"
# path = "/pnfs/desy.de/cms/tier2/store/user/beozek/TT01j1lCA_HT500_v2/TT01j1lCA_HT500_v2/231004_134141/0000"
# client_instance = client.FileSystem(url)
# status, listing = client_instance.dirlist(path, DirListFlags.STAT)
# root_files = [entry.name for entry in listing if entry.name.endswith('.root')]

# for root_file in root_files:
#     full_path = url + os.path.join(path, root_file)
#     analyze(full_path)
    
path = "/nfs/dust/cms/user/beozek/EFT/CMSSW_10_6_27/src/EFT_gen_old/EFT_samples/nanogen_folder/nano_files/"
root_files = [f for f in os.listdir(path) if f.endswith('.root')]
root_files = root_files[:5]

for root_file in root_files:
    full_path = os.path.join(path, root_file)
    analyze(full_path)
    
    
canvas = ROOT.TCanvas("canvas", "Analysis Plots", 4000, 4000)
canvas.Divide(4, 5)

canvas.cd(1)
h_leptonPt.Draw()
canvas.cd(2)
h_topPt.Draw()
canvas.cd(3)
h_antitopPt.Draw()
canvas.cd(4)
h_leptoneta.Draw()
canvas.cd(5)
h_leptonphi.Draw()
canvas.cd(6)
h_invariantMass.Draw()
canvas.cd(7)
h_decayChannel.Draw()
canvas.cd(8)
h_MET.Draw()
canvas.cd(9)
h_leptonFlavor.Draw()
canvas.cd(10)
h_bquark_pt.Draw()
canvas.cd(11)
h_bquark_eta.Draw()
canvas.cd(12)
h_angle_top_antitop.Draw()
canvas.cd(13)
h_partonMultiplicity.Draw()
canvas.cd(14)
h_nonTopMotherJets.Draw()
canvas.cd(15)
h_topMultiplicity.Draw()
canvas.cd(16)
h_jetMultiplicity.Draw()
canvas.cd(17)
h_topMother.Draw()
canvas.cd(18)
h_missingParticles.Draw()
canvas.cd(19)
h_motherPdgId.Draw()


canvas.SaveAs("allPlots.png")

c_leptonPt = ROOT.TCanvas("c_leptonPt", "Lepton pT Distribution", 800, 600)
h_leptonPt.Draw()
ROOT.gPad.SetLogy(1)
c_leptonPt.SaveAs("leptonPtDistribution.png")

c_topPt = ROOT.TCanvas("c_topPt", "Top Quark pT Distribution", 800, 600)
h_topPt.Draw()
ROOT.gPad.SetLogy(1)
c_topPt.SaveAs("topPtDistribution.png")

c_antitopPt = ROOT.TCanvas("c_antitopPt", "Anti-Top Quark pT Distribution", 800, 600)
h_antitopPt.Draw()
ROOT.gPad.SetLogy(1)
c_antitopPt.SaveAs("antitopPtDistribution.png")

c_eta = ROOT.TCanvas("c_eta", "Lepton Eta Distribution", 800, 600)
h_leptoneta.Draw()
c_eta.SaveAs("etaDistribution.png")

c_phi = ROOT.TCanvas("c_phi", "Lepton Azimuthal Angle Distribution", 800, 600)
h_leptonphi.Draw()
c_phi.SaveAs("phiDistribution.png")

c_invariantMass = ROOT.TCanvas("c_invariantMass", "Invariant Mass Distribution", 800, 600)
h_invariantMass.Draw()
ROOT.gPad.SetLogy(1)
c_invariantMass.SaveAs("invariantMassDistribution.png")

c_decay = ROOT.TCanvas("c_decay", "Decay Channel Canvas", 800, 600)
h_decayChannel.Draw()
c_decay.SaveAs("topDecayChannel.png")

c_MET = ROOT.TCanvas("cMET", "MET Distribution", 800, 600)
h_MET.Draw()
ROOT.gPad.SetLogy(1)
c_MET.SaveAs("METDistribution.png")

c_leptonFlavor = ROOT.TCanvas("c_leptonFlavor", "Lepton Flavor Distribution", 800, 600)
h_leptonFlavor.Draw()
c_leptonFlavor.SaveAs("leptonFlavorDistribution.png")

c_bquark_pt = ROOT.TCanvas("cbquarkPt", "b-quark pT Distribution", 800, 600)
h_bquark_pt.Draw()
ROOT.gPad.SetLogy(1)
c_bquark_pt.SaveAs("bquarkPtDistribution.png")

c_bquark_eta = ROOT.TCanvas("cbquarkEta", "b-quark Eta Distribution", 800, 600)
h_bquark_eta.Draw()
c_bquark_eta.SaveAs("bquarkEtaDistribution.png")

c_angle = ROOT.TCanvas("cangle", "Angle between top and antitop", 800, 600)
h_angle_top_antitop.Draw()
c_angle.SaveAs("angleTopAntitop.png")

c_partonMultiplicity = ROOT.TCanvas("c_partonMultiplicity", "Parton Multiplicity Distribution", 800, 600)
h_partonMultiplicity.SetFillColor(ROOT.kBlue - 10)
h_partonMultiplicity.SetLineColor(ROOT.kBlue)
h_partonMultiplicity.Draw()
ROOT.gPad.SetLogy(1)
c_partonMultiplicity.SaveAs("partonMultiplicityDistribution.png")

c_nonTopMotherJets = ROOT.TCanvas("c_nonTopMotherJets", "Jets without Top as Mother", 800, 600)
h_nonTopMotherJets.SetFillColor(ROOT.kBlue - 10)
h_nonTopMotherJets.SetLineColor(ROOT.kBlue)
h_nonTopMotherJets.Draw()
ROOT.gPad.SetLogy(1)
c_nonTopMotherJets.SaveAs("nonTopMotherJets.png")

# c_antitopMultiplicity = ROOT.TCanvas("c_antitopMultiplicity", "Anti-Top Multiplicity Distribution", 800, 600)
# h_antitopMultiplicity.SetFillColor(ROOT.kBlue - 10)
# h_antitopMultiplicity.SetLineColor(ROOT.kBlue)
# h_antitopMultiplicity.Draw()
# c_antitopMultiplicity.SaveAs("antitopMultiplicityDistribution.png")

c_topMultiplicity = ROOT.TCanvas("c_topMultiplicity", "Top Multiplicity Distribution", 800, 600)
h_topMultiplicity.SetFillColor(ROOT.kBlue - 10)
h_topMultiplicity.SetLineColor(ROOT.kBlue)
h_topMultiplicity.Draw()
ROOT.gPad.SetLogy(1)
c_topMultiplicity.SaveAs("topMultiplicityDistribution.png")

c_jetMultiplicity = ROOT.TCanvas("c_jetMultiplicity", "Number of Jets per Event", 800, 600)
h_jetMultiplicity.SetFillColor(ROOT.kBlue - 10)
h_jetMultiplicity.SetLineColor(ROOT.kBlue)
h_jetMultiplicity.Draw()
ROOT.gPad.SetLogy(1)
c_jetMultiplicity.SaveAs("jetMultiplicity.png")

c_topMother = ROOT.TCanvas("c_topMother", "Mothers of the top quark", 800, 600)
h_topMother.Draw()
c_topMother.SaveAs("topMother.png")

c_missingParticles = ROOT.TCanvas("c_missingParticles", "Missing Particles", 800, 600)
h_missingParticles.Draw()
ROOT.gPad.SetLogy(1)
c_missingParticles.SaveAs("missingpart.png")

c_motherPdgId = ROOT.TCanvas("c_motherPdgId", "PDG ID of Top's Mother", 800,600)
h_motherPdgId.Draw()
c_motherPdgId.SaveAs("motherPDG.png")



print("Total number of events:", totalEvents)



# analyze("/nfs/dust/cms/user/beozek/EFT/CMSSW_10_6_27/src/EFT_gen_old/root_files/0000/GEN_LO_01j_102X_14.root")