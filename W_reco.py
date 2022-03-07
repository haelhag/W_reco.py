import ROOT 
from ROOT import *
from numpy import *


#rootFiles = loadtxt('Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt', dtype=str, delimiter = "root:")
rootFiles = loadtxt('/nfs/dust/cms/user/sobhatta/work/TopTagPol/TreeMaker/CMSSW_10_5_0/src/sourceFiles/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt', dtype=str, delimiter = " ")

TopHist = ROOT.TH1D('TopHist', '', 100, 70, 500)
WmHist = ROOT.TH1D('WmHist', '', 100, 80.39, 80.41)
file_count = 1
events_count = 0

for rootFile in rootFiles:

    print(rootFile, str(file_count)+'/21 files')

    inFile = ROOT.TFile.Open(str(rootFile),"READ")
    #inFile = ROOT.TFile.Open('root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv7/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/F0D392FD-B5A1-F145-B4F6-C24F3ABC7898.root',"READ")

    tree = inFile.Get('Events')

    #TopHist = ROOT.TH1D('TopHist', '', 100, 70, 500)

    print(tree.GetEntries())

    entries = tree.GetEntries()
    for iEvent in range(0, entries):
        entry = tree.GetEntry(iEvent)

        Ele = ROOT.TLorentzVector()
        nu = ROOT.TLorentzVector()

        Ele_pt = getattr(tree, 'Electron_pt')
        Ele_eta = getattr(tree,'Electron_eta')
        Ele_phi = getattr(tree, 'Electron_phi')
        Ele_mass = getattr(tree, 'Electron_mass')
        nu_pt = getattr(tree,'MET_pt')
        nu_phi = getattr(tree,'MET_phi')
        Muon_pt = getattr(tree,'Muon_pt')
        Muon_eta = getattr(tree,'Muon_eta')
        Jet_pt = getattr(tree,'Jet_pt')
        Jet_eta = getattr(tree,'Jet_eta')
        Jet_phi = getattr(tree,'Jet_phi')
        Jet_mass = getattr(tree,'Jet_mass')
        JetBtag = getattr(tree, 'Jet_btagDeepFlavB')
        Ele_noIso = getattr(tree,'Electron_mvaFall17V2noIso_WP90')
        Muon_tightId = getattr(tree,'Muon_tightId')

        ele_cut_count = 0
        ele_maxpt_idx = 0
        for iEle in range(0,len(Ele_pt)):
            #print(Ele_pt[iEle])
            if Ele_pt[iEle] < 180 or abs(Ele_eta[iEle]) > 2.5 or Ele_noIso[iEle] != 1:
                continue
            ele_cut_count += 1
            if Ele_pt[iEle] > Ele_pt[ele_maxpt_idx]:
                ele_maxpt_idx = iEle        
        
        mu_cut_count = 0
        mu_maxpt_idx = 0 
        for iMu in range(0, len(Muon_pt)):
            if Muon_pt[iMu] < 180 or abs(Muon_eta[iMu]) > 2.1 or Muon_tightId[iMu] != 1:
                continue
            mu_cut_count += 1
            if Muon_pt[iMu] > Muon_pt[mu_maxpt_idx]:
                mu_maxpt_idx = iMu

        if mu_cut_count < 1 and ele_cut_count < 1:
            continue

        jet_cut_count = 0
        jet_Btag_count = 0
        jet_maxpt_idx = 0
        jet_sub_maxpt_idx = 0
        for iJet in range(0, len(Jet_pt)):
            if Jet_pt[iJet] < 30 or abs(Jet_eta[iJet]) > 2.4:
                continue
            jet_cut_count += 1
            if JetBtag[iJet] < 0.2270:
                continue
            jet_Btag_count += 1
            if Jet_pt[iJet] > Jet_pt[jet_maxpt_idx]:
                jet_maxpt_idx = iJet
            if abs(Jet_pt[jet_maxpt_idx] - Jet_pt[iJet]) < abs(Jet_pt[jet_maxpt_idx] - Jet_pt[jet_sub_maxpt_idx]):
                jet_sub_maxpt_idx = iJet
             
        if jet_cut_count < 2 or jet_Btag_count < 1 or Jet_pt[jet_maxpt_idx] < 350:
            continue

        # if nu_pt < 120:
        #     continue

        #if ele_cut_count >= 1:
            #print(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx] )

        if ele_cut_count >= 1:
            Ele = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
            #nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), sqrt(nu_pt**2 + nu.Pz()**2))
            nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), nu.E())
            #nu = ROOT.Math.PtEtaPhiEVector(nu_pt, nu.Eta(), nu_phi, nu.E())
        
        if jet_cut_count >= 2 and ele_cut_count >= 1:
            JetP1 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_sub_maxpt_idx], Jet_eta[jet_sub_maxpt_idx], Jet_phi[jet_sub_maxpt_idx], Jet_mass[jet_sub_maxpt_idx])
            JetP2 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_maxpt_idx], Jet_eta[jet_maxpt_idx], Jet_phi[jet_maxpt_idx], Jet_mass[jet_maxpt_idx])   

        Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))

        # 80.4 GeV is the W-boson pole mass
        if Wm > 80.4:
            k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
            if k < 0.0001:
                k = 0.0001
            scf = 1/2 * 80.4**2/k
            nu.SetPx(nu.Px()*scf)
            nu.SetPy(nu.Py()*scf)
            nu.SetE(nu.P())
            #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                
        lamda = (80.4)**2/2 + Ele.Px()*nu.Px() + Ele.Py()*nu.Py()
            
        if Ele.Pt() > 0 or Ele.Pt() < 0:
            discr = (lamda*Ele.Pz())**2/(Ele.Pt())**4 - ((Ele.E()*nu.Pt())**2 - lamda**2)/(Ele.Pt())**2 

            if discr < 0:
                s = (lamda*Ele.Pz())/(Ele.Pt())**2
                nu.SetPz(s)
                nu.SetE(nu.P())
                #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            else:
                s1 = (lamda*Ele.Pz())/(Ele.Pt())**2 + sqrt(discr)
                s2 = (lamda*Ele.Pz())/(Ele.Pt())**2 - sqrt(discr)

                nu.SetPz(s1)
                nu.SetE(nu.P())
                #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                Wm1 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))

                nu.SetPz(s2)
                nu.SetE(nu.P())
                #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                Wm2 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))

                if abs(Wm2 - 80.4) > abs(Wm1 - 80.4):    
                    nu.SetPz(s1)
                    nu.SetE(nu.P())
                    #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                    
                else:
                    nu.SetPz(s2)
                    nu.SetE(nu.P())
                    #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                        
    
        if ele_cut_count >= 1 and jet_cut_count >= 2:
            #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            Wm = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            print(Wm)
            WmHist.Fill(Wm)
            Pw = Ele + nu
            #Mw2 = Pw.Dot(Pw)
            #print(sqrt(Mw2))

        if (Jet_pt[jet_sub_maxpt_idx] > 25 and Jet_eta[jet_sub_maxpt_idx] < 2.4) and jet_cut_count >= 2 and ele_cut_count >= 1:
            Ptop1 = Pw + JetP1
            Mtop1 = sqrt(Ptop1.Dot(Ptop1))
            Ptop2 = Pw + JetP2
            Mtop2 = sqrt(Ptop2.Dot(Ptop2))
            Ptop12 = Ptop1 + Ptop2
            if abs(Mtop1-172.5) > abs(Mtop2-172.5):
                Mtop = Mtop2
            else:
                Mtop = Mtop1
            if (Ptop1.Pt() > 250 and Ptop12.Pt() > 350) or (Ptop2.Pt() > 250 and Ptop12.Pt() > 350):
                if Mtop > 120 and Mtop < 220:
                    TopHist.Fill(Mtop)
    file_count += 1
    events_count += entries
    if events_count >= 15000:
        break

c2 = ROOT.TCanvas( 'c2', 'W_mass', 1000, 875 )
c2.cd()
WmHist.GetYaxis().SetTitle('Number of events')
WmHist.GetXaxis().SetTitle('W reco mass [GeV]')
WmHist.SetStats(0)
WmHist.SetLineColor(kRed)
WmHist.SetLineWidth(2)
WmHist.Draw()
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c2.SaveAs('W_reco_mass.pdf')

c3 = ROOT.TCanvas( 'c3', 'Top_mass', 1000, 875 )
c3.cd()
TopHist.GetYaxis().SetTitle('Number of events')
TopHist.GetXaxis().SetTitle('Top reco mass [GeV]')
TopHist.SetStats(0)
TopHist.SetLineColor(kRed)
TopHist.SetLineWidth(2)
TopHist.Draw()
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c3.SaveAs('Top_reco_mass.pdf')





    



    


