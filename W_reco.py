import ROOT 
from ROOT import *
from numpy import *


#rootFiles = loadtxt('Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt', dtype=str, delimiter = "root:")
rootFiles = loadtxt('/nfs/dust/cms/user/sobhatta/work/TopTagPol/TreeMaker/CMSSW_10_5_0/src/sourceFiles/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt', dtype=str, delimiter = " ")

TopHist = ROOT.TH1D('TopHist', '', 100, 0, 500)
MaxPtJetHist = ROOT.TH1D('MaxPtJetHist', '', 100, 0, 1000)
SubMaxPtJetHist = ROOT.TH1D('SubMaxPtJetHist', '', 100, 0, 1000)
WmTHist = ROOT.TH1D('WmTHist', '', 100, 0, 200)
WmIHist = ROOT.TH1D('WmIHist', '', 200, 80.395, 80.405)
WmtHist = ROOT.TH1D('WmtHist', '', 100, 0, 1400)
LepNuMtHist = ROOT.TH1D('LepNuMtHist', '', 100, 0, 1400)
DrEleJet1Hist = ROOT.TH1D('DrEleJet1Hist', '', 100, 0, 10)
DrEleJet2Hist = ROOT.TH1D('DrEleJet2Hist', '', 100, 0, 10)
DrMuJet1Hist = ROOT.TH1D('DrMuJet1Hist', '', 100, 0, 10)
DrMuJet2Hist = ROOT.TH1D('DrMuJet2Hist', '', 100, 0, 10)
DrEleJet1TopHist = ROOT.TH1D('DrEleJet1TopHist', '', 100, 0, 10)
DrEleJet2TopHist = ROOT.TH1D('DrEleJet2TopHist', '', 100, 0, 10)
DrMuJet1TopHist = ROOT.TH1D('DrMuJet1TopHist', '', 100, 0, 10)
DrMuJet2TopHist = ROOT.TH1D('DrMuJet2TopHist', '', 100, 0, 10)
file_count = 1
events_count = 0

lep_channel = 0

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
        mu = ROOT.TLorentzVector()

        Ele_pt = getattr(tree, 'Electron_pt')
        Ele_eta = getattr(tree,'Electron_eta')
        Ele_phi = getattr(tree, 'Electron_phi')
        Ele_mass = getattr(tree, 'Electron_mass')
        Ele_noIso = getattr(tree,'Electron_mvaFall17V2noIso_WP90')
        Ele_miniIso = getattr(tree,'Electron_miniPFRelIso_all')

        nu_pt = getattr(tree,'MET_pt')
        nu_phi = getattr(tree,'MET_phi')

        Muon_pt = getattr(tree,'Muon_pt')
        Muon_eta = getattr(tree,'Muon_eta')
        Muon_phi = getattr(tree,'Muon_phi')
        Muon_mass = getattr(tree,'Muon_mass')
        Muon_tightId = getattr(tree,'Muon_tightId')
        Muon_miniIso = getattr(tree,'Muon_miniPFRelIso_all')

        Jet_pt = getattr(tree,'Jet_pt')
        Jet_eta = getattr(tree,'Jet_eta')
        Jet_phi = getattr(tree,'Jet_phi')
        Jet_mass = getattr(tree,'Jet_mass')
        JetBtag = getattr(tree, 'Jet_btagDeepFlavB')

        GenPart_pt = getattr(tree, 'GenPart_pt')
        GenPart_eta = getattr(tree,'GenPart_eta')
        GenPart_phi = getattr(tree, 'GenPart_phi')
        GenPart_mass = getattr(tree, 'GenPart_mass')
        GenPart_genPartIdxMother = getattr(tree,'GenPart_genPartIdxMother')
        GenPart_pdgId = getattr(tree,'GenPart_pdgId')
        
        

        ele_cut_count = 0
        ele_maxpt_idx = 0
        #ele_maxpt = max(Ele_pt)
        for iEle in range(0,len(Ele_pt)):
            if Ele_pt[iEle] < 100 or abs(Ele_eta[iEle]) > 2.5 or Ele_noIso[iEle] != 1 or Ele_miniIso[iEle] > 0.1:
                continue
            ele_cut_count += 1
            if Ele_pt[iEle] > Ele_pt[ele_maxpt_idx]:
            #if Ele_pt[iEle] == ele_maxpt:
                ele_maxpt_idx = iEle
        ele_gen_count = 0
        if ele_cut_count >= 1:  
            RecEle = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
            for iPart in range(0, len(GenPart_pt)):
                mother_idx = GenPart_genPartIdxMother[iPart]
                if (mother_idx >= 0 and abs(GenPart_pdgId[mother_idx]) == 24) and abs(GenPart_pdgId[iPart]) == 11:
                    if GenPart_pt[iPart] < 100 and abs(GenPart_eta[iPart]) > 2.5:
                        continue
                    GenEle = ROOT.Math.PtEtaPhiMVector(GenPart_pt[iPart], GenPart_eta[iPart], GenPart_phi[iPart], GenPart_mass[iPart])
                    DrGenRecEle = ROOT.Math.VectorUtil.DeltaR(RecEle, GenEle)
                    if DrGenRecEle < 0.3:
                        ele_gen_count += 1
                    print('GenPart_pdgId=',GenPart_pdgId[iPart])
                    print('GenPart_pt=',GenPart_pt[iPart])
                    print('Ele_pt',Ele_pt[ele_maxpt_idx])
            if ele_gen_count < 1:
                continue
                    
        
        mu_cut_count = 0
        mu_maxpt_idx = 0
        #mu_maxpt = max(Muon_pt) 
        for iMu in range(0, len(Muon_pt)):
            if Muon_pt[iMu] < 50 or abs(Muon_eta[iMu]) > 2.1 or Muon_tightId[iMu] != 1 or Muon_miniIso[iMu] > 0.1:
                continue
            mu_cut_count += 1
            if Muon_pt[iMu] > Muon_pt[mu_maxpt_idx]:
            #if Muon_pt[iMu] == mu_maxpt:
                mu_maxpt_idx = iMu
        mu_gen_count = 0
        if mu_cut_count >= 1:
            RecMu = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
            for iPart in range(0, len(GenPart_pt)):
                mother_idx = GenPart_genPartIdxMother[iPart]
                if (mother_idx >= 0 and abs(GenPart_pdgId[mother_idx]) == 24) and abs(GenPart_pdgId[iPart]) == 13:
                    if GenPart_pt[iPart] < 50 and abs(GenPart_eta[iPart]) > 2.4:
                        continue
                    GenEle = ROOT.Math.PtEtaPhiMVector(GenPart_pt[iPart], GenPart_eta[iPart], GenPart_phi[iPart], GenPart_mass[iPart])
                    DrGenRecEle = ROOT.Math.VectorUtil.DeltaR(RecMu, GenEle)
                    if DrGenRecEle < 0.3:
                        mu_gen_count += 1
                    print('GenPart_pdgId=',GenPart_pdgId[iPart])
                    print('GenPart_pt=',GenPart_pt[iPart])
                    print('Muon_pt',Muon_pt[mu_maxpt_idx])
            if mu_gen_count < 1:
                continue

        
        #Jet_pt_list = []       
        jet_cut_count = 0
        jet_Btag_count = 0
        jet_maxpt_idx = -1
        jet_sub_maxpt_idx = -1
        jet_maxpt = max(Jet_pt)
        for iJet in range(0, len(Jet_pt)):
            if Jet_pt[iJet] < 100 or abs(Jet_eta[iJet]) > 2.4:
                continue
            jet_cut_count += 1
            #Jet_pt_list.append(Jet_pt[iJet])
            if Jet_pt[iJet] == jet_maxpt:
                jet_maxpt_idx = iJet
            if iJet != jet_maxpt_idx:
                if Jet_pt[iJet] > Jet_pt[jet_sub_maxpt_idx]:
                    jet_sub_maxpt_idx = iJet
        
        if jet_cut_count >=2:
            #print('max jet:',jet_maxpt)
            MaxPtJetHist.Fill(jet_maxpt)
            #print('sub max jet:',Jet_pt[jet_sub_maxpt_idx])
            SubMaxPtJetHist.Fill(Jet_pt[jet_sub_maxpt_idx])
        
        if jet_cut_count < 2 or Jet_pt[jet_maxpt_idx] < 300 or Jet_pt[jet_sub_maxpt_idx] < 150:
            continue

        if JetBtag[jet_maxpt_idx] > 0.2270 or JetBtag[jet_sub_maxpt_idx] > 0.2270:
            jet_Btag_count += 1

        if jet_Btag_count < 1:
            continue

        #if nu_pt < 120:
        #    continue

        #if ele_cut_count >= 1:
            #print(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx] )

        if ((ele_cut_count >= 1 and ele_gen_count >= 1)  or (mu_cut_count >=1 and mu_gen_count >= 1)) and jet_cut_count >= 2:
            lep = ROOT.Math.PtEtaPhiMVector()
            nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), sqrt(nu_pt**2 + nu.Pz()**2))

            if ele_cut_count >= 1 and mu_cut_count >= 1:
                if Ele_pt[ele_maxpt_idx] > Muon_pt[mu_maxpt_idx]:
                    lep = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
                    lep_channel = 0
                else:
                    lep = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
                    lep_channel = 1
            elif ele_cut_count >= 1 and mu_cut_count < 1:
                lep = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
                lep_channel = 0
            elif ele_cut_count < 1 and mu_cut_count >=1:
                lep = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
                lep_channel = 1

            
            JetP1 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_sub_maxpt_idx], Jet_eta[jet_sub_maxpt_idx], Jet_phi[jet_sub_maxpt_idx], Jet_mass[jet_sub_maxpt_idx])
            JetP2 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_maxpt_idx], Jet_eta[jet_maxpt_idx], Jet_phi[jet_maxpt_idx], Jet_mass[jet_maxpt_idx])   
            
            # WmT2 = (lep.Pt()+nu.Et())**2-(lep.Px()+nu.Px())**2-(lep.Py()+nu.Py())**2

            # print('WmT^2=',WmT2)

            WmT =  sqrt((lep.Pt()+nu.Et())**2-(lep.Px()+nu.Px())**2-(lep.Py()+nu.Py())**2)

            # 80.4 GeV is the W-boson pole mass
            if WmT > 80.4:
                k = nu.Et() * lep.Pt() - nu.Px() * lep.Px() - nu.Py() * lep.Py()
                if k < 0.0001:
                    k = 0.0001
                scf = 1/2 * 80.4**2/k
                nu.SetPx(nu.Px()*scf)
                nu.SetPy(nu.Py()*scf)
                nu.SetE(nu.P())
         
            lamda = (80.4)**2/2 + lep.Px()*nu.Px() + lep.Py()*nu.Py()

            discr = (lamda*lep.Pz())**2/(lep.Pt())**4 - ((lep.E()*nu.Pt())**2 - lamda**2)/(lep.Pt())**2 

            if WmT > 80.4 or discr < 0:
                s = (lamda*lep.Pz())/(lep.Pt())**2
                nu.SetPz(s)
                nu.SetE(nu.P())
            else:
                s1 = (lamda*lep.Pz())/(lep.Pt())**2 + sqrt(discr)
                s2 = (lamda*lep.Pz())/(lep.Pt())**2 - sqrt(discr)


                nu.SetPz(s1)
                nu.SetE(nu.P())
                Wm1 =  sqrt((lep.Pt()+nu.Et())**2-(lep.Px()+nu.Px())**2-(lep.Py()+nu.Py())**2)

                nu.SetPz(s2)
                nu.SetE(nu.P())
                Wm2 =  sqrt((lep.Pt()+nu.Et())**2-(lep.Px()+nu.Px())**2-(lep.Py()+nu.Py())**2)

                if abs(Wm2 - 80.4) > abs(Wm1 - 80.4):    
                    nu.SetPz(s1)
                    nu.SetE(nu.P())
                    
                else:
                    nu.SetPz(s2)
                    nu.SetE(nu.P())
                        

            WmI = sqrt((lep.E()+nu.E())**2-(lep.Px()+nu.Px())**2-(lep.Py()+nu.Py())**2-(lep.Pz()+nu.Pz())**2)
            WmT = sqrt((lep.Pt()+nu.Et())**2-(lep.Px()+nu.Px())**2-(lep.Py()+nu.Py())**2)
            #change this when you want to plot for one channel only
            if lep_channel == 0 or lep_channel == 1:
                print('WmT=',WmT)
                WmTHist.Fill(WmT)
                Pw = lep + nu
                WmtHist.Fill(Pw.Mt())
                WmIHist.Fill(WmI)
                M12T2 = (lep.Et()+nu.Et())**2 - (lep.Pt()+nu.Pt())**2
                LepNuMtHist.Fill(sqrt(M12T2))

                Ptop1 = Pw + JetP1
                Mtop1 = Ptop1.M()
                Ptop2 = Pw + JetP2
                Mtop2 = Ptop2.M()
                if abs(Mtop1-172.5) > abs(Mtop2-172.5):
                    Mtop = Mtop2
                    if lep_channel == 1:
                        DrMuonJet2Top = ROOT.Math.VectorUtil.DeltaR(lep, JetP2)
                        DrMuJet2TopHist.Fill(DrMuonJet2Top)
                    elif lep_channel == 0:
                        DrEleJet2Top = ROOT.Math.VectorUtil.DeltaR(lep, JetP2)
                        DrEleJet2TopHist.Fill(DrEleJet2Top)

                else:
                    Mtop = Mtop1
                    if lep_channel == 1:
                        DrMuonJet1Top = ROOT.Math.VectorUtil.DeltaR(lep, JetP1)
                        DrMuJet1TopHist.Fill(DrMuonJet1Top)
                    elif lep_channel == 0:
                        DrEleJet1Top = ROOT.Math.VectorUtil.DeltaR(lep, JetP1)
                        DrEleJet1TopHist.Fill(DrEleJet1Top)
                TopHist.Fill(Mtop)
            
            # print('ele_cut_count:',ele_cut_count)
            # print('mu_cut_count:',mu_cut_count)
            # print('jet_cut_count:',jet_cut_count)
            # if lep == ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx]) :
            
            if lep_channel == 1:
                DrMuonJet1 = ROOT.Math.VectorUtil.DeltaR(lep, JetP1)
                DrMuJet1Hist.Fill(DrMuonJet1)

                DrMuonJet2 = ROOT.Math.VectorUtil.DeltaR(lep, JetP2)
                DrMuJet2Hist.Fill(DrMuonJet2)

            elif lep_channel == 0:
                DrEleJet1 = ROOT.Math.VectorUtil.DeltaR(lep, JetP1)
                DrEleJet1Hist.Fill(DrEleJet1)

                DrEleJet2 = ROOT.Math.VectorUtil.DeltaR(lep, JetP2)
                DrEleJet2Hist.Fill(DrEleJet2)
          
            # if ele_cut_count >= 1 and mu_cut_count >= 1:
            #     if Ele_pt[ele_maxpt_idx] > Muon_pt[mu_maxpt_idx]:
            #         Ele = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
            #         nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), sqrt(nu_pt**2 + nu.Pz()**2))
            #         #mu = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
                    
            #         JetP1 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_sub_maxpt_idx], Jet_eta[jet_sub_maxpt_idx], Jet_phi[jet_sub_maxpt_idx], Jet_mass[jet_sub_maxpt_idx])
            #         JetP2 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_maxpt_idx], Jet_eta[jet_maxpt_idx], Jet_phi[jet_maxpt_idx], Jet_mass[jet_maxpt_idx])   
                    
            #         WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

            #         print('WmT^2=',WmT2)

            #         #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #         WmT =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
            #         #WmT =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)
            #         #WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)

            #         # 80.4 GeV is the W-boson pole mass
            #         if WmT > 80.4:
            #             k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
            #             #k = nu.Et() * mu.Pt() - nu.Px() * mu.Px() - nu.Py() * mu.Py()
            #             if k < 0.0001:
            #                 k = 0.0001
            #             scf = 1/2 * 80.4**2/k
            #             nu.SetPx(nu.Px()*scf)
            #             nu.SetPy(nu.Py()*scf)
            #             nu.SetE(nu.P())


                            
            #         lamda = (80.4)**2/2 + Ele.Px()*nu.Px() + Ele.Py()*nu.Py()
            #         #lamda = (80.4)**2/2 + mu.Px()*nu.Px() + mu.Py()*nu.Py()
                        
            #         #if Ele.Pt() > 0 or Ele.Pt() < 0:
            #         discr = (lamda*Ele.Pz())**2/(Ele.Pt())**4 - ((Ele.E()*nu.Pt())**2 - lamda**2)/(Ele.Pt())**2 
            #         #discr = (lamda*mu.Pz())**2/(mu.Pt())**4 - ((mu.E()*nu.Pt())**2 - lamda**2)/(mu.Pt())**2

            #         if WmT > 80.4 or discr < 0:

            #             s = (lamda*Ele.Pz())/(Ele.Pt())**2
            #             #s = (lamda*mu.Pz())/(mu.Pt())**2
            #             nu.SetPz(s)
            #             nu.SetE(nu.P())
            #             #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            #         else:
            #             s1 = (lamda*Ele.Pz())/(Ele.Pt())**2 + sqrt(discr)
            #             s2 = (lamda*Ele.Pz())/(Ele.Pt())**2 - sqrt(discr)
            #             #s1 = (lamda*mu.Pz())/(mu.Pt())**2 + sqrt(discr)
            #             #s2 = (lamda*mu.Pz())/(mu.Pt())**2 - sqrt(discr)

            #             nu.SetPz(s1)
            #             nu.SetE(nu.P())
            #             #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            #             #Wm1 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            #             #Wm1 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #             Wm1 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
            #             #Wm1 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

            #             nu.SetPz(s2)
            #             nu.SetE(nu.P())
            #             #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            #             #Wm2 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            #             #Wm2 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #             Wm2 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
            #             #Wm2 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

            #             if abs(Wm2 - 80.4) > abs(Wm1 - 80.4):    
            #                 nu.SetPz(s1)
            #                 nu.SetE(nu.P())
            #                 #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                            
            #             else:
            #                 nu.SetPz(s2)
            #                 nu.SetE(nu.P())
            #                 #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                                
            
            #     #if ele_cut_count >= 1 and jet_cut_count >= 2:
            #         #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #         #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())

            #         WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            #         WmT = sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)

            #         # WmI2 = (Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2
            #         # WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

            #         # print('WmI^2=',WmI2,',WmT^2=',WmT2)


            #         #WmI = sqrt((mu.E()+nu.E())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2-(mu.Pz()+nu.Pz())**2)
            #         #WmT = sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

            #         print(WmT)
            #         WmTHist.Fill(WmT)
            #         Pw = Ele + nu
            #         #Pw = mu + nu
            #         WmtHist.Fill(Pw.Mt())
            #         WmIHist.Fill(WmI)
            #         #Mw2 = Pw.Dot(Pw)
            #         #print(sqrt(Mw2))

            #     #if (Jet_pt[jet_sub_maxpt_idx] > 25 and Jet_eta[jet_sub_maxpt_idx] < 2.4) and jet_cut_count >= 2 and ele_cut_count >= 1:
            #         Ptop1 = Pw + JetP1
            #         #Mtop1 = sqrt(Ptop1.Dot(Ptop1))
            #         Mtop1 = Ptop1.M()
            #         Ptop2 = Pw + JetP2
            #         #Mtop2 = sqrt(Ptop2.Dot(Ptop2))
            #         Mtop2 = Ptop2.M()
            #         #Ptop12 = Ptop1 + Ptop2
            #         if abs(Mtop1-172.5) > abs(Mtop2-172.5):
            #             Mtop = Mtop2
            #         else:
            #             Mtop = Mtop1
            #         #if (Ptop1.Pt() > 250 and Ptop12.Pt() > 350) or (Ptop2.Pt() > 250 and Ptop12.Pt() > 350):
            #             #if Mtop > 120 and Mtop < 220:
            #         TopHist.Fill(Mtop)

            #         # DrMuonJet1 = ROOT.Math.VectorUtil.DeltaR(mu, JetP1)
            #         # DrMuonJet1Hist.Fill(DrMuonJet1)

            #         # DrMuonJet2 = ROOT.Math.VectorUtil.DeltaR(mu, JetP2)
            #         # DrMuonJet2Hist.Fill(DrMuonJet2)

            #         DrEleJet1 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP1)
            #         DrEleJet1Hist.Fill(DrEleJet1)

            #         DrEleJet2 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP2)
            #         DrEleJet2Hist.Fill(DrEleJet2)
            # elif ele_cut_count >= 1 and mu_cut_count < 1:
            #     Ele = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
            #     nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), sqrt(nu_pt**2 + nu.Pz()**2))
            #     #mu = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
                
            #     #nu = ROOT.Math.PxPyPzMVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), 0)
            #     #nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), nu.E())
            #     #nu = ROOT.Math.PtEtaPhiEVector(nu_pt, nu.Eta(), nu_phi, nu.E())
            #     #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            
            # #if jet_cut_count >= 2 and ele_cut_count >= 1:
            #     JetP1 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_sub_maxpt_idx], Jet_eta[jet_sub_maxpt_idx], Jet_phi[jet_sub_maxpt_idx], Jet_mass[jet_sub_maxpt_idx])
            #     JetP2 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_maxpt_idx], Jet_eta[jet_maxpt_idx], Jet_phi[jet_maxpt_idx], Jet_mass[jet_maxpt_idx])   
                
            #     WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

            #     print('WmT^2=',WmT2)

            #     #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #     WmT =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
            #     #WmT =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)
            #     #WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)

            #     # 80.4 GeV is the W-boson pole mass
            #     if WmT > 80.4:
            #         k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
            #         #k = nu.Et() * mu.Pt() - nu.Px() * mu.Px() - nu.Py() * mu.Py()
            #         if k < 0.0001:
            #             k = 0.0001
            #         scf = 1/2 * 80.4**2/k
            #         nu.SetPx(nu.Px()*scf)
            #         nu.SetPy(nu.Py()*scf)
            #         nu.SetE(nu.P())


                        
            #     lamda = (80.4)**2/2 + Ele.Px()*nu.Px() + Ele.Py()*nu.Py()
            #     #lamda = (80.4)**2/2 + mu.Px()*nu.Px() + mu.Py()*nu.Py()
                    
            #     #if Ele.Pt() > 0 or Ele.Pt() < 0:
            #     discr = (lamda*Ele.Pz())**2/(Ele.Pt())**4 - ((Ele.E()*nu.Pt())**2 - lamda**2)/(Ele.Pt())**2 
            #     #discr = (lamda*mu.Pz())**2/(mu.Pt())**4 - ((mu.E()*nu.Pt())**2 - lamda**2)/(mu.Pt())**2

            #     if WmT > 80.4 or discr < 0:
            #         # k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
            #         # if k < 0.0001:
            #         #     k = 0.0001
            #         # scf = 1/2 * 80.4**2/k
            #         # nu.SetPx(nu.Px()*scf)
            #         # nu.SetPy(nu.Py()*scf)
            #         #nu.SetE(nu.P())

            #         s = (lamda*Ele.Pz())/(Ele.Pt())**2
            #         #s = (lamda*mu.Pz())/(mu.Pt())**2
            #         nu.SetPz(s)
            #         nu.SetE(nu.P())
            #         #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            #     else:
            #         s1 = (lamda*Ele.Pz())/(Ele.Pt())**2 + sqrt(discr)
            #         s2 = (lamda*Ele.Pz())/(Ele.Pt())**2 - sqrt(discr)
            #         #s1 = (lamda*mu.Pz())/(mu.Pt())**2 + sqrt(discr)
            #         #s2 = (lamda*mu.Pz())/(mu.Pt())**2 - sqrt(discr)

            #         nu.SetPz(s1)
            #         nu.SetE(nu.P())
            #         #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            #         #Wm1 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            #         #Wm1 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #         Wm1 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
            #         #Wm1 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

            #         nu.SetPz(s2)
            #         nu.SetE(nu.P())
            #         #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            #         #Wm2 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            #         #Wm2 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #         Wm2 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
            #         #Wm2 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

            #         if abs(Wm2 - 80.4) > abs(Wm1 - 80.4):    
            #             nu.SetPz(s1)
            #             nu.SetE(nu.P())
            #             #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                        
            #         else:
            #             nu.SetPz(s2)
            #             nu.SetE(nu.P())
            #             #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                            
        
            # #if ele_cut_count >= 1 and jet_cut_count >= 2:
            #     #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
            #     #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())

            #     WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
            #     WmT = sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)

            #     # WmI2 = (Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2
            #     # WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

            #     # print('WmI^2=',WmI2,',WmT^2=',WmT2)


            #     #WmI = sqrt((mu.E()+nu.E())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2-(mu.Pz()+nu.Pz())**2)
            #     #WmT = sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

            #     print(WmT)
            #     WmTHist.Fill(WmT)
            #     Pw = Ele + nu
            #     #Pw = mu + nu
            #     WmtHist.Fill(Pw.Mt())
            #     WmIHist.Fill(WmI)
            #     #Mw2 = Pw.Dot(Pw)
            #     #print(sqrt(Mw2))

            # #if (Jet_pt[jet_sub_maxpt_idx] > 25 and Jet_eta[jet_sub_maxpt_idx] < 2.4) and jet_cut_count >= 2 and ele_cut_count >= 1:
            #     Ptop1 = Pw + JetP1
            #     #Mtop1 = sqrt(Ptop1.Dot(Ptop1))
            #     Mtop1 = Ptop1.M()
            #     Ptop2 = Pw + JetP2
            #     #Mtop2 = sqrt(Ptop2.Dot(Ptop2))
            #     Mtop2 = Ptop2.M()
            #     #Ptop12 = Ptop1 + Ptop2
            #     if abs(Mtop1-172.5) > abs(Mtop2-172.5):
            #         Mtop = Mtop2
            #     else:
            #         Mtop = Mtop1
            #     #if (Ptop1.Pt() > 250 and Ptop12.Pt() > 350) or (Ptop2.Pt() > 250 and Ptop12.Pt() > 350):
            #         #if Mtop > 120 and Mtop < 220:
            #     TopHist.Fill(Mtop)

            #     # DrMuonJet1 = ROOT.Math.VectorUtil.DeltaR(mu, JetP1)
            #     # DrMuonJet1Hist.Fill(DrMuonJet1)

            #     # DrMuonJet2 = ROOT.Math.VectorUtil.DeltaR(mu, JetP2)
            #     # DrMuonJet2Hist.Fill(DrMuonJet2)

            #     DrEleJet1 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP1)
            #     DrEleJet1Hist.Fill(DrEleJet1)

            #     DrEleJet2 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP2)
            #     DrEleJet2Hist.Fill(DrEleJet2)

    file_count += 1
    events_count += entries
    print('events count =', events_count)
    if events_count >= 10000:
        break
    

c1 = ROOT.TCanvas( 'c1', 'W_transverse_mass', 1000, 875 )
c1.cd()
WmtHist.GetYaxis().SetTitle('Number of events')
WmtHist.GetXaxis().SetTitle('W transverse mass Mt() [GeV]')
WmtHist.SetStats(0)
WmtHist.SetLineColor(kRed)
WmtHist.SetLineWidth(2)
WmtHist.Draw()
c1.SaveAs('W_trans_mass.pdf')

c2 = ROOT.TCanvas( 'c2', 'W_mass', 1000, 875 )
c2.cd()
WmTHist.GetYaxis().SetTitle('Number of events')
WmTHist.GetXaxis().SetTitle('W reco mT [GeV]')
WmTHist.SetStats(0)
WmTHist.SetLineColor(kRed)
WmTHist.SetLineWidth(2)
WmTHist.Draw()
c2.SaveAs('W_reco_mT.pdf')

c3 = ROOT.TCanvas( 'c3', 'Top_mass', 1000, 875 )
c3.cd()
TopHist.GetYaxis().SetTitle('Number of events')
TopHist.GetXaxis().SetTitle('Top reco mass [GeV]')
TopHist.SetStats(0)
TopHist.SetLineColor(kRed)
TopHist.SetLineWidth(2)
TopHist.Draw()
c3.SaveAs('Top_reco_mass.pdf')

c4 = ROOT.TCanvas( 'c4', 'W_inv_mass', 1000, 875 )
c4.cd()
WmIHist.GetYaxis().SetTitle('Number of events')
WmIHist.GetXaxis().SetTitle('W mI [GeV]')
WmIHist.SetStats(0)
WmIHist.SetLineColor(kRed)
WmIHist.SetLineWidth(2)
WmIHist.Draw()
c4.SaveAs('W_mI.pdf')

c5 = ROOT.TCanvas( 'c5', 'DrEleJet1', 1000, 875 )
c5.cd()
DrEleJet1Hist.GetYaxis().SetTitle('Number of events')
DrEleJet1Hist.GetXaxis().SetTitle('dR(Electron,Jet1) [GeV]')
DrEleJet1Hist.SetStats(0)
DrEleJet1Hist.SetLineColor(kRed)
DrEleJet1Hist.SetLineWidth(2)
DrEleJet1Hist.Draw()
c5.SaveAs('DrEleJet1.pdf')

c10 = ROOT.TCanvas( 'c10', 'DrMuJet1', 1000, 875 )
c10.cd()
DrMuJet1Hist.GetYaxis().SetTitle('Number of events')
DrMuJet1Hist.GetXaxis().SetTitle('dR(Muon,Jet1) [GeV]')
DrMuJet1Hist.SetStats(0)
DrMuJet1Hist.SetLineColor(kRed)
DrMuJet1Hist.SetLineWidth(2)
DrMuJet1Hist.Draw()
c10.SaveAs('DrMuJet1.pdf')

c12 = ROOT.TCanvas( 'c12', 'DrEleJet1Top', 1000, 875 )
c12.cd()
DrEleJet1TopHist.GetYaxis().SetTitle('Number of events')
DrEleJet1TopHist.GetXaxis().SetTitle('dR(Electron,Jet1Top) [GeV]')
DrEleJet1TopHist.SetStats(0)
DrEleJet1TopHist.SetLineColor(kRed)
DrEleJet1TopHist.SetLineWidth(2)
DrEleJet1TopHist.Draw()
c12.SaveAs('DrEleJet1Top.pdf')

c13 = ROOT.TCanvas( 'c13', 'DrMuJet1Top', 1000, 875 )
c13.cd()
DrMuJet1TopHist.GetYaxis().SetTitle('Number of events')
DrMuJet1TopHist.GetXaxis().SetTitle('dR(Muon,Jet1Top) [GeV]')
DrMuJet1TopHist.SetStats(0)
DrMuJet1TopHist.SetLineColor(kRed)
DrMuJet1TopHist.SetLineWidth(2)
DrMuJet1TopHist.Draw()
c13.SaveAs('DrMuJet1Top.pdf')

c11 = ROOT.TCanvas( 'c11', 'DrMuJet2', 1000, 875 )
c11.cd()
DrMuJet2Hist.GetYaxis().SetTitle('Number of events')
DrMuJet2Hist.GetXaxis().SetTitle('dR(Muon,Jet2) [GeV]')
DrMuJet2Hist.SetStats(0)
DrMuJet2Hist.SetLineColor(kRed)
DrMuJet2Hist.SetLineWidth(2)
DrMuJet2Hist.Draw()
c11.SaveAs('DrMuJet2.pdf')

c6 = ROOT.TCanvas( 'c6', 'DrMuonJet2', 1000, 875 )
c6.cd()
DrEleJet2Hist.GetYaxis().SetTitle('Number of events')
DrEleJet2Hist.GetXaxis().SetTitle('dR(Ele,Jet2) [GeV]')
DrEleJet2Hist.SetStats(0)
DrEleJet2Hist.SetLineColor(kRed)
DrEleJet2Hist.SetLineWidth(2)
DrEleJet2Hist.Draw()
c6.SaveAs('DrEleJet2.pdf')

c14 = ROOT.TCanvas( 'c14', 'DrMuJet2Top', 1000, 875 )
c14.cd()
DrMuJet2TopHist.GetYaxis().SetTitle('Number of events')
DrMuJet2TopHist.GetXaxis().SetTitle('dR(Muon,Jet2Top) [GeV]')
DrMuJet2TopHist.SetStats(0)
DrMuJet2TopHist.SetLineColor(kRed)
DrMuJet2TopHist.SetLineWidth(2)
DrMuJet2TopHist.Draw()
c14.SaveAs('DrMuJet2Top.pdf')

c15 = ROOT.TCanvas( 'c15', 'DrMuonJet2Top', 1000, 875 )
c15.cd()
DrEleJet2TopHist.GetYaxis().SetTitle('Number of events')
DrEleJet2TopHist.GetXaxis().SetTitle('dR(Ele,Jet2Top) [GeV]')
DrEleJet2TopHist.SetStats(0)
DrEleJet2TopHist.SetLineColor(kRed)
DrEleJet2TopHist.SetLineWidth(2)
DrEleJet2TopHist.Draw()
c15.SaveAs('DrEleJet2Top.pdf')

c7 = ROOT.TCanvas( 'c7', 'MaxPtJetHist', 1000, 875 )
c7.cd()
MaxPtJetHist.GetYaxis().SetTitle('Number of events')
MaxPtJetHist.GetXaxis().SetTitle('Leading Jet Pt [GeV]')
MaxPtJetHist.SetStats(0)
MaxPtJetHist.SetLineColor(kRed)
MaxPtJetHist.SetLineWidth(2)
MaxPtJetHist.Draw()
c7.SaveAs('MaxPtJet.pdf')

c8 = ROOT.TCanvas( 'c8', 'SubMaxPtJetHist', 1000, 875 )
c8.cd()
SubMaxPtJetHist.GetYaxis().SetTitle('Number of events')
SubMaxPtJetHist.GetXaxis().SetTitle('Subleading Jet Pt [GeV]')
SubMaxPtJetHist.SetStats(0)
SubMaxPtJetHist.SetLineColor(kRed)
SubMaxPtJetHist.SetLineWidth(2)
SubMaxPtJetHist.Draw()
c8.SaveAs('SubMaxPtJet.pdf')

c9 = ROOT.TCanvas( 'c9', 'LepNu_transverse_mass', 1000, 875 )
c9.cd()
LepNuMtHist.GetYaxis().SetTitle('Number of events')
LepNuMtHist.GetXaxis().SetTitle('Lep+Nu Mt() [GeV]')
LepNuMtHist.SetStats(0)
LepNuMtHist.SetLineColor(kRed)
LepNuMtHist.SetLineWidth(2)
LepNuMtHist.Draw()
c9.SaveAs('LepNu_trans_mass.pdf')



    



    


