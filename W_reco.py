import ROOT 
from ROOT import *
from numpy import *


#rootFiles = loadtxt('Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt', dtype=str, delimiter = "root:")
rootFiles = loadtxt('/nfs/dust/cms/user/sobhatta/work/TopTagPol/TreeMaker/CMSSW_10_5_0/src/sourceFiles/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/Wprimetotb_M2000W20_LH_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt', dtype=str, delimiter = " ")

TopHist = ROOT.TH1D('TopHist', '', 50, 0, 500)
WmTHist = ROOT.TH1D('WmTHist', '', 100, 0, 200)
WmIHist = ROOT.TH1D('WmIHist', '', 200, 80.395, 80.405)
WmtHist = ROOT.TH1D('WmtHist', '', 100, 0, 1400)
DrEleJet1Hist = ROOT.TH1D('DrEleJet1Hist', '', 100, 0, 10)
DrEleJet2Hist = ROOT.TH1D('DrEleJet2Hist', '', 100, 0, 10)
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
        Muon_phi = getattr(tree,'Muon_phi')
        Muon_mass = getattr(tree,'Muon_mass')
        Jet_pt = getattr(tree,'Jet_pt')
        Jet_eta = getattr(tree,'Jet_eta')
        Jet_phi = getattr(tree,'Jet_phi')
        Jet_mass = getattr(tree,'Jet_mass')
        JetBtag = getattr(tree, 'Jet_btagDeepFlavB')
        Ele_noIso = getattr(tree,'Electron_mvaFall17V2noIso_WP90')
        Muon_tightId = getattr(tree,'Muon_tightId')
        Ele_miniIso = getattr(tree,'Electron_miniPFRelIso_all')
        Muon_miniIso = getattr(tree,'Muon_miniPFRelIso_all')

        ele_cut_count = 0
        ele_maxpt_idx = 0
        for iEle in range(0,len(Ele_pt)):
            #print(Ele_pt[iEle])
            if Ele_pt[iEle] < 100 or abs(Ele_eta[iEle]) > 2.5 or Ele_noIso[iEle] != 1 or Ele_miniIso[iEle] > 0.1:
                continue
            ele_cut_count += 1
            if Ele_pt[iEle] > Ele_pt[ele_maxpt_idx]:
                ele_maxpt_idx = iEle        
        
        mu_cut_count = 0
        mu_maxpt_idx = 0 
        for iMu in range(0, len(Muon_pt)):
            if Muon_pt[iMu] < 100 or abs(Muon_eta[iMu]) > 2.1 or Muon_tightId[iMu] != 1 or Muon_miniIso[iMu] > 0.1:
                continue
            mu_cut_count += 1
            if Muon_pt[iMu] > Muon_pt[mu_maxpt_idx]:
                mu_maxpt_idx = iMu

        if mu_cut_count < 1: 
            continue

        jet_cut_count = 0
        jet_Btag_count = 0
        jet_maxpt_idx = 0
        jet_sub_maxpt_idx = 0
        for iJet in range(0, len(Jet_pt)):
            if Jet_pt[iJet] < 100 or abs(Jet_eta[iJet]) > 2.4:
                continue
            jet_cut_count += 1
            # if JetBtag[iJet] < 0.2270:
            #     continue
            # jet_Btag_count += 1
            if Jet_pt[iJet] > Jet_pt[jet_maxpt_idx]:
                jet_maxpt_idx = iJet
            if abs(Jet_pt[jet_maxpt_idx] - Jet_pt[iJet]) < abs(Jet_pt[jet_maxpt_idx] - Jet_pt[jet_sub_maxpt_idx]):
                jet_sub_maxpt_idx = iJet
        
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

        if (ele_cut_count >= 1 or mu_cut_count >=1) and jet_cut_count >= 2:
            if ele_cut_count >= 1 and mu_cut_count >= 1:
                if Ele_pt[ele_maxpt_idx] > Muon_pt[mu_maxpt_idx]:
                    Ele = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
                    nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), sqrt(nu_pt**2 + nu.Pz()**2))
                    #mu = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
                    
                    JetP1 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_sub_maxpt_idx], Jet_eta[jet_sub_maxpt_idx], Jet_phi[jet_sub_maxpt_idx], Jet_mass[jet_sub_maxpt_idx])
                    JetP2 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_maxpt_idx], Jet_eta[jet_maxpt_idx], Jet_phi[jet_maxpt_idx], Jet_mass[jet_maxpt_idx])   
                    
                    WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

                    print('WmT^2=',WmT2)

                    #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                    WmT =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
                    #WmT =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)
                    #WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)

                    # 80.4 GeV is the W-boson pole mass
                    if WmT > 80.4:
                        k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
                        #k = nu.Et() * mu.Pt() - nu.Px() * mu.Px() - nu.Py() * mu.Py()
                        if k < 0.0001:
                            k = 0.0001
                        scf = 1/2 * 80.4**2/k
                        nu.SetPx(nu.Px()*scf)
                        nu.SetPy(nu.Py()*scf)
                        nu.SetE(nu.P())


                            
                    lamda = (80.4)**2/2 + Ele.Px()*nu.Px() + Ele.Py()*nu.Py()
                    #lamda = (80.4)**2/2 + mu.Px()*nu.Px() + mu.Py()*nu.Py()
                        
                    #if Ele.Pt() > 0 or Ele.Pt() < 0:
                    discr = (lamda*Ele.Pz())**2/(Ele.Pt())**4 - ((Ele.E()*nu.Pt())**2 - lamda**2)/(Ele.Pt())**2 
                    #discr = (lamda*mu.Pz())**2/(mu.Pt())**4 - ((mu.E()*nu.Pt())**2 - lamda**2)/(mu.Pt())**2

                    if WmT > 80.4 or discr < 0:
                        # k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
                        # if k < 0.0001:
                        #     k = 0.0001
                        # scf = 1/2 * 80.4**2/k
                        # nu.SetPx(nu.Px()*scf)
                        # nu.SetPy(nu.Py()*scf)
                        #nu.SetE(nu.P())

                        s = (lamda*Ele.Pz())/(Ele.Pt())**2
                        #s = (lamda*mu.Pz())/(mu.Pt())**2
                        nu.SetPz(s)
                        nu.SetE(nu.P())
                        #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                    else:
                        s1 = (lamda*Ele.Pz())/(Ele.Pt())**2 + sqrt(discr)
                        s2 = (lamda*Ele.Pz())/(Ele.Pt())**2 - sqrt(discr)
                        #s1 = (lamda*mu.Pz())/(mu.Pt())**2 + sqrt(discr)
                        #s2 = (lamda*mu.Pz())/(mu.Pt())**2 - sqrt(discr)

                        nu.SetPz(s1)
                        nu.SetE(nu.P())
                        #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                        #Wm1 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
                        #Wm1 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                        Wm1 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
                        #Wm1 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

                        nu.SetPz(s2)
                        nu.SetE(nu.P())
                        #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                        #Wm2 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
                        #Wm2 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                        Wm2 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
                        #Wm2 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

                        if abs(Wm2 - 80.4) > abs(Wm1 - 80.4):    
                            nu.SetPz(s1)
                            nu.SetE(nu.P())
                            #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                            
                        else:
                            nu.SetPz(s2)
                            nu.SetE(nu.P())
                            #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                                
            
                #if ele_cut_count >= 1 and jet_cut_count >= 2:
                    #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                    #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())

                    WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
                    WmT = sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)

                    # WmI2 = (Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2
                    # WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

                    # print('WmI^2=',WmI2,',WmT^2=',WmT2)


                    #WmI = sqrt((mu.E()+nu.E())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2-(mu.Pz()+nu.Pz())**2)
                    #WmT = sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

                    print(WmT)
                    WmTHist.Fill(WmT)
                    Pw = Ele + nu
                    #Pw = mu + nu
                    WmtHist.Fill(Pw.Mt())
                    WmIHist.Fill(WmI)
                    #Mw2 = Pw.Dot(Pw)
                    #print(sqrt(Mw2))

                #if (Jet_pt[jet_sub_maxpt_idx] > 25 and Jet_eta[jet_sub_maxpt_idx] < 2.4) and jet_cut_count >= 2 and ele_cut_count >= 1:
                    Ptop1 = Pw + JetP1
                    #Mtop1 = sqrt(Ptop1.Dot(Ptop1))
                    Mtop1 = Ptop1.M()
                    Ptop2 = Pw + JetP2
                    #Mtop2 = sqrt(Ptop2.Dot(Ptop2))
                    Mtop2 = Ptop2.M()
                    #Ptop12 = Ptop1 + Ptop2
                    if abs(Mtop1-172.5) > abs(Mtop2-172.5):
                        Mtop = Mtop2
                    else:
                        Mtop = Mtop1
                    #if (Ptop1.Pt() > 250 and Ptop12.Pt() > 350) or (Ptop2.Pt() > 250 and Ptop12.Pt() > 350):
                        #if Mtop > 120 and Mtop < 220:
                    TopHist.Fill(Mtop)

                    # DrMuonJet1 = ROOT.Math.VectorUtil.DeltaR(mu, JetP1)
                    # DrMuonJet1Hist.Fill(DrMuonJet1)

                    # DrMuonJet2 = ROOT.Math.VectorUtil.DeltaR(mu, JetP2)
                    # DrMuonJet2Hist.Fill(DrMuonJet2)

                    DrEleJet1 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP1)
                    DrEleJet1Hist.Fill(DrEleJet1)

                    DrEleJet2 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP2)
                    DrEleJet2Hist.Fill(DrEleJet2)
            elif ele_cut_count >= 1 and mu_cut_count < 1:
                Ele = ROOT.Math.PtEtaPhiMVector(Ele_pt[ele_maxpt_idx], Ele_eta[ele_maxpt_idx], Ele_phi[ele_maxpt_idx], Ele_mass[ele_maxpt_idx])
                nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), sqrt(nu_pt**2 + nu.Pz()**2))
                #mu = ROOT.Math.PtEtaPhiMVector(Muon_pt[mu_maxpt_idx], Muon_eta[mu_maxpt_idx], Muon_phi[mu_maxpt_idx], Muon_mass[mu_maxpt_idx])
                
                #nu = ROOT.Math.PxPyPzMVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), 0)
                #nu = ROOT.Math.PxPyPzEVector(nu_pt*cos(nu_phi), nu_pt*sin(nu_phi), nu.Pz(), nu.E())
                #nu = ROOT.Math.PtEtaPhiEVector(nu_pt, nu.Eta(), nu_phi, nu.E())
                #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
            
            #if jet_cut_count >= 2 and ele_cut_count >= 1:
                JetP1 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_sub_maxpt_idx], Jet_eta[jet_sub_maxpt_idx], Jet_phi[jet_sub_maxpt_idx], Jet_mass[jet_sub_maxpt_idx])
                JetP2 = ROOT.Math.PtEtaPhiMVector(Jet_pt[jet_maxpt_idx], Jet_eta[jet_maxpt_idx], Jet_phi[jet_maxpt_idx], Jet_mass[jet_maxpt_idx])   
                
                WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

                print('WmT^2=',WmT2)

                #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                WmT =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
                #WmT =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)
                #WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)

                # 80.4 GeV is the W-boson pole mass
                if WmT > 80.4:
                    k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
                    #k = nu.Et() * mu.Pt() - nu.Px() * mu.Px() - nu.Py() * mu.Py()
                    if k < 0.0001:
                        k = 0.0001
                    scf = 1/2 * 80.4**2/k
                    nu.SetPx(nu.Px()*scf)
                    nu.SetPy(nu.Py()*scf)
                    nu.SetE(nu.P())


                        
                lamda = (80.4)**2/2 + Ele.Px()*nu.Px() + Ele.Py()*nu.Py()
                #lamda = (80.4)**2/2 + mu.Px()*nu.Px() + mu.Py()*nu.Py()
                    
                #if Ele.Pt() > 0 or Ele.Pt() < 0:
                discr = (lamda*Ele.Pz())**2/(Ele.Pt())**4 - ((Ele.E()*nu.Pt())**2 - lamda**2)/(Ele.Pt())**2 
                #discr = (lamda*mu.Pz())**2/(mu.Pt())**4 - ((mu.E()*nu.Pt())**2 - lamda**2)/(mu.Pt())**2

                if WmT > 80.4 or discr < 0:
                    # k = nu.Et() * Ele.Pt() - nu.Px() * Ele.Px() - nu.Py() * Ele.Py()
                    # if k < 0.0001:
                    #     k = 0.0001
                    # scf = 1/2 * 80.4**2/k
                    # nu.SetPx(nu.Px()*scf)
                    # nu.SetPy(nu.Py()*scf)
                    #nu.SetE(nu.P())

                    s = (lamda*Ele.Pz())/(Ele.Pt())**2
                    #s = (lamda*mu.Pz())/(mu.Pt())**2
                    nu.SetPz(s)
                    nu.SetE(nu.P())
                    #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                else:
                    s1 = (lamda*Ele.Pz())/(Ele.Pt())**2 + sqrt(discr)
                    s2 = (lamda*Ele.Pz())/(Ele.Pt())**2 - sqrt(discr)
                    #s1 = (lamda*mu.Pz())/(mu.Pt())**2 + sqrt(discr)
                    #s2 = (lamda*mu.Pz())/(mu.Pt())**2 - sqrt(discr)

                    nu.SetPz(s1)
                    nu.SetE(nu.P())
                    #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                    #Wm1 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
                    #Wm1 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                    Wm1 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
                    #Wm1 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

                    nu.SetPz(s2)
                    nu.SetE(nu.P())
                    #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                    #Wm2 = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
                    #Wm2 =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                    Wm2 =  sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)
                    #Wm2 =  sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

                    if abs(Wm2 - 80.4) > abs(Wm1 - 80.4):    
                        nu.SetPz(s1)
                        nu.SetE(nu.P())
                        #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                        
                    else:
                        nu.SetPz(s2)
                        nu.SetE(nu.P())
                        #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())
                            
        
            #if ele_cut_count >= 1 and jet_cut_count >= 2:
                #Wm =  sqrt(2*(Ele.E()*nu.E()-Ele.Px()*nu.Px()-Ele.Py()*nu.Py()-Ele.Pz()*nu.Pz()))
                #nu = ROOT.Math.PxPyPzEVector(nu.Px(), nu.Py(), nu.Pz(), nu.E())

                WmI = sqrt((Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2)
                WmT = sqrt((Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2)

                # WmI2 = (Ele.E()+nu.E())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2-(Ele.Pz()+nu.Pz())**2
                # WmT2 = (Ele.Pt()+nu.Et())**2-(Ele.Px()+nu.Px())**2-(Ele.Py()+nu.Py())**2

                # print('WmI^2=',WmI2,',WmT^2=',WmT2)


                #WmI = sqrt((mu.E()+nu.E())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2-(mu.Pz()+nu.Pz())**2)
                #WmT = sqrt((mu.Pt()+nu.Et())**2-(mu.Px()+nu.Px())**2-(mu.Py()+nu.Py())**2)

                print(WmT)
                WmTHist.Fill(WmT)
                Pw = Ele + nu
                #Pw = mu + nu
                WmtHist.Fill(Pw.Mt())
                WmIHist.Fill(WmI)
                #Mw2 = Pw.Dot(Pw)
                #print(sqrt(Mw2))

            #if (Jet_pt[jet_sub_maxpt_idx] > 25 and Jet_eta[jet_sub_maxpt_idx] < 2.4) and jet_cut_count >= 2 and ele_cut_count >= 1:
                Ptop1 = Pw + JetP1
                #Mtop1 = sqrt(Ptop1.Dot(Ptop1))
                Mtop1 = Ptop1.M()
                Ptop2 = Pw + JetP2
                #Mtop2 = sqrt(Ptop2.Dot(Ptop2))
                Mtop2 = Ptop2.M()
                #Ptop12 = Ptop1 + Ptop2
                if abs(Mtop1-172.5) > abs(Mtop2-172.5):
                    Mtop = Mtop2
                else:
                    Mtop = Mtop1
                #if (Ptop1.Pt() > 250 and Ptop12.Pt() > 350) or (Ptop2.Pt() > 250 and Ptop12.Pt() > 350):
                    #if Mtop > 120 and Mtop < 220:
                TopHist.Fill(Mtop)

                # DrMuonJet1 = ROOT.Math.VectorUtil.DeltaR(mu, JetP1)
                # DrMuonJet1Hist.Fill(DrMuonJet1)

                # DrMuonJet2 = ROOT.Math.VectorUtil.DeltaR(mu, JetP2)
                # DrMuonJet2Hist.Fill(DrMuonJet2)

                DrEleJet1 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP1)
                DrEleJet1Hist.Fill(DrEleJet1)

                DrEleJet2 = ROOT.Math.VectorUtil.DeltaR(Ele, JetP2)
                DrEleJet2Hist.Fill(DrEleJet2)
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
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c1.SaveAs('W_trans_mass.pdf')

c2 = ROOT.TCanvas( 'c2', 'W_mass', 1000, 875 )
c2.cd()
WmTHist.GetYaxis().SetTitle('Number of events')
WmTHist.GetXaxis().SetTitle('W reco mT [GeV]')
WmTHist.SetStats(0)
WmTHist.SetLineColor(kRed)
WmTHist.SetLineWidth(2)
WmTHist.Draw()
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c2.SaveAs('W_reco_mT.pdf')

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

c4 = ROOT.TCanvas( 'c4', 'W_inv_mass', 1000, 875 )
c4.cd()
WmIHist.GetYaxis().SetTitle('Number of events')
WmIHist.GetXaxis().SetTitle('W mI [GeV]')
WmIHist.SetStats(0)
WmIHist.SetLineColor(kRed)
WmIHist.SetLineWidth(2)
WmIHist.Draw()
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c4.SaveAs('W_mI.pdf')

c5 = ROOT.TCanvas( 'c5', 'DrEleJet1', 1000, 875 )
c5.cd()
DrEleJet1Hist.GetYaxis().SetTitle('Number of events')
DrEleJet1Hist.GetXaxis().SetTitle('dR(Electron,Jet1) [GeV]')
DrEleJet1Hist.SetStats(0)
DrEleJet1Hist.SetLineColor(kRed)
DrEleJet1Hist.SetLineWidth(2)
DrEleJet1Hist.Draw()
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c5.SaveAs('DrEleJet1.pdf')

c6 = ROOT.TCanvas( 'c6', 'DrMuonJet2', 1000, 875 )
c6.cd()
DrEleJet2Hist.GetYaxis().SetTitle('Number of events')
DrEleJet2Hist.GetXaxis().SetTitle('dR(Ele,Jet2) [GeV]')
DrEleJet2Hist.SetStats(0)
DrEleJet2Hist.SetLineColor(kRed)
DrEleJet2Hist.SetLineWidth(2)
DrEleJet2Hist.Draw()
#legend = ROOT.TLegend(0.15 ,0.7 ,0.45 ,0.8)
#legend.AddEntry(JetHist, "p_{T} > 100 GeV and #eta < 2.5")
#legend.SetLineWidth (0)
#legend.Draw("same")
c6.SaveAs('DrEleJet2.pdf')




    



    


