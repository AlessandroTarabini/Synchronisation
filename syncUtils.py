import math, os, sys#, ctypes
from ROOT import *
#from import c_int, c_float
from ctypes import *

gConstant_g4 = TFile("gConstant_HZZ2e2mu_g4.root")
spline_g4 = gConstant_g4.Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4")

gConstant_g2 = TFile("gConstant_HZZ2e2mu_g2.root")
spline_g2 = gConstant_g2.Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2")

gConstant_L1 = TFile("gConstant_HZZ2e2mu_L1.root")
spline_L1 = gConstant_L1.Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1")

gConstant_L1Zgs = TFile("gConstant_HZZ2e2mu_L1Zgs.root")
spline_L1Zgs = gConstant_L1Zgs.Get("sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs")

def sign(x):
    if x > 0:
        return 1.
    elif x < 0:
        return -1.
    elif x == 0:
        return 0.
    else:
        return x


class Candidate:

    def __init__(self, treeEntry) :
        self.lib = CDLL('libZZAnalysisAnalysisStep.so')

        isMC = False
        isReco = False
        if treeEntry.GetBranch('genHEPMCweight') :
            isMC = True
        if treeEntry.GetBranch('ZZMass'):
            isReco = True

        self.run                             = treeEntry.RunNumber
        self.lumi                            = treeEntry.LumiNumber
        self.event                           = treeEntry.EventNumber
        self.passedFiducial                  = treeEntry.passedFiducialSelection_bbf

        self.GENmass4l                       = treeEntry.GENmass4l
        self.GENpT4l                         = treeEntry.GENpT4l
        self.GENeta4l                        = treeEntry.GENeta4l
        self.GENphi4l                        = treeEntry.GENphi4l

        self.GENlep_id                       = treeEntry.GENlep_id
        self.GENlep_pt                       = list(treeEntry.GENlep_pt)
        self.GENlep_eta                      = treeEntry.GENlep_eta
        self.GENlep_phi                      = treeEntry.GENlep_phi
        self.GENlep_mass                     = treeEntry.GENlep_mass
        self.GENlep_Hindex                   = list(treeEntry.GENlep_Hindex)

        self.GENlep_pt_ordered = sorted(self.GENlep_pt, reverse = True)
        self.GENlep1_pt = self.GENlep_pt_ordered[0]
        self.GENlep1_eta = self.GENlep_eta[self.GENlep_pt.index(self.GENlep_pt_ordered[0])]
        self.GENlep1_mass = self.GENlep_mass[self.GENlep_pt.index(self.GENlep_pt_ordered[0])]
        self.GENlep1_phi = self.GENlep_phi[self.GENlep_pt.index(self.GENlep_pt_ordered[0])]
        self.GENlep2_pt = self.GENlep_pt_ordered[1]
        self.GENlep2_eta = self.GENlep_eta[self.GENlep_pt.index(self.GENlep_pt_ordered[1])]
        self.GENlep2_mass = self.GENlep_mass[self.GENlep_pt.index(self.GENlep_pt_ordered[1])]
        self.GENlep2_phi = self.GENlep_phi[self.GENlep_pt.index(self.GENlep_pt_ordered[1])]
        self.GENlep3_pt = self.GENlep_pt_ordered[2]
        self.GENlep3_eta = self.GENlep_eta[self.GENlep_pt.index(self.GENlep_pt_ordered[2])]
        self.GENlep3_mass = self.GENlep_mass[self.GENlep_pt.index(self.GENlep_pt_ordered[2])]
        self.GENlep3_phi = self.GENlep_phi[self.GENlep_pt.index(self.GENlep_pt_ordered[2])]
        self.GENlep4_pt = self.GENlep_pt_ordered[3]
        self.GENlep4_eta = self.GENlep_eta[self.GENlep_pt.index(self.GENlep_pt_ordered[3])]
        self.GENlep4_mass = self.GENlep_mass[self.GENlep_pt.index(self.GENlep_pt_ordered[3])]
        self.GENlep4_phi = self.GENlep_phi[self.GENlep_pt.index(self.GENlep_pt_ordered[3])]

        self.GENmassZ1                       = treeEntry.GENmassZ1
        self.GENmassZ2                       = treeEntry.GENmassZ2

        self.GENnjets_pt30_eta4p7            = len(list(treeEntry.GENjetsPt_pt30_eta4p7))
        self.GENjetsPt_pt30_eta4p7           = treeEntry.GENjetsPt_pt30_eta4p7
        self.GENjetsEta_pt30_eta4p7          = treeEntry.GENjetsEta_pt30_eta4p7
        self.GENjetsPhi_pt30_eta4p7          = treeEntry.GENjetsPhi_pt30_eta4p7
        self.GENjetsMass_pt30_eta4p7         = treeEntry.GENjetsMass_pt30_eta4p7

        # self.GENnjets_pt30_eta2p5            = treeEntry.GENnjets_pt30_eta2p5
        # self.GENjetsPt_pt30_eta2p5           = treeEntry.GENjetsPt_pt30_eta2p5
        # self.GENjetsEta_pt30_eta2p5          = treeEntry.GENjetsEta_pt30_eta2p5
        # self.GENjetsPhi_pt30_eta2p5          = treeEntry.GENjetsPhi_pt30_eta2p5
        # self.GENjetsMass_pt30_eta2p5         = treeEntry.GENjetsMass_pt30_eta2p5

        self.GENpTHj = -1.
        self.GENmHj = -1.
        self.GENpTHjj = -1.
        self.GENmHjj = -1.
        self.GENdetajj = -1.
        self.GENabsdetajj = -1.
        self.GENmjj = -1.
        self.GENdphijj = -1.
        self.GENabsdphijj = -1.

        self.GENpTj1                         = -1.
        self.GENmj1                          = -1.
        self.GENetaj1                        = -1.
        self.GENphij1                        = -1.
        for i in range(len(self.GENjetsPt_pt30_eta4p7)):
            if self.GENjetsPt_pt30_eta4p7[i] > self.GENpTj1:
                self.GENpTj1 = self.GENjetsPt_pt30_eta4p7[i]
                self.GENmj1 = self.GENjetsMass_pt30_eta4p7[i]
                self.GENetaj1 = self.GENjetsEta_pt30_eta4p7[i]
                self.GENphij1 = self.GENjetsPhi_pt30_eta4p7[i]
        if self.GENpTj1 > 0 and self.passedFiducial:
            H = TLorentzVector()
            H.SetPtEtaPhiM(self.GENpT4l, self.GENeta4l, self.GENphi4l, self.GENmass4l)
            j1 = TLorentzVector()
            j1.SetPtEtaPhiM(self.GENpTj1, self.GENetaj1, self.GENphij1, self.GENmj1)
            self.GENpTHj = (H+j1).Pt();
            self.GENmHj = (H+j1).M();

        self.GENpTj2                         = -1.
        self.GENmj2                          = -1.
        self.GENetaj2                        = -1.
        self.GENphij2                        = -1.
        for i in range(len(self.GENjetsPt_pt30_eta4p7)):
            if (self.GENjetsPt_pt30_eta4p7[i] > self.GENpTj2) and (self.GENpTj1 != self.GENjetsPt_pt30_eta4p7[i]):
                self.GENpTj2 = self.GENjetsPt_pt30_eta4p7[i]
                self.GENmj2 = self.GENjetsMass_pt30_eta4p7[i]
                self.GENetaj2 = self.GENjetsEta_pt30_eta4p7[i]
                self.GENphij2 = self.GENjetsPhi_pt30_eta4p7[i]

        if self.GENpTj2 > 0 and self.passedFiducial:
            H = TLorentzVector()
            H.SetPtEtaPhiM(self.GENpT4l, self.GENeta4l, self.GENphi4l, self.GENmass4l)
            j1 = TLorentzVector()
            j1.SetPtEtaPhiM(self.GENpTj1, self.GENetaj1, self.GENphij1, self.GENmj1)
            j2 = TLorentzVector()
            j2.SetPtEtaPhiM(self.GENpTj2, self.GENetaj2, self.GENphij2, self.GENmj2)
            self.GENpTHj = (H+j1).Pt();
            self.GENmHj = (H+j1).M();
            self.GENpTHjj = (H+j1+j2).Pt();
            self.GENmHjj = (H+j1+j2).M();
            self.GENdetajj = j1.Eta()-j2.Eta();
            self.GENabsdetajj = abs(self.GENdetajj);
            self.GENmjj = (j1+j2).M();
            self.GENdphijj = math.atan2(math.sin(j1.Phi()-j2.Phi()), math.cos(j1.Phi()-j2.Phi()));
            self.GENabsdphijj = abs(self.GENdphijj);


        if self.passedFiducial:
            self.GEND0m           = treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + (treeEntry.p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen * pow(spline_g4.Eval(self.GENmass4l),2)))
            self.GENDcp           = treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen / (2 * math.sqrt(treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen * treeEntry.p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen))
            self.GEND0hp          = treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + (treeEntry.p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen * pow(spline_g2.Eval(self.GENmass4l),2)))
            self.GENDint          = treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * math.sqrt(treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen * treeEntry.p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen))
            self.GENDL1           = treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(spline_L1.Eval(self.GENmass4l),2)))
            self.GENDL1Zg         = treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((treeEntry.p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(spline_L1Zgs.Eval(self.GENmass4l),2)))
        else:
            self.GEND0m           = -1
            self.GENDcp           = -1
            self.GEND0hp          = -1
            self.GENDint          = -1
            self.GENDL1           = -1
            self.GENDL1Zg         = -1

        self.GENTCjmax        = -1.
        self.GENTBjmax        = -1.
        for index in range(len(self.GENjetsPt_pt30_eta4p7)):
            theJet = TLorentzVector()
            theJet.SetPtEtaPhiM(self.GENjetsPt_pt30_eta4p7[index], self.GENjetsEta_pt30_eta4p7[index], self.GENjetsPhi_pt30_eta4p7[index], self.GENjetsMass_pt30_eta4p7[index])
            H = TLorentzVector()
            H.SetPtEtaPhiM(self.GENpT4l, self.GENeta4l, self.GENphi4l, self.GENmass4l)
            GENTCj = math.sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*math.cosh(theJet.Rapidity() - H.Rapidity()))
            GENTBj = math.sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*math.exp(-1*abs(theJet.Rapidity() - H.Rapidity()))
            if GENTCj > self.GENTCjmax: self.GENTCjmax = GENTCj
            if GENTBj > self.GENTBjmax: self.GENTBjmax = GENTBj

        self.mass4l = -1.
        self.eta4l = -1.
        self.phi4l = -1.
        self.pt4l = -1.
        self.Z1Mass = -1.
        self.Z2Mass = -1.
        self.lepBDT = [-1., -1., -1., -1.]
        self.njets30 = -1.
        self.jet1pt = -1.
        self.jet1eta = -1.
        self.jet1phi = -1.
        self.jet1mass = -1.
        self.jet2pt = -1.
        self.jet2eta = -1.
        self.jet2phi = -1.
        self.jet2mass = -1.
        self.njets30eta4p7 = -1.
        self.weight = -1.
        self.weight_sign = -1.
        self.D0m = -1.
        self.Dcp = -1.
        self.D0hp = -1.
        self.Dint = -1.
        self.DL1 = -1.
        self.DL1Zg = -1.
        self.TCjmax = -1.
        self.TBjmax = -1.
        self.mjj = -1.
        self.detajj = -1.
        self.pTHj = -1.
        self.mHj = -1.
        self.pTHjj = -1.
        self.mHjj = -1.
        self.absdetajj = -1.
        self.mjj = -1.
        self.dphijj = -1.
        self.absdphijj = -1.
        self.LepPt = [-1., -1., -1., -1.]
        self.LepEta = [-1., -1., -1., -1.]
        self.LepPhi = [-1., -1., -1., -1.]
        self.LepMass = [-1., -1., -1., -1.]
        self.LepPtOrdered = [-1., -1., -1., -1.]
        self.LepEtaOrdered = [-1., -1., -1., -1.]
        self.LepPhiOrdered = [-1., -1., -1., -1.]
        self.LepMassOrdered = [-1., -1., -1., -1.]


        if isReco:
            self.Z1Flav           = treeEntry.Z1Flav
            self.Z2Flav           = treeEntry.Z2Flav
            self.ZZFlav           = self.Z1Flav * self.Z2Flav
            self.Z1Mass           = treeEntry.Z1Mass
            self.Z2Mass           = treeEntry.Z2Mass

            self.mass4l           = treeEntry.ZZMass
            self.eta4l            = treeEntry.ZZEta
            self.phi4l            = treeEntry.ZZPhi
            self.pt4l             = treeEntry.ZZPt

            self.LepPt = []
            self.LepEta = []
            self.LepPhi = []
            self.LepMass = []
            for i in range(len(treeEntry.LepPt)):
                self.LepPt.append(treeEntry.LepPt[i])
                self.LepEta.append(treeEntry.LepEta[i])
                self.LepPhi.append(treeEntry.LepPhi[i])
                self.LepMass.append(0)

            #Order leptons in decreasing order in each pair
            self.LepPtOrdered = []
            self.LepEtaOrdered = []
            self.LepPhiOrdered = []
            self.LepMassOrdered = []
            lepts = [(self.LepPt[0],self.LepEta[0],self.LepPhi[0],self.LepMass[0]), (self.LepPt[1],self.LepEta[1],self.LepPhi[1],self.LepMass[1])]
            lepts = sorted(lepts, key=lambda lept: lept[0], reverse=True)
            for i in lepts:
                self.LepPtOrdered.append(i[0])
                self.LepEtaOrdered.append(i[1])
                self.LepPhiOrdered.append(i[2])
                self.LepMassOrdered.append(i[3])
            lepts = [(self.LepPt[2],self.LepEta[2],self.LepPhi[2],self.LepMass[2]), (self.LepPt[3],self.LepEta[3],self.LepPhi[3],self.LepMass[3])]
            lepts = sorted(lepts, key=lambda lept: lept[0], reverse=True)
            for i in lepts:
                self.LepPtOrdered.append(i[0])
                self.LepEtaOrdered.append(i[1])
                self.LepPhiOrdered.append(i[2])
                self.LepMassOrdered.append(i[3])


            self.nExtraLep        = treeEntry.nExtraLep
            self.nExtraZ          = treeEntry.nExtraZ
            self.jetpt            = treeEntry.JetPt
            self.jeteta           = treeEntry.JetEta
            self.jetphi           = treeEntry.JetPhi
            self.jetmass          = treeEntry.JetMass
            self.njets30          = treeEntry.nCleanedJetsPt30
            self.njets30Btag      = treeEntry.nCleanedJetsPt30BTagged
            self.mjj              = treeEntry.DiJetMass
            self.detajj           = treeEntry.DiJetDEta
            self.pfMet            = treeEntry.PFMET

            self.genHEPMC         = treeEntry.genHEPMCweight
            self.PUweight         = treeEntry.PUWeight
            self.dataMC           = treeEntry.dataMCWeight
            self.L1pref           = treeEntry.L1prefiringWeight
            self.trigEff          = treeEntry.trigEffWeight
            self.weight           = 1
            if (isMC):
                self.weight_sign = sign(treeEntry.genHEPMCweight) * treeEntry.PUWeight * treeEntry.dataMCWeight * treeEntry.L1prefiringWeight * treeEntry.trigEffWeight
                self.weight = treeEntry.genHEPMCweight * treeEntry.PUWeight * treeEntry.dataMCWeight * treeEntry.L1prefiringWeight * treeEntry.trigEffWeight

            self.jets30pt             = []
            self.jets30eta            = []
            self.jets30phi            = []
            self.jets30mass           = []
            self.jets30bTag           = []
            jet_counter = 0
            for i in range(len(treeEntry.JetPt)):
                if treeEntry.JetPt[i]>30. and abs(treeEntry.JetEta[i])<4.7: #Double-check on pT and cut on eta
                    jet_counter += 1
                    self.jets30pt.append(treeEntry.JetPt[i])
                    self.jets30eta.append(treeEntry.JetEta[i])
                    self.jets30phi.append(treeEntry.JetPhi[i])
                    self.jets30mass.append(treeEntry.JetMass[i])
                    self.jets30bTag.append(treeEntry.JetBTagger[i])
            self.njets30eta4p7 = len(self.jets30pt)


            self.lepBDT = []
            for i in range(len(treeEntry.LepBDT)):
                if (treeEntry.LepLepId[i] == 13 or treeEntry.LepLepId[i] == -13):
                    self.lepBDT.append(-1.)
                else:
                    self.lepBDT.append(treeEntry.LepBDT[i])


            # self.jet1pt           = -1.
            # self.jet2pt           = -1.
            self.jet1bTag         = -2000.
            self.jet2bTag         = -2000.
            self.fillJetInfo()

            self.D0m           = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + (treeEntry.p_GG_SIG_ghg2_1_ghz4_1_JHUGen * pow(spline_g4.Eval(self.mass4l),2)))
            self.Dcp           = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen / (2 * math.sqrt(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen * treeEntry.p_GG_SIG_ghg2_1_ghz4_1_JHUGen))
            self.D0hp          = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + (treeEntry.p_GG_SIG_ghg2_1_ghz2_1_JHUGen * pow(spline_g2.Eval(self.mass4l),2)))
            self.Dint          = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * math.sqrt(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen * treeEntry.p_GG_SIG_ghg2_1_ghz2_1_JHUGen))
            self.DL1           = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((treeEntry.p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(spline_L1.Eval(self.mass4l),2)))
            self.DL1Zg         = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((treeEntry.p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(spline_L1Zgs.Eval(self.mass4l),2)))

            self.TCjmax       = -1.
            self.TBjmax       = -1.
            for index in range(len(self.jets30pt)):
                if abs(self.jets30eta[index]) < 4.7:
                    theJet = TLorentzVector()
                    theJet.SetPtEtaPhiM(self.jets30pt[index], self.jets30eta[index], self.jets30phi[index], self.jets30mass[index])
                    H = TLorentzVector()
                    H.SetPtEtaPhiM(self.pt4l, self.eta4l, self.phi4l, self.mass4l)
                    TCj = math.sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*math.cosh(theJet.Rapidity() - H.Rapidity()))
                    TBj = math.sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*math.exp(-1*abs(theJet.Rapidity() - H.Rapidity()))
                    if TCj > self.TCjmax: self.TCjmax = TCj
                    if TBj > self.TBjmax: self.TBjmax = TBj


    def fillJetInfo(self):

        if self.njets30eta4p7 == 0: self.njets30eta4p7 = -99

        if self.njets30==1:
            self.jet1pt = self.jets30pt[0]
            self.jet1eta = self.jets30eta[0]
            self.jet1phi = self.jets30phi[0]
            self.jet1mass = self.jets30mass[0]
            self.jet1bTag = self.jets30bTag[0]
            self.jet2pt = -99
            self.jet2eta = -99
            self.jet2phi = -99
            self.jet2mass = -99
            self.mjj = -1.
            self.detajj = -1.

            H = TLorentzVector()
            H.SetPtEtaPhiM(self.pt4l, self.eta4l, self.phi4l, self.mass4l)
            j1 = TLorentzVector()
            j1.SetPtEtaPhiM(self.jet1pt, self.jet1eta, self.jet1phi, self.jet1mass)
            self.pTHj = (H+j1).Pt();
            self.mHj = (H+j1).M();
            self.pTHjj = -1.
            self.mHjj = -1.
            self.detajj = -1.
            self.absdetajj = -1.
            self.mjj = -1.
            self.dphijj = -1.
            self.absdphijj = -1.

        elif self.njets30>=2:
            self.jet1pt = self.jets30pt[0]
            self.jet1eta = self.jets30eta[0]
            self.jet1phi = self.jets30phi[0]
            self.jet1mass = self.jets30mass[0]
            self.jet2pt = self.jets30pt[1]
            self.jet2eta = self.jets30eta[1]
            self.jet2phi = self.jets30phi[1]
            self.jet2mass = self.jets30mass[1]
            self.jet1bTag = self.jets30bTag[0]
            self.jet2bTag = self.jets30bTag[1]

            H = TLorentzVector()
            H.SetPtEtaPhiM(self.pt4l, self.eta4l, self.phi4l, self.mass4l)
            j1 = TLorentzVector()
            j1.SetPtEtaPhiM(self.jet1pt, self.jet1eta, self.jet1phi, self.jet1mass)
            j2 = TLorentzVector()
            j2.SetPtEtaPhiM(self.jet2pt, self.jet2eta, self.jet2phi, self.jet2mass)
            self.pTHj = (H+j1).Pt();
            self.mHj = (H+j1).M();
            self.pTHjj = (H+j1+j2).Pt();
            self.mHjj = (H+j1+j2).M();
            self.detajj = j1.Eta()-j2.Eta();
            self.absdetajj = abs(self.detajj);
            self.mjj = (j1+j2).M();
            self.dphijj = math.atan2(math.sin(j1.Phi()-j2.Phi()), math.cos(j1.Phi()-j2.Phi()));
            self.absdphijj = abs(self.dphijj);

        else:
            self.jet1pt = -99
            self.jet1eta = -99
            self.jet1phi = -99
            self.jet1mass = -99
            self.jet2pt = -99
            self.jet2eta = -99
            self.jet2phi = -99
            self.jet2mass = -99

            self.mjj = -1.
            self.detajj = -1.

            self.pTHj = -1.
            self.mHj = -1.
            self.pTHjj = -1.
            self.mHjj = -1.
            self.detajj = -1.
            self.absdetajj = -1.
            self.mjj = -1.
            self.dphijj = -1.
            self.absdphijj = -1.



    def printOut(self):
        # line = ""
        # line  += str(int(self.run))
        # line  += ":" + str(int(self.lumi))
        # line  += ":" + str(int(self.event))
        # line  += ":{0:.2f}".format(self.mass4l)
        # line  += ":" + "{0:.2f}".format(self.Z1Mass)
        # line  += ":" + "{0:.2f}".format(self.Z2Mass)
        # line  += ":" + "{0:.6f}".format(self.lepBDT[0])
        # line  += ":" + "{0:.6f}".format(self.lepBDT[1])
        # line  += ":" + "{0:.6f}".format(self.lepBDT[2])
        # line  += ":" + "{0:.6f}".format(self.lepBDT[3])
        # line  += ":" + "{0:.4f}".format(self.D0m)
        # line  += ":" + "{0:.4f}".format(self.Dcp)
        # line  += ":" + "{0:.4f}".format(self.D0hp)
        # line  += ":" + "{0:.4f}".format(self.Dint)
        # line  += ":" + "{0:.4f}".format(self.DL1)
        # line  += ":" + "{0:.4f}".format(self.DL1Zg)
        # line  += ":" + "{0:.4f}".format(self.TCjmax)
        # line  += ":" + "{0:.4f}".format(self.TBjmax)
        # line  += ":" + "{0:f}".format(self.njets30)
        # line  += ":" + "{0:.2f}".format(self.jet1pt)
        # line  += ":" + "{0:.2f}".format(self.jet2pt)
        # line  += ":" + "{0:.2f}".format(self.mjj)
        # line  += ":" + "{0:.2f}".format(self.detajj)
        # line  += ":" + "{0:.3f}".format(self.weight_sign)
        # #GEN variables
        # line  += ":{0:.2f}".format(self.GENmass4l)
        # line  += ":" + "{0:.2f}".format(self.GENpT4l)
        # line  += ":" + "{0:.2f}".format(self.GENeta4l)
        # line  += ":" + "{0:.2f}".format(self.GENphi4l)
        # line  += ":" + "{0:.2f}".format(self.GENmassZ1)
        # line  += ":" + "{0:.2f}".format(self.GENmassZ2)
        # line  += ":" + "{0:.4f}".format(self.GEND0m)
        # line  += ":" + "{0:.4f}".format(self.GENDcp)
        # line  += ":" + "{0:.4f}".format(self.GEND0hp)
        # line  += ":" + "{0:.4f}".format(self.GENDint)
        # line  += ":" + "{0:.4f}".format(self.GENDL1)
        # line  += ":" + "{0:.4f}".format(self.GENDL1Zg)
        # line  += ":" + "{0:.4f}".format(self.GENTCjmax)
        # line  += ":" + "{0:.4f}".format(self.GENTBjmax)
        # line  += ":" + "{0:.2f}".format(self.GENpTj1)
        # line  += ":" + "{0:.2f}".format(self.GENmj1)
        # line  += ":" + "{0:.2f}".format(self.GENetaj1)
        # line  += ":" + "{0:.2f}".format(self.GENphij1)
        # line  += ":" + "{0:.2f}".format(self.GENpTj2)
        # line  += ":" + "{0:.2f}".format(self.GENmj2)
        # line  += ":" + "{0:.2f}".format(self.GENetaj2)
        # line  += ":" + "{0:.2f}".format(self.GENphij2)
        # line  += ":" + "{0:f}".format(self.passedFiducial)

        # # Lepton kinematics
        # line = ""
        # line  += str(int(self.run))
        # line  += ":" + str(int(self.lumi))
        # line  += ":" + str(int(self.event))
        # line  += ":" + "{0:.4f}".format(self.LepPtOrdered[0])
        # line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[0])
        # line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[0])
        # line  += ":" + "{0:.4f}".format(self.LepMassOrdered[0])
        # line  += ":" + "{0:.4f}".format(self.LepPtOrdered[1])
        # line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[1])
        # line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[1])
        # line  += ":" + "{0:.4f}".format(self.LepMassOrdered[1])
        # line  += ":" + "{0:.4f}".format(self.LepPtOrdered[2])
        # line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[2])
        # line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[2])
        # line  += ":" + "{0:.4f}".format(self.LepMassOrdered[2])
        # line  += ":" + "{0:.4f}".format(self.LepPtOrdered[3])
        # line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[3])
        # line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[3])
        # line  += ":" + "{0:.4f}".format(self.LepMassOrdered[3])

        # Jets kinematics
        # line = ""
        # line  += str(int(self.run))
        # line  += ":" + str(int(self.lumi))
        # line  += ":" + str(int(self.event))
        # line  += ":" + "{0:.4f}".format(self.jet1pt)
        # line  += ":" + "{0:.4f}".format(self.jet1eta)
        # line  += ":" + "{0:.4f}".format(self.jet1phi)
        # line  += ":" + "{0:.4f}".format(self.jet1mass)
        # line  += ":" + "{0:.4f}".format(self.jet2pt)
        # line  += ":" + "{0:.4f}".format(self.jet2eta)
        # line  += ":" + "{0:.4f}".format(self.jet2phi)
        # line  += ":" + "{0:.4f}".format(self.jet2mass)
        # line  += ":" + "{0:.4f}".format(self.njets30eta4p7)

        # # GEN kinematics
        # if self.passedFiducial == 0.:
        #     self.GENmassZ1 = -1.
        #     self.GENmassZ2 = -1.
        #     self.GENmass4l = -1.
        # line = ""
        # line  += str(int(self.run))
        # line  += ":" + str(int(self.lumi))
        # line  += ":" + str(int(self.event))
        # line  += ":" + "{0:.4f}".format(self.GENmassZ1)
        # line  += ":" + "{0:.4f}".format(self.GENmassZ2)
        # line  += ":" + "{0:.4f}".format(self.GENmass4l)
        # line  += ":" + "{0:.4f}".format(self.GENlep1_pt)
        # line  += ":" + "{0:.4f}".format(self.GENlep1_eta)
        # line  += ":" + "{0:.4f}".format(self.GENlep1_phi)
        # line  += ":" + "{0:.4f}".format(self.GENlep1_mass)
        # line  += ":" + "{0:.4f}".format(self.GENlep2_pt)
        # line  += ":" + "{0:.4f}".format(self.GENlep2_eta)
        # line  += ":" + "{0:.4f}".format(self.GENlep2_phi)
        # line  += ":" + "{0:.4f}".format(self.GENlep2_mass)
        # line  += ":" + "{0:.4f}".format(self.GENlep3_pt)
        # line  += ":" + "{0:.4f}".format(self.GENlep3_eta)
        # line  += ":" + "{0:.4f}".format(self.GENlep3_phi)
        # line  += ":" + "{0:.4f}".format(self.GENlep3_mass)
        # line  += ":" + "{0:.4f}".format(self.GENlep4_pt)
        # line  += ":" + "{0:.4f}".format(self.GENlep4_eta)
        # line  += ":" + "{0:.4f}".format(self.GENlep4_phi)
        # line  += ":" + "{0:.4f}".format(self.GENlep4_mass)
        # line  += ":" + "{0:.4f}".format(self.GENpTj1)
        # line  += ":" + "{0:.4f}".format(self.GENetaj1)
        # line  += ":" + "{0:.4f}".format(self.GENphij1)
        # line  += ":" + "{0:.4f}".format(self.GENmj1)
        # line  += ":" + "{0:.4f}".format(self.GENpTj2)
        # line  += ":" + "{0:.4f}".format(self.GENetaj2)
        # line  += ":" + "{0:.4f}".format(self.GENphij2)
        # line  += ":" + "{0:.4f}".format(self.GENmj2)
        # line  += ":" + str(int(self.GENnjets_pt30_eta4p7))
        # line  += ":" + str(int(self.passedFiducial))

        # # Higher level
        # line = ""
        # line  += str(int(self.run))
        # line  += ":" + str(int(self.lumi))
        # line  += ":" + str(int(self.event))
        # line  += ":" + "{0:.3f}".format(self.D0m)
        # line  += ":" + "{0:.3f}".format(self.Dcp)
        # line  += ":" + "{0:.3f}".format(self.D0hp)
        # line  += ":" + "{0:.3f}".format(self.Dint)
        # line  += ":" + "{0:.3f}".format(self.DL1)
        # line  += ":" + "{0:.3f}".format(self.DL1Zg)
        # line  += ":" + "{0:.3f}".format(self.TCjmax)
        # line  += ":" + "{0:.3f}".format(self.TBjmax)
        # line  += ":" + "{0:.3f}".format(self.mjj)
        # line  += ":" + "{0:.3f}".format(self.detajj)
        # line  += ":" + "{0:.3f}".format(self.dphijj)
        # line  += ":" + "{0:.3f}".format(self.pTHj)
        # line  += ":" + "{0:.3f}".format(self.pTHjj)
        # line  += ":" + "{0:.3f}".format(self.mHj)
        # line  += ":" + "{0:.3f}".format(self.mHjj)
        # line  += ":" + "{0:.4f}".format(self.weight_sign)
        # line  += ":" + "{0:.3f}".format(self.GEND0m)
        # line  += ":" + "{0:.3f}".format(self.GENDcp)
        # line  += ":" + "{0:.3f}".format(self.GEND0hp)
        # line  += ":" + "{0:.3f}".format(self.GENDint)
        # line  += ":" + "{0:.3f}".format(self.GENDL1)
        # line  += ":" + "{0:.3f}".format(self.GENDL1Zg)
        # line  += ":" + "{0:.3f}".format(self.GENTCjmax)
        # line  += ":" + "{0:.3f}".format(self.GENTBjmax)
        # line  += ":" + "{0:.3f}".format(self.GENmjj)
        # line  += ":" + "{0:.3f}".format(self.GENdetajj)
        # line  += ":" + "{0:.3f}".format(self.GENdphijj)
        # line  += ":" + "{0:.3f}".format(self.GENpTHj)
        # line  += ":" + "{0:.3f}".format(self.GENpTHjj)
        # line  += ":" + "{0:.3f}".format(self.GENmHj)
        # line  += ":" + "{0:.3f}".format(self.GENmHjj)


        # general format
        line = ""
        line  += str(int(self.run))
        line  += ":" + str(int(self.lumi))
        line  += ":" + str(int(self.event))
        line  += ":" + "{0:.4f}".format(self.mass4l)
        line  += ":" + "{0:.4f}".format(self.eta4l)
        line  += ":" + "{0:.4f}".format(self.phi4l)
        line  += ":" + "{0:.4f}".format(self.pt4l)
        line  += ":" + "{0:.4f}".format(self.LepPtOrdered[0])
        line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[0])
        line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[0])
        line  += ":" + "{0:.4f}".format(self.LepPtOrdered[1])
        line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[1])
        line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[1])
        line  += ":" + "{0:.4f}".format(self.LepPtOrdered[2])
        line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[2])
        line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[2])
        line  += ":" + "{0:.4f}".format(self.LepPtOrdered[3])
        line  += ":" + "{0:.4f}".format(self.LepEtaOrdered[3])
        line  += ":" + "{0:.4f}".format(self.LepPhiOrdered[3])
        line  += ":" + "{0:.4f}".format(self.jet1pt)
        line  += ":" + "{0:.4f}".format(self.jet1eta)
        line  += ":" + "{0:.4f}".format(self.jet1phi)
        line  += ":" + "{0:.4f}".format(self.jet1mass)
        line  += ":" + "{0:.4f}".format(self.jet2pt)
        line  += ":" + "{0:.4f}".format(self.jet2eta)
        line  += ":" + "{0:.4f}".format(self.jet2phi)
        line  += ":" + "{0:.4f}".format(self.jet2mass)
        line  += ":" + "{0:.4f}".format(self.njets30eta4p7)
        line  += ":" + "{0:.4f}".format(self.GENmassZ1)
        line  += ":" + "{0:.4f}".format(self.GENmassZ2)
        line  += ":" + "{0:.4f}".format(self.GENmass4l)
        line  += ":" + "{0:.4f}".format(self.GENlep1_pt)
        line  += ":" + "{0:.4f}".format(self.GENlep1_eta)
        line  += ":" + "{0:.4f}".format(self.GENlep1_phi)
        line  += ":" + "{0:.4f}".format(self.GENlep1_mass)
        line  += ":" + "{0:.4f}".format(self.GENlep2_pt)
        line  += ":" + "{0:.4f}".format(self.GENlep2_eta)
        line  += ":" + "{0:.4f}".format(self.GENlep2_phi)
        line  += ":" + "{0:.4f}".format(self.GENlep2_mass)
        line  += ":" + "{0:.4f}".format(self.GENlep3_pt)
        line  += ":" + "{0:.4f}".format(self.GENlep3_eta)
        line  += ":" + "{0:.4f}".format(self.GENlep3_phi)
        line  += ":" + "{0:.4f}".format(self.GENlep3_mass)
        line  += ":" + "{0:.4f}".format(self.GENlep4_pt)
        line  += ":" + "{0:.4f}".format(self.GENlep4_eta)
        line  += ":" + "{0:.4f}".format(self.GENlep4_phi)
        line  += ":" + "{0:.4f}".format(self.GENlep4_mass)
        line  += ":" + "{0:.4f}".format(self.GENpTj1)
        line  += ":" + "{0:.4f}".format(self.GENetaj1)
        line  += ":" + "{0:.4f}".format(self.GENphij1)
        line  += ":" + "{0:.4f}".format(self.GENmj1)
        line  += ":" + "{0:.4f}".format(self.GENpTj2)
        line  += ":" + "{0:.4f}".format(self.GENetaj2)
        line  += ":" + "{0:.4f}".format(self.GENphij2)
        line  += ":" + "{0:.4f}".format(self.GENmj2)
        line  += ":" + str(int(self.GENnjets_pt30_eta4p7))
        line  += ":" + str(int(self.passedFiducial))
        line  += ":" + "{0:.4f}".format(self.weight_sign)
        line  += ":" + "{0:.3f}".format(self.D0m)
        line  += ":" + "{0:.3f}".format(self.Dcp)
        line  += ":" + "{0:.3f}".format(self.D0hp)
        line  += ":" + "{0:.3f}".format(self.Dint)
        line  += ":" + "{0:.3f}".format(self.DL1)
        line  += ":" + "{0:.3f}".format(self.DL1Zg)
        line  += ":" + "{0:.3f}".format(self.TCjmax)
        line  += ":" + "{0:.3f}".format(self.TBjmax)
        line  += ":" + "{0:.3f}".format(self.mjj)
        line  += ":" + "{0:.3f}".format(self.detajj)
        line  += ":" + "{0:.3f}".format(self.dphijj)
        line  += ":" + "{0:.3f}".format(self.pTHj)
        line  += ":" + "{0:.3f}".format(self.pTHjj)
        line  += ":" + "{0:.3f}".format(self.mHj)
        line  += ":" + "{0:.3f}".format(self.mHjj)
        line  += ":" + "{0:.3f}".format(self.GEND0m)
        line  += ":" + "{0:.3f}".format(self.GENDcp)
        line  += ":" + "{0:.3f}".format(self.GEND0hp)
        line  += ":" + "{0:.3f}".format(self.GENDint)
        line  += ":" + "{0:.3f}".format(self.GENDL1)
        line  += ":" + "{0:.3f}".format(self.GENDL1Zg)
        line  += ":" + "{0:.3f}".format(self.GENTCjmax)
        line  += ":" + "{0:.3f}".format(self.GENTBjmax)
        line  += ":" + "{0:.3f}".format(self.GENmjj)
        line  += ":" + "{0:.3f}".format(self.GENdetajj)
        line  += ":" + "{0:.3f}".format(self.GENdphijj)
        line  += ":" + "{0:.3f}".format(self.GENpTHj)
        line  += ":" + "{0:.3f}".format(self.GENpTHjj)
        line  += ":" + "{0:.3f}".format(self.GENmHj)
        line  += ":" + "{0:.3f}".format(self.GENmHjj)


        return line
