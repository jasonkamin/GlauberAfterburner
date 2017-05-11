// An example for adding a ZDC signal to 
// the glauber output.
//
// -J.Kamin
// 8 May 2017
//


void UpdateZDCvalues(char TreeFileName[], char TreeObjectName[])
{

  TRandom3 randy;
  randy.SetSeed(5281978);
  char saythis[500];
  cout << TreeFileName << endl;

  // back up the input file... 
  sprintf(saythis,"cp %s %s.original_zdc",TreeFileName,TreeFileName);
  gSystem->Exec(saythis);

  TFile *f = new TFile(TreeFileName,"update");
  TTree *T = (TTree*)f->Get(TreeObjectName);
  Int_t   Ngray;
  Int_t   Nblack;
  Int_t   NgrayHit;
  Int_t   NblackHit;
  Float_t SNenergy;
  Float_t ZDCresponse;
  Float_t Npart;
  Float_t Ncoll;

  Float_t  phiGray    = 0.0;
  Float_t  theGray    = 0.0;
  Float_t  phiBlack   = 0.0;
  Float_t  theBlack   = 0.0;
  Float_t  meanNslow  = 0.0;
  Float_t  meanNgray  = 0.0;
  Float_t  meanNblack = 0.0;
  // handy variables... 
  Float_t meanNgrayp  = 0.0;
  Float_t meanNblackp = 0.0;
  Float_t meanNslowp  = 0.0;
  Float_t N_LCF      = 0.0;


  // my new ZDC branches... 
  TBranch *bmeanNgrayp    = T->Branch("meanNgrayp",    &meanNgrayp,     "meanNgrayp/F");
  TBranch *bmeanNblackp   = T->Branch("meanNblackp",   &meanNblackp,    "meanNblackp/F");
  TBranch *bmeanNslowp    = T->Branch("meanNslowp",    &meanNslowp,     "meanNslowp/F");
  TBranch *bN_LCF         = T->Branch("N_LCF",         &N_LCF,          "N_LCF/F");
  TBranch *bmeanNgray     = T->Branch("meanNgray",     &meanNgray,      "meanNgray/F");
  TBranch *bmeanNblack    = T->Branch("meanNblack",    &meanNblack,     "meanNblack/F");
  TBranch *bmeanNslow     = T->Branch("meanNslow",     &meanNslow,      "meanNslow/F");

  TBranch *bphiGray       = T->Branch("phiGray",     &phiGray,      "phiGray/F");
  TBranch *btheGray       = T->Branch("theGray",     &theGray,      "theGray/F");
  TBranch *bphiBlack      = T->Branch("phiBlack",    &phiBlack,     "phiBlack/F");
  TBranch *btheBlack      = T->Branch("theBlack",    &theBlack,     "theBlack/F");

  TBranch *bNgray     = T->Branch("Ngray",     &Ngray,      "Ngray/I");
  TBranch *bNblack    = T->Branch("Nblack",    &Nblack,     "Nblack/I");
  TBranch *bNgrayHit  = T->Branch("NgrayHit",  &NgrayHit,   "NgrayHit/I");
  TBranch *bNblackHit = T->Branch("NblackHit", &NblackHit,  "NblackHit/I");
  TBranch *bSNenergy  = T->Branch("SNenergy",  &SNenergy,   "SNenergy/F");
  TBranch *bZDC       = T->Branch("ZDC",       &ZDCresponse,"ZDC/F");
  T->SetBranchAddress("Npart",&Npart);
  T->SetBranchAddress("Ncoll",&Ncoll);

  // some input variables to the toy MC...
  // stolen from Ferenc's paper: 
  // https://arxiv.org/pdf/hep-ph/0304065.pdf
  Double_t betaGray  = 0.05;
  Double_t betaBlack = 0.0;
  Double_t E0Gray  = 0.050;// GeV
  Double_t E0Black = 0.005;// GeV
  Double_t pAveGray  = 0.750;// GeV/c
  Double_t pAveBlack = 0.200;// GeV/c
  Double_t m_neutron = 0.939565;//GeV


  // dN/dtheta functional forms: 
  TF1 *dN_dThe_gray  = new TF1("dN_dThe_gray", "exp( ([0]*[1]/[2]) * cos(x))",0.0,TMath::Pi());
  TF1 *dN_dThe_black = new TF1("dN_dThe_black","exp( ([0]*[1]/[2]) * cos(x))",0.0,TMath::Pi());
  dN_dThe_gray ->SetParameters(betaGray, pAveGray, E0Gray);
  dN_dThe_black->SetParameters(betaBlack,pAveBlack,E0Black);

  TF1 *dN_dE_black  = new TF1("dN_dE_black", "exp(-x/[0])",0.0,10.0);//GeV !
  dN_dE_black ->SetParameter(0,E0Black);
  TF1 *dN_dE_gray  = new TF1("dN_dE_gray", "exp(-x/[0])",0.0,10.0);//GeV !
  dN_dE_gray ->SetParameter(0,E0Gray);


  // fit results from neutron peaks in 8.16TeV pPb data. 
  Double_t mu  = 67.43;
  Double_t sig = 15.27;

  // just a guess... (consider it a placeholder).
  Double_t ZdcAccTheta = 0.4e-3;//0.4e-3
  //Double_t ZdcAccTheta = 0.4;//dumb guess
  //Double_t ZdcAccTheta = 1e6;//inf


  // rough apprxs for relationship btwn NColl and Ngray/black.
  // taken from xyscan of plots in short ALICE Internal Note
  TF1 *f_p1 = new TF1("f_p1","[0]+[1]*x",0,25);
  f_p1->SetParameters(0.241111, 1.28119);
  TF1 *f_sat = new TF1("f_sat","([0]-[3])/(1+exp(-[2]*(x-[1])))+[3]",0,25);
  f_sat->SetParameters(18.1686, 2.09637, 0.518023, -6.13282);

  // taken from long ALICE centrality paper
  TF1 *f_Ngrayp = new TF1("f_Ngrayp","[0]+[1]*x+[2]*x*x",0,100);
  f_Ngrayp->SetParameters(-0.24,0.55,0.0007);
  TF1 *f_Nslown = new TF1("f_Nslown","[0]*x + ([1] - [2]/([3]+x))",0,100);
  f_Nslown->SetParameters(0.48,50,230,4.2);


  // Loop over the tree ! 
  Long64_t nentries = T->GetEntries();
  for(Long64_t i=0; i<nentries; i++){
    T->GetEntry(i);

    // rough apprx for relationship btwn NColl and Ngray/black.
    //meanNgray  = 2.0*Ncoll;
    //meanNblack = 4.0*Ncoll;
    //meanNgray  = f_p1->Eval(Ncoll);
    //meanNblack = f_sat->Eval(Ncoll);

    //complicated model from long ALICE centrality paper:
    meanNgrayp  = f_Ngrayp->Eval(Ncoll);
    meanNblackp = 0.65*meanNgrayp;
    meanNslowp  = meanNgrayp + meanNblackp;
    N_LCF       = 1.71*meanNslowp;
    meanNslow   = f_Nslown->Eval(N_LCF);
    meanNgray   = 0.1*meanNslow;
    meanNblack  = 0.9*meanNslow;

    // trying out gaus vs poisson vs binomial
    //Ngray  = randy.Gaus(meanNgray ,sqrt(meanNgray) );  if(Ngray <0) Ngray =0;
    //Nblack = randy.Gaus(meanNblack,sqrt(meanNblack));  if(Nblack<0) Nblack=0;
    //Ngray  = randy.Poisson(meanNgray);
    //Nblack = randy.Poisson(meanNblack);
    Ngray  = randy.Binomial(126, meanNgray/126.0);
    Nblack = randy.Binomial(126, meanNblack/126.0);
    //Ngray  = 0;
    //Nblack = randy.Binomial(126, meanNslow/126.0);

    NgrayHit = 0;
    NblackHit = 0;
    SNenergy = 0.0;
    ZDCresponse = 0.0;

    // Check which slow gray neutrons hit the ZDC:
    for(Long64_t j=0; j<Ngray; j++){
      TLorentzVector vecGray;
      Double_t mom =  dN_dE_gray->GetRandom();
      //dN_dThe_gray->SetParameter(1,mom);//commented because it's quite slow (and not so important).
      Double_t the = dN_dThe_gray->GetRandom();
      Double_t phi = TMath::Pi()*randy.Rndm()-0.5*TMath::Pi();
      vecGray.SetPtEtaPhiM(mom*sin(the), -1.0*TMath::Log(TMath::Tan(0.5*the)), phi, m_neutron);
      vecGray.Boost(0,0,TMath::Sqrt(1.0 - m_neutron*m_neutron/(2560.0*2650.0)));
      phiGray = vecGray.Phi();
      theGray = vecGray.Theta();
      if(theGray<ZdcAccTheta){
        NgrayHit++;
        SNenergy += vecGray.E();
      }
    }

    // Check which slow black neutrons hit the ZDC:
    for(Long64_t j=0; j<Nblack; j++){
      TLorentzVector vecBlack;
      Double_t mom =  dN_dE_black->GetRandom();
      //dN_dThe_black->SetParameter(1,mom);//commented because it's quite slow (and not so important).
      Double_t the = dN_dThe_black->GetRandom();
      Double_t phi = TMath::Pi()*randy.Rndm()-0.5*TMath::Pi();
      vecBlack.SetPtEtaPhiM(mom*sin(the), -1.0*TMath::Log(TMath::Tan(0.5*the)), phi, m_neutron);
      vecBlack.Boost(0,0,TMath::Sqrt(1.0 - m_neutron*m_neutron/(2560.0*2650.0)));
      phiBlack = vecBlack.Phi();
      theBlack = vecBlack.Theta();
      if(theBlack<ZdcAccTheta){
        NblackHit++;
        SNenergy += vecBlack.E();
      }
    }

    // Toy-Sim the ZDC response for a given slow neutron
    for(Long64_t j =0; j<NgrayHit; j++){
      Double_t resp = gRandom->Gaus(mu,sig);
      if(resp>0)
        ZDCresponse += resp;
    }
    for(Long64_t j =0; j<NblackHit; j++){
      Double_t resp = gRandom->Gaus(mu,sig);
      if(resp>0)
        ZDCresponse += resp;
    }


    // have a look at a couple values...
    if(i%int(0.1*nentries)==0){
      cout << "tree entry " << i << "/" << nentries << "    NColl: " << Ncoll << endl;
      cout << "NColl       " << Ncoll       << endl;
      cout << "meanNgrayp  " << meanNgrayp  << endl;
      cout << "meanNblackp " << meanNblackp << endl;
      cout << "meanNslowp  " << meanNslowp  << endl;
      cout << "N_LCF       " << N_LCF       << endl;

      cout << " <Ngray>: " << meanNgray << " <Nblack>: " << meanNblack << endl;
      cout << "   Ngray: " << Ngray << "   Nblack: " << Nblack << endl;
      cout << "   NGhit: " << NgrayHit << "   NBlhit: " << NblackHit << endl;
      cout << "   SNenergy: " << SNenergy << "    ZDC:" << ZDCresponse << endl;
      cout << endl;
    }

    bphiGray ->Fill();
    btheGray ->Fill();
    bphiBlack->Fill();
    btheBlack->Fill();

    bmeanNgrayp ->Fill();
    bmeanNblackp->Fill();
    bmeanNslowp ->Fill();
    bN_LCF      ->Fill();
    bmeanNgray ->Fill();
    bmeanNblack->Fill();
    bNgray     ->Fill();
    bNblack    ->Fill();
    bNgrayHit  ->Fill();
    bNblackHit ->Fill();
    bSNenergy  ->Fill();
    bZDC       ->Fill();

  }
  //T->Print();
  T->Write();
  delete f;

  sprintf(saythis,"%s.original_zdc",TreeFileName);
  cout << "-- Added ZDC response branch as 'ZDC' (plus other SNM variables)" << endl;
  cout << "-- Wrote " << TreeFileName << endl;
  cout << "-- and moved the original file to " << saythis << endl;

  return;
}

