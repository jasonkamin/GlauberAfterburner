// An example for adding a ZDC signal to 
// the glauber output.
//
// -J.Kamin
// 8 May 2017
//

void UpdateZDCvalues(char TreeFileName[], char TreeObjectName[], 
    Int_t algo   = 1    ,
    Double_t sp  = 0.0  , 
    Double_t bog = 0.65 , 
    Double_t gm  = 1.71 , 
    Double_t al  = 0.48 , 
    Double_t a   = 50.0 , 
    Double_t b   = 230.0, 
    Double_t c   = 4.2  ,
    Double_t fGn = 0.1  ,
    Double_t c1  = 1.47 ,
    Double_t sat = 30.4 )
{

  sigma_p   = sp;
  Nboverg   = bog;
  gm_par    = gm;
  alpha     = al;
  a_par     = a;
  b_par     = b;
  c_par     = c;
  fracGneut = fGn;
  c1_par    = c1;
  sat_par   = sat;

  cout << endl << endl;
  cout << "###################################" << endl;
  cout << "algo     : " << algo      << endl;
  cout << "sigma_p  : " << sigma_p   << endl;
  cout << "Nboverg  : " << Nboverg   << endl;
  cout << "gm_par   : " << gm_par    << endl;
  cout << "alpha    : " << alpha     << endl;
  cout << "a_par    : " << a_par     << endl;
  cout << "b_par    : " << b_par     << endl;
  cout << "c_par    : " << c_par     << endl;
  cout << "fracGneut: " << fracGneut << endl;
  cout << "c1_par   : " << c1_par    << endl;
  cout << "sat_par  : " << sat_par   << endl;
  cout << "###################################" << endl;
  cout << endl;

  TH1D *MYPARS = new TH1D("MYPARS","MYPARS",20,0.5,20.5);//store parameters and algo
  MYPARS->SetBinContent(1, sigma_p  );  MYPARS->GetXaxis()->SetBinLabel(1, "sigma_p  ");
  MYPARS->SetBinContent(2, Nboverg  );  MYPARS->GetXaxis()->SetBinLabel(2, "Nboverg  ");
  MYPARS->SetBinContent(3, gm_par   );  MYPARS->GetXaxis()->SetBinLabel(3, "gm_par   ");
  MYPARS->SetBinContent(4, alpha    );  MYPARS->GetXaxis()->SetBinLabel(4, "alpha    ");
  MYPARS->SetBinContent(5, a_par    );  MYPARS->GetXaxis()->SetBinLabel(5, "a_par    ");
  MYPARS->SetBinContent(6, b_par    );  MYPARS->GetXaxis()->SetBinLabel(6, "b_par    ");
  MYPARS->SetBinContent(7, c_par    );  MYPARS->GetXaxis()->SetBinLabel(7, "c_par    ");
  MYPARS->SetBinContent(8, fracGneut);  MYPARS->GetXaxis()->SetBinLabel(8, "fracGneut");
  MYPARS->SetBinContent(9, c1_par   );  MYPARS->GetXaxis()->SetBinLabel(9, "c1_par   ");
  MYPARS->SetBinContent(10,sat_par  );  MYPARS->GetXaxis()->SetBinLabel(10,"sat_par  ");
  MYPARS->SetBinContent(11,algo     );  MYPARS->GetXaxis()->SetBinLabel(11,"algo     ");

  TFile *f_NCVB = TFile::Open("NColl_vs_b.root");
  Ncollb_pfx = (TH1D*)f_NCVB->Get("Ncollb_pfx")->Clone("Ncollb_pfx");

  randy.SetSeed(5281978);
  char saythis[500];
  char saythis1[500];
  cout << TreeFileName << endl;

  // back up the input file... 
  sprintf(saythis,"cp %s %s.original_zdc",TreeFileName,TreeFileName);
  gSystem->Exec(saythis);

  TFile *f = new TFile(TreeFileName,"update");
  TTree *T = (TTree*)f->Get(TreeObjectName);
  Float_t Npart;
  Float_t Ncoll;
  Float_t B;
  T->SetBranchAddress("Npart",&Npart);
  T->SetBranchAddress("Ncoll",&Ncoll);
  T->SetBranchAddress("B",&B);


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



  dN_dThe_black->SetParameters(betaBlack,pAveBlack,E0Black);
  dN_dThe_gray ->SetParameters(betaGray, pAveGray, E0Gray);
  dN_dE_black  ->SetParameter(0,E0Black);
  dN_dE_gray   ->SetParameter(0,E0Gray);
  f_Ngrayp     ->SetParameters(-0.24,0.55,0.0007);
  f_Nslown     ->SetParameters(alpha,a_par,b_par,c_par);
  f_p1         ->SetParameter(0,c1_par);
  f_sat        ->SetParameters(sat_par, Nboverg, c1_par);
  f_DensPath   ->SetParameters(sat_par/6.62,6.62,0.546);

  double parsN[3] = {15.0,6.62,0.546};
  f_rho  ->SetParameters(parsN);

  // Loop over the tree ! 
  Long64_t nentries = T->GetEntries();
  for(Long64_t i=0; i<nentries; i++){
    T->GetEntry(i);


    int nSlown = 0;
    if(algo==1)      nSlown = GetNSlow_ALICE (Ncoll);
    else if(algo==2) nSlown = GetNSlow_ALICEp(Ncoll);
    else if(algo==3) nSlown = GetNSlow_LinSat(Ncoll);
    else if(algo==4) nSlown = GetNSlow_AveNColl_vB(Ncoll,B);
    else if(algo==5) nSlown = GetNSlow_Ferenc(Ncoll);
    else if(algo==6) nSlown = GetNSlow_ALICE_noFerencElse(Ncoll);
    else if(algo==7) nSlown = GetNSlow_ALICE_OnlyFerencElse(Ncoll);
    else if(algo==8) nSlown = GetNSlow_ALICE_OnlyElse(Ncoll);
    else if(algo==9) nSlown = GetNSlow_Black_b(Ncoll,B);


    NgrayHit = 0;
    NblackHit = 0;
    SNenergy = 0.0;
    ZDCresponse = 0.0;

    // Check which slow gray neutrons hit the ZDC:
    for(Long64_t j=0; j<Ngray; j++){
      TLorentzVector vecGray;
      Double_t mom =  TMath::Sqrt(2.0*m_neutron*dN_dE_gray->GetRandom());//go from KE to p. 
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
      Double_t mom =  TMath::Sqrt(2.0*m_neutron*dN_dE_black->GetRandom());//go from KE to p. 
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
    for(Long64_t j=0; j<NgrayHit; j++){
      Double_t resp = gRandom->Gaus(mu,sig)+ped;
      if(resp>0)
        ZDCresponse += resp;
    }
    for(Long64_t j=0; j<NblackHit; j++){
      Double_t resp = gRandom->Gaus(mu,sig)+ped;
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
    bmeanNgray  ->Fill();
    bmeanNblack ->Fill();
    bNgray      ->Fill();
    bNblack     ->Fill();
    bNgrayHit   ->Fill();
    bNblackHit  ->Fill();
    bSNenergy   ->Fill();
    bZDC        ->Fill();

  }
  //T->Print();
  T->Write();
  MYPARS->Write();
  delete f;

  // rename output file...
  sprintf(saythis1,"algo%d_pars_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f__%s",algo, sigma_p, Nboverg, gm_par , alpha  , a_par  , b_par  , c_par, c1_par, sat_par, TreeFileName);
  sprintf(saythis,"mv %s %s",TreeFileName,saythis1);
  gSystem->Exec(saythis);
  // recopy original file to here. 
  sprintf(saythis,"cp ../%s .",TreeFileName);
  gSystem->Exec(saythis);


  sprintf(saythis,"%s.original_zdc",TreeFileName);
  cout << "-- Added ZDC response branch as 'ZDC' (plus other SNM variables)" << endl;
  cout << "-- Wrote " << TreeFileName << endl;
  cout << "-- and moved the original file to " << saythis << endl;

  return;
}

