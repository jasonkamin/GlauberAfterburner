// An example for adding a ZDC signal to 
// the glauber output.
//
// -J.Kamin
// 8 May 2017
//

void UpdateZDCvalues(char TreeFileName[], char TreeObjectName[])
{

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

  // my new ZDC branches... 
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
  Double_t meanNgray  = 0.0;
  Double_t meanNblack = 0.0;
  Double_t betaGray  = 0.05;
  Double_t betaBlack = 0.0;
  Double_t E0Gray  = 0.050;// GeV
  Double_t E0Black = 0.005;// GeV
  Double_t pAveGray  = 0.750;// GeV/c
  Double_t pAveBlack = 0.200;// GeV/c

  // dN/dtheta functional forms: 
  TF1 *dN_dThe_gray  = new TF1("dN_dThe_gray", "exp( ([0]*[1]/[2]) * cos(x))",0.0,TMath::Pi());
  TF1 *dN_dThe_black = new TF1("dN_dThe_black","exp( ([0]*[1]/[2]) * cos(x))",0.0,TMath::Pi());
  dN_dThe_gray ->SetParameters(betaGray, pAveGray, E0Gray);
  dN_dThe_black->SetParameters(betaBlack,pAveBlack,E0Black);

  // fit results from neutron peaks in 8.16TeV pPb data. 
  Double_t mu  = 67.43;
  Double_t sig = 15.27;

  // just a guess... (consider it a placeholder).
  Double_t ZdcAccTheta = 0.4;

  // Loop over the tree ! 
  Long64_t nentries = T->GetEntries();
  for(Long64_t i=0; i<nentries; i++){
    T->GetEntry(i);

    // rough apprx for relationship btwn NColl and Ngray/black.
    meanNgray  = 2.0*Ncoll;
    meanNblack = 4.0*Ncoll;

    // these should be binomial but poisson is close enough for the moment.
    Ngray  = gRandom->Poisson(meanNgray);
    Nblack = gRandom->Poisson(meanNblack);

    NgrayHit = 0;
    SNenergy = 0.0;
    ZDCresponse = 0.0;

    // Check which slow neutrons hit the ZDC:
    for(Long64_t j=0; j<Ngray; j++){
      Double_t the = dN_dThe_gray->GetRandom();
      //Double_t phi = gRandom->Rndm(-0.5*TMath::Pi(),0.5*TMath::Pi());//not yet used.
      if(the<ZdcAccTheta){
        NgrayHit++;
        SNenergy += E0Gray;
      }
    }

    // Check which slow neutrons hit the ZDC:
    NblackHit = 0;
    for(Long64_t j=0; j<Nblack; j++){
      Double_t the = dN_dThe_black->GetRandom();
      //Double_t phi = gRandom->Rndm(-0.5*TMath::Pi(),0.5*TMath::Pi());//not yet used.
      if(the<ZdcAccTheta){
        NblackHit++;
        SNenergy += E0Black;
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
        ZDCresponse += 0.1*resp;
    }


    // have a look at a couple values...
    if(i%int(0.1*nentries)==0){
      cout << "tree entry " << i << "/" << nentries << "    NColl: " << Ncoll << endl;
      cout << "   Ngray: " << Ngray << "   Nblack: " << Nblack << endl;
      cout << "   NGhit: " << NgrayHit << "   NBlhit: " << NblackHit << endl;
      cout << "   SNenergy: " << SNenergy << "    ZDC:" << ZDCresponse << endl;
    }

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

