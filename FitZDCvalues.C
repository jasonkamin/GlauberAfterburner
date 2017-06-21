// An example for adding a ZDC signal to 
// the glauber output.
//
// -J.Kamin
// 8 May 2017
//

TH1D *h_KoloHist[1000];
#include "/Users/jason/Dropbox/work/code/ZJet/CompatibilityTests/CompatibilityHelperFunctions.C"
// https://github.com/jasonkamin/CompatibilityTests/blob/master/CompatibilityHelperFunctions.C


//Double_t init_sigma_p =0.00  /*0.0  */  /*0.00 */;  Double_t step_sigma_p = 0.1 ;  int nStep_sigma_p =  1;//5;
//Double_t init_Nboverg =2.00  /*7.4  */  /*2.00 */;  Double_t step_Nboverg = 0.05 ;  int nStep_Nboverg = 1;//3;
//Double_t init_gm_par  =1.71  /*1.71 */  /*1.71 */;  Double_t step_gm_par  = 0.1  ;  int nStep_gm_par  = 1;//3;
//Double_t init_alpha   =0.48  /*0.48 */  /*0.48 */;  Double_t step_alpha   = 0.05 ;  int nStep_alpha   = 1;//4;//1;
//Double_t init_a_par   =50.0  /*50.0 */  /*50.0 */;  Double_t step_a_par   = 2.5  ;  int nStep_a_par   = 1;//4;//1;
//Double_t init_b_par   =230.0 /*230.0*/  /*230.0*/;  Double_t step_b_par   = 50.0 ;  int nStep_b_par   = 1;//1;
//Double_t init_c_par   =4.2   /*4.2  */  /*4.2  */;  Double_t step_c_par   = 0.8  ;  int nStep_c_par   = 1;//5;
//Double_t init_fGn     =0.1   /*0.1  */  /*0.1  */;  Double_t step_fGn     = 0.1  ;  int nStep_fGn     = 1; 
//Double_t init_c1      =0.50  /*1.50 */  /*0.40 */;  Double_t step_c1      = 0.1  ;  int nStep_c1      = 50; 
//Double_t init_sat     =30.4  /*30.4 */  /*30.4 */;  Double_t step_sat     = 30.4 ;  int nStep_sat     = 1; 

Double_t init_sigma_p =0.00  /*0.0  */  /*0.00 */;  Double_t step_sigma_p = 0.5 ;  int nStep_sigma_p =  4;//5;
Double_t init_Nboverg =0.00  /*7.4  */  /*2.00 */;  Double_t step_Nboverg = 0.05 ;  int nStep_Nboverg = 1;//3;
Double_t init_gm_par  =0.00  /*1.71 */  /*1.71 */;  Double_t step_gm_par  = 0.1  ;  int nStep_gm_par  = 1;//3;
Double_t init_alpha   =0.00  /*0.48 */  /*0.48 */;  Double_t step_alpha   = 0.05 ;  int nStep_alpha   = 1;//4;//1;
Double_t init_a_par   =0.00  /*50.0 */  /*50.0 */;  Double_t step_a_par   = 2.5  ;  int nStep_a_par   = 1;//4;//1;
Double_t init_b_par   =0.00  /*230.0*/  /*230.0*/;  Double_t step_b_par   = 50.0 ;  int nStep_b_par   = 1;//1;
Double_t init_c_par   =0.00  /*4.2  */  /*4.2  */;  Double_t step_c_par   = 0.8  ;  int nStep_c_par   = 1;//5;
Double_t init_fGn     =0.00  /*0.1  */  /*0.1  */;  Double_t step_fGn     = 0.1  ;  int nStep_fGn     = 1; 
Double_t init_c1      =0.50  /*1.50 */  /*0.40 */;  Double_t step_c1      = 0.3  ;  int nStep_c1      = 10; 
Double_t init_sat     =12.0  /*30.4 */  /*30.4 */;  Double_t step_sat     = 3.0  ;  int nStep_sat     = 5; 

TH1D *h_zdc;
TH1D *h_zdctest;
Double_t Best_chi2_sigma_p = 1.5;//0.0;
Double_t Best_chi2_Nboverg = 0.65;
Double_t Best_chi2_gm_par  = 1.71;
Double_t Best_chi2_alpha   = 0.48;
Double_t Best_chi2_a_par   = 50.0;
Double_t Best_chi2_b_par   = 230.0;
Double_t Best_chi2_c_par   = 4.6;
Double_t Best_chi2_fGn     = 0.1  ;
Double_t Best_chi2_c1      = 1.47 ;
Double_t Best_chi2_sat     = 30.4 ;

Double_t Best_kolo_sigma_p = 0.0;
Double_t Best_kolo_Nboverg = 0.65;
Double_t Best_kolo_gm_par  = 1.71;
Double_t Best_kolo_alpha   = 0.48;
Double_t Best_kolo_a_par   = 50.0;
Double_t Best_kolo_b_par   = 230.0;
Double_t Best_kolo_c_par   = 4.6;
Double_t Best_kolo_fGn     = 0.1  ;
Double_t Best_kolo_c1      = 1.47 ;
Double_t Best_kolo_sat     = 30.4 ;

Double_t mychi2 = 1.1e10;
Double_t mykolo = 1.1e10;
Double_t Best_mychi2 = 1.0e10;
Double_t Best_mykolo = 1.0e10;

TH1D *h_chi2;
TH1D *h_kolo;

TH2D *h2_chi2;
TH2D *h2_kolo;
TH2D *h2_chi2kolo;
TH2D *h2_zdc;

TH1D *h1MYPARS;
TH2D *h2MYPARS;


void FitZDCvalues(char TreeFileName[], char TreeObjectName[], Int_t algo=1)
{

  randy.SetSeed(5281978);
  char saythis[500];
  char saythis1[500];
  cout << TreeFileName << endl;

  TFile *filed = TFile::Open("ZDCDATA.root");
  h_zdc = (TH1D*)filed->Get("h_zdc")->Clone("h_zdc");
  h_zdctest = (TH1D*)h_zdc->Clone("h_zdctest");

  TFile *f_NCVB = TFile::Open("NColl_vs_b.root");
  Ncollb_pfx = (TH1D*)f_NCVB->Get("Ncollb_pfx")->Clone("Ncollb_pfx");


  h_chi2 = new TH1D("h_chi2","h_chi2",1000,0,500000);
  h_kolo = new TH1D("h_kolo","h_kolo",1000,0,5);

  h2_chi2 = new TH2D("h2_chi2","h2_chi2",1000,0,500000, 5001,-0.5,5000.5);
  h2_kolo = new TH2D("h2_kolo","h2_kolo",1000,0,5,      5001,-0.5,5000.5);
  h2_chi2kolo = new TH2D("h2_chi2kolo","h2_chi2kolo", 1000,0,500000, 1000,0,5);
  h2_zdc  = new TH2D("h2_zdc","h2_zdc",h_zdc->GetNbinsX(),h_zdc->GetXaxis()->GetBinLowEdge(1),h_zdc->GetXaxis()->GetBinUpEdge(h_zdc->GetNbinsX()), 5001,-0.5,5000.5);

  h1MYPARS = new TH1D("h1MYPARS","h1MYPARS",20,0.5,20.5);//npars + 3. skip a bin and last 2 bins are chi2, kolo.
  h2MYPARS = new TH2D("h2MYPARS","h2MYPARS",20,0.5,20.5, 5001,-0.5,5000.5);//npars + 3. skip a bin and last 2 bins are chi2, kolo.

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

  int nTrials = 0;
  for(int iV_sigma_p=0; iV_sigma_p<nStep_sigma_p; iV_sigma_p++){
    for(int iV_Nboverg=0; iV_Nboverg<nStep_Nboverg; iV_Nboverg++){
      for(int iV_gm_par =0; iV_gm_par <nStep_gm_par ; iV_gm_par ++){
        for(int iV_alpha  =0; iV_alpha  <nStep_alpha  ; iV_alpha  ++){
          for(int iV_a_par  =0; iV_a_par  <nStep_a_par  ; iV_a_par  ++){
            for(int iV_b_par  =0; iV_b_par  <nStep_b_par  ; iV_b_par  ++){
              for(int iV_c_par  =0; iV_c_par  <nStep_c_par  ; iV_c_par  ++){
                for(int iV_fGn    =0; iV_fGn    <nStep_fGn    ; iV_fGn  ++){
                  for(int iV_c1     =0; iV_c1     <nStep_c1     ; iV_c1  ++){
                    for(int iV_sat    =0; iV_sat    <nStep_sat    ; iV_sat  ++){


                      if(iV_sigma_p==0) sigma_p  = init_sigma_p;
                      if(iV_Nboverg==0) Nboverg  = init_Nboverg;
                      if(iV_gm_par ==0) gm_par   = init_gm_par ;
                      if(iV_alpha  ==0) alpha    = init_alpha  ;
                      if(iV_a_par  ==0) a_par    = init_a_par  ;
                      if(iV_b_par  ==0) b_par    = init_b_par  ;
                      if(iV_c_par  ==0) c_par    = init_c_par  ;
                      if(iV_fGn    ==0) fracGneut= init_fGn    ;
                      if(iV_c1     ==0) c1_par   = init_c1     ;
                      if(iV_sat    ==0) sat_par  = init_sat    ;


                      cout << nTrials << " " << "####- Attempting Parameters -########" << endl;
                      cout << "    sigma_p " << sigma_p   << endl;
                      cout << "    Nboverg " << Nboverg   << endl;
                      cout << "    gm_par  " << gm_par    << endl;
                      cout << "    alpha   " << alpha     << endl;
                      cout << "    a_par   " << a_par     << endl;
                      cout << "    b_par   " << b_par     << endl;
                      cout << "    c_par   " << c_par     << endl;
                      cout << "    fGn     " << fracGneut << endl;
                      cout << "    c1      " << c1_par    << endl;
                      cout << "    sat     " << sat_par   << endl;
                      cout << "#########################################" << endl;


                      h_zdctest->Reset();

                      f_Nslown     ->SetParameters(alpha,a_par,b_par,c_par);
                      f_p1         ->SetParameter(0,c1_par);
                      f_sat        ->SetParameters(sat_par, Nboverg, c1_par);
                      f_DensPath   ->SetParameters(sat_par/6.62,6.62,0.546);


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
                        for(Long64_t j =0; j<NgrayHit; j++){
                          Double_t resp = gRandom->Gaus(mu,sig)+ped;
                          if(resp>0)
                            ZDCresponse += resp;
                        }
                        for(Long64_t j =0; j<NblackHit; j++){
                          Double_t resp = gRandom->Gaus(mu,sig)+ped;
                          if(resp>0)
                            ZDCresponse += resp;
                        }

                        h_zdctest->Fill(ZDCresponse);
                        h2_zdc   ->Fill(ZDCresponse,nTrials);

                        // have a look at a couple values...
                        if(i%int(0.2*nentries)==0){
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

                        //bphiGray ->Fill();
                        //btheGray ->Fill();
                        //bphiBlack->Fill();
                        //btheBlack->Fill();

                        //bmeanNgrayp ->Fill();
                        //bmeanNblackp->Fill();
                        //bmeanNslowp ->Fill();
                        //bN_LCF      ->Fill();
                        //bmeanNgray ->Fill();
                        //bmeanNblack->Fill();
                        //bNgray     ->Fill();
                        //bNblack    ->Fill();
                        //bNgrayHit  ->Fill();
                        //bNblackHit ->Fill();
                        //bSNenergy  ->Fill();
                        //bZDC       ->Fill();

                      }
                      if(nStep_sigma_p>1)
                        madepoissonslices = 0;

                      for(int i=1; i<h_zdctest->GetNbinsX()+1; i++){
                        h_zdctest->SetBinError(i,TMath::Sqrt(h_zdctest->GetBinContent(i)));
                      }
                      h_zdctest->Scale(h_zdc->Integral(-1,-1)/h_zdctest->Integral(-1,-1));

                      mychi2 = CalcChiSq (h_zdc, h_zdctest, -1, -1);
                      mykolo = CalcKoloSm(h_zdc, h_zdctest, -1, -1, 1, -1);
                      if(mychi2<Best_mychi2){
                        Best_mychi2       = mychi2;

                        Best_chi2_sigma_p = sigma_p;
                        Best_chi2_Nboverg = Nboverg;
                        Best_chi2_gm_par  = gm_par ;
                        Best_chi2_alpha   = alpha  ;
                        Best_chi2_a_par   = a_par  ;
                        Best_chi2_b_par   = b_par  ;
                        Best_chi2_c_par   = c_par  ;
                        Best_chi2_fGn     = fracGneut;
                        Best_chi2_c1      = c1_par   ;
                        Best_chi2_sat     = sat_par  ;
                      }
                      if(mykolo<Best_mykolo){
                        Best_mykolo       = mykolo;

                        Best_kolo_sigma_p = sigma_p;
                        Best_kolo_Nboverg = Nboverg;
                        Best_kolo_gm_par  = gm_par ;
                        Best_kolo_alpha   = alpha  ;
                        Best_kolo_a_par   = a_par  ;
                        Best_kolo_b_par   = b_par  ;
                        Best_kolo_c_par   = c_par  ;
                        Best_kolo_fGn     = fracGneut;
                        Best_kolo_c1      = c1_par   ;
                        Best_kolo_sat     = sat_par  ;
                      }
                      cout << endl;
                      cout << "chi2: " << mychi2 << endl;
                      cout << "kolo: " << mykolo << endl;
                      cout << endl;

                      h1MYPARS->SetBinContent(1, sigma_p);    h1MYPARS->GetXaxis()->SetBinLabel(1, "sigma_p  ");
                      h1MYPARS->SetBinContent(2, Nboverg);    h1MYPARS->GetXaxis()->SetBinLabel(2, "Nboverg  ");
                      h1MYPARS->SetBinContent(3, gm_par );    h1MYPARS->GetXaxis()->SetBinLabel(3, "gm_par   ");
                      h1MYPARS->SetBinContent(4, alpha  );    h1MYPARS->GetXaxis()->SetBinLabel(4, "alpha    ");
                      h1MYPARS->SetBinContent(5, a_par  );    h1MYPARS->GetXaxis()->SetBinLabel(5, "a_par    ");
                      h1MYPARS->SetBinContent(6, b_par  );    h1MYPARS->GetXaxis()->SetBinLabel(6, "b_par    ");
                      h1MYPARS->SetBinContent(7, c_par  );    h1MYPARS->GetXaxis()->SetBinLabel(7, "c_par    ");
                      h1MYPARS->SetBinContent(8, fracGneut);  h1MYPARS->GetXaxis()->SetBinLabel(8, "fracGneut");
                      h1MYPARS->SetBinContent(9, c1_par   );  h1MYPARS->GetXaxis()->SetBinLabel(9, "c1_par   ");
                      h1MYPARS->SetBinContent(10,sat_par  );  h1MYPARS->GetXaxis()->SetBinLabel(10,"sat_par  ");
                      h1MYPARS->SetBinContent(11,algo     );  h1MYPARS->GetXaxis()->SetBinLabel(11,"algo     ");
                      h1MYPARS->SetBinContent(19, mychi2);
                      h1MYPARS->SetBinContent(20,mykolo);

                      h2MYPARS->SetBinContent(1,  nTrials, sigma_p  );
                      h2MYPARS->SetBinContent(2,  nTrials, Nboverg  );
                      h2MYPARS->SetBinContent(3,  nTrials, gm_par   );
                      h2MYPARS->SetBinContent(4,  nTrials, alpha    );
                      h2MYPARS->SetBinContent(5,  nTrials, a_par    );
                      h2MYPARS->SetBinContent(6,  nTrials, b_par    );
                      h2MYPARS->SetBinContent(7,  nTrials, c_par    );
                      h2MYPARS->SetBinContent(8,  nTrials, fracGneut);
                      h2MYPARS->SetBinContent(9,  nTrials, c1_par   );
                      h2MYPARS->SetBinContent(10, nTrials, sat_par  );
                      h2MYPARS->SetBinContent(11, nTrials, algo     );
                      h2MYPARS->SetBinContent(19, nTrials, mychi2);
                      h2MYPARS->SetBinContent(20, nTrials, mykolo);

                      h_chi2  ->Fill(mychi2);
                      h_kolo  ->Fill(mykolo);
                      h2_chi2 ->Fill(mychi2, nTrials);
                      h2_kolo ->Fill(mykolo, nTrials);
                      h2_chi2kolo ->Fill(mychi2,mykolo);
                      nTrials++;

                      //T->Print();
                      //sprintf(saythis1,"pars_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%s",sigma_p, Nboverg, gm_par , alpha  , a_par  , b_par  , c_par, TreeFileName);
                      //TFile *fout = new TFile(saythis1,"NEW");
                      //T->Write();
                      //h1MYPARS->Write();
                      //fout->Close();
                      //delete f;



                      // back up the input file... 
                      //sprintf(saythis1,"pars_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f__%s",sigma_p, Nboverg, gm_par , alpha  , a_par  , b_par  , c_par, TreeFileName);
                      //sprintf(saythis,"mv %s %s",TreeFileName,saythis1);
                      //gSystem->Exec(saythis);
                      //// back up the input file... 
                      //sprintf(saythis,"cp ../%s .",TreeFileName);
                      //gSystem->Exec(saythis);



                      sat_par   += step_sat  ;
                    }
                    c1_par    += step_c1  ;
                  }
                  fracGneut += step_fGn  ;
                }
                c_par   += step_c_par  ;
              }
              b_par   += step_b_par  ;
            }
            a_par   += step_a_par  ;
          }
          alpha   += step_alpha  ;
        }
        gm_par  += step_gm_par ;
      }
      Nboverg += step_Nboverg;
    }
    sigma_p += step_sigma_p;
  }

  cout << "Best_mychi2       " << Best_mychi2       << endl;
  cout << "Best_chi2_sigma_p " << Best_chi2_sigma_p << endl;
  cout << "Best_chi2_Nboverg " << Best_chi2_Nboverg << endl;
  cout << "Best_chi2_gm_par  " << Best_chi2_gm_par  << endl;
  cout << "Best_chi2_alpha   " << Best_chi2_alpha   << endl;
  cout << "Best_chi2_a_par   " << Best_chi2_a_par   << endl;
  cout << "Best_chi2_b_par   " << Best_chi2_b_par   << endl;
  cout << "Best_chi2_c_par   " << Best_chi2_c_par   << endl;
  cout << "Best_chi2_fGn     " << Best_chi2_fGn     << endl;
  cout << "Best_chi2_c1      " << Best_chi2_c1      << endl;
  cout << "Best_chi2_sat     " << Best_chi2_sat     << endl << endl;

  cout << "Best_mykolo       " << Best_mykolo       << endl;
  cout << "Best_kolo_sigma_p " << Best_kolo_sigma_p << endl;
  cout << "Best_kolo_Nboverg " << Best_kolo_Nboverg << endl;
  cout << "Best_kolo_gm_par  " << Best_kolo_gm_par  << endl;
  cout << "Best_kolo_alpha   " << Best_kolo_alpha   << endl;
  cout << "Best_kolo_a_par   " << Best_kolo_a_par   << endl;
  cout << "Best_kolo_b_par   " << Best_kolo_b_par   << endl;
  cout << "Best_kolo_c_par   " << Best_kolo_c_par   << endl;
  cout << "Best_kolo_fGn     " << Best_kolo_fGn     << endl;
  cout << "Best_kolo_c1      " << Best_kolo_c1      << endl;
  cout << "Best_kolo_sat     " << Best_kolo_sat     << endl << endl;

  TFile *fcompat = new TFile("compat.root","NEW");
  h_chi2      ->Write();
  h_kolo      ->Write();
  h2_chi2     ->Write();
  h2_kolo     ->Write();
  h2_chi2kolo ->Write();
  h2MYPARS    ->Write();
  h2_zdc      ->Write();
  fcompat->Close();

  sprintf(saythis,"%s.original_zdc",TreeFileName);
  cout << "-- Added ZDC response branch as 'ZDC' (plus other SNM variables)" << endl;
  cout << "-- Wrote " << TreeFileName << endl;
  cout << "-- and moved the original file to " << saythis << endl;

  return;
}


