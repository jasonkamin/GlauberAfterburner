double funcDensPath(double *x, double *par);

int GetNSlow_ALICE (int nColl); // algo = 1
int GetNSlow_ALICEp(int nColl); // algo = 2
int GetNSlow_LinSat(int nColl); // algo = 3
int GetNSlow_Ferenc(int nColl); // algo = 4
int GetNSlow_AveNColl_vB(int nColl, double b); // algo = 5
int GetNSlow_ALICE_noFerencElse(int nColl); // algo = 6
int GetNSlow_ALICE_OnlyFerencElse(int nColl); // algo = 7
int GetNSlow_ALICE_OnlyElse(int nColl); // algo = 8
int GetNSlow_Black_b(int nColl, double b); // algo = 9

void SetupPoissonSlices(int nc, TF1 *f1);
bool madepoissonslices = 0;

Int_t   Ngray;
Int_t   Nblack;
Int_t   NgrayHit;
Int_t   NblackHit;
Float_t SNenergy;
Float_t ZDCresponse;

Float_t meanNgrayp  = 0.0;
Float_t meanNblackp = 0.0;
Float_t meanNslowp  = 0.0;
Float_t N_LCF       = 0.0;
Float_t  meanNgray  = 0.0;
Float_t  meanNblack = 0.0;
Float_t  meanNslow  = 0.0;

Float_t  phiGray    = 0.0;
Float_t  theGray    = 0.0;
Float_t  phiBlack   = 0.0;
Float_t  theBlack   = 0.0;

TRandom3 randy;

static const int nSlices = 50;
TH1D *myPslices[nSlices];

//fittable parameters...
Double_t sigma_p   = 0.5;//0.0;
Double_t Nboverg   = 0.65;//0.65;
Double_t gm_par    = 1.71;//1.71;
Double_t alpha     = 0.48;
Double_t a_par     = 50.0;
Double_t b_par     = 230.0;
Double_t c_par     = 4.6;// alice uses 4.2 but we changed to 4.6 to enforce f(0)=0.
Double_t fracGneut = 0.1;// alice uses 0.1 but Ferenc thinks this should be more like 0.65
                         // which makes sense given that they use 0.65 for protons too. 
Double_t c1_par    = 1.47;
Double_t sat_par   = 30.4;


// some input variables to the toy MC...
// stolen from Ferenc's paper: 
// https://arxiv.org/pdf/hep-ph/0304065.pdf
Double_t betaGray  = 0.05;
Double_t betaBlack = 0.0;
Double_t E0Gray  = 0.050;// GeV
Double_t E0Black = 0.005;// GeV
Double_t pAveGray  = 0.271;// GeV/c
Double_t pAveBlack = 0.086;// GeV/c
Double_t m_neutron = 0.939565;//GeV

// fit results from neutron peaks in 8.16TeV pPb data. 
Double_t mu  = 67.43;
Double_t sig = 15.27;
Double_t ped = 0.6171;

// ZDC accpetance (roughly)
Double_t ZdcAccTheta = 0.4e-3;//0.4e-3
//Double_t ZdcAccTheta = 1e6;//inf


// dN/dtheta functional forms: 
TF1 *dN_dThe_black = new TF1("dN_dThe_black","exp( ([0]*[1]/[2]) * cos(x))",0.0,TMath::Pi());
TF1 *dN_dThe_gray  = new TF1("dN_dThe_gray", "exp( ([0]*[1]/[2]) * cos(x))",0.0,TMath::Pi());

TF1 *dN_dE_black  = new TF1("dN_dE_black", "exp(-x/[0])",0.0,20.0*E0Black);//GeV !
TF1 *dN_dE_gray  = new TF1("dN_dE_gray", "exp(-x/[0])",0.0,20.0*E0Gray);//GeV !

// rough apprxs for relationship btwn NColl and Ngray/black.
// taken from xyscan of plots in short ALICE Internal Note
TF1 *f_p1 = new TF1("f_p1","[0]*x",0,100);//[0] is the same parameter as f_sat par[2]
TF1 *f_sat = new TF1("f_sat","[0]*(2.0/(1+exp(-2*[1]*[2]*x/[0])) - 1)",0,100);//[2] is the same parameter as f_p1 par[0]

// taken from long ALICE centrality paper
TF1 *f_Ngrayp = new TF1("f_Ngrayp","[0]+[1]*x+[2]*x*x",0,100);//bascially linear.
TF1 *f_Nslown = new TF1("f_Nslown","[0]*x + ([1] - [2]/([3]+x))",0,200);//as a fcn of NCF

TF1 *f_rho  = new TF1("f_rho", "[0]* ( 1/(1+exp((x-[1])/[2])) )",0,20);

TH1D *Ncollb_pfx;

TF1 *f_DensPath = new TF1("f_DensPath", funcDensPath,0,20,3);

//##########################################################################
int GetNSlow_ALICE(int nColl) // algo = 1
{
  //complicated model from long ALICE centrality paper:

  meanNgrayp  = f_Ngrayp->Eval(nColl);
  meanNgrayp  = randy.Gaus(meanNgrayp,sigma_p);
  if(meanNgrayp<0.)  meanNgrayp = 0;

  
  meanNblackp = Nboverg*meanNgrayp;
  meanNslowp  = meanNgrayp + meanNblackp;
  N_LCF       = gm_par*meanNslowp;//in the code alice divides by a par called "gamma" but it may be ok anyway. 
  if(N_LCF>0.){
    meanNslow   = f_Nslown->Eval(N_LCF);
    meanNgray   = fracGneut*meanNslow;
    meanNblack  = meanNslow-meanNgray;
  }
  else{
    meanNblack  = 126.0/208.0 * 4.0 * nColl;
    meanNblack  = randy.Gaus(meanNblack,sigma_p);
    meanNgray   = meanNblack/9.0;
  }

  if(meanNblack<0.) meanNblack = 0;
  if(meanNgray <0.) meanNgray  = 0;

  Nblack = randy.Binomial(126, meanNblack/126.0);
  Ngray  = randy.Binomial(126-Nblack, meanNgray/(126.0-Nblack));
  if(Nblack<0.) Nblack = 0;
  if(Ngray <0.) Ngray  = 0;

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_ALICEp(int nColl) // algo = 2
{
  //complicated model from long ALICE centrality paper (w/ poiss smearing):

  meanNgrayp  = f_Ngrayp->Eval(nColl);

  if(sigma_p>0){
    if(madepoissonslices==0)  SetupPoissonSlices(nSlices, f_Ngrayp);
    if(nColl>=nSlices){
      cout << "you don't have enough poisson slices !!  breaking !! " << endl;
      return 0;
    }
    double tempmu = (meanNgrayp/sigma_p)*(meanNgrayp/sigma_p);
    meanNgrayp  = myPslices[nColl]->GetRandom()*(meanNgrayp/tempmu);
  }
  if(meanNgrayp<0.)  meanNgrayp = 0;
  
  meanNblackp = Nboverg*meanNgrayp;
  meanNslowp  = meanNgrayp + meanNblackp;
  N_LCF       = gm_par*meanNslowp;
  meanNslow   = f_Nslown->Eval(N_LCF);
  if(meanNslow<0.) meanNslow = 0;
  meanNgray   = fracGneut*meanNslow;
  meanNblack  = meanNslow-meanNgray;

  Ngray  = randy.Binomial(126, meanNgray/126.0);
  Nblack = randy.Binomial(126-Ngray, meanNblack/(126.0-Ngray));

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_LinSat(int nColl) // algo = 3
{

  meanNgrayp  = f_p1 ->Eval(nColl);
  meanNblackp = f_sat->Eval(nColl);
  //meanNgrayp  = randy.Gaus(meanNgrayp, sigma_p);
  //if(meanNgrayp<0.)   meanNgrayp  = 0;
  //if(meanNblackp<0.)  meanNblackp = 0;
  //meanNblackp = f_sat->Eval(meanNgrayp/f_p1->GetParameter(0));
  //meanNblackp = randy.Gaus(meanNblackp,sigma_p);
  //if(meanNgrayp<0.)   meanNgrayp  = 0;
  //if(meanNblackp<0.)  meanNblackp = 0;


  if(sigma_p>0){
    if(madepoissonslices==0)  SetupPoissonSlices(nSlices, f_p1);
    if(nColl>=nSlices){
      cout << "you don't have enough poisson slices !!  breaking !! " << endl;
      return 0;
    }

    double tempmu = (meanNgrayp/sigma_p)*(meanNgrayp/sigma_p);
    meanNgrayp  = myPslices[nColl]->GetRandom()*(meanNgrayp/tempmu);

    meanNblackp = myPslices[nColl]->GetRandom()*(f_Ngrayp->Eval(nColl)/tempmu);
    meanNblackp = f_sat->Eval(meanNblackp/f_p1->GetParameter(0));
  }

  //meanNblackp = Nboverg*meanNgrayp;
  if(meanNgrayp<0.)   meanNgrayp  = 0;
  if(meanNblackp<0.)  meanNblackp = 0;
  meanNslowp  = meanNgrayp + meanNblackp;

  meanNgray   = meanNgrayp;
  meanNblack  = meanNblackp;
  meanNslow   = meanNgray+meanNblack;

  Ngray  = randy.Binomial(126, meanNgray/126.0);
  Nblack = randy.Binomial(126-Ngray, meanNblack/(126.0-Ngray));

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_AveNColl_vB(int nColl, double b) // algo = 4
{


  int binnum = Ncollb_pfx->GetXaxis()->FindBin(b);
  nColl = Ncollb_pfx->GetBinContent(binnum);

  meanNgray  = f_p1 ->Eval(nColl);// 1.4
  meanNblack = f_sat->Eval(nColl);// 12.5, 1.375, 1.47


  //Ngray  = randy.Binomial(126, meanNgray/126.0);
  //Nblack = randy.Binomial(126-Ngray, meanNblack/(126.0-Ngray));
  Ngray  = randy.Binomial(126, meanNgray/126.0);
  Nblack = randy.Binomial(126, meanNblack/126.0);

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_Ferenc(int nColl) // algo = 5
{
  //Double_t sigma_p = 2.0;


  // rough apprx for relationship btwn NColl and Ngray/black.
  meanNgray  = 2.0*nColl;
  meanNblack = 4.0*nColl;

  // trying out gaus vs poisson vs binomial
  //Ngray  = randy.Gaus(meanNgray ,sqrt(meanNgray) );  if(Ngray <0) Ngray =0;
  //Nblack = randy.Gaus(meanNblack,sqrt(meanNblack));  if(Nblack<0) Nblack=0;
  //Ngray  = randy.Poisson(meanNgray);
  //Nblack = randy.Poisson(meanNblack);
  Ngray  = randy.Binomial(126, meanNgray/126.0);
  Nblack = randy.Binomial(126, meanNblack/126.0);
  //Ngray  = 0;
  //Nblack = randy.Binomial(126, meanNslow/126.0);


  return Ngray+Nblack;
}//#########################################################################

//##########################################################################
int GetNSlow_ALICE_noFerencElse(int nColl) // algo = 6
{
  //complicated model from long ALICE centrality paper:

  meanNgrayp  = f_Ngrayp->Eval(nColl);
  meanNgrayp  = randy.Gaus(meanNgrayp,sigma_p);
  if(meanNgrayp<0)  meanNgrayp = 0;

  
  meanNblackp = Nboverg*meanNgrayp;
  meanNslowp  = meanNgrayp + meanNblackp;
  N_LCF       = gm_par*meanNslowp;//in the code alice divides by gamma but it may be ok anyway. 
  if(N_LCF>0.){
    meanNslow   = f_Nslown->Eval(N_LCF);
    meanNgray   = fracGneut*meanNslow;
    meanNblack  = meanNslow-meanNgray;
  }
  else{
    meanNblack  = 0;
    meanNgray   = 0;
  }

  if(meanNblack<0.) meanNblack = 0;
  if(meanNgray <0.) meanNgray  = 0;

  Nblack = randy.Binomial(126, meanNblack/126.0);
  Ngray  = randy.Binomial(126-Nblack, meanNgray/(126.0-Nblack));
  if(Nblack<0.) Nblack = 0;
  if(Ngray <0.) Ngray  = 0;

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_ALICE_OnlyFerencElse(int nColl) // algo = 7
{
  //complicated model from long ALICE centrality paper:

  meanNgrayp  = f_Ngrayp->Eval(nColl);
  meanNgrayp  = randy.Gaus(meanNgrayp,sigma_p);
  if(meanNgrayp<0)  meanNgrayp = 0;

  
  meanNblackp = Nboverg*meanNgrayp;
  meanNslowp  = meanNgrayp + meanNblackp;
  N_LCF       = gm_par*meanNslowp;//in the code alice divides by gamma but it may be ok anyway. 
  if(N_LCF>0.){
    meanNblack  = 0;
    meanNgray   = 0;
  }
  else{
    meanNblack  = 126.0/208.0 * 4.0 * nColl;
    meanNblack  = randy.Gaus(meanNblack,sigma_p);
    meanNgray   = meanNblack/9.0;
  }

  if(meanNblack<0.) meanNblack = 0;
  if(meanNgray <0.) meanNgray  = 0;

  Nblack = randy.Binomial(126, meanNblack/126.0);
  Ngray  = randy.Binomial(126-Nblack, meanNgray/(126.0-Nblack));
  if(Nblack<0.) Nblack = 0;
  if(Ngray <0.) Ngray  = 0;

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_ALICE_OnlyElse(int nColl) // algo = 8
{
  // the infamous 'else' statement from the model mixing
  // which kicks in when the smearing is less than 0. 

  meanNblack  = 126.0/208.0 * 4.0 * nColl;
  meanNblack  = randy.Gaus(meanNblack,sigma_p);
  meanNgray   = meanNblack/9.0;

  if(meanNblack<0.) meanNblack = 0;
  if(meanNgray <0.) meanNgray  = 0;

  Nblack = randy.Binomial(126, meanNblack/126.0);
  Ngray  = randy.Binomial(126-Nblack, meanNgray/(126.0-Nblack));
  if(Nblack<0.) Nblack = 0;
  if(Ngray <0.) Ngray  = 0;

  return Ngray+Nblack;
}//#########################################################################


//##########################################################################
int GetNSlow_Black_b(int nColl, double b) // algo = 9
{

  meanNgrayp  = f_p1 ->Eval(nColl);
  //meanNblackp = f_sat->Eval(nColl);
  //meanNblackp = f_rho->Eval(b);
  meanNblackp = f_DensPath->Eval(b);
  //meanNblackp = -1.25*b + 12.5;// didn't seem to work.

  if(sigma_p>0){
    if(madepoissonslices==0)  SetupPoissonSlices(nSlices, f_p1);
    if(nColl>=nSlices){
      cout << "you don't have enough poisson slices !!  breaking !! " << endl;
      return 0;
    }

    double tempmu = (meanNgrayp/sigma_p)*(meanNgrayp/sigma_p);
    meanNgrayp  = myPslices[nColl]->GetRandom()*(meanNgrayp/tempmu);

    //meanNblackp = myPslices[nColl]->GetRandom()*(f_Ngrayp->Eval(nColl)/tempmu);
    //meanNblackp = f_sat->Eval(meanNblackp/f_p1->GetParameter(0));
  }

  //meanNblackp = Nboverg*meanNgrayp;
  if(meanNgrayp<0.)   meanNgrayp  = 0;
  if(meanNblackp<0.)  meanNblackp = 0;
  meanNslowp  = meanNgrayp + meanNblackp;

  meanNgray   = meanNgrayp;
  meanNblack  = meanNblackp;
  meanNslow   = meanNgray+meanNblack;

  Ngray  = randy.Binomial(126, meanNgray/126.0);
  Nblack = randy.Binomial(126-Ngray, meanNblack/(126.0-Ngray));

  return Ngray+Nblack;
}//#########################################################################




//##########################################################################
void SetupPoissonSlices(int nS=10, TF1 *f1=0)
{

  cout << "... Setting up Poisson distributions for [a,b,c] = [" << a_par <<","<< b_par <<","<< c_par << "] and sigma_p = " << sigma_p << "  ..." << endl;
  char saythis[500];

  //TF1 *f_mypois = new TF1("f_mypois","[0]*TMath::Poisson(x,[1])",0,(0.6*nS/sigma_p)*(0.6*nS/sigma_p)+5.0*(0.6*nS/sigma_p));
  TF1 *f_mypois = new TF1("f_mypois","[0]*TMath::Poisson(x,[1])",0,(f1->Eval(nS)/sigma_p)*(f1->Eval(nS)/sigma_p)+5.0*(f1->Eval(nS)/sigma_p));
  f_mypois->SetParameter(0,100);
  f_mypois->SetParameter(1,5);
  
  Float_t meanNgrayp  = 0.0;

  for(int nc=0; nc<nS; nc++){

    delete myPslices[nc];

    //meanNgrayp  = f_Ngrayp->Eval(nc);
    meanNgrayp  = f1->Eval(nc);
    if(meanNgrayp<0) meanNgrayp = 0;

    double tempmu = (meanNgrayp/sigma_p)*(meanNgrayp/sigma_p);
    cout << nc << "  " << meanNgrayp << "  " << tempmu << endl;
    f_mypois->SetParameter(1,tempmu);

    //sprintf(saythis,"myPslices_%d",nc);
    //cout << "histo limits: central(" << (0.6*nS/sigma_p)*(0.6*nS/sigma_p) << ") + 5sig(" << (0.6*nS/sigma_p)*(0.6*nS/sigma_p)+5.0*(0.6*nS/sigma_p) << ")" << endl;
    //myPslices[nc] = new TH1D(saythis,saythis,50000,0,(0.6*nS/sigma_p)*(0.6*nS/sigma_p)+5.0*(0.6*nS/sigma_p));
    //myPslices[nc]->FillRandom("f_mypois",1000000);

    sprintf(saythis,"myPslices_%d",nc);
    cout << "histo limits: central(" << (f1->Eval(nc)/sigma_p)*(f1->Eval(nc)/sigma_p) << ") + 5sig(" << 5.0*(f1->Eval(nc)/sigma_p) << ") = " << (f1->Eval(nc)/sigma_p)*(f1->Eval(nc)/sigma_p)+5.0*(f1->Eval(nc)/sigma_p) << endl;
    myPslices[nc] = new TH1D(saythis,saythis,50000,0,(f1->Eval(nc)/sigma_p)*(f1->Eval(nc)/sigma_p)+5.0*(f1->Eval(nc)/sigma_p));
    myPslices[nc]->FillRandom("f_mypois",1000000);

  }

  cout << "... done (nSlices=" << nS << ") ! " << endl << endl;
  madepoissonslices = 1;

  return;
}//#########################################################################


double funcDensPath(double *x, double *par)
{

  double value = 0.0;

  if(x[0]<par[1])
    value = par[0]* ( 1/(1+exp((x[0]-par[1])/par[2])) ) *(TMath::Sqrt(par[1]*par[1]-x[0]*x[0]));
  else
    value = par[0]* ( 1/(1+exp((x[0]-par[1])/par[2])) );

  return value;
}

