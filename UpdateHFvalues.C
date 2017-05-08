// This glauber afterburner module adds the HF 
// signal to glauber tree. 
//
// -J.Kamin
// 8 May 2017
//

void UpdateHFvalues(char TreeFileName[],char TreeObjectName[])
{

  char saythis[500];
  //sprintf(TreeFileName,"~/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  cout << TreeFileName << endl;

  TH1D *h_HFresponse[500];

  TFile *fincms = TFile::Open("./HFslices.root");
  for(int j=0; j<500; j++){
    sprintf(saythis,"h_HFresponse_%d",j);
    h_HFresponse[j] = (TH1D*)fincms->Get(saythis)->Clone(saythis);
  }

  sprintf(saythis,"cp %s %s.original_hf",TreeFileName,TreeFileName);
  gSystem->Exec(saythis);

  TFile *f = new TFile(TreeFileName,"update");
  TTree *T = (TTree*)f->Get(TreeObjectName);
  float HFresponse;
  float Npart;

  TBranch *bHF = T->Branch("HF",&HFresponse,"HF/F");
  T->SetBranchAddress("Npart",&Npart);
  Long64_t nentries = T->GetEntries();
  for(Long64_t i=0; i<nentries; i++){
    T->GetEntry(i);
    int iNpart = (int)Npart;
    HFresponse = h_HFresponse[iNpart]->GetRandom();
    if(i%int(0.1*nentries)==0)
      cout << "tree entry " << i << "/" << nentries << "    NPart: " << Npart << "    HF: " << HFresponse << endl;
    bHF->Fill();
  }
  //T->Print();
  T->Write();
  delete f;

  sprintf(saythis,"%s.original_hf",TreeFileName);
  cout << "-- Added HF response branch as 'HF'" << endl;
  cout << "-- Wrote " << TreeFileName << endl;
  cout << "-- and moved the original file to " << saythis << endl;

  return;
}

