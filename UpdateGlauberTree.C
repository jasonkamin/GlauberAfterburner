// This is simply a controller for 

#include "UpdateHFvalues.C"
#include "UpdateZDCvalues.C"

void UpdateGlauberTree(char treeinfile[],char treeobjname[])
{

  //sprintf(treeinfile,"./Jason_sigNN70.0_dave0.8_dsig0.00_dmin0.40_normaldensitypars_HCP_10000.root");

  UpdateHFvalues(treeinfile,treeobjname);
  UpdateZDCvalues(treeinfile,treeobjname);



}

