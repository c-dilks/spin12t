// prints run number list for pi0s using phiset/all.root
// (uses s0 g0 p0 e0 TObjArray)

void PrintRunNumberList(const char * infile = "../phiset/all.root")
{
  char outfile[64];
  strcpy(outfile,"RUNLIST");
  TFile * tf = new TFile(infile,"READ");
  //tf->cd("pi0");
  TObjArray * oa = (TObjArray*) tf->Get("pi0/phi_dist_pi0_s0_g0_p0_e0");
  Int_t runnum;
  for(Int_t i=0; i<oa->GetEntries(); i++)
  {
    sscanf(oa->At(i)->GetName(),"phi_pi0_s0_g0_p0_e0_r%d",&runnum);
    if(i==0) gSystem->RedirectOutput(outfile,"w");
    else gSystem->RedirectOutput(outfile,"a");
    printf("%d\n",runnum);
    gSystem->RedirectOutput(0);
  };
  printf("%s created\n",outfile);
}
