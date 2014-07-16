// prints list of runs which have consistent rellum measurements

void PrintConsistentRuns()
{
  TFile * infile = new TFile("rtree.root","READ");
  TTree * rellum = (TTree*) infile->Get("rellum");

  Int_t runnum;
  Bool_t isConsistent;

  rellum->SetBranchAddress("runnum",&runnum);
  rellum->SetBranchAddress("isConsistent",&isConsistent);

  gROOT->ProcessLine(".!touch runnum_list; rm runnum_list; touch runnum_list");
  for(Int_t i=0; i<rellum->GetEntries(); i++)
  {
    rellum->GetEntry(i);
    if(isConsistent)
    {
      gSystem->RedirectOutput("runnum_list");
      printf("%d\n",runnum);
      gSystem->RedirectOutput(0);
    };
  };
  gROOT->ProcessLine(".!cat runnum_list | uniq > runnum_list_tmp");
  gROOT->ProcessLine(".!mv runnum_list_tmp runnum_list");
  printf("runnum_list built\n");
};
