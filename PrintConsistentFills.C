// prints list of fills which have at least one run which has
// consistent rellum measurements

void PrintConsistentFills()
{
  TFile * infile = new TFile("rtree.root","READ");
  TTree * rellum = (TTree*) infile->Get("rellum");

  Int_t fill;
  Bool_t isConsistent;

  rellum->SetBranchAddress("fill",&fill);
  rellum->SetBranchAddress("isConsistent",&isConsistent);

  gROOT->ProcessLine(".!touch fill_list; rm fill_list; touch fill_list");
  for(Int_t i=0; i<rellum->GetEntries(); i++)
  {
    rellum->GetEntry(i);
    if(isConsistent)
    {
      gSystem->RedirectOutput("fill_list");
      printf("%d\n",fill);
      gSystem->RedirectOutput(0);
    };
  };
  gROOT->ProcessLine(".!cat fill_list | uniq > fill_list_tmp");
  gROOT->ProcessLine(".!mv fill_list_tmp fill_list");
  printf("fill_list built\n");
};
