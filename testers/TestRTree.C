void TestRTree()
{
  TFile * infile = new TFile("../rtree.root","READ");
  TTree * tr = (TTree*) infile->Get("rellum");

  Int_t i,runnum,fill;
  Float_t R3,R3err;
  Bool_t isConsistent;
  Int_t pattern;
  tr->SetBranchAddress("i",&i);
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("R3_vpd_mean",&R3);
  tr->SetBranchAddress("R3_vpd_mean_err",&R3err);
  tr->SetBranchAddress("isConsistent",&isConsistent);
  tr->SetBranchAddress("pattern",&pattern);

  for(Int_t q=0; q<tr->GetEntries(); q++)
  {
    tr->GetEntry(q);
    gSystem->RedirectOutput("out_rtree");
    printf("%d %d %d %f %f %d %d\n",q,runnum,fill,R3,R3err,(Int_t)isConsistent,pattern);
    gSystem->RedirectOutput(0);
  };
};
