void TestPol()
{
  TFile * infile = new TFile("../pol.root","READ");
  TTree * tr = (TTree*) infile->Get("pol");

  Int_t fill;
  Float_t b_pol,b_pol_e,y_pol,y_pol_e;
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("b_pol",&b_pol);
  tr->SetBranchAddress("b_pol_e",&b_pol_e);
  tr->SetBranchAddress("y_pol",&y_pol);
  tr->SetBranchAddress("y_pol_e",&y_pol_e);

  for(Int_t q=0; q<tr->GetEntries(); q++)
  {
    tr->GetEntry(q);
    gSystem->RedirectOutput("out_pol");
    printf("%d %f %f %f %f\n",fill,b_pol,b_pol_e,y_pol,y_pol_e);
    gSystem->RedirectOutput(0);
  };
};
