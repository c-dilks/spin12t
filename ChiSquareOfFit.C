// analyses pT 1-5 GeV region, fits to constant line at A_LL=0 and 
// prints the chi2 & ndf to an output file
//
// run using asym_call_{fill,run} !!!

void ChiSquareOfFit(Int_t num=0,
                    const char * filename="spin_sy_sum.root",
                    const char * graph_name="pt_dep_e0")
{
  TFile * infile = new TFile(filename,"READ");
  TGraphErrors * gr = (TGraphErrors*)infile->Get(graph_name);

  Float_t pt_low=1;
  Float_t pt_high=5;

  TF1 * ff = new TF1("ff","pol0",pt_low,pt_high);
  ff->FixParameter(0,0);
  gr->Fit(ff,"","",pt_low,pt_high);

  Float_t chisquare = ff->GetChisquare();
  Float_t ndf = ff->GetNDF();

  gSystem->RedirectOutput("chisq","a");
  printf("%d %f %f\n",num,chisquare,ndf);
  gSystem->RedirectOutput(0);
};
