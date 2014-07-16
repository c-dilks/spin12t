// draws error(x) vs. x for various quantities
// -- rellum
// -- polarization
// -- A_LL
//
// CURRENTLY BROKEN, SINCE *_WIDTH VARIABLES NEED TO BE REDEFINED AFTER BINNING UPDATE
// -- see log 14.05.28 for details about the update

void DrawErrors(Int_t NBINS=50, const char * filename="spin.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("asy");

  // set binning
  // -- *_div = lower limit of each bin; last one is the upper limit
  // -- *_width = bin width
  Int_t phi_bins0, eta_bins0, pt_bins0, en_bins0;
  if(gSystem->Getenv("PHI")==NULL){fprintf(stderr,"ERROR: source env vars\n"); return;};
  sscanf(gSystem->Getenv("PHI"),"%d",&phi_bins0);
  sscanf(gSystem->Getenv("ETA"),"%d",&eta_bins0);
  sscanf(gSystem->Getenv("PT"),"%d",&pt_bins0);
  sscanf(gSystem->Getenv("EN"),"%d",&en_bins0);
  const Double_t pi=3.1415;
  const Double_t phi_bins=phi_bins0; 
  const Double_t eta_bins=eta_bins0;
    const Double_t eta_low=2.6; const Double_t eta_high=4.2;
  const Double_t pt_bins=pt_bins0; 
    const Double_t pt_low=0; const Double_t pt_high=10;
  const Double_t en_bins=en_bins0; 
    const Double_t en_low=0; const Double_t en_high=100;
  const Double_t phi_low = (-1*pi)-0.1;
  const Double_t phi_high = pi+0.1;
  Double_t phi_div[phi_bins+1];
  Double_t phi_width = (phi_high - phi_low)/phi_bins;
  for(Int_t i=0; i<phi_bins; i++) phi_div[i] = phi_low + i * phi_width;
  Double_t eta_div[eta_bins+1];
  Double_t eta_width = (eta_high - eta_low)/eta_bins;
  for(Int_t i=0; i<eta_bins; i++) eta_div[i] = eta_low + i * eta_width;
  Double_t pt_div[pt_bins+1];
  Double_t pt_width = (pt_high - pt_low)/pt_bins;
  for(Int_t i=0; i<pt_bins; i++) pt_div[i] = pt_low + i * pt_width;
  Double_t en_div[en_bins+1];
  Double_t en_width = (en_high - en_low)/en_bins;
  for(Int_t i=0; i<en_bins; i++) en_div[i] = en_low + i * en_width;


  eta_div[eta_bins] = eta_high;
  pt_div[pt_bins] = pt_high;
  en_div[en_bins] = en_high;

  Float_t R3_min, R3_max;
  Float_t R3_err_min, R3_err_max;
  Float_t PB_min, PB_max;
  Float_t PB_err_min, PB_err_max;
  Float_t PY_min, PY_max;
  Float_t PY_err_min, PY_err_max;
  Float_t A_LL_min, A_LL_max;
  Float_t A_LL_err_min, A_LL_err_max;

  R3_min = tr->GetMinimum("R3");
  R3_max = tr->GetMaximum("R3");
  R3_err_min = tr->GetMinimum("R3_err");
  R3_err_max = tr->GetMaximum("R3_err");
  PB_min = tr->GetMinimum("PB");
  PB_max = tr->GetMaximum("PB");
  PB_err_min = tr->GetMinimum("PB_err");
  PB_err_max = tr->GetMaximum("PB_err");
  PY_min = tr->GetMinimum("PY");
  PY_max = tr->GetMaximum("PY");
  PY_err_min = tr->GetMinimum("PY_err");
  PY_err_max = tr->GetMaximum("PY_err");
  A_LL_min = tr->GetMinimum("A_LL");
  A_LL_max = tr->GetMaximum("A_LL");
  A_LL_err_min = tr->GetMinimum("A_LL_err");
  A_LL_err_max = tr->GetMaximum("A_LL_err");

  TH2F * R3_dist = new TH2F("R3_dist","#sigma(R_{3}) vs. R_{3}",
    NBINS,R3_min,R3_max,NBINS,R3_err_min,R3_err_max);
  TH2F * PB_dist = new TH2F("PB_dist","#sigma(P_{B}) vs. P_{B}",
    NBINS,PB_min,PB_max,NBINS,PB_err_min,PB_err_max);
  TH2F * PY_dist = new TH2F("PY_dist","#sigma(P_{Y}) vs. P_{Y}",
    NBINS,PY_min,PY_max,NBINS,PY_err_min,PY_err_max);
  TH2F * A_LL_dist[eta_bins][pt_bins][en_bins];
  TH2F * A_LL_dist_all  = new TH2F("A_LL_dist_all","#sigma(A_{LL}) vs. A_{LL} :: all bins",
    NBINS,A_LL_min,A_LL_max,NBINS,A_LL_err_min,A_LL_err_max);

  char zero_bin[128];
  strcpy(zero_bin,"phi_bin==1 && eta_bin==0 && pt_bin==0 && pt_bin==0");

  tr->Project("R3_dist","R3_err:R3",zero_bin);
  tr->Project("PB_dist","PB_err:PB",zero_bin);
  tr->Project("PY_dist","PY_err:PY",zero_bin);
  tr->Project("A_LL_dist_all","A_LL_err:A_LL");

  TFile * outfile = new TFile("error.root","RECREATE");
  R3_dist->Write();
  PB_dist->Write();
  PY_dist->Write();
  A_LL_dist_all->Write();

  char A_LL_dist_n[eta_bins][pt_bins][en_bins][256];
  char A_LL_dist_t[eta_bins][pt_bins][en_bins][256];
  char A_LL_cut[eta_bins][pt_bins][en_bins][256];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        sprintf(A_LL_cut[g][p][e],"eta_bin==%d && pt_bin==%d && en_bin==%d",g,p,e);
        sprintf(A_LL_dist_n[g][p][e],"A_LL_dist_g%d_p%d_e%d",g,p,e);
        sprintf(A_LL_dist_t[g][p][e],
          "#sigma(A_{LL}) vs. A_{LL} :: #eta#in[%.2f,%.2f), p_{T}#in[%.2f,%.2f), E#in[%.2f,%.2f)",
          eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1]);
        A_LL_dist[g][p][e] = new TH2F(A_LL_dist_n[g][p][e],A_LL_dist_t[g][p][e],
          NBINS,A_LL_min,A_LL_max,NBINS,A_LL_err_min,A_LL_err_max);
        tr->Project(A_LL_dist_n[g][p][e],"A_LL_err:A_LL",A_LL_cut[g][p][e]);
        A_LL_dist[g][p][e]->Write();
      };
    };
  };
}
