// computes average A_LL over runs for all bins
// -- draws kinematic dependence plot if eta_bins==1 AND ( pt_bins==1 OR en_bins==1 )
//
// CURRENTLY BROKEN, SINCE *_WIDTH VARIABLES NEED TO BE REDEFINED AFTER BINNING UPDATE
// -- see log 14.05.28 for details about the update

void DrawAverages(Bool_t eliminateGTone=0, const char * filename="spin.root")
{
  TFile * infile = new TFile(filename,"UPDATE");
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

  if(phi_bins!=1)
  {
    fprintf(stderr,"ERROR: phi_bins!=1.... I haven't implemented code for phi_bins>1 yet....\n");
    return;
  };

  Float_t A_LL,A_LL_err;
  Int_t eta_bin,pt_bin,en_bin;
  Int_t runnum;
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("A_LL",&A_LL);
  tr->SetBranchAddress("A_LL_err",&A_LL_err);
  tr->SetBranchAddress("eta_bin",&eta_bin);
  tr->SetBranchAddress("pt_bin",&pt_bin);
  tr->SetBranchAddress("en_bin",&en_bin);

  Float_t Ave_A_LL[eta_bins][pt_bins][en_bins];
  Float_t Ave_A_LL_err[eta_bins][pt_bins][en_bins];
  Int_t A_LL_cnt[eta_bins][pt_bins][en_bins];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        Ave_A_LL[g][p][e]=0;
        Ave_A_LL_err[g][p][e]=0;
        A_LL_cnt[g][p][e]=0;
      };
    };
  };
  printf("\n");
  if(eliminateGTone==0) printf("not eliminating cases with |A_LL|>1\n");
  for(Int_t q=0; q<tr->GetEntries(); q++)
  {
    tr->GetEntry(q);
    if(fabs(A_LL)<=1 || eliminateGTone==0)
    {
      Ave_A_LL[eta_bin][pt_bin][en_bin] += A_LL;
      Ave_A_LL_err[eta_bin][pt_bin][en_bin] += pow(A_LL_err,2);
      A_LL_cnt[eta_bin][pt_bin][en_bin] += 1;
    }
    else printf("g=%d p=%d e=%d r=%d omitted from average computation since |A_LL| > 1\n",g,p,e,runnum);
  };
  printf("\n");

  Int_t kin_dep_switch=0; // 0 = no plot;  1 = pt dep.  2 = en dep.
  if(eta_bins==1 && en_bins==1 && pt_bins>1) kin_dep_switch=1;
  else if(eta_bins==1 && pt_bins==1 && en_bins>1) kin_dep_switch=2;

  TGraphErrors * kin_dep = new TGraphErrors();
  kin_dep->SetName("kin_dep");
  if(kin_dep_switch==1) kin_dep->SetTitle("A_{LL} vs. p_{T}");
  else if(kin_dep_switch==2) kin_dep->SetTitle("A_{LL} vs. E_{#gamma#gamma}");
  kin_dep->SetMarkerStyle(kFullCircle);
  kin_dep->SetMarkerColor(kRed);
  kin_dep->SetMarkerSize(1.5);
  Int_t counter=0;

  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        if(A_LL_cnt[g][p][e]>0)
        {
          // average is sum divided by num
          Ave_A_LL[g][p][e] /= ((Float_t)A_LL_cnt[g][p][e]);
          // propagated error is sqrt of sum of squares of errors, divided by num
          Ave_A_LL_err[g][p][e] = 1/((Float_t)A_LL_cnt[g][p][e]) * sqrt(Ave_A_LL_err[g][p][e]);

          if(kin_dep_switch==1)
          {
            kin_dep->SetPoint(counter,(pt_div[p]+pt_div[p+1])/2.0,Ave_A_LL[g][p][e]);
            kin_dep->SetPointError(counter,pt_width/2.0,Ave_A_LL_err[g][p][e]);
            counter++;
          }
          else if(kin_dep_switch==2)
          {
            kin_dep->SetPoint(counter,(en_div[e]+en_div[e+1])/2.0,Ave_A_LL[g][p][e]);
            kin_dep->SetPointError(counter,en_width/2.0,Ave_A_LL_err[g][p][e]);
            counter++;
          };
        };
      };
    };
  };

  if(kin_dep_switch>0) 
  {
    kin_dep->Draw("ape");
    kin_dep->Write("kin_dep");
    printf("kin_dep added to %s\n",filename);
  }
  else printf("(pt_bins==1 || en_bins==1) = false  ==>  nothing drawn.\n");
};
