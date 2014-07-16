// DEPRECATED
// computes (currently only double) spin asymmetry
// - for phi dists or sin(phi) dists (reads file in subdirectory phiset/)

void Spin(const char * filename="all.root", Int_t phi_bins0=16,
          Int_t eta_bins0=1, Int_t pt_bins0=3, Int_t en_bins0=3)
{
  // set output file name
  char setname[32];
  if(!strcmp(filename,"all.root"))
    sprintf(setname,"_all.root");
  else
    sscanf(filename,"phi%s",setname);
  char outname[128];
  sprintf(outname,"spinset/spin%s",setname);
  // open input file
  sprintf(filename,"phiset/%s",filename);
  TFile * infile = new TFile(filename,"READ");

  // set ranges
  const Double_t pi=3.1415;
  const Double_t phi_bins=phi_bins0; 
  const Double_t eta_bins=eta_bins0;
    const Double_t eta_low=2.6; const Double_t eta_high=4.2;
  const Double_t pt_bins=pt_bins0; 
    const Double_t pt_low=0; const Double_t pt_high=10;
  const Double_t en_bins=en_bins0; 
    const Double_t en_low=0; const Double_t en_high=100;

  // set binning
  // -- *_div = lower limit of each bin; last one is the upper limit
  // -- *_width = bin width
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

  // read phi distributions, called "dist"
  // and polarized weighted phi distributions, called dist_P
  TH1D * dist[4][eta_bins][pt_bins][en_bins];
  char dist_name[4][eta_bins][pt_bins][en_bins][32];
  TH1D * dist_P[4][eta_bins][pt_bins][en_bins];
  char dist_P_name[4][eta_bins][pt_bins][en_bins][32];
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(dist_name[s][g][p][e],"phi_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(dist_P_name[s][g][p][e],"phi_P_s%d_g%d_p%d_e%d",s,g,p,e);
          dist[s][g][p][e] = (TH1D*) infile->Get(dist_name[s][g][p][e]);
          dist_P[s][g][p][e] = (TH1D*) infile->Get(dist_P_name[s][g][p][e]);
        };
      };
    };
  };

  // get phi_low and phi_high (i.e. determine if distribution is
  // over phi or sin(phi) )
  Double_t low_edge = dist[0][0][0][0]->GetBinLowEdge(1);
  Double_t phi_low0,phi_high0;
  char var_str[16];
  if(fabs(low_edge+1)<0.2)
  {
    phi_low0=-1;
    phi_high0=1;
    strcpy(var_str,"sin(#phi)");
  }
  else if(fabs(low_edge+pi)<0.2)
  {
    phi_low0=-1*pi;
    phi_high0=pi;
    strcpy(var_str,"#phi");
  }
  else
  {
    fprintf(stderr,"ERROR: unfamiliar binning of distribution (binning recently changed?)\n");
    return;
  };
  const Double_t phi_low = phi_low0;
  const Double_t phi_high = phi_high0;


  // compute asymmetry
  TH1D * dist_same[eta_bins][pt_bins][en_bins]; // S = N++ + N--
  TH1D * dist_diff[eta_bins][pt_bins][en_bins]; // D = R*(N+- + N-+)
  TH1D * dist_same_P[eta_bins][pt_bins][en_bins]; // S_P = 1/P * (N++ + N--)
  TH1D * dist_diff_P[eta_bins][pt_bins][en_bins]; // D_P = 1/P * ( R*(N+- + N-+) )
  TH1D * numer[eta_bins][pt_bins][en_bins]; // S_P - D_P
  TH1D * denom[eta_bins][pt_bins][en_bins]; // S + D
  TH1D * asym[eta_bins][pt_bins][en_bins]; // numer / denom
  char dist_same_n[eta_bins][pt_bins][en_bins][64];
  char dist_diff_n[eta_bins][pt_bins][en_bins][64];
  char dist_same_P_n[eta_bins][pt_bins][en_bins][64];
  char dist_diff_P_n[eta_bins][pt_bins][en_bins][64];
  char asym_n[eta_bins][pt_bins][en_bins][64];
  char numer_n[eta_bins][pt_bins][en_bins][64];
  char denom_n[eta_bins][pt_bins][en_bins][64];
  char asym_t[eta_bins][pt_bins][en_bins][256];
  TF1 * cos_fit[eta_bins][pt_bins][en_bins];
  char cos_fit_n[eta_bins][pt_bins][en_bins][32];
  Float_t p0,p0e,chi2,ndf;
  gROOT->ProcessLine(".! touch fit_data.txt; rm fit_data.txt; touch fit_data.txt");
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        sprintf(dist_same_n[g][p][e],"dist_same_g%d_p%d_e%d",g,p,e);
        sprintf(dist_diff_n[g][p][e],"dist_diff_g%d_p%d_e%d",g,p,e);
        sprintf(dist_same_P_n[g][p][e],"dist_same_P_g%d_p%d_e%d",g,p,e);
        sprintf(dist_diff_P_n[g][p][e],"dist_diff_P_g%d_p%d_e%d",g,p,e);
        sprintf(numer_n[g][p][e],"numer_g%d_p%d_e%d",g,p,e);
        sprintf(denom_n[g][p][e],"denom_g%d_p%d_e%d",g,p,e);
        sprintf(asym_n[g][p][e],"asym_g%d_p%d_e%d",g,p,e);
        sprintf(asym_t[g][p][e],
          "A_{LL} vs. %s :: #eta#in[%.2f,%.2f), p_{T}#in[%.2f,%.2f), E#in[%.2f,%.2f)",
          var_str,eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1]);
        dist_same[g][p][e] = new TH1D(dist_same_n[g][p][e],
            dist_same_n[g][p][e],phi_bins,phi_low,phi_high);
        dist_diff[g][p][e] = new TH1D(dist_diff_n[g][p][e],
            dist_diff_n[g][p][e],phi_bins,phi_low,phi_high);
        dist_same_P[g][p][e] = new TH1D(dist_same_P_n[g][p][e],
            dist_same_P_n[g][p][e],phi_bins,phi_low,phi_high);
        dist_diff_P[g][p][e] = new TH1D(dist_diff_P_n[g][p][e],
            dist_diff_P_n[g][p][e],phi_bins,phi_low,phi_high);
        numer[g][p][e] = new TH1D(numer_n[g][p][e],
            numer_n[g][p][e],phi_bins,phi_low,phi_high);
        denom[g][p][e] = new TH1D(denom_n[g][p][e],
            denom_n[g][p][e],phi_bins,phi_low,phi_high);
        asym[g][p][e] = new TH1D(asym_n[g][p][e],
            asym_t[g][p][e],phi_bins,phi_low,phi_high);

        dist_same[g][p][e]->Add(dist[0][g][p][e],dist[3][g][p][e],1.0,1.0);
        dist_diff[g][p][e]->Add(dist[1][g][p][e],dist[2][g][p][e],1.0,1.0);

        dist_same_P[g][p][e]->Add(dist_P[0][g][p][e],dist_P[3][g][p][e],1.0,1.0);
        dist_diff_P[g][p][e]->Add(dist_P[1][g][p][e],dist_P[2][g][p][e],1.0,1.0);
        
        numer[g][p][e]->Add(dist_same_P[g][p][e],dist_diff_P[g][p][e],1.0,-1.0);
        denom[g][p][e]->Add(dist_same[g][p][e],dist_diff[g][p][e],1.0,1.0);

        asym[g][p][e]->Divide(numer[g][p][e],denom[g][p][e],1.0,1.0);

        printf("\n");
        printf(asym[g][p][e]->GetTitle());

        /*
        sprintf(cos_fit_n[g][p][e],"cos_fit_g%d_p%d_e%d",g,p,e);
        cos_fit[g][p][e] = new TF1(cos_fit_n[g][p][e],"[0]*cos([1]*x)+[2]",-1*pi,pi/3);
        cos_fit[g][p][e]->SetParLimits(0,0.00001,100); // amplitude
        cos_fit[g][p][e]->SetParLimits(1,0.00001,2); // frequency
        cos_fit[g][p][e]->SetParLimits(2,-1,1); // offset
        asym[g][p][e]->Fit(cos_fit_n[g][p][e],"","",-1*pi,pi/3);
        */
        asym[g][p][e]->Fit("pol0","","",phi_low,phi_high);

        // output fit parameters (for binning in Pt only)
        p0 = asym[g][p][e]->GetFunction("pol0")->GetParameter(0);
        p0e = asym[g][p][e]->GetFunction("pol0")->GetParError(0);
        chi2 = asym[g][p][e]->GetFunction("pol0")->GetChisquare();
        ndf = asym[g][p][e]->GetFunction("pol0")->GetNDF();
        gSystem->RedirectOutput("fit_data.txt","a");
        printf("%f %f %f %f %f %f\n",pt_div[p],pt_div[p+1],p0,p0e,chi2,ndf);
        gSystem->RedirectOutput(0);
      };
    };
  };


  // plot pt or en dependence
  Double_t aLL_pt[pt_bins];
  Double_t aLL_en[en_bins];
  Double_t err_pt[pt_bins];
  Double_t err_en[en_bins];
  Double_t pt_cent[pt_bins];
  Double_t en_cent[en_bins];
  Double_t pt_binerr[pt_bins];
  Double_t en_binerr[en_bins];
  TGraphErrors * kin_dep;
  if(eta_bins==1)
  {
    if(en_bins==1 && pt_bins>1)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        if(asym[0][p][0]->GetFunction("pol0"))
        {
          aLL_pt[p] = asym[0][p][0]->GetFunction("pol0")->GetParameter(0);
          err_pt[p] = asym[0][p][0]->GetFunction("pol0")->GetParError(0);
          pt_cent[p] = pt_div[p] + (pt_width/2.0);
          pt_binerr[p] = pt_width/2.0;
        };
      };
      kin_dep = new TGraphErrors(pt_bins,pt_cent,aLL_pt,pt_binerr,err_pt);
      kin_dep->SetTitle("A_{LL} vs. P_{T}");
      kin_dep->GetXaxis()->SetTitle("P_{T} (GeV)");
    }
    else if(pt_bins==1 && en_bins>1)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        if(asym[0][0][e]->GetFunction("pol0"))
        {
          aLL_en[e] = asym[0][0][e]->GetFunction("pol0")->GetParameter(0);
          err_en[e] = asym[0][0][e]->GetFunction("pol0")->GetParError(0);
          en_cent[e] = en_div[e] + (en_width/2.0);
          en_binerr[e] = en_width/2.0;
        };
      };
      kin_dep = new TGraphErrors(en_bins,en_cent,aLL_en,en_binerr,err_en);
      kin_dep->SetTitle("A_{LL} vs. E_{#gamma#gamma}");
      kin_dep->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
    };
  };
  kin_dep->GetYaxis()->SetTitle("A_{LL}");
  kin_dep->SetMarkerStyle(kFullCircle);
  kin_dep->SetMarkerColor(kRed);
  kin_dep->GetYaxis()->SetTitleOffset(1.5);


  // write output
  TFile * outfile = new TFile(outname,"RECREATE");
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        dist_same[g][p][e]->Write();
        dist_diff[g][p][e]->Write();
        asym[g][p][e]->Write();
      };
    };
  };
  if(en_bins==1 || pt_bins==1) kin_dep->Write("kin_dep");
  printf("\n%s written\n",outname);
};
