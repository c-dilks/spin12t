// computes asymmetry between fully summed (i.e. over all runs) phi distributions
// -- single helicity asymmetry for yellow beam A_L^Y
//    (see staszak thesis, eq. 7.13; polarization weighted as in A_LL)


void Asym_sy_sum(const char * filter_type="all",Int_t filter_low=0, Int_t filter_high=0)
{
  const Float_t pi=3.1415;
  TFile * infile = new TFile("phiset/all.root","READ");
  gSystem->Load("src/RunData.so");
  RunData * RD = new RunData();


  // get bins from environment
  Int_t phi_bins0, eta_bins0, pt_bins0, en_bins0;
  if(gSystem->Getenv("PHI_BINS")==NULL){fprintf(stderr,"ERROR: source env vars\n"); return;};
  sscanf(gSystem->Getenv("PHI_BINS"),"%d",&phi_bins0);
  sscanf(gSystem->Getenv("ETA_BINS"),"%d",&eta_bins0);
  sscanf(gSystem->Getenv("PT_BINS"),"%d",&pt_bins0);
  sscanf(gSystem->Getenv("EN_BINS"),"%d",&en_bins0);
  const Int_t phi_bins = phi_bins0;
  const Int_t eta_bins = eta_bins0;
  const Int_t pt_bins = pt_bins0;
  const Int_t en_bins = en_bins0;
  Float_t phi_div[phi_bins+1];
  Float_t eta_div[eta_bins+1];
  Float_t pt_div[pt_bins+1];
  Float_t en_div[en_bins+1];
  char phi_div_env[phi_bins+1][20];
  char eta_div_env[eta_bins+1][20];
  char pt_div_env[pt_bins+1][20];
  char en_div_env[en_bins+1][20];
  for(Int_t i=0; i<=phi_bins; i++)
  {
    sprintf(phi_div_env[i],"PHI_DIV_%d",i);
    sscanf(gSystem->Getenv(phi_div_env[i]),"%f",&(phi_div[i]));
    printf("%s %f\n",phi_div_env[i],phi_div[i]);
  };
  for(Int_t i=0; i<=eta_bins; i++)
  {
    sprintf(eta_div_env[i],"ETA_DIV_%d",i);
    sscanf(gSystem->Getenv(eta_div_env[i]),"%f",&(eta_div[i]));
    printf("%s %f\n",eta_div_env[i],eta_div[i]);
  };
  for(Int_t i=0; i<=pt_bins; i++)
  {
    sprintf(pt_div_env[i],"PT_DIV_%d",i);
    sscanf(gSystem->Getenv(pt_div_env[i]),"%f",&(pt_div[i]));
    printf("%s %f\n",pt_div_env[i],pt_div[i]);
  };
  for(Int_t i=0; i<=en_bins; i++)
  {
    sprintf(en_div_env[i],"EN_DIV_%d",i);
    sscanf(gSystem->Getenv(en_div_env[i]),"%f",&(en_div[i]));
    printf("%s %f\n",en_div_env[i],en_div[i]);
  };
  Float_t phi_low; sscanf(gSystem->Getenv("PHI_LOW"),"%f",&phi_low);
  Float_t phi_high; sscanf(gSystem->Getenv("PHI_HIGH"),"%f",&phi_high);
  Float_t eta_low; sscanf(gSystem->Getenv("ETA_LOW"),"%f",&eta_low);
  Float_t eta_high; sscanf(gSystem->Getenv("ETA_HIGH"),"%f",&eta_high);
  Float_t pt_low; sscanf(gSystem->Getenv("PT_LOW"),"%f",&pt_low);
  Float_t pt_high; sscanf(gSystem->Getenv("PT_HIGH"),"%f",&pt_high);
  Float_t en_low; sscanf(gSystem->Getenv("EN_LOW"),"%f",&en_low);
  Float_t en_high; sscanf(gSystem->Getenv("EN_HIGH"),"%f",&en_high);

  // read TObjArrays
  TObjArray * phi_dist_arr[4][eta_bins][pt_bins][en_bins];
  char phi_dist_arr_n[4][eta_bins][pt_bins][en_bins][64];
  Int_t NRUNS_tmp;
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(phi_dist_arr_n[s][g][p][e],"phi_dist_s%d_g%d_p%d_e%d",s,g,p,e);
          phi_dist_arr[s][g][p][e] = (TObjArray*) infile->Get(phi_dist_arr_n[s][g][p][e]);
          printf("phi_dist_arr[%d][%d][%d][%d] @ %p\n",s,g,p,e,(void*)phi_dist_arr[s][g][p][e]);
          if(s==0 && g==0 && p==0 && e==0)
          {
            NRUNS_tmp=phi_dist_arr[s][g][p][e]->GetEntries();
          }
          else
          {
            if(phi_dist_arr[s][g][p][e]->GetEntries() != NRUNS_tmp)
            {
              fprintf(stderr,"ERROR: TObjArrays have different sizes\n");
              return;
            };
          };
        };
      };
    };
  };
  const Int_t NRUNS = NRUNS_tmp;
  printf("NRUNS=%d\n",NRUNS);


  char var_str[16]; strcpy(var_str,"#phi");
  /*
  // get phi_low and phi_high (i.e. determine if distribution is
  // over phi or sin(phi) )
  //phi_dist_arr[0][0][0][0]->Print();
  Double_t low_edge = ((TH1D*)phi_dist_arr[0][0][0][0]->At(0))->GetBinLowEdge(1);
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
  //const Double_t phi_low = phi_low0;
  //const Double_t phi_high = phi_high0;
  const Double_t phi_low = -1*pi;
  const Double_t phi_high = pi;
  */

  


  // build summed phi distributions with the following weights / corrections
  // phi_dist_num:  spinbit 3 & 1, weight P;    spinbit 2 & 0, weight R*P
  // phi_dist_den:  spinbit 3 & 1, weight P^2;  spinbit 2 & 0, weight R*P^2
  TH1D * phi_dist_num[4][eta_bins][pt_bins][en_bins];
  TH1D * phi_dist_den[4][eta_bins][pt_bins][en_bins];
  char phi_dist_num_n[4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_den_n[4][eta_bins][pt_bins][en_bins][128];
  Int_t runnum;
  Float_t rellum,polar_b,polar_y,weight_num,weight_den;
  Float_t rellum_3;
  Int_t fill,pattern;
  Bool_t isConsistent;
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(phi_dist_num_n[s][g][p][e],"phi_num_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_den_n[s][g][p][e],"phi_den_s%d_g%d_p%d_e%d",s,g,p,e);
          phi_dist_num[s][g][p][e] = new TH1D(phi_dist_num_n[s][g][p][e],phi_dist_num_n[s][g][p][e],
            phi_bins,phi_low,phi_high);
          phi_dist_den[s][g][p][e] = new TH1D(phi_dist_den_n[s][g][p][e],phi_dist_den_n[s][g][p][e],
            phi_bins,phi_low,phi_high);

          for(Int_t r=0; r<NRUNS; r++)
          {
            sscanf(phi_dist_arr[0][g][p][e]->At(r)->GetName(),
              "phi_s%*d_g%*d_p%*d_e%*d_r%d",&runnum);

            rellum = RD->Rellum(runnum,1,"zdc"); // R1 for yellow single spin asym
            //rellum=1;
            //rellum=0.9986;
            polar_b = RD->BluePol(runnum);
            polar_y = RD->YellPol(runnum);
            fill = RD->GetFill(runnum);
            isConsistent = RD->RellumConsistent(runnum);
            pattern = RD->Pattern(runnum);

            weight_num = polar_y;
            weight_den = pow(polar_y, 2);
            if(s==2 || s==0) 
            {
              weight_num *= rellum;
              weight_den *= rellum;
            };

            rellum_3 = RD->Rellum(runnum,3,"zdc");

            if(isConsistent)
            {
              if( ( !strcmp(filter_type,"fill") && (fill>=filter_low && fill<=filter_high) ) ||
                  ( !strcmp(filter_type,"run") && (runnum>=filter_low && runnum<=filter_high) ) ||
                  ( !strcmp(filter_type,"runout") && !(runnum>=filter_low && runnum<=filter_high) ) ||
                  !strcmp(filter_type,"all"))
              {
                phi_dist_num[s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_num);
                phi_dist_den[s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_den);
              };
            };
          };
        };
      };
    };
  };
        

  // compute asymmetry; polarization and rellum have been corrected for above
  TH1D * dist_same_num[eta_bins][pt_bins][en_bins]; // S_num = P * (N++ + N-+)
  TH1D * dist_diff_num[eta_bins][pt_bins][en_bins]; // D_num = P * R * (N+- + N--)
  TH1D * dist_same_den[eta_bins][pt_bins][en_bins]; // S_den = P^2 * (N++ + N-+)
  TH1D * dist_diff_den[eta_bins][pt_bins][en_bins]; // D_den = P^2 * R * (N+- + N--)
  TH1D * numer[eta_bins][pt_bins][en_bins]; // S_num - D_num
  TH1D * denom[eta_bins][pt_bins][en_bins]; // S_num + D_num
  TH1D * asym[eta_bins][pt_bins][en_bins]; // numer / denom
  Int_t asym_pts[eta_bins][pt_bins][en_bins]; // number of points in asym
  char dist_same_num_n[eta_bins][pt_bins][en_bins][128];
  char dist_diff_num_n[eta_bins][pt_bins][en_bins][128];
  char dist_same_den_n[eta_bins][pt_bins][en_bins][128];
  char dist_diff_den_n[eta_bins][pt_bins][en_bins][128];
  char numer_n[eta_bins][pt_bins][en_bins][128];
  char denom_n[eta_bins][pt_bins][en_bins][128];
  char asym_n[eta_bins][pt_bins][en_bins][128];
  char asym_t[eta_bins][pt_bins][en_bins][256];
  Float_t p0,p0e,chi2,ndf;
  Float_t bc[4];
  Float_t bcent;
  Int_t runnum_0;
  Float_t A_LL[eta_bins][pt_bins][en_bins];
  Int_t A_LL_cnt[eta_bins][pt_bins][en_bins];
  TF1 * asym_fit[eta_bins][pt_bins][en_bins];
  Float_t asym_tmp;
  Float_t asym_max[eta_bins][pt_bins][en_bins];
  Float_t asym_min[eta_bins][pt_bins][en_bins];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        A_LL[g][p][e]=0.0;
        A_LL_cnt[g][p][e]=0;
        asym_max[g][p][e]=0;
        asym_min[g][p][e]=0;
      };
    };
  };
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        /*
        sscanf(phi_dist_arr[0][g][p][e]->At(r)->GetName(),
          "phi_s%*d_g%*d_p%*d_e%*d_r%d",&runnum);
        if(e==0 && g==0 && p==0) 
        {
          printf("%d %p %p r%d\n",r,(void*)phi_dist_arr[0][g][p][e],
              (void*)phi_dist_arr[0][g][p][e]->At(0),runnum);
          runnum_0 = runnum;
        }
        else
        {
          if(runnum != runnum_0)
          {
            fprintf(stderr,"Error: TObjArrays not synced\n");
            return;
          };
        };
        */

        sprintf(dist_same_num_n[g][p][e],"dist_same_num_g%d_p%d_e%d",g,p,e);
        sprintf(dist_diff_num_n[g][p][e],"dist_diff_num_g%d_p%d_e%d",g,p,e);
        sprintf(dist_same_den_n[g][p][e],"dist_same_den_g%d_p%d_e%d",g,p,e);
        sprintf(dist_diff_den_n[g][p][e],"dist_diff_den_g%d_p%d_e%d",g,p,e);
        sprintf(numer_n[g][p][e],"numer_g%d_p%d_e%d",g,p,e);
        sprintf(denom_n[g][p][e],"denom_g%d_p%d_e%d",g,p,e);
        sprintf(asym_n[g][p][e],"asym_g%d_p%d_e%d",g,p,e);
        sprintf(asym_t[g][p][e],
          "A_{L}^{Y} vs. %s :: #eta#in[%.2f,%.2f), p_{T}#in[%.2f,%.2f), E#in[%.2f,%.2f) (runsum)",
          var_str,eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1]);

        dist_same_num[g][p][e] = new TH1D(dist_same_num_n[g][p][e],
            dist_same_num_n[g][p][e],phi_bins,phi_low,phi_high);
        dist_diff_num[g][p][e] = new TH1D(dist_diff_num_n[g][p][e],
            dist_diff_num_n[g][p][e],phi_bins,phi_low,phi_high);
        dist_same_den[g][p][e] = new TH1D(dist_same_den_n[g][p][e],
            dist_same_den_n[g][p][e],phi_bins,phi_low,phi_high);
        dist_diff_den[g][p][e] = new TH1D(dist_diff_den_n[g][p][e],
            dist_diff_den_n[g][p][e],phi_bins,phi_low,phi_high);
        numer[g][p][e] = new TH1D(numer_n[g][p][e],
            numer_n[g][p][e],phi_bins,phi_low,phi_high);
        denom[g][p][e] = new TH1D(denom_n[g][p][e],
            denom_n[g][p][e],phi_bins,phi_low,phi_high);
        asym[g][p][e] = new TH1D(asym_n[g][p][e],
            asym_t[g][p][e],phi_bins,phi_low,phi_high);

        dist_same_num[g][p][e]->Add(phi_dist_num[3][g][p][e],phi_dist_num[1][g][p][e],1.0,1.0);
        dist_diff_num[g][p][e]->Add(phi_dist_num[2][g][p][e],phi_dist_num[0][g][p][e],1.0,1.0);
        dist_same_den[g][p][e]->Add(phi_dist_den[3][g][p][e],phi_dist_den[1][g][p][e],1.0,1.0);
        dist_diff_den[g][p][e]->Add(phi_dist_den[2][g][p][e],phi_dist_den[0][g][p][e],1.0,1.0);

        numer[g][p][e]->Add(dist_same_num[g][p][e],dist_diff_num[g][p][e],1.0,-1.0);
        denom[g][p][e]->Add(dist_same_den[g][p][e],dist_diff_den[g][p][e],1.0,1.0);

        asym[g][p][e]->Divide(numer[g][p][e],denom[g][p][e],1.0,1.0);


        printf(asym[g][p][e]->GetTitle());
        printf("\n");

        
        // n.b. for one phi bin, constant fit & error matches the bin & its error
        asym[g][p][e]->Fit("pol0","Q","",phi_low,phi_high);
        asym_fit[g][p][e] = asym[g][p][e]->GetFunction("pol0");
        if(asym_fit[g][p][e]!=NULL)
        {
          A_LL[g][p][e]+=asym_fit[g][p][e]->GetParameter(0);
          A_LL_cnt[g][p][e]++;
        };
      };
    };
  };

  // set plot ranges
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        //asym[g][p][e]->GetYaxis()->SetRangeUser(2*asym_min[g][p][e],2*asym_max[g][p][e]);
        asym[g][p][e]->GetXaxis()->SetRangeUser(phi_low,phi_high);
      };
    };
  };


  // average A_LL for each bin
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        A_LL[g][p][e] /= ((Float_t)A_LL_cnt[g][p][e]);
        printf("g%d p%d e%d <A_L>=%f\n",g,p,e,A_LL[g][p][e]);
      };
    };
  };


  // kinematic dependence plots

  TGraphErrors * en_dep[pt_bins]; // en dependent plots, one for each pt bin
  TGraphErrors * pt_dep[en_bins]; // pt dependent plots, one for each en bin
  char en_dep_t[pt_bins][128];
  char pt_dep_t[en_bins][128];
  Int_t en_dep_cnt[pt_bins]; // en dependent plots point counter
  Int_t pt_dep_cnt[en_bins]; // pt dependent plots point counter
  for(Int_t p=0; p<pt_bins; p++) en_dep_cnt[p]=0;
  for(Int_t e=0; e<en_bins; e++) pt_dep_cnt[e]=0;

  Double_t aLL_en[pt_bins][en_bins];     // arrays for en dependent plots, one 
  Double_t err_en[pt_bins][en_bins];     // for each pt bin
  Double_t cent_en[pt_bins][en_bins];
  Double_t width_en[pt_bins][en_bins];

  Double_t aLL_pt[en_bins][pt_bins];     // arrays for pt dependent plots, one
  Double_t err_pt[en_bins][pt_bins];     // for each en bin
  Double_t cent_pt[en_bins][pt_bins];
  Double_t width_pt[en_bins][pt_bins];

  if(eta_bins==1)
  {
    // en dependent points for each pt bin
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        if(asym[0][p][e]->GetFunction("pol0"))
        {
          aLL_en[p][en_dep_cnt[p]] = asym[0][p][e]->GetFunction("pol0")->GetParameter(0);
          err_en[p][en_dep_cnt[p]] = asym[0][p][e]->GetFunction("pol0")->GetParError(0);
          cent_en[p][en_dep_cnt[p]] = en_div[e] + ((en_div[e+1]-en_div[e])/2.0);
          width_en[p][en_dep_cnt[p]] = (en_div[e+1]-en_div[e])/2.0;
          en_dep_cnt[p]++;
        };
      };
      en_dep[p] = new TGraphErrors(en_dep_cnt[p],cent_en[p],aLL_en[p],width_en[p],err_en[p]);
      sprintf(en_dep_t[p],"A_{L}^{Y} vs. E_{#gamma#gamma} for p_{T}#in[%.2f,%.2f)",pt_div[p],pt_div[p+1]);
      en_dep[p]->SetTitle(en_dep_t[p]);
      en_dep[p]->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
    };

    // pt dependent points for each en bin
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        if(asym[0][p][e]->GetFunction("pol0"))
        {
          aLL_pt[e][pt_dep_cnt[e]] = asym[0][p][e]->GetFunction("pol0")->GetParameter(0);
          err_pt[e][pt_dep_cnt[e]] = asym[0][p][e]->GetFunction("pol0")->GetParError(0);
          cent_pt[e][pt_dep_cnt[e]] = pt_div[p] + ((pt_div[p+1]-pt_div[p])/2.0);
          width_pt[e][pt_dep_cnt[e]] = (pt_div[p+1]-pt_div[p])/2.0;
          pt_dep_cnt[e]++;
        };
      };
      pt_dep[e] = new TGraphErrors(pt_dep_cnt[e],cent_pt[e],aLL_pt[e],width_pt[e],err_pt[e]);
      sprintf(pt_dep_t[e],"A_{L}^{Y} vs. p_{T} for E_{#gamma#gamma}#in[%.2f,%.2f)",en_div[e],en_div[e+1]);
      pt_dep[e]->SetTitle(pt_dep_t[e]);
      pt_dep[e]->GetXaxis()->SetTitle("p_{T} (GeV)");
    };
  };

  for(Int_t p=0; p<pt_bins; p++)
  {
    en_dep[p]->GetYaxis()->SetTitle("A_{L}^{Y}");
    en_dep[p]->SetMarkerStyle(kFullCircle);
    en_dep[p]->SetMarkerColor(kRed);
    en_dep[p]->GetYaxis()->SetTitleOffset(1.5);
  };
  for(Int_t e=0; e<en_bins; e++)
  {
    pt_dep[e]->GetYaxis()->SetTitle("A_{L}^{Y}");
    pt_dep[e]->SetMarkerStyle(kFullCircle);
    pt_dep[e]->SetMarkerColor(kRed);
    pt_dep[e]->GetYaxis()->SetTitleOffset(1.5);
  };


  // write phi dists
  printf("writing spin_sy_sum.root...\n");
  TFile * outfile = new TFile("spin_sy_sum.root","RECREATE");
  char en_dep_n[pt_bins][32];
  char pt_dep_n[en_bins][32];
  for(Int_t p=0; p<pt_bins; p++)
  {
    sprintf(en_dep_n[p],"en_dep_p%d",p);
    en_dep[p]->Write(en_dep_n[p]);
  };
  for(Int_t e=0; e<en_bins; e++)
  {
    sprintf(pt_dep_n[e],"pt_dep_e%d",e);
    pt_dep[e]->Write(pt_dep_n[e]);
  };
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_num[s][g][p][e]->Write(phi_dist_num_n[s][g][p][e]);
          phi_dist_den[s][g][p][e]->Write(phi_dist_den_n[s][g][p][e]);
        };
      };
    };
  };
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        asym[g][p][e]->Write(asym_n[g][p][e]);
      };
    };
  };
  printf("written\n");
};
