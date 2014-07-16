// computes run-by run asymmetries

void Asym_run_by_run()
{
  // open phi dists file and construct RunData12 RD
  TFile * infile = new TFile("phiset/all.root","READ");
  gSystem->Load("src/RunData12.so");
  RunData12 * RD = new RunData12();


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
  TObjArray * phi_dist_arr_dP[4][eta_bins][pt_bins][en_bins];
  char phi_dist_arr_n[4][eta_bins][pt_bins][en_bins][64];
  char phi_dist_arr_dP_n[4][eta_bins][pt_bins][en_bins][64];
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
          sprintf(phi_dist_arr_dP_n[s][g][p][e],"phi_dist_dP_s%d_g%d_p%d_e%d",s,g,p,e);
          phi_dist_arr[s][g][p][e] = (TObjArray*) infile->Get(phi_dist_arr_n[s][g][p][e]);
          phi_dist_arr_dP[s][g][p][e] = (TObjArray*) infile->Get(phi_dist_arr_dP_n[s][g][p][e]);
          printf("phi_dist_arr[%d][%d][%d][%d] @ %p\n",s,g,p,e,(void*)phi_dist_arr[s][g][p][e]);
          if(s==0 && g==0 && p==0 && e==0)
          {
            NRUNS_tmp=phi_dist_arr[s][g][p][e]->GetEntries();
          }
          else
          {
            if(phi_dist_arr[s][g][p][e]->GetEntries() != NRUNS_tmp ||
               phi_dist_arr_dP[s][g][p][e]->GetEntries() != NRUNS_tmp)
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
  */
  char var_str[16]; // sin(phi) vs. phi is deprecated!
  strcpy(var_str,"#phi");


  // compute asymmetry
  Int_t runnum;
  Int_t runnum_0;
  Int_t index=0;
  Double_t DD,MM,SS; // D, M, and S definitions from uncertainty.pdf
  Float_t R3,PB,PY; // rellum & polarizations
  Float_t R3_e,PB_e,PY_e; // rellum & polarization uncertainties
  Double_t T1,T2,T3,T4; // "terms" (see below)
  Float_t asym_calc; // calculated asymmetry
  Float_t asym_unc; // sigma_{A_LL} = sqrt( (T1 + T2 + T3) / T4 )
  //Float_t phi_p,phi_w; // _p = point, _w = width (half)
  Float_t eta_p,eta_w;
  Float_t pt_p,pt_w;
  Float_t en_p,en_w;
  Int_t phi_b,eta_b,pt_b,en_b; // bin no.'s 
  
  //phi_w = phi_width / 2.0;
  //eta_w = eta_width / 2.0;
  //pt_w = pt_width / 2.0;
  //en_w = en_width / 2.0;

  TTree * asy = new TTree();
  asy->Branch("i",&index,"i/I"); // run index
  asy->Branch("runnum",&runnum,"runnum/I");
  asy->Branch("R3",&R3,"R3/F"); // rellum 
  asy->Branch("R3_err",&R3_e,"R3_err/F"); // rellum error
  asy->Branch("PB",&PB,"PB/F"); // blue polarization
  asy->Branch("PY",&PY,"PY/F"); // yellow polarization
  asy->Branch("PB_err",&PB_e,"PB_err/F"); // blue polarization error
  asy->Branch("PY_err",&PY_e,"PY_err/F"); // yellow polarization error
  asy->Branch("phi",&phi_p,"phi/F"); // phi point
  //asy->Branch("phi_width",&phi_w,"phi_width/F"); // phi width
  asy->Branch("phi_bin",&phi_b,"phi_bin/I"); // phi bin no.
  asy->Branch("eta",&eta_p,"eta/F"); // eta point
  //asy->Branch("eta_width",&eta_w,"eta_width/F"); // eta width
  asy->Branch("eta_bin",&eta_b,"eta_bin/I"); // eta bin no.
  asy->Branch("pt",&pt_p,"pt/F"); // pt point
  //asy->Branch("pt_width",&pt_w,"pt_width/F"); // pt width
  asy->Branch("pt_bin",&pt_b,"pt_bin/I"); // pt bin no.
  asy->Branch("en",&en_p,"en/F"); // en point
  //asy->Branch("en_width",&en_w,"en_width/F"); // en width
  asy->Branch("en_bin",&en_b,"en_bin/I"); // en bin no.
  asy->Branch("A_LL",&asym_calc,"A_LL/F"); 
  asy->Branch("A_LL_err",&asym_unc,"A_LL_err/F");

  TH1D * dist_same[eta_bins][pt_bins][en_bins][NRUNS]; // S = N++ + N--
  TH1D * dist_diff[eta_bins][pt_bins][en_bins][NRUNS]; // D = R*(N+- + N-+)
  TH1D * dist_same_dP[eta_bins][pt_bins][en_bins][NRUNS]; // S_dP = 1/P * (N++ + N--)
  TH1D * dist_diff_dP[eta_bins][pt_bins][en_bins][NRUNS]; // D_dP = 1/P * ( R*(N+- + N-+) )
  TH1D * numer[eta_bins][pt_bins][en_bins][NRUNS]; // S_dP - D_dP
  TH1D * denom[eta_bins][pt_bins][en_bins][NRUNS]; // S + D
  TGraphErrors * asym[eta_bins][pt_bins][en_bins][NRUNS]; // numer / denom
  Int_t asym_pts[eta_bins][pt_bins][en_bins][NRUNS]; // number of points in asym
  char dist_same_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char dist_diff_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char dist_same_dP_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char dist_diff_dP_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char numer_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char denom_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char asym_n[eta_bins][pt_bins][en_bins][NRUNS][128];
  char asym_t[eta_bins][pt_bins][en_bins][NRUNS][256];
  Float_t p0,p0e,chi2,ndf;
  Float_t bc[4];
  Float_t bc_dP[4];
  Float_t A_LL[eta_bins][pt_bins][en_bins];
  Int_t A_LL_cnt[eta_bins][pt_bins][en_bins];
  TF1 * asym_fit[eta_bins][pt_bins][en_bins];
  Float_t asym_max[eta_bins][pt_bins][en_bins][NRUNS];
  Float_t asym_min[eta_bins][pt_bins][en_bins][NRUNS];
  Int_t case_counter=0;
  Int_t case_fail=0;
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        A_LL[g][p][e]=0.0;
        A_LL_cnt[g][p][e]=0;
        for(Int_t r=0; r<NRUNS; r++)
        {
          asym_max[g][p][e][r]=0;
          asym_min[g][p][e][r]=0;
        };
      };
    };
  };
  for(Int_t r=0; r<NRUNS; r++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          // get bin numbers and positions
          eta_b = g;
          pt_b = p;
          en_b = e;
          eta_p = (eta_div[g+1] + eta_div[g]) / 2.0;
          pt_p = (pt_div[p+1] + pt_div[p]) / 2.0;
          en_p = (en_div[e+1] + en_div[e]) / 2.0;

          // get run number
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
          index=r+1;

          // define distributions
          sprintf(dist_same_n[g][p][e][r],"dist_same_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(dist_diff_n[g][p][e][r],"dist_diff_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(dist_same_dP_n[g][p][e][r],"dist_same_dP_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(dist_diff_dP_n[g][p][e][r],"dist_diff_dP_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(numer_n[g][p][e][r],"numer_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(denom_n[g][p][e][r],"denom_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(asym_n[g][p][e][r],"asym_g%d_p%d_e%d_r%d",g,p,e,runnum);
          sprintf(asym_t[g][p][e][r],
            "A_{LL} vs. %s :: #eta#in[%.2f,%.2f), p_{T}#in[%.2f,%.2f), E#in[%.2f,%.2f), run=%d",
            var_str,eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum);
          dist_same[g][p][e][r] = new TH1D(dist_same_n[g][p][e][r],
              dist_same_n[g][p][e][r],phi_bins,phi_low,phi_high);
          dist_diff[g][p][e][r] = new TH1D(dist_diff_n[g][p][e][r],
              dist_diff_n[g][p][e][r],phi_bins,phi_low,phi_high);
          dist_same_dP[g][p][e][r] = new TH1D(dist_same_dP_n[g][p][e][r],
              dist_same_dP_n[g][p][e][r],phi_bins,phi_low,phi_high);
          dist_diff_dP[g][p][e][r] = new TH1D(dist_diff_dP_n[g][p][e][r],
              dist_diff_dP_n[g][p][e][r],phi_bins,phi_low,phi_high);
          numer[g][p][e][r] = new TH1D(numer_n[g][p][e][r],
              numer_n[g][p][e][r],phi_bins,phi_low,phi_high);
          denom[g][p][e][r] = new TH1D(denom_n[g][p][e][r],
              denom_n[g][p][e][r],phi_bins,phi_low,phi_high);

          asym[g][p][e][r] = new TGraphErrors();
          asym[g][p][e][r]->SetName(asym_n[g][p][e][r]);
          asym[g][p][e][r]->SetTitle(asym_t[g][p][e][r]);
          asym[g][p][e][r]->SetMarkerStyle(kFullCircle);
          asym[g][p][e][r]->SetMarkerColor(kRed);
          asym[g][p][e][r]->SetMarkerSize(1.5);
          asym_pts[g][p][e][r] = 0;

          // loop through phi
          for(Int_t b=1; b<=phi_bins; b++)
          {
            phi_b = b;
            phi_p = ((TH1D*)(phi_dist_arr[0][g][p][e]->At(r)))->GetBinCenter(b);

            for(Int_t s=0; s<4; s++)
            {
              bc[s] = ((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)))->GetBinContent(b);
              bc_dP[s] = ((TH1D*)(phi_dist_arr_dP[s][g][p][e]->At(r)))->GetBinContent(b);
            };

            // see 14.03.14 log entry for details about choosing this filter
            //if(bc[0]*bc[1]*bc[2]*bc[3]>0)
            //if(!(bc[0]==0 && bc[1]==0 && bc[2]==0 && bc[3]==0))
            if( (bc[0]!=0 || bc[3]!=0) && (bc[1]!=0 || bc[2]!=0) )
            {
              case_counter++;

              // get rellum and polarization
              R3 = RD->Rellum(runnum,3,"zdc");
              R3_e = RD->RellumErr(runnum,3,"zdc");
              PB = RD->BluePol(runnum);
              PY = RD->YellPol(runnum);
              PB_e = RD->BluePolErr(runnum);
              PY_e = RD->YellPolErr(runnum);

              // FORMULA -- A_LL
              asym_calc = 
               ( (bc_dP[0]+bc_dP[3]) - R3*(bc_dP[1] + bc_dP[2]) ) / ( bc[0] + bc[3] + R3*(bc[1] + bc[2]) );
              asym_max[g][p][e][r] = (asym_calc > asym_max[g][p][e][r]) ? asym_calc:asym_max[g][p][e][r];
              asym_min[g][p][e][r] = (asym_calc < asym_min[g][p][e][r]) ? asym_calc:asym_min[g][p][e][r];

              // FORMULA - compute A_LL uncertainty
              DD = pow(bc[0] + bc[3], 2) - pow(R3,2) * pow(bc[1] + bc[2], 2);
              MM = (bc[1] + bc[2]) * (bc[0] + bc[3]);
              SS = bc[0] + bc[1] + bc[2] + bc[3];
              T1 = pow(PY_e,2) * pow(PB,2) * pow(DD,2);
              T2 = pow(PB_e,2) * pow(PY,2) * pow(DD,2);
              T3 = 4 * pow(PB,2) * pow(PY,2) * MM * ( pow(R3_e,2) * MM + pow(R3,2) * SS );
              T4 = pow(PB,4) * pow(PY,4) * pow( ((bc[0] + bc[3]) + R3*(bc[1] + bc[2])), 4);
              asym_unc = sqrt((T1+T2+T3)/T4);

              // set A_LL point with uncertainty
              asym[g][p][e][r]->SetPoint(asym_pts[g][p][e][r],phi_p,asym_calc);
              asym[g][p][e][r]->SetPointError(asym_pts[g][p][e][r],phi_p,asym_unc);
              asym_pts[g][p][e][r]++;

              asy->Fill();
            }
            else case_fail++;
          };

          // asym binned in phi distribution fits
          asym[g][p][e][r]->Fit("pol0","Q","",phi_low,phi_high);
          asym_fit[g][p][e] = asym[g][p][e][r]->GetFunction("pol0");
          if(asym_fit[g][p][e]!=NULL)
          {
            A_LL[g][p][e]+=asym_fit[g][p][e]->GetParameter(0);
            A_LL_cnt[g][p][e]++;
          };
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
        for(Int_t r=0; r<NRUNS; r++)
        {
          asym[g][p][e][r]->GetYaxis()->SetLimits(2*asym_min[g][p][e][r],2*asym_max[g][p][e][r]);
          asym[g][p][e][r]->GetXaxis()->SetLimits(phi_low,phi_high);
        };
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
        printf("g%d p%d e%d <A_LL>=%f\n",g,p,e,A_LL[g][p][e]);
      };
    };
  };



  // write arrays
  printf("writing spin.root...\n");
  TFile * outfile = new TFile("spin.root","RECREATE");
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_arr[s][g][p][e]->Write(phi_dist_arr_n[s][g][p][e],
              TObject::kSingleKey);
          phi_dist_arr_dP[s][g][p][e]->Write(phi_dist_arr_dP_n[s][g][p][e],
              TObject::kSingleKey);
        };
      };
    };
  };
  TObjArray * asym_arr[eta_bins][pt_bins][en_bins];
  char asym_arr_n[eta_bins][pt_bins][en_bins][64];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        asym_arr[g][p][e] = new TObjArray();
        for(Int_t r=0; r<NRUNS; r++) 
        {
          asym_arr[g][p][e]->AddLast(asym[g][p][e][r]);
          //printf("%p  %d\n",(void*)asym[g][p][e][r],asym_arr[g][p][e]->GetEntries());
        };
        sprintf(asym_arr_n[g][p][e],"asym_g%d_p%d_e%d",g,p,e);
        asym_arr[g][p][e]->Write(asym_arr_n[g][p][e],TObject::kSingleKey);
      };
    };
  };
  asy->Write("asy");
  printf("written\n");
  printf("cases passed: %d   cases failed: %d\n",case_counter,case_fail);
};
