// computes asymmetry between fully summed (i.e. over all runs) phi distributions
//    (see staszak thesis, eq. 7.5)
//
//  - this script does all asymmetries; currently A_LL and A_L for both beams, each classified
//    by an "asymmetry number", defined as the Rellum number used to compute the asymmetry:
//    -- asym=1 :: A_L yellow
//    -- asym=2 :: A_L blue
//    -- asym=3 :: A_LL double helicity asymmetry
//  
//  - filter type: defines a filter for the data, useful for consistency checks
//    -- all: entire data set
//    -- run: only runs from filter_low to filter_high
//      -- runout: only runs excluding those from filter_low to filter_high
//      -- runeven: only even run numbers
//      -- runodd: only odd run numbers
//    -- fill: only fills from filter_low to filter_high
//  --> use asym_call* scripts to call different filters
//

void Asym3(const char * jtype="pi0", const char * filter_type="all",Int_t filter_low=0, Int_t filter_high=0)
{
  const Float_t pi=3.1415;
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


  // read TObjArrays of phi distributions
  TObjArray * phi_dist_arr[4][eta_bins][pt_bins][en_bins];
  char phi_dist_arr_n[4][eta_bins][pt_bins][en_bins][64];
  Int_t NRUNS_tmp=0;
  Int_t ARR_size;
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(phi_dist_arr_n[s][g][p][e],"%s/phi_dist_%s_s%d_g%d_p%d_e%d",jtype,jtype,s,g,p,e);
          phi_dist_arr[s][g][p][e] = (TObjArray*) infile->Get(phi_dist_arr_n[s][g][p][e]);
          printf("phi_dist_arr[%d][%d][%d][%d] @ %p\n",s,g,p,e,(void*)phi_dist_arr[s][g][p][e]);
          if(s==0 && g==0 && p==0 && e==0)
          {
            ARR_size=phi_dist_arr[s][g][p][e]->GetEntries();
            for(Int_t kk=0; kk<ARR_size; kk++)
            {
              if(((TH1D*)(phi_dist_arr[s][g][p][e]->At(kk)))->GetEntries() > 0) NRUNS_tmp++;
            }
          }
          else
          {
            if(phi_dist_arr[s][g][p][e]->GetEntries() != ARR_size)
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
  printf("ARR_size=%d\n",ARR_size);
  printf("NRUNS=%d\n",NRUNS);



  // build summed phi distributions with appropriate polarization and rellum weights
  // *_num = numerator in MLM calculation
  // *_den = denominator in MLM calculation
  // *_num_e = numerator in statistical error
  // *_den_e = denominator in statistical error
  const Int_t asym_bins=4;
  TH1D * phi_dist_num[asym_bins][4][eta_bins][pt_bins][en_bins]; // [asymmetry number] [spinbit] [eta] [pt] [en] 
  TH1D * phi_dist_den[asym_bins][4][eta_bins][pt_bins][en_bins];
  TH1D * phi_dist_num_e[asym_bins][4][eta_bins][pt_bins][en_bins]; 
  TH1D * phi_dist_den_e[asym_bins][4][eta_bins][pt_bins][en_bins];
  char phi_dist_num_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_den_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_num_e_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_den_e_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  Double_t yield[4][eta_bins][pt_bins][en_bins]; // counts yields (only for a==3); used for error analysis
  Double_t scarat_num[asym_bins][4][eta_bins][pt_bins][en_bins]; // zdc/vpd scaler ratio; weighted by P
  Double_t scarat_den[asym_bins][4][eta_bins][pt_bins][en_bins]; // zdc/vpd scaler ratio; weighted by P^2
  Int_t runnum;
  Float_t rellum,polar_b,polar_y;
  Float_t weight_num,weight_den; // weight for asymmetry (MLM for A_LL)
  Float_t weight_num_e,weight_den_e; // weight for statistical error (by-product of MLM for A_LL);
  Float_t weight_num_s,weight_den_s; // weight for systematic error (MLM for R_LL)
  Float_t scarat;
  Int_t fill,pattern;
  Bool_t isConsistent;
  gROOT->ProcessLine(".! touch runlist; rm runlist; touch runlist");
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          // initialise arrays
          yield[s][g][p][e] = 0;
          for(Int_t a=1; a<asym_bins; a++)
          {
            scarat_num[a][s][g][p][e] = 0;
            scarat_den[a][s][g][p][e] = 0;
          };
        };
      };
    };
  };
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            if(a==3) yield[s][g][p][e] = 0; // initialise yield counter

            sprintf(phi_dist_num_n[a][s][g][p][e],"phi_num_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            sprintf(phi_dist_den_n[a][s][g][p][e],"phi_den_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            sprintf(phi_dist_num_e_n[a][s][g][p][e],"phi_num_e_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            sprintf(phi_dist_den_e_n[a][s][g][p][e],"phi_den_e_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            phi_dist_num[a][s][g][p][e] = new TH1D(phi_dist_num_n[a][s][g][p][e],phi_dist_num_n[a][s][g][p][e],
              phi_bins,phi_low,phi_high);
            phi_dist_den[a][s][g][p][e] = new TH1D(phi_dist_den_n[a][s][g][p][e],phi_dist_den_n[a][s][g][p][e],
              phi_bins,phi_low,phi_high);
            phi_dist_num_e[a][s][g][p][e] = new TH1D(phi_dist_num_e_n[a][s][g][p][e],phi_dist_num_e_n[a][s][g][p][e],
              phi_bins,phi_low,phi_high);
            phi_dist_den_e[a][s][g][p][e] = new TH1D(phi_dist_den_e_n[a][s][g][p][e],phi_dist_den_e_n[a][s][g][p][e],
              phi_bins,phi_low,phi_high);

            for(Int_t r=0; r<ARR_size; r++)
            {
              if(!strcmp(jtype,"sph"))
                sscanf(phi_dist_arr[0][g][p][e]->At(r)->GetName(), "phi_sph_s%*d_g%*d_p%*d_e%*d_r%d",&runnum);
              else if(!strcmp(jtype,"pi0"))
                sscanf(phi_dist_arr[0][g][p][e]->At(r)->GetName(), "phi_pi0_s%*d_g%*d_p%*d_e%*d_r%d",&runnum);
              else if(!strcmp(jtype,"thr"))
                sscanf(phi_dist_arr[0][g][p][e]->At(r)->GetName(), "phi_thr_s%*d_g%*d_p%*d_e%*d_r%d",&runnum);

              //printf("phi_dist_name=%s\n",phi_dist_arr[0][g][p][e]->At(r)->GetName());
              //printf("jtype=%s runnum=%d\n",jtype,runnum);

              rellum = RD->Rellum(runnum,a,"zdc"); // note that asym no. = rellum no. needed for this asymmetry
              scarat = RD->Scarat(runnum,"x",s);
              //rellum=1;     // for testing
              polar_b = RD->BluePol(runnum);
              polar_y = RD->YellPol(runnum);
              fill = RD->GetFill(runnum);
              isConsistent = RD->RellumConsistent(runnum);
              pattern = RD->Pattern(runnum);


              // polarization and rellum weighting; depends on which asymmetry
              if(a==3)
              {
                weight_num = polar_b * polar_y;
                weight_den = pow(polar_b * polar_y, 2);
                weight_num_e = pow(polar_b * polar_y, 2);
                weight_den_e = pow(polar_b * polar_y, 2);
                weight_num_s = polar_b * polar_y;
                weight_den_s = pow(polar_b * polar_y, 2);
                if(s==1 || s==2) 
                {
                  weight_num *= rellum;
                  weight_den *= rellum;
                  weight_num_e *= pow(rellum, 2);
                  weight_den_e *= rellum;
                  //weight_num_s *= rellum; // should we do this?
                  //weight_den_s *= rellum;
                };
              }
              else if(a==1)
              {
                weight_num = polar_y;
                weight_den = pow(polar_y, 2);
                weight_num_e = pow(polar_y, 2);
                weight_den_e = pow(polar_y, 2);
                weight_num_s = polar_y;
                weight_den_s = pow(polar_y, 2);
                if(s==0 || s==2) 
                {
                  weight_num *= rellum;
                  weight_den *= rellum;
                  weight_num_e *= pow(rellum, 2);
                  weight_den_e *= rellum;
                  //weight_num_s *= rellum; // should we do this?
                  //weight_den_s *= rellum;
                };
              }
              else if(a==2)
              {
                weight_num = polar_b;
                weight_den = pow(polar_b, 2);
                weight_num_e = pow(polar_b, 2);
                weight_den_e = pow(polar_b, 2);
                weight_num_s = polar_b;
                weight_den_s = pow(polar_b, 2);
                if(s==0 || s==1) 
                {
                  weight_num *= rellum;
                  weight_den *= rellum;
                  weight_num_e *= pow(rellum, 2);
                  weight_den_e *= rellum;
                  //weight_num_s *= rellum; // should we do this?
                  //weight_den_s *= rellum;
                };
              };

              
              // print out runlist with fill no. and R3
              if(a==3 && s==0 && g==0 && p==0 && e==0 && ((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)))->GetEntries()>0)
              {
                gSystem->RedirectOutput("runlist");
                printf("%d %d %d %f\n",r,runnum,fill,rellum);
                gSystem->RedirectOutput(0);
              };

              // add weighted dists if we have consistent rellum measurement and filter_type 
              // defined above is passed
              if(isConsistent)
              {
                if( ( !strcmp(filter_type,"fill") && (fill>=filter_low && fill<=filter_high) ) ||
                    ( !strcmp(filter_type,"run") && (runnum>=filter_low && runnum<=filter_high) ) ||
                    ( !strcmp(filter_type,"runout") && !(runnum>=filter_low && runnum<=filter_high) ) ||
                    ( !strcmp(filter_type,"runeven") && (runnum % 2)==0 ) ||
                    ( !strcmp(filter_type,"runodd") && (runnum % 2)==1 ) ||
                      !strcmp(filter_type,"all"))
                {
                  phi_dist_num[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_num);
                  phi_dist_den[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_den);
                  phi_dist_num_e[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_num_e);
                  phi_dist_den_e[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_den_e);
                  if(a==3) yield[s][g][p][e] += ((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)))->GetEntries(); // increment yield counter
                  scarat_num[a][s][g][p][e] += scarat * weight_num_s;
                  scarat_den[a][s][g][p][e] += scarat * weight_den_s;
                };
              };
            };
          };
        };
      };
    };
  };


  // compute asymmetry; polarization and rellum have been corrected for above
  // asym = (ll_num - rr_num) / (ll_den + rr_den)
  // -- a=1 :: ll=s1+s3 :: rr=s0+s2
  // -- a=2 :: ll=s2+s3 :: rr=s0+s1
  // -- a=3 :: ll=s0+s3 :: rr=s1+s2
  // also compute statistical error bars
  // asym_e = sqrt(ll_num_e + rr_num_e) / (ll_den_e + rr_den_e)
  TH1D * dist_ll_num[asym_bins][eta_bins][pt_bins][en_bins]; // left & right asym terms
  TH1D * dist_rr_num[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_ll_den[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_rr_den[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_ll_num_e[asym_bins][eta_bins][pt_bins][en_bins]; // left & right error terms
  TH1D * dist_rr_num_e[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_ll_den_e[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_rr_den_e[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * numer[asym_bins][eta_bins][pt_bins][en_bins]; // numer = ll_num - rr_num
  TH1D * denom[asym_bins][eta_bins][pt_bins][en_bins]; // denom = ll_den + rr_den
  TH1D * numer_e[asym_bins][eta_bins][pt_bins][en_bins]; // numer_e = ll_num_e + rr_num_e
  TH1D * denom_e[asym_bins][eta_bins][pt_bins][en_bins]; // denom_e = ll_den_e + rr_den_e
  TH1D * numer_e_sqrt[asym_bins][eta_bins][pt_bins][en_bins]; // numer_e_sqrt = sqrt(numer_e)
  TH1D * asym[asym_bins][eta_bins][pt_bins][en_bins]; // asym = numer / denom
  TH1D * asym_e[asym_bins][eta_bins][pt_bins][en_bins]; // asym_e = numer_e_sqrt / denom_e
  Int_t asym_pts[asym_bins][eta_bins][pt_bins][en_bins]; // number of points in asym

  char dist_ll_num_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_num_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_ll_den_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_den_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_ll_num_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_num_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_ll_den_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_den_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char numer_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char denom_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char numer_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char denom_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char numer_e_sqrt_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char asym_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char asym_t[asym_bins][eta_bins][pt_bins][en_bins][256];
  char asym_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];

  Float_t p0,p0e,chi2,ndf;
  Float_t bc[4];
  Float_t bcent;
  Double_t bc_e,bc_ee;
  Int_t runnum_0;
  Float_t asym_value[asym_bins][eta_bins][pt_bins][en_bins];
  Int_t asym_value_cnt[asym_bins][eta_bins][pt_bins][en_bins];
  TF1 * asym_fit[asym_bins][eta_bins][pt_bins][en_bins];
  Float_t asym_tmp;
  Float_t asym_max[asym_bins][eta_bins][pt_bins][en_bins];
  Float_t asym_min[asym_bins][eta_bins][pt_bins][en_bins];
  char var_str[16]; strcpy(var_str,"#phi");

  Double_t RLL_ll_num[asym_bins][eta_bins][pt_bins][en_bins];
  Double_t RLL_rr_num[asym_bins][eta_bins][pt_bins][en_bins];
  Double_t RLL_ll_den[asym_bins][eta_bins][pt_bins][en_bins];
  Double_t RLL_rr_den[asym_bins][eta_bins][pt_bins][en_bins];
  Double_t RLL_numer[asym_bins][eta_bins][pt_bins][en_bins];
  Double_t RLL_denom[asym_bins][eta_bins][pt_bins][en_bins];
  Double_t RLL_asym[asym_bins][eta_bins][pt_bins][en_bins];

  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          asym_value[a][g][p][e]=0.0;
          asym_value_cnt[a][g][p][e]=0;
          asym_max[a][g][p][e]=0;
          asym_min[a][g][p][e]=0;
          RLL_ll_num[a][g][p][e]=0;
          RLL_rr_num[a][g][p][e]=0;
          RLL_ll_den[a][g][p][e]=0;
          RLL_rr_den[a][g][p][e]=0;
          RLL_numer[a][g][p][e]=0;
          RLL_denom[a][g][p][e]=0;
          RLL_asym[a][g][p][e]=0;
        };
      };
    };
  };
  char asym_title[asym_bins][16];
  strcpy(asym_title[1],"A_{L}^{Y}");
  strcpy(asym_title[2],"A_{L}^{B}");
  strcpy(asym_title[3],"A_{LL}");
  char asym_title_scarat[asym_bins][16];
  strcpy(asym_title_scarat[1],"R_{L}^{Y}");
  strcpy(asym_title_scarat[2],"R_{L}^{B}");
  strcpy(asym_title_scarat[3],"R_{LL}");
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          // initialise ll,rr,numer,denom,asym histograms
          sprintf(dist_ll_num_n[a][g][p][e],"dist_ll_num_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_num_n[a][g][p][e],"dist_rr_num_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_ll_den_n[a][g][p][e],"dist_ll_den_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_den_n[a][g][p][e],"dist_rr_den_a%d_g%d_p%d_e%d",a,g,p,e);

          sprintf(dist_ll_num_e_n[a][g][p][e],"dist_ll_num_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_num_e_n[a][g][p][e],"dist_rr_num_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_ll_den_e_n[a][g][p][e],"dist_ll_den_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_den_e_n[a][g][p][e],"dist_rr_den_e_a%d_g%d_p%d_e%d",a,g,p,e);

          sprintf(numer_n[a][g][p][e],"numer_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(denom_n[a][g][p][e],"denom_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(numer_e_n[a][g][p][e],"numer_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(denom_e_n[a][g][p][e],"denom_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(numer_e_sqrt_n[a][g][p][e],"numer_e_sqrt_a%d_g%d_p%d_e%d",a,g,p,e);

          sprintf(asym_n[a][g][p][e],"asym_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(asym_e_n[a][g][p][e],"asym_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(asym_t[a][g][p][e],
             "%s vs. %s :: #eta#in[%.2f,%.2f), p_{T}#in[%.2f,%.2f), E#in[%.2f,%.2f) (runsum)",
             asym_title[a],var_str,eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1]);


          dist_ll_num[a][g][p][e] = new TH1D(dist_ll_num_n[a][g][p][e],
              dist_ll_num_n[a][g][p][e],phi_bins,phi_low,phi_high);
          dist_rr_num[a][g][p][e] = new TH1D(dist_rr_num_n[a][g][p][e],
              dist_rr_num_n[a][g][p][e],phi_bins,phi_low,phi_high);
          dist_ll_den[a][g][p][e] = new TH1D(dist_ll_den_n[a][g][p][e],
              dist_ll_den_n[a][g][p][e],phi_bins,phi_low,phi_high);
          dist_rr_den[a][g][p][e] = new TH1D(dist_rr_den_n[a][g][p][e],
              dist_rr_den_n[a][g][p][e],phi_bins,phi_low,phi_high);
          
          dist_ll_num_e[a][g][p][e] = new TH1D(dist_ll_num_e_n[a][g][p][e],
              dist_ll_num_e_n[a][g][p][e],phi_bins,phi_low,phi_high);
          dist_rr_num_e[a][g][p][e] = new TH1D(dist_rr_num_e_n[a][g][p][e],
              dist_rr_num_e_n[a][g][p][e],phi_bins,phi_low,phi_high);
          dist_ll_den_e[a][g][p][e] = new TH1D(dist_ll_den_e_n[a][g][p][e],
              dist_ll_den_e_n[a][g][p][e],phi_bins,phi_low,phi_high);
          dist_rr_den_e[a][g][p][e] = new TH1D(dist_rr_den_e_n[a][g][p][e],
              dist_rr_den_e_n[a][g][p][e],phi_bins,phi_low,phi_high);

          numer[a][g][p][e] = new TH1D(numer_n[a][g][p][e],
              numer_n[a][g][p][e],phi_bins,phi_low,phi_high);
          denom[a][g][p][e] = new TH1D(denom_n[a][g][p][e],
              denom_n[a][g][p][e],phi_bins,phi_low,phi_high);
          numer_e[a][g][p][e] = new TH1D(numer_e_n[a][g][p][e],
              numer_e_n[a][g][p][e],phi_bins,phi_low,phi_high);
          denom_e[a][g][p][e] = new TH1D(denom_e_n[a][g][p][e],
              denom_e_n[a][g][p][e],phi_bins,phi_low,phi_high);

          numer_e_sqrt[a][g][p][e] = new TH1D(numer_e_sqrt_n[a][g][p][e],
              numer_e_sqrt_n[a][g][p][e],phi_bins,phi_low,phi_high);

          asym[a][g][p][e] = new TH1D(asym_n[a][g][p][e],
              asym_t[a][g][p][e],phi_bins,phi_low,phi_high);
          asym_e[a][g][p][e] = new TH1D(asym_e_n[a][g][p][e],
              asym_e_n[a][g][p][e],phi_bins,phi_low,phi_high);

          // build left & right terms
          if(a==3)
          {
            dist_ll_num[a][g][p][e]->Add(phi_dist_num[a][0][g][p][e],phi_dist_num[a][3][g][p][e],1.0,1.0);
            dist_rr_num[a][g][p][e]->Add(phi_dist_num[a][1][g][p][e],phi_dist_num[a][2][g][p][e],1.0,1.0);
            dist_ll_den[a][g][p][e]->Add(phi_dist_den[a][0][g][p][e],phi_dist_den[a][3][g][p][e],1.0,1.0);
            dist_rr_den[a][g][p][e]->Add(phi_dist_den[a][1][g][p][e],phi_dist_den[a][2][g][p][e],1.0,1.0);
            dist_ll_num_e[a][g][p][e]->Add(phi_dist_num_e[a][0][g][p][e],phi_dist_num_e[a][3][g][p][e],1.0,1.0);
            dist_rr_num_e[a][g][p][e]->Add(phi_dist_num_e[a][1][g][p][e],phi_dist_num_e[a][2][g][p][e],1.0,1.0);
            dist_ll_den_e[a][g][p][e]->Add(phi_dist_den_e[a][0][g][p][e],phi_dist_den_e[a][3][g][p][e],1.0,1.0);
            dist_rr_den_e[a][g][p][e]->Add(phi_dist_den_e[a][1][g][p][e],phi_dist_den_e[a][2][g][p][e],1.0,1.0);
            RLL_ll_num[a][g][p][e] = scarat_num[a][0][g][p][e] + scarat_num[a][3][g][p][e];
            RLL_rr_num[a][g][p][e] = scarat_num[a][1][g][p][e] + scarat_num[a][2][g][p][e];
            RLL_ll_den[a][g][p][e] = scarat_den[a][0][g][p][e] + scarat_den[a][3][g][p][e];
            RLL_rr_den[a][g][p][e] = scarat_den[a][1][g][p][e] + scarat_den[a][2][g][p][e];
          }
          else if(a==1)
          {
            dist_ll_num[a][g][p][e]->Add(phi_dist_num[a][1][g][p][e],phi_dist_num[a][3][g][p][e],1.0,1.0);
            dist_rr_num[a][g][p][e]->Add(phi_dist_num[a][0][g][p][e],phi_dist_num[a][2][g][p][e],1.0,1.0);
            dist_ll_den[a][g][p][e]->Add(phi_dist_den[a][1][g][p][e],phi_dist_den[a][3][g][p][e],1.0,1.0);
            dist_rr_den[a][g][p][e]->Add(phi_dist_den[a][0][g][p][e],phi_dist_den[a][2][g][p][e],1.0,1.0);
            dist_ll_num_e[a][g][p][e]->Add(phi_dist_num_e[a][1][g][p][e],phi_dist_num_e[a][3][g][p][e],1.0,1.0);
            dist_rr_num_e[a][g][p][e]->Add(phi_dist_num_e[a][0][g][p][e],phi_dist_num_e[a][2][g][p][e],1.0,1.0);
            dist_ll_den_e[a][g][p][e]->Add(phi_dist_den_e[a][1][g][p][e],phi_dist_den_e[a][3][g][p][e],1.0,1.0);
            dist_rr_den_e[a][g][p][e]->Add(phi_dist_den_e[a][0][g][p][e],phi_dist_den_e[a][2][g][p][e],1.0,1.0);
            RLL_ll_num[a][g][p][e] = scarat_num[a][1][g][p][e] + scarat_num[a][3][g][p][e];
            RLL_rr_num[a][g][p][e] = scarat_num[a][0][g][p][e] + scarat_num[a][2][g][p][e];
            RLL_ll_den[a][g][p][e] = scarat_den[a][1][g][p][e] + scarat_den[a][3][g][p][e];
            RLL_rr_den[a][g][p][e] = scarat_den[a][0][g][p][e] + scarat_den[a][2][g][p][e];
          }
          else if(a==2)
          {
            dist_ll_num[a][g][p][e]->Add(phi_dist_num[a][2][g][p][e],phi_dist_num[a][3][g][p][e],1.0,1.0);
            dist_rr_num[a][g][p][e]->Add(phi_dist_num[a][0][g][p][e],phi_dist_num[a][1][g][p][e],1.0,1.0);
            dist_ll_den[a][g][p][e]->Add(phi_dist_den[a][2][g][p][e],phi_dist_den[a][3][g][p][e],1.0,1.0);
            dist_rr_den[a][g][p][e]->Add(phi_dist_den[a][0][g][p][e],phi_dist_den[a][1][g][p][e],1.0,1.0);
            dist_ll_num_e[a][g][p][e]->Add(phi_dist_num_e[a][2][g][p][e],phi_dist_num_e[a][3][g][p][e],1.0,1.0);
            dist_rr_num_e[a][g][p][e]->Add(phi_dist_num_e[a][0][g][p][e],phi_dist_num_e[a][1][g][p][e],1.0,1.0);
            dist_ll_den_e[a][g][p][e]->Add(phi_dist_den_e[a][2][g][p][e],phi_dist_den_e[a][3][g][p][e],1.0,1.0);
            dist_rr_den_e[a][g][p][e]->Add(phi_dist_den_e[a][0][g][p][e],phi_dist_den_e[a][1][g][p][e],1.0,1.0);
            RLL_ll_num[a][g][p][e] = scarat_num[a][2][g][p][e] + scarat_num[a][3][g][p][e];
            RLL_rr_num[a][g][p][e] = scarat_num[a][0][g][p][e] + scarat_num[a][1][g][p][e];
            RLL_ll_den[a][g][p][e] = scarat_den[a][2][g][p][e] + scarat_den[a][3][g][p][e];
            RLL_rr_den[a][g][p][e] = scarat_den[a][0][g][p][e] + scarat_den[a][1][g][p][e];
          }

          // build numer and denom
          numer[a][g][p][e]->Add(dist_ll_num[a][g][p][e],dist_rr_num[a][g][p][e],1.0,-1.0);
          denom[a][g][p][e]->Add(dist_ll_den[a][g][p][e],dist_rr_den[a][g][p][e],1.0,1.0);
          numer_e[a][g][p][e]->Add(dist_ll_num_e[a][g][p][e],dist_rr_num_e[a][g][p][e],1.0,1.0); // n.b.: ll & rr summed
          denom_e[a][g][p][e]->Add(dist_ll_den_e[a][g][p][e],dist_rr_den_e[a][g][p][e],1.0,1.0);
          RLL_numer[a][g][p][e] = RLL_ll_num[a][g][p][e] - RLL_rr_num[a][g][p][e];
          RLL_denom[a][g][p][e] = RLL_ll_den[a][g][p][e] + RLL_rr_den[a][g][p][e];
          
          // numer_e_sqrt = sqrt(numer_e)
          for(Int_t b=1; b<=numer_e[a][g][p][e]->GetNbinsX(); b++)
          {
            bc_e = numer_e[a][g][p][e]->GetBinContent(b);
            bc_e = sqrt(bc_e);
            numer_e_sqrt[a][g][p][e]->SetBinContent(b,bc_e);
          };

          asym[a][g][p][e]->Divide(numer[a][g][p][e],denom[a][g][p][e],1.0,1.0);
          asym_e[a][g][p][e]->Divide(numer_e_sqrt[a][g][p][e],denom_e[a][g][p][e],1.0,1.0);
          RLL_asym[a][g][p][e] = RLL_numer[a][g][p][e] / RLL_denom[a][g][p][e];

          printf(asym[a][g][p][e]->GetTitle());
          printf("\n");

          
          // set statistical errors of asym bins using values computed in asym_e
          for(Int_t b=1; b<=asym[a][g][p][e]->GetNbinsX(); b++)
          {
            bc_e = asym_e[a][g][p][e]->GetBinContent(b);
            // compare computed stat error with ROOT's error propagation
            //bc_ee = asym[a][g][p][e]->GetBinError(b);
            //gSystem->RedirectOutput("error_check.dat","a");
            //printf("%.20f\n",bc_ee-bc_e); 
            //gSystem->RedirectOutput(0);
            asym[a][g][p][e]->SetBinError(b,bc_e);
          };

          // n.b. for one phi bin, constant fit & error matches the bin & its error
          asym[a][g][p][e]->Fit("pol0","Q","",phi_low,phi_high);
          asym_fit[a][g][p][e] = asym[a][g][p][e]->GetFunction("pol0");
          if(asym_fit[a][g][p][e]!=NULL)
          {
            asym_value[a][g][p][e]+=asym_fit[a][g][p][e]->GetParameter(0);
            asym_value_cnt[a][g][p][e]++;
          };
        };
      };
    };
  };


  // set plot ranges
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          //asym[a][g][p][e]->GetYaxis()->SetRangeUser(2*asym_min[a][g][p][e],2*asym_max[a][g][p][e]);
          asym[a][g][p][e]->GetXaxis()->SetRangeUser(phi_low,phi_high);
        };
      };
    };
  };


  // average asym_value for each bin
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          asym_value[a][g][p][e] /= ((Float_t)asym_value_cnt[a][g][p][e]);
          printf("g%d p%d e%d <%s>=%f\n",g,p,e,asym_title[a],asym_value[a][g][p][e]);
        };
      };
    };
  };


  // kinematic dependence plots
  TGraphErrors * en_dep[asym_bins][eta_bins][pt_bins]; // en dependent plots, one for each pt bin (statistical errors)
  TGraphErrors * pt_dep[asym_bins][eta_bins][en_bins]; // pt dependent plots, one for each en bin
  TGraphErrors * en_sys[asym_bins][eta_bins][pt_bins]; // en dependent plots, one for each pt bin (systematic errors)
  TGraphErrors * pt_sys[asym_bins][eta_bins][en_bins]; // pt dependent plots, one for each en bin
  char en_dep_t[asym_bins][eta_bins][pt_bins][256];
  char pt_dep_t[asym_bins][eta_bins][en_bins][256];
  char en_sys_t[asym_bins][eta_bins][pt_bins][256];
  char pt_sys_t[asym_bins][eta_bins][en_bins][256];
  Int_t en_dep_cnt[asym_bins][eta_bins][pt_bins]; // point counter
  Int_t pt_dep_cnt[asym_bins][eta_bins][en_bins];
  Int_t en_sys_cnt[asym_bins][eta_bins][pt_bins];
  Int_t pt_sys_cnt[asym_bins][eta_bins][en_bins];
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        en_dep_cnt[a][g][p]=0;
        en_sys_cnt[a][g][p]=0;
      };
      for(Int_t e=0; e<en_bins; e++) 
      {
        pt_dep_cnt[a][g][e]=0;
        pt_sys_cnt[a][g][e]=0;
      };
    };
  };

  Double_t val_en[asym_bins][eta_bins][pt_bins][en_bins];     // arrays for en dependent plots, one for each pt bin
  Double_t err_en[asym_bins][eta_bins][pt_bins][en_bins];     // (statistical error)
  Double_t sys_en[asym_bins][eta_bins][pt_bins][en_bins];     // (systematic error)
  Double_t cent_en[asym_bins][eta_bins][pt_bins][en_bins];    // (energy bin center)
  Double_t width_en[asym_bins][eta_bins][pt_bins][en_bins];   // (energy bin width)
  Double_t zeroz_en[asym_bins][eta_bins][pt_bins][en_bins];   // array of zeros

  Double_t val_pt[asym_bins][eta_bins][en_bins][pt_bins];     // arrays for pt dependent plots, one
  Double_t err_pt[asym_bins][eta_bins][en_bins][pt_bins];     // (statistical error)
  Double_t sys_pt[asym_bins][eta_bins][en_bins][pt_bins];     // (systematic error)
  Double_t cent_pt[asym_bins][eta_bins][en_bins][pt_bins];    // (pt bin center)
  Double_t width_pt[asym_bins][eta_bins][en_bins][pt_bins];   // (pt bin width)
  Double_t zeroz_pt[asym_bins][eta_bins][en_bins][pt_bins];   // array of zeros

  // en dependent points for each pt bin
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          zeroz_en[a][g][p][e]=0;
          if(asym[a][g][p][e]->GetFunction("pol0"))
          {
            val_en[a][g][p][en_dep_cnt[a][g][p]] = asym[a][g][p][e]->GetFunction("pol0")->GetParameter(0);
            // error of fit to asym; asym statistical errors set via asym_e
              err_en[a][g][p][en_dep_cnt[a][g][p]] = asym[a][g][p][e]->GetFunction("pol0")->GetParError(0);
            // estimated statistical error
              //err_en[a][g][p][en_dep_cnt[a][g][p]]=1/(0.55*0.55)*1/sqrt(yield[0][g][p][e]+yield[1][g][p][e]+yield[2][g][p][e]+yield[3][g][p][e]);
            sys_en[a][g][p][en_dep_cnt[a][g][p]] = RLL_asym[a][g][p][e]; // MLM R_LL, R_Lb, R_Ly
            cent_en[a][g][p][en_dep_cnt[a][g][p]] = en_div[e] + ((en_div[e+1]-en_div[e])/2.0);
            width_en[a][g][p][en_dep_cnt[a][g][p]] = (en_div[e+1]-en_div[e])/2.0;
            en_dep_cnt[a][g][p]++;
          };
        };
        // asymmetry vs. en with statistical error bars
          en_dep[a][g][p] = new TGraphErrors(en_dep_cnt[a][g][p],cent_en[a][g][p],val_en[a][g][p],width_en[a][g][p],err_en[a][g][p]);
          sprintf(en_dep_t[a][g][p],"%s #pm #sigma %s vs. E_{#gamma#gamma} for p_{T}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",asym_title[a],
                  asym_title[a],pt_div[p],pt_div[p+1],eta_div[g],eta_div[g+1]);
          
        // --o--

        // set systematic errors on asymmetry nominal vals:
          /*
          for(Int_t e=0; e<en_bins; e++) sys_en[a][g][p][e] = fabs(sys_en[a][g][p][e]);
          en_sys[a][g][p] = new TGraphErrors(en_dep_cnt[a][g][p],cent_en[a][g][p],val_en[a][g][p],width_en[a][g][p],sys_en[a][g][p]);
          sprintf(en_sys_t[a][g][p],"%s #pm %s vs. E_{#gamma#gamma} for p_{T}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",asym_title[a],
                  asym_title_scarat[a],pt_div[p],pt_div[p+1],eta_div[g],eta_div[g+1]);
          */
        // nominal value is R_L* (scarat asymmetry); error bars currently set to zero
          en_sys[a][g][p] = new TGraphErrors(en_dep_cnt[a][g][p],cent_en[a][g][p],sys_en[a][g][p],width_en[a][g][p],zeroz_en[a][g][p]);
          sprintf(en_sys_t[a][g][p],"%s vs. E_{#gamma#gamma} for p_{T}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",asym_title_scarat[a],
                  pt_div[p],pt_div[p+1],eta_div[g],eta_div[g+1]);
        
        en_dep[a][g][p]->SetTitle(en_dep_t[a][g][p]);
        en_sys[a][g][p]->SetTitle(en_sys_t[a][g][p]);
        en_dep[a][g][p]->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
        en_sys[a][g][p]->GetXaxis()->SetTitle("E_{#gamma#gamma} (GeV)");
      };
    };
  };

  // pt dependent points for each en bin
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          zeroz_pt[a][g][e][p]=0;
          if(asym[a][g][p][e]->GetFunction("pol0"))
          {
            val_pt[a][g][e][pt_dep_cnt[a][g][e]] = asym[a][g][p][e]->GetFunction("pol0")->GetParameter(0);
            // error of fit to asym; asym statistical errors set via asym_e
              err_pt[a][g][e][pt_dep_cnt[a][g][e]] = asym[a][g][p][e]->GetFunction("pol0")->GetParError(0);
            // estimated statistical error
              //err_pt[a][g][e][pt_dep_cnt[a][g][e]]=1/(0.55*0.55)*1/sqrt(yield[0][g][p][e]+yield[1][g][p][e]+yield[2][g][p][e]+yield[3][g][p][e]);
            sys_pt[a][g][e][pt_dep_cnt[a][g][e]] = RLL_asym[a][g][p][e]; // MLM R_LL, R_Lb, R_Ly
            cent_pt[a][g][e][pt_dep_cnt[a][g][e]] = pt_div[p] + ((pt_div[p+1]-pt_div[p])/2.0);
            width_pt[a][g][e][pt_dep_cnt[a][g][e]] = (pt_div[p+1]-pt_div[p])/2.0;
            pt_dep_cnt[a][g][e]++;
          };
        };
        // asymmetry vs. pt with statistical error bars
          pt_dep[a][g][e] = new TGraphErrors(pt_dep_cnt[a][g][e],cent_pt[a][g][e],val_pt[a][g][e],width_pt[a][g][e],err_pt[a][g][e]);
          sprintf(pt_dep_t[a][g][e],"%s #pm #sigma %s vs. p_{T} for E_{#gamma#gamma}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",asym_title[a],
                  asym_title[a],en_div[e],en_div[e+1],eta_div[g],eta_div[g+1]);

        // --o--

        // set systematic errors on asymmetry nominal vals:
          /*
          for(Int_t p=0; p<pt_bins; p++) sys_pt[a][g][e][p] = fabs(sys_pt[a][g][e][p]);
          pt_sys[a][g][e] = new TGraphErrors(pt_dep_cnt[a][g][e],cent_pt[a][g][e],val_pt[a][g][e],width_pt[a][g][e],sys_pt[a][g][e]);
          sprintf(pt_sys_t[a][g][e],"%s #pm %s vs. p_{T} for E_{#gamma#gamma}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",asym_title[a],
                  asym_title_scarat[a],en_div[e],en_div[e+1],eta_div[g],eta_div[g+1]);
          */
        // nominal value is R_L* (scarat asymmetry); error bars currently set to zero
          pt_sys[a][g][e] = new TGraphErrors(pt_dep_cnt[a][g][e],cent_pt[a][g][e],sys_pt[a][g][e],width_pt[a][g][e],zeroz_pt[a][g][e]);
          sprintf(pt_sys_t[a][g][e],"%s vs. p_{T} for E_{#gamma#gamma}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",asym_title_scarat[a],
                  en_div[e],en_div[e+1],eta_div[g],eta_div[g+1]);

        pt_dep[a][g][e]->SetTitle(pt_dep_t[a][g][e]);
        pt_sys[a][g][e]->SetTitle(pt_sys_t[a][g][e]);
        pt_dep[a][g][e]->GetXaxis()->SetTitle("p_{T} (GeV)");
        pt_sys[a][g][e]->GetXaxis()->SetTitle("p_{T} (GeV)");
      };
    };
  };

  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        en_dep[a][g][p]->GetYaxis()->SetTitle(asym_title[a]);
        en_dep[a][g][p]->SetMarkerStyle(kFullCircle);
        en_dep[a][g][p]->SetMarkerColor(kRed);
        en_dep[a][g][p]->GetYaxis()->SetTitleOffset(1.5);
        en_sys[a][g][p]->GetYaxis()->SetTitle(asym_title[a]);
        en_sys[a][g][p]->SetMarkerStyle(kFullCircle);
        en_sys[a][g][p]->SetMarkerColor(kRed);
        en_sys[a][g][p]->GetYaxis()->SetTitleOffset(1.5);
      };
      for(Int_t e=0; e<en_bins; e++)
      {
        pt_dep[a][g][e]->GetYaxis()->SetTitle(asym_title[a]);
        pt_dep[a][g][e]->SetMarkerStyle(kFullCircle);
        pt_dep[a][g][e]->SetMarkerColor(kRed);
        pt_dep[a][g][e]->GetYaxis()->SetTitleOffset(1.5);
        pt_sys[a][g][e]->GetYaxis()->SetTitle(asym_title[a]);
        pt_sys[a][g][e]->SetMarkerStyle(kFullCircle);
        pt_sys[a][g][e]->SetMarkerColor(kRed);
        pt_sys[a][g][e]->GetYaxis()->SetTitleOffset(1.5);
      };
    };
  };


  // write phi dists
  printf("writing spin.root...\n");
  TFile * outfile = new TFile("spin.root","RECREATE");
  outfile->mkdir("A_LL");
  outfile->mkdir("A_L_blue");
  outfile->mkdir("A_L_yellow");
  char en_dep_n[asym_bins][eta_bins][pt_bins][32];
  char pt_dep_n[asym_bins][eta_bins][en_bins][32];
  char en_sys_n[asym_bins][eta_bins][pt_bins][32];
  char pt_sys_n[asym_bins][eta_bins][en_bins][32];
  for(Int_t a=1; a<asym_bins; a++)
  {
    if(a==3) outfile->cd("/A_LL");
    else if(a==1) outfile->cd("/A_L_yellow");
    else if(a==2) outfile->cd("/A_L_blue");
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        sprintf(en_dep_n[a][g][p],"en_dep_a%d_g%d_p%d",a,g,p);
        en_dep[a][g][p]->Write(en_dep_n[a][g][p]);
      };
      for(Int_t e=0; e<en_bins; e++)
      {
        sprintf(pt_dep_n[a][g][e],"pt_dep_a%d_g%d_e%d",a,g,e);
        pt_dep[a][g][e]->Write(pt_dep_n[a][g][e]);
      };
    };
  };
  for(Int_t a=1; a<asym_bins; a++)
  {
    if(a==3) outfile->cd("/A_LL");
    else if(a==1) outfile->cd("/A_L_yellow");
    else if(a==2) outfile->cd("/A_L_blue");
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        sprintf(en_sys_n[a][g][p],"en_sys_a%d_g%d_p%d",a,g,p);
        en_sys[a][g][p]->Write(en_sys_n[a][g][p]);
      };
      for(Int_t e=0; e<en_bins; e++)
      {
        sprintf(pt_sys_n[a][g][e],"pt_sys_a%d_g%d_e%d",a,g,e);
        pt_sys[a][g][e]->Write(pt_sys_n[a][g][e]);
      };
    };
  };

  for(Int_t a=1; a<asym_bins; a++)
  {
    if(a==3) outfile->cd("/A_LL");
    else if(a==1) outfile->cd("/A_L_yellow");
    else if(a==2) outfile->cd("/A_L_blue");
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          for(Int_t s=0; s<4; s++)
          {
            phi_dist_num[a][s][g][p][e]->Write(phi_dist_num_n[a][s][g][p][e]);
            phi_dist_den[a][s][g][p][e]->Write(phi_dist_den_n[a][s][g][p][e]);
          };
        };
      };
    };
  };
  for(Int_t a=1; a<asym_bins; a++)
  {
    if(a==3) outfile->cd("/A_LL");
    else if(a==1) outfile->cd("/A_L_yellow");
    else if(a==2) outfile->cd("/A_L_blue");
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          asym[a][g][p][e]->Write(asym_n[a][g][p][e]);
        };
      };
    };
  };
  printf("written\n");
};
