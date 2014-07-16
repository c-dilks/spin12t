// -- create phi distributions with various kinematic cuts (eta,pt,E) for each spinbit
//    (kinematic cuts are set as Double_t's below)
// -- phi distributions for "+ -" and "- +" are weighted by relative luminosity
//    A_LL = 1/(P_b*P_y)*(N++ + N-- - RN+- -RN-+)/(N++ + N-- + RN+- + RN-+)
// -- phi distributions are written to phiset/ directory with similar name; 
//    they are named: phi_s[spinbit]_g[eta bin]_p[pt bin]_e[en bin]

void PhiDists(const char * filename="RedOutputset132ha.root",Int_t phi_bins0=16,
              Int_t eta_bins0=1, Int_t pt_bins0=3, Int_t en_bins0=3)
{
  // set output file name and read input tree
  char setname[32];
  sscanf(filename,"RedOutputset%s",setname);
  sprintf(filename,"redset/%s",filename);
  TFile * infile = new TFile(filename,"READ");
  TTree * tree = (TTree*) infile->Get("str");
  char outname[256];
  sprintf(outname,"phiset/phi%s",setname);


  // set ranges (set # bins = 1 to average over all values from low to high)
  const Double_t pi=3.1415;
  const Double_t phi_bins=phi_bins0; 
    const Double_t phi_low=-1*pi; const Double_t phi_high=pi;
  const Double_t eta_bins=eta_bins0;
    const Double_t eta_low=2.6; const Double_t eta_high=4.2;
  const Double_t pt_bins=pt_bins0; 
    const Double_t pt_low=0; const Double_t pt_high=10;
  const Double_t en_bins=en_bins0; 
    const Double_t en_low=0; const Double_t en_high=100;


  // define spinbit strings
  char spinbit_t[4][4];
  char spinbit_n[4][4];
  strcpy(spinbit_t[0],"--"); strcpy(spinbit_n[0],"nn");
  strcpy(spinbit_t[1],"-+"); strcpy(spinbit_n[1],"np");
  strcpy(spinbit_t[2],"+-"); strcpy(spinbit_n[2],"pn");
  strcpy(spinbit_t[3],"++"); strcpy(spinbit_n[3],"pp");


  // set binning
  // -- *_div = lower limit of each bin; last extra one is the upper limit of last bin
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

  // print ranges
  /*
  for(Int_t g=0; g<eta_bins+1; g++) printf("eta_div[%d] = %.2f\n",g,eta_div[g]);
  for(Int_t p=0; p<pt_bins+1; p++) printf("pt_div[%d] = %.2f\n",p,pt_div[p]);
  for(Int_t e=0; e<en_bins+1; e++) printf("en_div[%d] = %.2f\n",e,en_div[e]);
  */

  // get number of runs in the redset file and build array of run numbers
  Int_t runnum;
  Int_t runnum_arr[30]; // (n.b. maximum size arbitrarily defined)
  Int_t runnum_tmp=0;
  Int_t NRUNS_tmp = 0;
  tree->SetBranchAddress("runnum",&runnum);
  if(tree->GetEntries() == 0)
  {
    fprintf(stderr,"ERROR: no str entries for %s\n --> phi file not produced!\n",filename);
    return;
  };
  for(Int_t i=0; i<tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    if(runnum != runnum_tmp)
    {
      runnum_tmp = runnum;
      if(NRUNS_tmp<30) runnum_arr[NRUNS_tmp] = runnum;
      else 
      {
        fprintf(stderr,"ERROR: more than 30 runs in root file; increase arbitrarily defined max\n");
        return;
      };
      NRUNS_tmp++;
    };
  };
  const Int_t NRUNS = NRUNS_tmp;
  for(Int_t r=0; r<NRUNS; r++) printf("%d\n",runnum_arr[r]);
  printf("NRUNS=%d\n",NRUNS);


  // define phi distributions 
  // - phi_dist -- phi distribution for each spin bit and kinematic bin
  // - phi_dist_dP -- phi distribution divided by polarisation (fill-by-fill basis)
  // - phi_dist_mP -- phi distribution multiplied by polarisation
  // - phi_dist_mP2 -- phi distribution multiplied by square of polarization
  TH1D * phi_dist[4][eta_bins][pt_bins][en_bins][NRUNS]; // [spinbit] [eta] [pt] [en] [run number]
  char phi_dist_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];
  char phi_dist_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];
  TH1D * phi_dist_dP[4][eta_bins][pt_bins][en_bins][NRUNS]; // [spinbit] [eta] [pt] [en] [run number]
  char phi_dist_dP_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];
  char phi_dist_dP_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];
  TH1D * phi_dist_mP[4][eta_bins][pt_bins][en_bins][NRUNS]; // [spinbit] [eta] [pt] [en] [run number]
  char phi_dist_mP_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];
  char phi_dist_mP_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];
  TH1D * phi_dist_mP2[4][eta_bins][pt_bins][en_bins][NRUNS]; // [spinbit] [eta] [pt] [en] [run number]
  char phi_dist_mP2_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];
  char phi_dist_mP2_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];
  for(Int_t r=0; r<NRUNS; r++)
  {
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            sprintf(phi_dist_n[s][g][p][e][r],"phi_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);
            sprintf(phi_dist_t[s][g][p][e][r],
             "#phi distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);
            phi_dist[s][g][p][e][r] = new TH1D(phi_dist_n[s][g][p][e][r],phi_dist_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);

            sprintf(phi_dist_dP_n[s][g][p][e][r],"phi_dP_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);
            sprintf(phi_dist_dP_t[s][g][p][e][r],
             "#phi/(P^{B}P^{Y}) distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);
            phi_dist_dP[s][g][p][e][r] = new TH1D(phi_dist_dP_n[s][g][p][e][r],phi_dist_dP_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);

            sprintf(phi_dist_mP_n[s][g][p][e][r],"phi_mP_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);
            sprintf(phi_dist_mP_t[s][g][p][e][r],
             "#phi/(P^{B}P^{Y}) distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);
            phi_dist_mP[s][g][p][e][r] = new TH1D(phi_dist_mP_n[s][g][p][e][r],phi_dist_mP_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);

            sprintf(phi_dist_mP2_n[s][g][p][e][r],"phi_mP2_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);
            sprintf(phi_dist_mP2_t[s][g][p][e][r],
             "#phi/(P^{B}P^{Y}) distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);
            phi_dist_mP2[s][g][p][e][r] = new TH1D(phi_dist_mP2_n[s][g][p][e][r],phi_dist_mP2_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);
          };
        }; 
      };
    };
  };


  // define cuts
  char kinematic_cut[4][eta_bins][pt_bins][en_bins][NRUNS][512];
  char kinematic_cut_dP[4][eta_bins][pt_bins][en_bins][NRUNS][512]; // divided by polarisation
  char kinematic_cut_mP[4][eta_bins][pt_bins][en_bins][NRUNS][512]; // mulitplied by polarisation
  char kinematic_cut_mP2[4][eta_bins][pt_bins][en_bins][NRUNS][512]; // multiplied by polarisation squared
  char spin_cut[4][16];
  char eta_cut[eta_bins][32];
  char pt_cut[pt_bins][32];
  char en_cut[en_bins][32];
  char run_cut[NRUNS][32];
  char trig_cut[32]; sprintf(trig_cut,"(TrigBits&0x200)"); // other cuts
  char z_cut[16]; sprintf(z_cut,"Z<0.8");
  char mass_cut[32]; sprintf(mass_cut,"abs(M12-0.135)<0.1");
  // run 13 patterns: 13, 14, 23, 24, 31, 32, 41, 42
  char pattern_cut[32]; 
    //strcpy(pattern_cut,"pattern==42");
    //strcpy(pattern_cut,"(pattern==14 || pattern==24)");
    strcpy(pattern_cut,"1"); // (no pattern cut)
  for(Int_t r=0; r<NRUNS; r++)
  {
    sprintf(run_cut[r],"runnum==%d",runnum_arr[r]);
    for(Int_t s=0; s<4; s++)
    {
      sprintf(spin_cut[s],"spin==%d",s);
      for(Int_t g=0; g<eta_bins; g++)
      {
        sprintf(eta_cut[g],"Eta>=%.2f && Eta<%.2f",eta_div[g],eta_div[g+1]);
        for(Int_t p=0; p<pt_bins; p++)
        {
          sprintf(pt_cut[p],"Pt>=%.2f && Pt<%.2f",pt_div[p],pt_div[p+1]);
          for(Int_t e=0; e<en_bins; e++)
          {
            sprintf(en_cut[e],"E12>=%.2f && E12<%.2f",en_div[e],en_div[e+1]);
            sprintf(kinematic_cut[s][g][p][e][r],"%s && %s && %s && %s && %s && %s && %s && %s && isConsistent && b_pol*y_pol!=0 && %s",
              spin_cut[s],eta_cut[g],pt_cut[p],en_cut[e],trig_cut,z_cut,mass_cut,run_cut[r],pattern_cut);
            // weight + - and - + distributions by rellum -- DEPRECATED; moved to Asym script
            //if(s==1||s==2) sprintf(kinematic_cut[s][g][p][e][r],"R3_zdc*(%s)",kinematic_cut[s][g][p][e][r]);

            sprintf(kinematic_cut_dP[s][g][p][e][r],"1/(b_pol*y_pol)*(%s)",kinematic_cut[s][g][p][e][r]);
            sprintf(kinematic_cut_mP[s][g][p][e][r],"(b_pol*y_pol)*(%s)",kinematic_cut[s][g][p][e][r]);
            sprintf(kinematic_cut_mP2[s][g][p][e][r],"(b_pol*b_pol*y_pol*y_pol)*(%s)",kinematic_cut[s][g][p][e][r]);

            printf("kinematic_cut[%d][%d][%d][%d][%d] = %s\n\n",s,g,p,e,r,kinematic_cut[s][g][p][e][r]);
            //printf("kinematic_cut_dP[%d][%d][%d][%d][%d] = %s\n\n",s,g,p,e,r,kinematic_cut_dP[s][g][p][e][r]);
            //printf("kinematic_cut_mP[%d][%d][%d][%d][%d] = %s\n\n",s,g,p,e,r,kinematic_cut_mP[s][g][p][e][r]);
            //printf("kinematic_cut_mP2[%d][%d][%d][%d][%d] = %s\n\n",s,g,p,e,r,kinematic_cut_mP2[s][g][p][e][r]);
          };
        };
      };
    };
  };


  // projections
  for(Int_t r=0; r<NRUNS; r++)
  {
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            phi_dist[s][g][p][e][r]->Sumw2();
            phi_dist_dP[s][g][p][e][r]->Sumw2();
            phi_dist_mP[s][g][p][e][r]->Sumw2();
            phi_dist_mP2[s][g][p][e][r]->Sumw2();
            tree->Project(phi_dist_n[s][g][p][e][r],"Phi",kinematic_cut[s][g][p][e][r]);
            tree->Project(phi_dist_dP_n[s][g][p][e][r],"Phi",kinematic_cut_dP[s][g][p][e][r]);
            tree->Project(phi_dist_mP_n[s][g][p][e][r],"Phi",kinematic_cut_mP[s][g][p][e][r]);
            tree->Project(phi_dist_mP2_n[s][g][p][e][r],"Phi",kinematic_cut_mP2[s][g][p][e][r]);
            printf("phi_dist[%d][%d][%d][%d][%d] entries projected: %d\n",s,g,p,e,r,
              phi_dist[s][g][p][e][r]->GetEntries());
          };
        };
      };
    };
  };


  // make object arrays
  TFile * outfile = new TFile(outname,"RECREATE");
  TObjArray * phi_dist_arr[4][eta_bins][pt_bins][en_bins];
  TObjArray * phi_dist_dP_arr[4][eta_bins][pt_bins][en_bins];
  TObjArray * phi_dist_mP_arr[4][eta_bins][pt_bins][en_bins];
  TObjArray * phi_dist_mP2_arr[4][eta_bins][pt_bins][en_bins];
  char phi_dist_arr_name[4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_dP_arr_name[4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_mP_arr_name[4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_mP2_arr_name[4][eta_bins][pt_bins][en_bins][128];
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_arr[s][g][p][e] = new TObjArray();
          phi_dist_dP_arr[s][g][p][e] = new TObjArray();
          phi_dist_mP_arr[s][g][p][e] = new TObjArray();
          phi_dist_mP2_arr[s][g][p][e] = new TObjArray();

          sprintf(phi_dist_arr_name[s][g][p][e],"phi_dist_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_dP_arr_name[s][g][p][e],"phi_dist_dP_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_mP_arr_name[s][g][p][e],"phi_dist_mP_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_mP2_arr_name[s][g][p][e],"phi_dist_mP2_s%d_g%d_p%d_e%d",s,g,p,e);
          
          for(Int_t r=0; r<NRUNS; r++)
          {
            phi_dist_arr[s][g][p][e]->AddLast(phi_dist[s][g][p][e][r]);
            phi_dist_dP_arr[s][g][p][e]->AddLast(phi_dist_dP[s][g][p][e][r]);
            phi_dist_mP_arr[s][g][p][e]->AddLast(phi_dist_mP[s][g][p][e][r]);
            phi_dist_mP2_arr[s][g][p][e]->AddLast(phi_dist_mP2[s][g][p][e][r]);
          };

        };
      };
    };
  };

  
  // write object arrays
  for(Int_t c=0; c<4; c++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t s=0; s<4; s++)
          {
            if(c==0)
              phi_dist_arr[s][g][p][e]->Write(phi_dist_arr_name[s][g][p][e],TObject::kSingleKey);
            else if(c==1)
              phi_dist_dP_arr[s][g][p][e]->Write(phi_dist_dP_arr_name[s][g][p][e],TObject::kSingleKey);
            else if(c==2)
              phi_dist_mP_arr[s][g][p][e]->Write(phi_dist_mP_arr_name[s][g][p][e],TObject::kSingleKey);
            else if(c==3)
              phi_dist_mP2_arr[s][g][p][e]->Write(phi_dist_mP2_arr_name[s][g][p][e],TObject::kSingleKey);
          };
        };
      };
    };
  };
};
