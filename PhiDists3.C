// -- create phi distributions with various kinematic cuts (eta,pt,E) for each spinbit
//    (kinematic cuts are set as Double_t's below)
// -- phi distributions for "+ -" and "- +" are weighted by relative luminosity
//    A_LL = 1/(P_b*P_y)*(N++ + N-- - RN+- -RN-+)/(N++ + N-- + RN+- + RN-+)
// -- phi distributions are written to phiset/ directory with similar name; 
//    they are named: phi_s[spinbit]_g[eta bin]_p[pt bin]_e[en bin]

void PhiDists3(const char * filename="RedOutputset132ha.root")
{
  // load polarization and rellum data
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


  // read redset tree and set output file name
  Int_t runnum,bx,blue,yell,TrigBits,pattern;
  Float_t M12,N12,E12,Z,Phi,Eta,Pt,b_pol,y_pol;
  Bool_t kicked,isConsistent;
  char setname[32];
  sscanf(filename,"RedOutputset%s",setname);
  sprintf(filename,"redset/%s",filename);
  TFile * infile = new TFile(filename,"READ");
  char outname[256];
  sprintf(outname,"phiset/phi%s",setname);
  TTree * tree = (TTree*) infile->Get("str");
  tree->SetBranchAddress("runnum",&runnum);
  tree->SetBranchAddress("Bunchid7bit",&bx);
  tree->SetBranchAddress("TrigBits",&TrigBits);
  tree->SetBranchAddress("M12",&M12);
  tree->SetBranchAddress("N12",&N12);
  tree->SetBranchAddress("E12",&E12);
  tree->SetBranchAddress("Z",&Z);
  tree->SetBranchAddress("Phi",&Phi);
  tree->SetBranchAddress("Eta",&Eta);
  tree->SetBranchAddress("Pt",&Pt);


  // define spinbit strings
  char spinbit_t[4][4];
  char spinbit_n[4][4];
  strcpy(spinbit_t[0],"--"); strcpy(spinbit_n[0],"nn");
  strcpy(spinbit_t[1],"-+"); strcpy(spinbit_n[1],"np");
  strcpy(spinbit_t[2],"+-"); strcpy(spinbit_n[2],"pn");
  strcpy(spinbit_t[3],"++"); strcpy(spinbit_n[3],"pp");


  // get number of runs in the redset file and build array of run numbers
  Int_t runnum_arr[30]; // (n.b. maximum size arbitrarily defined)
  Int_t runnum_tmp=0;
  Int_t NRUNS_tmp = 0;
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


  // define phi distributions  for each spinbit and kinematic bin
  // _sph = single photon
  // _pi0 = neutral pions
  // _thr = >=3 photon jet
  // ---> n.b. I can't use higher dimensional arrays
  TH1D * phi_dist_sph[4][eta_bins][pt_bins][en_bins][NRUNS]; // [spinbit] [eta] [pt] [en] [run number]
  TH1D * phi_dist_pi0[4][eta_bins][pt_bins][en_bins][NRUNS]; 
  TH1D * phi_dist_thr[4][eta_bins][pt_bins][en_bins][NRUNS]; 

  char phi_dist_sph_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];
  char phi_dist_pi0_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];
  char phi_dist_thr_n[4][eta_bins][pt_bins][en_bins][NRUNS][256];

  char phi_dist_sph_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];
  char phi_dist_pi0_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];
  char phi_dist_thr_t[4][eta_bins][pt_bins][en_bins][NRUNS][400];

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
            sprintf(phi_dist_sph_n[s][g][p][e][r],"phi_sph_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);
            sprintf(phi_dist_pi0_n[s][g][p][e][r],"phi_pi0_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);
            sprintf(phi_dist_thr_n[s][g][p][e][r],"phi_thr_s%d_g%d_p%d_e%d_r%d",s,g,p,e,runnum_arr[r]);

            sprintf(phi_dist_sph_t[s][g][p][e][r],
             "single-#gamma #phi distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);
            sprintf(phi_dist_pi0_t[s][g][p][e][r],
             "#pi^{0} #phi distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);
            sprintf(phi_dist_thr_t[s][g][p][e][r],
             "N_{#gamma}>3 jet #phi distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
             spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1],runnum_arr[r]);

            phi_dist_sph[s][g][p][e][r] = new TH1D(phi_dist_sph_n[s][g][p][e][r],phi_dist_sph_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);
            phi_dist_pi0[s][g][p][e][r] = new TH1D(phi_dist_pi0_n[s][g][p][e][r],phi_dist_pi0_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);
            phi_dist_thr[s][g][p][e][r] = new TH1D(phi_dist_thr_n[s][g][p][e][r],phi_dist_thr_t[s][g][p][e][r],
              phi_bins,phi_low,phi_high);

            phi_dist_sph[s][g][p][e][r]->Sumw2();
            phi_dist_pi0[s][g][p][e][r]->Sumw2();
            phi_dist_thr[s][g][p][e][r]->Sumw2();
          };
        }; 
      };
    };
  };


  // fill phi distributions
  Int_t ss,gg,pp,ee,rr;
  rr=-1; runnum_tmp=0;
  printf("fill phi dists...\n");
  for(Int_t x=0; x<tree->GetEntries(); x++)
  {
    if((x%10000)==0) printf("%.2f%%\n",100*((Float_t)x)/((Float_t)tree->GetEntries()));
    ss=gg=pp=ee=-1; // reset 
    tree->GetEntry(x);
    blue = RD->BlueSpin(runnum,bx);
    yell = RD->YellSpin(runnum,bx);
    kicked = RD->Kicked(runnum,bx);

    // spin --> array index
    if(blue==-1 && yell==-1) ss=0;
    else if(blue==-1 && yell==1) ss=1;
    else if(blue==1 && yell==-1) ss=2;
    else if(blue==1 && yell==1) ss=3;

    // kinematic bins --> array indices
    for(Int_t g=0; g<eta_bins; g++) { if(Eta>=eta_div[g] && Eta<eta_div[g+1]) gg=g; };
    for(Int_t p=0; p<pt_bins;  p++) { if(Pt>=pt_div[p]   && Pt<pt_div[p+1]  ) pp=p; };
    for(Int_t e=0; e<en_bins;  e++) { if(E12>=en_div[e]  && E12<en_div[e+1] ) ee=e; };

    // run number --> array index
    if(runnum != runnum_tmp)
    {
      rr=-1;
      for(Int_t r=0; r<NRUNS; r++) { if(runnum_arr[r] == runnum) rr=r; }; 
      runnum_tmp = runnum;
      isConsistent = RD->RellumConsistent(runnum);
      pattern = RD->Pattern(runnum);
      b_pol = RD->BluePol(runnum);
      y_pol = RD->YellPol(runnum);
    };

    // check for valid array indices (filters out events outside kinematic boundaries)
    //printf("%d %d %d %d %d\n",ss,gg,pp,ee,rr);
    if(ss>=0 && gg>=0 && pp>=0 && ee>=0 && rr>=0)
    {
      // rellum consistency & polarization cut
      if( kicked==0 && isConsistent==1 && b_pol>0 && y_pol>0)
      {
        // n photon cut
        if( fabs(N12-1)<0.01 ) phi_dist_sph[ss][gg][pp][ee][rr]->Fill(Phi);

        // pi0 cut
        else if( (TrigBits&0x200) && 
                  fabs(N12-2)<0.01 &&
                  Z<0.8 &&
                  fabs(M12-0.135)<0.1) phi_dist_pi0[ss][gg][pp][ee][rr]->Fill(Phi);

        // >=3 photon jet cut
        else if( N12>2 ) phi_dist_thr[ss][gg][pp][ee][rr]->Fill(Phi);
      };
    };
  };
        

  // make object arrays
  TFile * outfile = new TFile(outname,"RECREATE");
  TObjArray * phi_dist_sph_arr[4][eta_bins][pt_bins][en_bins];
  TObjArray * phi_dist_pi0_arr[4][eta_bins][pt_bins][en_bins];
  TObjArray * phi_dist_thr_arr[4][eta_bins][pt_bins][en_bins];
  char phi_dist_sph_arr_name[4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_pi0_arr_name[4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_thr_arr_name[4][eta_bins][pt_bins][en_bins][128];
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_sph_arr[s][g][p][e] = new TObjArray();
          phi_dist_pi0_arr[s][g][p][e] = new TObjArray();
          phi_dist_thr_arr[s][g][p][e] = new TObjArray();

          sprintf(phi_dist_sph_arr_name[s][g][p][e],"phi_dist_sph_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_pi0_arr_name[s][g][p][e],"phi_dist_pi0_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_thr_arr_name[s][g][p][e],"phi_dist_thr_s%d_g%d_p%d_e%d",s,g,p,e);
          
          for(Int_t r=0; r<NRUNS; r++)
          {
            phi_dist_sph_arr[s][g][p][e]->AddLast(phi_dist_sph[s][g][p][e][r]);
            phi_dist_pi0_arr[s][g][p][e]->AddLast(phi_dist_pi0[s][g][p][e][r]);
            phi_dist_thr_arr[s][g][p][e]->AddLast(phi_dist_thr[s][g][p][e][r]);
          };
        };
      };
    };
  };

  
  // write object arrays
  outfile->mkdir("sph");
  outfile->cd("sph");
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_sph_arr[s][g][p][e]->Write(phi_dist_sph_arr_name[s][g][p][e],TObject::kSingleKey);
        };
      };
    };
  };
  outfile->cd();
  outfile->mkdir("pi0");
  outfile->cd("pi0");
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_pi0_arr[s][g][p][e]->Write(phi_dist_pi0_arr_name[s][g][p][e],TObject::kSingleKey);
        };
      };
    };
  };
  outfile->cd();
  outfile->mkdir("thr");
  outfile->cd("thr");
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist_thr_arr[s][g][p][e]->Write(phi_dist_thr_arr_name[s][g][p][e],TObject::kSingleKey);
        };
      };
    };
  };
  outfile->cd();
};
