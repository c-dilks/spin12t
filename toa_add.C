// TObjArray Adder
// 
// combines and organises phi distributions from TObjArrays
// of phiset files
//
// -- reads rtree.root


void toa_add()
{
  gROOT->Reset();
  // build array of phiset/*.root TFile pointers
  const Int_t MAX_NUM_FILES=200;
  TFile * phi_file[MAX_NUM_FILES]; 
  Int_t phi_file_cnt=0;
  gROOT->ProcessLine(".! ls phiset/phi*.root > toa_files.txt");
  const Int_t filename_buffer=64;
  char filename[MAX_NUM_FILES][filename_buffer];
  char temp[filename_buffer];
  FILE * toa_files;
  toa_files = fopen("toa_files.txt","r");
  if(toa_files==NULL)
  {
    fprintf(stderr,"Error opening toa_files.txt\n");
    return;
  }
  else
  {
    while(!feof(toa_files))
    {
      fgets(filename[phi_file_cnt],filename_buffer,toa_files);

      // fgets reads in "returns"; this hack gets rid of them
      sscanf(filename[phi_file_cnt],"%s",filename[phi_file_cnt]);

      if(strcmp(filename[phi_file_cnt],""))
      {
        printf("%d: %s\n",phi_file_cnt,filename[phi_file_cnt]);
        phi_file[phi_file_cnt] = new TFile(filename[phi_file_cnt],"READ");
        phi_file_cnt++;
      };
    };
  };
  const Int_t NFILES=phi_file_cnt;
  gROOT->ProcessLine(".! rm toa_files.txt");



  // build run number hash table
  TFile * rtree_file = new TFile("rtree.root","READ");
  TTree * rtree = (TTree*) rtree_file->Get("rellum");
  Int_t index,runnum;
  Int_t rtree_ent_tmp = rtree->GetEntries();
  const Int_t rtree_ent = rtree_ent_tmp;
  Int_t runnum_arr[rtree_ent];
  rtree->SetBranchAddress("i",&index);
  rtree->SetBranchAddress("runnum",&runnum);
  for(Int_t i=0; i<rtree->GetEntries(); i++)
  {
    rtree->GetEntry(i);
    if(i+1 == index) runnum_arr[i] = runnum;
    else 
    {
      fprintf(stderr,"ERROR: rtree file problem\n");
      return;
    };
  };

  
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


  // phi_dist [spin] [eta] [pt] [energy] [run index - 1]
  TH1D * phi_dist_sph[4][eta_bins][pt_bins][en_bins][rtree_ent];
  TH1D * phi_dist_pi0[4][eta_bins][pt_bins][en_bins][rtree_ent];
  TH1D * phi_dist_thr[4][eta_bins][pt_bins][en_bins][rtree_ent];
  
  // infile_arr [spin] [eta] [pt] [energy] [phi file]
  TObjArray * infile_sph_arr[4][eta_bins][pt_bins][en_bins][NFILES];
  TObjArray * infile_pi0_arr[4][eta_bins][pt_bins][en_bins][NFILES];
  TObjArray * infile_thr_arr[4][eta_bins][pt_bins][en_bins][NFILES];
  char infile_sph_arr_n[4][eta_bins][pt_bins][en_bins][200];
  char infile_pi0_arr_n[4][eta_bins][pt_bins][en_bins][200];
  char infile_thr_arr_n[4][eta_bins][pt_bins][en_bins][200];

  Int_t inrun;
  Bool_t filter_sph[rtree_ent];
  Bool_t filter_pi0[rtree_ent];
  Bool_t filter_thr[rtree_ent];
  for(Int_t rr=0; rr<rtree_ent; rr++) 
  {
    filter_sph[rr]=false;
    filter_pi0[rr]=false;
    filter_thr[rr]=false;
  };

  for(Int_t f=0; f<NFILES; f++)
  {
    phi_file[f]->cd(); // focus on next TFile
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            // set up TObjArray names to be read (only needs one execution)
            if(f==0)
            {
              sprintf(infile_sph_arr_n[s][g][p][e],"/sph/phi_dist_sph_s%d_g%d_p%d_e%d",s,g,p,e);
              sprintf(infile_pi0_arr_n[s][g][p][e],"/pi0/phi_dist_pi0_s%d_g%d_p%d_e%d",s,g,p,e);
              sprintf(infile_thr_arr_n[s][g][p][e],"/thr/phi_dist_thr_s%d_g%d_p%d_e%d",s,g,p,e);
            };
            // read TObjArrays
            infile_sph_arr[s][g][p][e][f] = (TObjArray*)phi_file[f]->Get(infile_sph_arr_n[s][g][p][e]);
            infile_pi0_arr[s][g][p][e][f] = (TObjArray*)phi_file[f]->Get(infile_pi0_arr_n[s][g][p][e]);
            infile_thr_arr[s][g][p][e][f] = (TObjArray*)phi_file[f]->Get(infile_thr_arr_n[s][g][p][e]);
            printf("%s: %p %p %p\n",filename[f],(void*)infile_sph_arr[s][g][p][e][f],
                                                (void*)infile_pi0_arr[s][g][p][e][f],
                                                (void*)infile_thr_arr[s][g][p][e][f]);

            // loop through sph TObjArrays
            for(Int_t o=0; o<infile_sph_arr[s][g][p][e][f]->GetEntries(); o++)
            {
              // get run number "inrun"
              sscanf(infile_sph_arr[s][g][p][e][f]->At(o)->GetName(),
                "phi_sph_s%*d_g%*d_p%*d_e%*d_r%d",&inrun);
              // linear hash --> typcast phi_dist's
              for(Int_t h=0; h<rtree_ent; h++)
              {
                if(inrun == runnum_arr[h])
                {
                  phi_dist_sph[s][g][p][e][h] = (TH1D*) infile_sph_arr[s][g][p][e][f]->At(o);
                  filter_sph[h]=true;
                };
              };
            };

            // loop through pi0 TObjArrays
            for(Int_t o=0; o<infile_pi0_arr[s][g][p][e][f]->GetEntries(); o++)
            {
              // get run number "inrun"
              sscanf(infile_pi0_arr[s][g][p][e][f]->At(o)->GetName(),
                "phi_pi0_s%*d_g%*d_p%*d_e%*d_r%d",&inrun);
              // linear hash --> typcast phi_dist's
              for(Int_t h=0; h<rtree_ent; h++)
              {
                if(inrun == runnum_arr[h])
                {
                  phi_dist_pi0[s][g][p][e][h] = (TH1D*) infile_pi0_arr[s][g][p][e][f]->At(o);
                  filter_pi0[h]=true;
                };
              };
            };

            // loop through thr TObjArrays
            for(Int_t o=0; o<infile_thr_arr[s][g][p][e][f]->GetEntries(); o++)
            {
              // get run number "inrun"
              sscanf(infile_thr_arr[s][g][p][e][f]->At(o)->GetName(),
                "phi_thr_s%*d_g%*d_p%*d_e%*d_r%d",&inrun);
              // linear hash --> typcast phi_dist's
              for(Int_t h=0; h<rtree_ent; h++)
              {
                if(inrun == runnum_arr[h])
                {
                  phi_dist_thr[s][g][p][e][h] = (TH1D*) infile_thr_arr[s][g][p][e][f]->At(o);
                  filter_thr[h]=true;
                };
              };
            };

          };
        };
      };
    };
  };


  // build final TObjArrays, one for each kinematic/geometric bin
  TFile * outfile = new TFile("phiset/all.root","RECREATE");
  outfile->cd();
  outfile->mkdir("sph");
  outfile->mkdir("pi0");
  outfile->mkdir("thr");
  TObjArray * combined_sph_array[4][eta_bins][pt_bins][en_bins];
  TObjArray * combined_pi0_array[4][eta_bins][pt_bins][en_bins];
  TObjArray * combined_thr_array[4][eta_bins][pt_bins][en_bins];
  char combined_sph_array_n[4][eta_bins][pt_bins][en_bins][200];
  char combined_pi0_array_n[4][eta_bins][pt_bins][en_bins][200];
  char combined_thr_array_n[4][eta_bins][pt_bins][en_bins][200];

  printf("--------------------------------------------------------\n");

  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          combined_sph_array[s][g][p][e] = new TObjArray();
          combined_pi0_array[s][g][p][e] = new TObjArray();
          combined_thr_array[s][g][p][e] = new TObjArray();

          sprintf(combined_sph_array_n[s][g][p][e],"phi_dist_sph_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(combined_pi0_array_n[s][g][p][e],"phi_dist_pi0_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(combined_thr_array_n[s][g][p][e],"phi_dist_thr_s%d_g%d_p%d_e%d",s,g,p,e);

          for(Int_t r=0; r<rtree_ent; r++)
          {
            printf("%d.%p %d.%p %d.%p\n",filter_sph[r],(void*)phi_dist_sph[s][g][p][e][r],
                                         filter_pi0[r],(void*)phi_dist_pi0[s][g][p][e][r],
                                         filter_thr[r],(void*)phi_dist_thr[s][g][p][e][r]);
            if(phi_dist_sph[s][g][p][e][r]!=NULL && filter_sph[r]==1)
              combined_sph_array[s][g][p][e]->AddLast(phi_dist_sph[s][g][p][e][r]);
            if(phi_dist_pi0[s][g][p][e][r]!=NULL && filter_pi0[r]==1)
              combined_pi0_array[s][g][p][e]->AddLast(phi_dist_pi0[s][g][p][e][r]);
            if(phi_dist_thr[s][g][p][e][r]!=NULL && filter_thr[r]==1)
              combined_thr_array[s][g][p][e]->AddLast(phi_dist_thr[s][g][p][e][r]);
          };
          outfile->cd("/sph"); combined_sph_array[s][g][p][e]->Write(combined_sph_array_n[s][g][p][e],TObject::kSingleKey);
          outfile->cd("/pi0"); combined_pi0_array[s][g][p][e]->Write(combined_pi0_array_n[s][g][p][e],TObject::kSingleKey);
          outfile->cd("/thr"); combined_thr_array[s][g][p][e]->Write(combined_thr_array_n[s][g][p][e],TObject::kSingleKey);
        };
      };
    };
  };
  printf("phiset/all.root written\n");
}
