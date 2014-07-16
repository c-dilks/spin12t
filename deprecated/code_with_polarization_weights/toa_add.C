// TObjArray Adder
// 
// combines and organises phi distributions from TObjArrays
// of phiset files
//
// -- reads rtree.root


void toa_add()
{
  // build array of phiset/*.root TFile pointers
  const Int_t MAX_NUM_FILES=200;
  TFile * phi_file[MAX_NUM_FILES]; 
  Int_t phi_file_cnt=0;
  gROOT->ProcessLine(".! ls phiset/phi*.root > toa_files.txt");
  const Int_t filename_buffer=32;
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


  // set bins
  Int_t phi_bins0, eta_bins0, pt_bins0, en_bins0;
  if(gSystem->Getenv("PHI")==NULL){fprintf(stderr,"ERROR: source env vars\n"); return;};
  sscanf(gSystem->Getenv("PHI"),"%d",&phi_bins0);
  sscanf(gSystem->Getenv("ETA"),"%d",&eta_bins0);
  sscanf(gSystem->Getenv("PT"),"%d",&pt_bins0);
  sscanf(gSystem->Getenv("EN"),"%d",&en_bins0);
  const Double_t pi=3.1415;
  const Double_t phi_bins=phi_bins0; 
  const Double_t eta_bins=eta_bins0;
  const Double_t pt_bins=pt_bins0; 
  const Double_t en_bins=en_bins0; 

  // phi_dist [spin] [eta] [pt] [energy] [run index - 1]
  TH1D * phi_dist[4][eta_bins][pt_bins][en_bins][rtree_ent];
  TH1D * phi_dist_dP[4][eta_bins][pt_bins][en_bins][rtree_ent];
  TH1D * phi_dist_mP[4][eta_bins][pt_bins][en_bins][rtree_ent];
  TH1D * phi_dist_mP2[4][eta_bins][pt_bins][en_bins][rtree_ent];
  
  // infile_arr [spin] [eta] [pt] [energy] [phi file]
  TObjArray * infile_arr[4][eta_bins][pt_bins][en_bins][NFILES];
  TObjArray * infile_arr_dP[4][eta_bins][pt_bins][en_bins][NFILES];
  TObjArray * infile_arr_mP[4][eta_bins][pt_bins][en_bins][NFILES];
  TObjArray * infile_arr_mP2[4][eta_bins][pt_bins][en_bins][NFILES];
  char infile_arr_n[4][eta_bins][pt_bins][en_bins][128];
  char infile_arr_dP_n[4][eta_bins][pt_bins][en_bins][128];
  char infile_arr_mP_n[4][eta_bins][pt_bins][en_bins][128];
  char infile_arr_mP2_n[4][eta_bins][pt_bins][en_bins][128];

  Int_t inrun;

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
              sprintf(infile_arr_n[s][g][p][e],"phi_dist_s%d_g%d_p%d_e%d",s,g,p,e);
              sprintf(infile_arr_dP_n[s][g][p][e],"phi_dist_dP_s%d_g%d_p%d_e%d",s,g,p,e);
              sprintf(infile_arr_mP_n[s][g][p][e],"phi_dist_mP_s%d_g%d_p%d_e%d",s,g,p,e);
              sprintf(infile_arr_mP2_n[s][g][p][e],"phi_dist_mP2_s%d_g%d_p%d_e%d",s,g,p,e);
            };
            // read TObjArrays
            infile_arr[s][g][p][e][f] = (TObjArray*)
              phi_file[f]->Get(infile_arr_n[s][g][p][e]);
            infile_arr_dP[s][g][p][e][f] = (TObjArray*)
              phi_file[f]->Get(infile_arr_dP_n[s][g][p][e]);
            infile_arr_mP[s][g][p][e][f] = (TObjArray*)
              phi_file[f]->Get(infile_arr_mP_n[s][g][p][e]);
            infile_arr_mP2[s][g][p][e][f] = (TObjArray*)
              phi_file[f]->Get(infile_arr_mP2_n[s][g][p][e]);
            printf("%s: %s @ %p\n",filename[f],infile_arr_n[s][g][p][e],(void*)infile_arr[s][g][p][e][f]);

            // loop through TObjArrays
            for(Int_t o=0; o<infile_arr[s][g][p][e][f]->GetEntries(); o++)
            {
              // get run number "inrun"
              sscanf(infile_arr[s][g][p][e][f]->At(o)->GetName(),
                "phi_s%*d_g%*d_p%*d_e%*d_r%d",&inrun);
              // linear hash --> typcast phi_dist's
              for(Int_t h=0; h<rtree_ent; h++)
              {
                if(inrun == runnum_arr[h])
                {
                  phi_dist[s][g][p][e][h] = (TH1D*) infile_arr[s][g][p][e][f]->At(o);
                  phi_dist_dP[s][g][p][e][h] = (TH1D*) infile_arr_dP[s][g][p][e][f]->At(o);
                  phi_dist_mP[s][g][p][e][h] = (TH1D*) infile_arr_mP[s][g][p][e][f]->At(o);
                  phi_dist_mP2[s][g][p][e][h] = (TH1D*) infile_arr_mP2[s][g][p][e][f]->At(o);
                };
              };
              //printf("%d/%d\n",o,infile_arr[s][g][p][e][f]->GetEntries());
            };
          };
        };
      };
    };
  };


  // build final TObjArrays, one for each kinematic/geometric bin
  TFile * outfile = new TFile("phiset/all.root","RECREATE");
  TObjArray * combined_array[4][eta_bins][pt_bins][en_bins];
  TObjArray * combined_array_dP[4][eta_bins][pt_bins][en_bins];
  TObjArray * combined_array_mP[4][eta_bins][pt_bins][en_bins];
  TObjArray * combined_array_mP2[4][eta_bins][pt_bins][en_bins];
  char combined_array_n[4][eta_bins][pt_bins][en_bins][128];
  char combined_array_dP_n[4][eta_bins][pt_bins][en_bins][128];
  char combined_array_mP_n[4][eta_bins][pt_bins][en_bins][128];
  char combined_array_mP2_n[4][eta_bins][pt_bins][en_bins][128];
  for(Int_t c=0; c<4; c++)
  {
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            combined_array[s][g][p][e] = new TObjArray();
            combined_array_dP[s][g][p][e] = new TObjArray();
            combined_array_mP[s][g][p][e] = new TObjArray();
            combined_array_mP2[s][g][p][e] = new TObjArray();
            sprintf(combined_array_n[s][g][p][e],"phi_dist_s%d_g%d_p%d_e%d",s,g,p,e);
            sprintf(combined_array_dP_n[s][g][p][e],"phi_dist_dP_s%d_g%d_p%d_e%d",s,g,p,e);
            sprintf(combined_array_mP_n[s][g][p][e],"phi_dist_mP_s%d_g%d_p%d_e%d",s,g,p,e);
            sprintf(combined_array_mP2_n[s][g][p][e],"phi_dist_mP2_s%d_g%d_p%d_e%d",s,g,p,e);
            for(Int_t r=0; r<rtree_ent; r++)
            {
              if(phi_dist[s][g][p][e][r]!=NULL)
              {
                combined_array[s][g][p][e]->AddLast(phi_dist[s][g][p][e][r]);
                combined_array_dP[s][g][p][e]->AddLast(phi_dist_dP[s][g][p][e][r]);
                combined_array_mP[s][g][p][e]->AddLast(phi_dist_mP[s][g][p][e][r]);
                combined_array_mP2[s][g][p][e]->AddLast(phi_dist_mP2[s][g][p][e][r]);
                //printf("phi_dist[%d][%d][%d][%d][%d] entries = %d\n",s,g,p,e,r,
                  //phi_dist[s][g][p][e][r]->GetEntries());
              };
            };
            if(c==0)
              combined_array[s][g][p][e]->Write(combined_array_n[s][g][p][e],TObject::kSingleKey);
            else if(c==1)
              combined_array_dP[s][g][p][e]->Write(combined_array_dP_n[s][g][p][e],TObject::kSingleKey);
            else if(c==2)
              combined_array_mP[s][g][p][e]->Write(combined_array_mP_n[s][g][p][e],TObject::kSingleKey);
            else if(c==3)
              combined_array_mP2[s][g][p][e]->Write(combined_array_mP2_n[s][g][p][e],TObject::kSingleKey);
          };
        };
      };
    };
  };
  printf("phiset/all.root written\n");
}
