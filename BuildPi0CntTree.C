// builds Pi0Cnt tree, which includes 
// - runnum
// - fill
// - number of pi0s
// - number of pi0s in low pT (1-5)
// - number of pi0s in mid pT (6-7)
// - number of pi0s in high pT (8-10)
// - rellum
// - polarisation
// - pattern

void BuildPi0CntTree(const char * filename="RedOutputset107ha.root")
{
  char outfilename[128];
  sscanf(filename,"RedOutputset%s",outfilename);
  sprintf(outfilename,"pi0cntset/pi0cnt%s",outfilename);
  sprintf(filename,"redset/%s",filename);
  printf("infile = %s\n",filename);
  printf("outfile = %s\n",outfilename);
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("str");
  tr->Print();
  printf("tr @ %p\n",(void*)tr);

  Int_t runnum_tmp=0;
  Int_t runnum,fill;
  Int_t count,count_l,count_m,count_h;
  Float_t R_bbc[10];
  Float_t R_zdc[10];
  Float_t R_vpd[10];
  char R_bbc_str[10][16];
  char R_zdc_str[10][16];
  char R_vpd_str[10][16];
  char R_bbc_typ[10][16];
  char R_zdc_typ[10][16];
  char R_vpd_typ[10][16];
  Float_t b_pol,y_pol;
  Int_t pattern;

  Int_t TrigBits;
  Float_t N12;
  Float_t Z,M12,E12,Pt;
  Bool_t kicked;

  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("b_pol",&b_pol);
  tr->SetBranchAddress("y_pol",&y_pol);
  tr->SetBranchAddress("pattern",&pattern);
  tr->SetBranchAddress("N12",&N12);
  for(Int_t r=1; r<10; r++)
  {
    sprintf(R_bbc_str[r],"R%d_bbc",r);
    sprintf(R_zdc_str[r],"R%d_zdc",r);
    sprintf(R_vpd_str[r],"R%d_vpd",r);
    sprintf(R_bbc_typ[r],"R%d_bbc/F",r);
    sprintf(R_zdc_typ[r],"R%d_zdc/F",r);
    sprintf(R_vpd_typ[r],"R%d_vpd/F",r);
    tr->SetBranchAddress(R_bbc_str[r],&(R_bbc[r]));
    tr->SetBranchAddress(R_zdc_str[r],&(R_zdc[r]));
    tr->SetBranchAddress(R_vpd_str[r],&(R_vpd[r]));
  };
  tr->SetBranchAddress("TrigBits",&TrigBits);
  tr->SetBranchAddress("Z",&Z);
  tr->SetBranchAddress("M12",&M12);
  tr->SetBranchAddress("E12",&E12);
  tr->SetBranchAddress("Pt",&Pt);
  tr->SetBranchAddress("kicked",&kicked);
  printf("branch addresses set\n");


  TFile * outfile = new TFile(outfilename,"RECREATE");
  TTree * ptr = new TTree();
  ptr->Branch("runnum",&runnum,"runnum/I");
  ptr->Branch("fill",&fill,"fill/I");
  ptr->Branch("count",&count,"count/I"); // no. pi0s per run
  ptr->Branch("count_l",&count_l,"count_l/I"); // no. pi0s per run pT 1-5 GeV
  ptr->Branch("count_m",&count_m,"count_m/I"); // no. pi0s per run pT 6-7 GeV
  ptr->Branch("count_h",&count_h,"count_h/I"); // no. pi0s per run pT 8-10 GeV
  for(Int_t r=1; r<10; r++)
  {
    ptr->Branch(R_bbc_str[r],&(R_bbc[r]),R_bbc_typ[r]);
    ptr->Branch(R_zdc_str[r],&(R_zdc[r]),R_zdc_typ[r]);
    ptr->Branch(R_vpd_str[r],&(R_vpd[r]),R_vpd_typ[r]);
  };
  ptr->Branch("b_pol",&b_pol,"b_pol/F");
  ptr->Branch("y_pol",&b_pol,"y_pol/F");
  ptr->Branch("pattern",&pattern,"pattern/I");
  printf("ptr defined\n");


  count = count_l = count_m = count_h = 0;
  tr->GetEntry(0);
  runnum_tmp = runnum;
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);

    if(runnum!=runnum_tmp || i+1==tr->GetEntries())
    {
      runnum_tmp = runnum;
      tr->GetEntry(i-1);
      ptr->Fill();
      count = count_l = count_m = count_h = 0;
      tr->GetEntry(i);
    };

    // pi0 cut
    if((TrigBits&0x200) && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && N12==2)
    {
      count++;
      if(Pt>=1 && Pt<=5) count_l++;
      if(Pt>=6 && Pt<=7) count_m++;
      if(Pt>=8 && Pt<=10) count_h++;
    };
  };
  ptr->Write("ptr");
};




