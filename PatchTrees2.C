// adds relative luminosity and polarization data to reduced data trees
// --> outputs modified redset files to subdirectory "patchset"
// --> this was used to patch old redset files with more complete rellum info and will
//     be deprecated when the ReduceData script is modified to make use of RunData12 methods

void PatchTrees2(const char * filename="RedOutputset156ha.root")
{
  char filename_in[512];
  sprintf(filename_in,"redset/%s",filename);
  TFile * redfile = new TFile(filename_in,"READ");
  TTree * red_tr = (TTree*) redfile->Get("str");

  Float_t M12;
  Float_t N12;
  Float_t E12;
  Float_t Z;
  Float_t Phi;
  Float_t Eta;
  Float_t Pt;
  Int_t spin;
  Int_t blue;
  Int_t yell;
  Bool_t kicked;
  Int_t Bunchid7bit;
  Int_t TrigBits;
  Int_t runnum;
  Int_t fill;
  Int_t pattern;

  red_tr->SetBranchAddress("M12",&M12);
  red_tr->SetBranchAddress("N12",&N12);
  red_tr->SetBranchAddress("E12",&E12);
  red_tr->SetBranchAddress("Z",&Z);
  red_tr->SetBranchAddress("Phi",&Phi);
  red_tr->SetBranchAddress("Eta",&Eta);
  red_tr->SetBranchAddress("Pt",&Pt);
  red_tr->SetBranchAddress("spin",&spin);
  red_tr->SetBranchAddress("blue",&blue);
  red_tr->SetBranchAddress("yell",&yell);
  //red_tr->SetBranchAddress("kicked",&kicked);
  red_tr->SetBranchAddress("Bunchid7bit",&Bunchid7bit);
  red_tr->SetBranchAddress("TrigBits",&TrigBits);
  red_tr->SetBranchAddress("runnum",&runnum);
  red_tr->SetBranchAddress("fill",&fill);
  red_tr->SetBranchAddress("pattern",&pattern);

  gSystem->Load("src/RunData12.so");
  RunData12 * RD = new RunData12();

  char patchname[512];
  sprintf(patchname,"patchset/%s",filename);
  TFile * patchfile = new TFile(patchname,"RECREATE");
  TTree * patch_tr = new TTree("str","str");

  Float_t R_bbc[10];
  Float_t R_zdc[10];
  Float_t R_vpd[10];
  Float_t R_bbc_err[10];
  Float_t R_zdc_err[10];
  Float_t R_vpd_err[10];
  Bool_t isConsistent;
  Float_t b_pol;
  Float_t y_pol;
  Float_t b_pol_err;
  Float_t y_pol_err;

  patch_tr->Branch("M12",&M12,"M12/F");
  patch_tr->Branch("N12",&N12,"N12/F");
  patch_tr->Branch("E12",&E12,"E12/F");
  patch_tr->Branch("Z",&Z,"Z/F");
  patch_tr->Branch("Phi",&Phi,"Phi/F");
  patch_tr->Branch("Eta",&Eta,"Eta/F");
  patch_tr->Branch("Pt",&Pt,"Pt/F");
  patch_tr->Branch("spin",&spin,"spin/I");
  patch_tr->Branch("blue",&blue,"blue/I");
  patch_tr->Branch("yell",&yell,"yell/I");
  patch_tr->Branch("kicked",&kicked,"kicked/O");
  patch_tr->Branch("Bunchid7bit",&Bunchid7bit,"Bunchid7bit/I");
  patch_tr->Branch("TrigBits",&TrigBits,"TrigBits/I");
  patch_tr->Branch("runnum",&runnum,"runnum/I");
  patch_tr->Branch("fill",&fill,"fill/I");

  char R_bbc_n[10][32];
  char R_zdc_n[10][32];
  char R_vpd_n[10][32];
  char R_bbc_err_n[10][32];
  char R_zdc_err_n[10][32];
  char R_vpd_err_n[10][32];
  char R_bbc_d[10][32];
  char R_zdc_d[10][32];
  char R_vpd_d[10][32];
  char R_bbc_err_d[10][32];
  char R_zdc_err_d[10][32];
  char R_vpd_err_d[10][32];
  for(Int_t r=1; r<10; r++)
  {
    sprintf(R_bbc_n[r],"R%d_bbc",r);
    sprintf(R_zdc_n[r],"R%d_zdc",r);
    sprintf(R_vpd_n[r],"R%d_vpd",r);
    sprintf(R_bbc_err_n[r],"R%d_bbc_err",r);
    sprintf(R_zdc_err_n[r],"R%d_zdc_err",r);
    sprintf(R_vpd_err_n[r],"R%d_vpd_err",r);

    sprintf(R_bbc_d[r],"R%d_bbc/F",r);
    sprintf(R_zdc_d[r],"R%d_zdc/F",r);
    sprintf(R_vpd_d[r],"R%d_vpd/F",r);
    sprintf(R_bbc_err_d[r],"R%d_bbc_err/F",r);
    sprintf(R_zdc_err_d[r],"R%d_zdc_err/F",r);
    sprintf(R_vpd_err_d[r],"R%d_vpd_err/F",r);

    patch_tr->Branch(R_bbc_n[r],&(R_bbc[r]),R_bbc_d[r]);
    patch_tr->Branch(R_bbc_err_n[r],&(R_bbc_err[r]),R_bbc_err_d[r]);
    patch_tr->Branch(R_zdc_n[r],&(R_zdc[r]),R_zdc_d[r]);
    patch_tr->Branch(R_zdc_err_n[r],&(R_zdc_err[r]),R_zdc_err_d[r]);
    patch_tr->Branch(R_vpd_n[r],&(R_vpd[r]),R_vpd_d[r]);
    patch_tr->Branch(R_vpd_err_n[r],&(R_vpd_err[r]),R_vpd_err_d[r]);
  };

  patch_tr->Branch("isConsistent",&isConsistent,"isConsistent/I");
  patch_tr->Branch("b_pol",&b_pol,"b_pol/F");
  patch_tr->Branch("y_pol",&y_pol,"y_pol/F");
  patch_tr->Branch("b_pol_err",&b_pol_err,"b_pol_err/F");
  patch_tr->Branch("y_pol_err",&y_pol_err,"y_pol_err/F");
  patch_tr->Branch("pattern",&pattern,"pattern/I");
  

  for(Int_t q=0; q<red_tr->GetEntries(); q++)
  {
    red_tr->GetEntry(q);

    if(q%1000 == 0) printf("%d/%d\n",q+1,red_tr->GetEntries());

    for(Int_t r=1; r<10; r++)
    {
      R_bbc[r] = RD->Rellum(runnum,r,"bbc");
      R_zdc[r] = RD->Rellum(runnum,r,"zdc");
      R_vpd[r] = RD->Rellum(runnum,r,"vpd");

      R_bbc_err[r] = RD->RellumErr(runnum,r,"bbc");
      R_zdc_err[r] = RD->RellumErr(runnum,r,"zdc");
      R_vpd_err[r] = RD->RellumErr(runnum,r,"vpd");
    };

    isConsistent = RD->RellumConsistent(runnum);
    b_pol = RD->BluePol(runnum);
    y_pol = RD->YellPol(runnum);
    b_pol_err = RD->BluePolErr(runnum);
    y_pol_err = RD->YellPolErr(runnum);
    kicked = RD->Kicked(runnum,Bunchid7bit);

    patch_tr->Fill();
  };

  patch_tr->Write("str");
}
