// tests RunData class

void TestRunData()
{
  TFile * infile = new TFile("../rtree.root","READ");
  TTree * tr = (TTree*) infile->Get("rellum");

  gSystem->Load("../src/RunData.so");
  RunData * rd = new RunData();

  Int_t runnum;
  tr->SetBranchAddress("runnum",&runnum);

  printf("i run fill R3 R3_e b_pol b_pol_e y_pol y_pol_e pattern\n");
  for(Int_t q=0; q<tr->GetEntries(); q++)
  {
    tr->GetEntry(q);
    gSystem->RedirectOutput("out_rundata");
    printf("%d %d %d %f %f %d %d %f %f %f %f\n",
      q,
      runnum,
      rd->GetFill(runnum),
      rd->Rellum(runnum,3,"vpd"),
      rd->RellumErr(runnum,3,"vpd"),
      (Int_t) rd->RellumConsistent(runnum),
      rd->Pattern(runnum),
      rd->BluePol(runnum),
      rd->BluePolErr(runnum),
      rd->YellPol(runnum),
      rd->YellPolErr(runnum));
    /*
      printf("%d\n",runnum);
      printf("%d\n",rd->GetFill(runnum));
      printf("%f\n",rd->Rellum(runnum,3,"zdc"));
      printf("%f\n",rd->RellumErr(runnum,3,"zdc"));
      printf("%f\n",rd->BluePol(runnum));
      printf("%f\n",rd->BluePolErr(runnum));
      printf("%f\n",rd->YellPol(runnum));
      printf("%f\n",rd->YellPolErr(runnum));
    */
    gSystem->RedirectOutput(0);
  };
};

