// test efficiency of "array map" algo

void test_array_map(Int_t runnum0=14121013)
{
  TFile * c_file = new TFile("counts.root","READ");
  TFile * r_file = new TFile("rtree.root","READ");

  TTree * ctr = (TTree*) c_file->Get("sca");
  TTree * rtr = (TTree*) r_file->Get("rellum");

  // build array map
  printf("building array map...\n");
  Int_t RUNS_tmp = ctr->GetMaximum("i");
  const Int_t RUNS = RUNS_tmp;
  const Int_t BXS = 120;
  const Int_t CONTENT = 3;

  Int_t array_map[RUNS][BXS][CONTENT]; // [run index-1] [bXing] [content number "j"]
  // "j" definition:
  //  j == 0 --> blue spin
  //  j == 1 --> yell spin
  //  j == 2 --> kicked typecasted to int)

  Int_t index_c,bx,blue,yell;
  Bool_t kicked;
  ctr->SetBranchAddress("i",&index_c);
  ctr->SetBranchAddress("bx",&bx);
  ctr->SetBranchAddress("blue",&blue);
  ctr->SetBranchAddress("yell",&yell);
  ctr->SetBranchAddress("kicked",&kicked);
  for(Int_t q=0; q<ctr->GetEntries(); q++)
  {
    ctr->GetEntry(q);
    array_map[index_c-1][bx][0] = blue;
    array_map[index_c-1][bx][1] = yell;
    if(kicked) array_map[index_c-1][bx][2] = 1;
    else array_map[index_c-1][bx][2] = 0;
  };


  // run number hashing (hash = run index - 1)
  printf("hashing run number %d\n",runnum0);
  Int_t index_r,runnum,hash;
  rtr->SetBranchAddress("i",&index_r);
  rtr->SetBranchAddress("runnum",&runnum);
  for(Int_t q=0; q<rtr->GetEntries(); q++)
  {
    rtr->GetEntry(q);
    if(runnum0 == runnum && q==index_r-1) hash=q;
  };

  
  // print out spin pattern from array map
  for(Int_t b=0; b<120; b++)
  {
    printf("%d %d %d %d\n",b,array_map[hash][b][0],array_map[hash][b][1],array_map[hash][b][2]);
  };
}
