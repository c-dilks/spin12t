// counts total number of runs in RedSet files

void count_runs()
{
  TChain * tr = new TChain("str");
  tr->Add("./redset/Red*.root");
  printf("ent=%d\n",tr->GetEntries());
  Int_t r;
  tr->SetBranchAddress("runnum",&r);
  Int_t r_tmp=0;
  Int_t cnt=0;
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    if(i%1000000==0) printf("%.1f%%\n",100*(Float_t)i/tr->GetEntries());
    tr->GetEntry(i);
    if(r!=r_tmp) 
    {
      cnt++;
      r_tmp=r;
    };
  };
  printf("cnt=%d\n",cnt);
};
