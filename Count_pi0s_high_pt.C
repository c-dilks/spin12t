// counts number of pi0's for each run in redset files
// INCLUDES PT>=8 cutoff

void Count_pi0s_high_pt()
{
  TChain * tr = new TChain("str");
  tr->Add("./redset/Red*.root");
  printf("ent=%d\n",tr->GetEntries());
  Int_t runnum_tmp=0;
  Int_t runnum;
  Int_t count;
  tr->SetBranchAddress("runnum",&runnum);

  TGraph * gr = new TGraph();
  gr->SetTitle("no. #pi^{0}s vs. run");
  gr->SetMarkerStyle(kFullCircle);
  gr->SetMarkerColor(kBlue);
  Int_t pt=0;

  char pi0_cut[256];
  char cut[512];
  strcpy(pi0_cut,"(TrigBits&0x200) && N12==2 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && Pt>=8");

  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    if(runnum!=runnum_tmp)
    {
      //printf("runnum_tmp=%d runnum=%d\n",runnum_tmp,runnum);
      sprintf(cut,"%s && runnum==%d",pi0_cut,runnum);
      count = tr->GetEntries(cut);
      tr->GetEntry(i);
      gr->SetPoint(pt,runnum,count);
      pt++;
      printf("i=%d runnum=%d count=%d pt=%d\n",i,runnum,count,pt);
      runnum_tmp=runnum;
    };
  };
  TCanvas * cc = new TCanvas("cc","cc",500,700);
  gr->Draw("APE");

  TFile * outfile = new TFile("pi0cnt_high_pt.root","RECREATE");
  gr->Write("gr");
};
