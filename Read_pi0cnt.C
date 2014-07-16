// reads pi0cnt.root and prints the data points to text file

void Read_pi0cnt(Int_t div=2, const char * filename="pi0cnt.root")
{
  TFile * infile = new TFile(filename,"READ");
  TGraph * gr = (TGraph*) infile->Get("gr");
  Double_t run,cnt;
  Int_t total;
  total=0;
  char tp[8];

  // count total # pi0s
  for(Int_t i=0; i<gr->GetN(); i++)
  {
    if(i>0) strcpy(tp,"a");
    else strcpy(tp,"w");
    gr->GetPoint(i,run,cnt);
    //gSystem->RedirectOutput("pi0_cnt.dat",tp);
    printf("%d %f\n",(Int_t)run,cnt);
    //gSystem->RedirectOutput(0);
    total+=((Int_t)cnt);
  };
  printf("\ntotal # pi0's: %d\n",(Int_t)total);

  // determine divisions
  Double_t cnt_div = total/((Double_t)div);
  Int_t ccc=0;
  Int_t ii=0;
  Double_t run_tmp;
  for(Int_t i=0; i<gr->GetN(); i++)
  {
    gr->GetPoint(i,run,cnt);
    if(ccc==0) run_tmp = run;
    ccc+=((Int_t)cnt);
    if(ccc>=cnt_div || i+1==gr->GetN())
    {
      if(ii>0) strcpy(tp,"a");
      else strcpy(tp,"w");
      printf("run_low=%d run_high=%d pi0s=%d\n",(Int_t)run_tmp,(Int_t)run,ccc);
      gSystem->RedirectOutput("div.dat",tp);
      printf("%d %d %d\n",(Int_t)run_tmp,(Int_t)run,ccc);
      gSystem->RedirectOutput(0);
      ii++;
      ccc=0;
    };
  };
  printf("div.dat created\n");
};
