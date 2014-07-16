// tests spin pattern / kicked status from RunData class

void Spin_Test(Int_t runnum=14096078)
{
  gSystem->Load("../src/RunData.so");
  RunData * RD = new RunData();

  printf("\nrun=%d\nfill=%d\n",runnum,RD->GetFill(runnum));
  for(Int_t bx=0; bx<120; bx++)
    printf("%d\t%d\t%d\n",
      RD->BlueSpin(runnum,bx),
      RD->YellSpin(runnum,bx),
      (Int_t)(RD->Kicked(runnum,bx)));
};
