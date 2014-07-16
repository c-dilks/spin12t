// counts number of pi0's for each run in redset files
// for different spinbits

void Count_pi0s_bx()
{
  TChain * tr = new TChain("str");
  tr->Add("./redset/Red*.root");
  printf("ent=%d\n",tr->GetEntries());
  Int_t runnum_tmp=0;
  Int_t runnum;
  Int_t count[4];
  tr->SetBranchAddress("runnum",&runnum);


  char cut[4][512];
  /*
  strcpy(cut[0],"Bunchid7bit<60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==-1 && yell==-1");
  strcpy(cut[1],"Bunchid7bit<60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==-1 && yell==1");
  strcpy(cut[2],"Bunchid7bit<60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==1 && yell==-1");
  strcpy(cut[3],"Bunchid7bit<60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==1 && yell==1");
  */
  ///*
  strcpy(cut[0],"Bunchid7bit>=60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==-1 && yell==-1");
  strcpy(cut[1],"Bunchid7bit>=60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==-1 && yell==1");
  strcpy(cut[2],"Bunchid7bit>=60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==1 && yell==-1");
  strcpy(cut[3],"Bunchid7bit>=60 && (TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==1 && yell==1");
  //*/
  /*
  strcpy(cut[0],"(TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==-1 && yell==-1");
  strcpy(cut[1],"(TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==-1 && yell==1");
  strcpy(cut[2],"(TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==1 && yell==-1");
  strcpy(cut[3],"(TrigBits&0x200) && Pt<15 && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0 && blue==1 && yell==1");
  */

  for(Int_t j=0; j<4; j++)
  {
    count[j] = tr->GetEntries(cut[j]);
    printf("spinbit=%d count=%d\n",j,count[j]);
  };
};
