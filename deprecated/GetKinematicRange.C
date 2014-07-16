// gets kinematic ranges, used for theory calculations
// uses hadded redset file, redset/all.root

void GetKinematicRange()
{
  TChain * tr = new TChain("str");
  tr->Add("./redset/Red*.root");
  Float_t E12,Pt,Eta,Phi,M12,Z,b_pol,y_pol;
  Bool_t kicked,isConsistent;
  Int_t TrigBits;

  Float_t E12_min,Pt_min,Eta_min;
  Float_t E12_max,Pt_max,Eta_max;
  E12_min=Pt_min=Eta_min=1000;
  E12_max=Pt_max=Eta_max=0;

  str->SetBranchAddress("E12",&E12);
  str->SetBranchAddress("Pt",&Pt);
  str->SetBranchAddress("Eta",&Eta);
  str->SetBranchAddress("Phi",&Phi);
  str->SetBranchAddress("M12",&M12);
  str->SetBranchAddress("Z",&Z);
  str->SetBranchAddress("b_pol",&b_pol);
  str->SetBranchAddress("y_pol",&y_pol);
  str->SetBranchAddress("kicked",&kicked);
  str->SetBranchAddress("isConsistent",&isConsistent);
  str->SetBranchAddress("TrigBits",&TrigBits);

  
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    if(fabs(M12-0.135)<0.1 && Z<0.8 && (TrigBits&0x200) && kicked==0 && isConsistent==1 && b_pol*y_pol!=0)
    {
      E12_min = (E12 < E12_min) ? E12:E12_min;
      Pt_min  = (Pt  < Pt_min)  ? Pt:Pt_min;
      Eta_min = (Eta < Eta_min) ? Eta:Eta_min;

      E12_max = (E12 > E12_max) ? E12:E12_max;
      Pt_max  = (Pt  > Pt_max)  ? Pt:Pt_max;
      Eta_max = (Eta > Eta_max) ? Eta:Eta_max;
    };
  };

  printf("E12 range: %f -- %f\n",E12_min,E12_max);
  printf("Pt range: %f -- %f\n",Pt_min,Pt_max);
  printf("Eta range: %f -- %f\n",Eta_min,Eta_max);
};


