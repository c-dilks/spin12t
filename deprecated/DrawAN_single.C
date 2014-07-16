// computes AN for a single Pt:E bin
// defined by CUTS below

void DrawAN_single(Bool_t useBlue=1, const char * filename="RedOutputset107ha.root")
{
  // open Output file
  char root_file[256];
  sprintf(root_file,"redset/%s",filename);
  TFile * infile = new TFile(root_file,"READ");

  //TTree * tr = (TTree*) infile->Get("str");
  TChain * tr = new TChain("str"); tr->Add("redset/*.root");

  // declare cos(phi) distributions
  const Int_t cos_phi_bins=32;
  TH1F * cos_phi_up = new TH1F("cos_phi_up","cos_phi_up",cos_phi_bins,-1,1);
  TH1F * cos_phi_dn = new TH1F("cos_phi_dn","cos_phi_dn",cos_phi_bins,-1,1);


  // cuts
  char spin_cut_up[32];
  char spin_cut_dn[32];
  char pt_cut[32];
  char e12_cut[32];
  char trig_cut[32];
  char m12_cut[32];
  char z_cut[32];
  char eta_cut[32];
  //char n12_cut[32]; // done in ReduceData.C

  if(useBlue)
  {
    // blue cuts
    strcpy(spin_cut_up,"(spin==2 || spin==3)");
    strcpy(spin_cut_dn,"(spin==0 || spin==1)");
  }
  else
  {
    // yellow cuts
    strcpy(spin_cut_up,"(spin==0 || spin==2)");
    strcpy(spin_cut_dn,"(spin==1 || spin==3)");
  };
  strcpy(trig_cut,"(TrigBits&0x200)");
  //strcpy(n12_cut,"N12==2");
  
  /* CUTS */
  // -- from Run11An_1.ppt
  sprintf(pt_cut,"Pt>0 && Pt<10");
  sprintf(e12_cut,"E12>45 && E12<95");
  sprintf(m12_cut,"abs(M12-%f)<%f",0.135,0.1);
  sprintf(z_cut,"Z<%f",0.8);
  sprintf(eta_cut,"abs(Eta-%f)<%f",3.7,0.2);

  char cut_up[512];
  char cut_dn[512];
  sprintf(cut_up,"%s && %s && %s && %s && %s && %s && %s",
          spin_cut_up,pt_cut,e12_cut,trig_cut,z_cut,m12_cut,eta_cut);
  sprintf(cut_dn,"%s && %s && %s && %s && %s && %s && %s",
          spin_cut_dn,pt_cut,e12_cut,trig_cut,z_cut,m12_cut,eta_cut);

  printf("cut_up=%s\n",cut_up);
  printf("cut_dn=%s\n",cut_dn);

  
  // tree projections
  printf("projecting cos_phi_up\n");
  tr->Project("cos_phi_up","cos(Phi)",cut_up);
  printf("projecting cos_phi_dn\n");
  tr->Project("cos_phi_dn","cos(Phi)",cut_dn);

  
  // analysing power
  //TH1F * an_cos_phi = new TH1F((*cos_phi_up-*cos_phi_dn)/(*cos_phi_up+*cos_phi_dn));
  TH1F * an_cos_phi = new TH1F("an_cos_phi","an_cos_phi",cos_phi_bins,-1,1);
  //cos_phi_up->Sumw2();  // not needed since computed in GetAsymmetry
  //cos_phi_dn->Sumw2();
  an_cos_phi = (TH1F*)cos_phi_up->GetAsymmetry(cos_phi_dn);
  if(useBlue) an_cos_phi->SetTitle("Blue A_{N}");
  else an_cos_phi->SetTitle("Yellow A_{N}");

  /*
  Float_t bcc_up,bcc_dn;
  for(Int_t i=1; i<cos_phi_up->GetNbinsX(); i++)
  {
    bcc_up = cos_phi_up->GetBinContent(i);
    bcc_dn = cos_phi_dn->GetBinContent(i);
    if(bcc_up>0 && bcc_dn>0) 
      an_cos_phi->SetBinContent(i,(bcc_up-bcc_dn)/(bcc_up+bcc_dn));
  };
  */

  gStyle->SetOptFit(1);
  an_cos_phi->Fit("pol1","b","",-1,1);
  c1->Close();

  
  // draw output
  TCanvas * cc = new TCanvas("cc","cc",700,500);
  an_cos_phi->Draw("e");
  
  // write output
  TFile * outfile = new TFile("single_AN.root","RECREATE");
  an_cos_phi->Write("an_cos_phi");
  cos_phi_up->Write("cos_phi_up");
  cos_phi_dn->Write("cos_phi_dn");

}
