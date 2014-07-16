// computes ALL for a single Pt:E bin
// defined by CUTS below

void DrawALL_single(const char * filename="RedOutputset107ha.root")
{
  // open Output file
  char root_file[256];
  sprintf(root_file,"redset/%s",filename);
  TFile * infile = new TFile(root_file,"READ");

  //TTree * tr = (TTree*) infile->Get("str");
  TChain * tr = new TChain("str"); tr->Add("redset/*.root");

  // declare cos(phi) distributions
  const Int_t cos_phi_bins=32;
  TH1F * cos_phi_Bup_Yup = new TH1F("cos_phi_Bup_Yup","cos_phi_Bup_Yup",cos_phi_bins,-1,1);
  TH1F * cos_phi_Bup_Ydn = new TH1F("cos_phi_Bup_Ydn","cos_phi_Bup_Ydn",cos_phi_bins,-1,1);
  TH1F * cos_phi_Bdn_Yup = new TH1F("cos_phi_Bdn_Yup","cos_phi_Bdn_Yup",cos_phi_bins,-1,1);
  TH1F * cos_phi_Bdn_Ydn = new TH1F("cos_phi_Bdn_Ydn","cos_phi_Bdn_Ydn",cos_phi_bins,-1,1);


  // cuts
  char spin_cut_Bup_Yup[32];
  char spin_cut_Bup_Ydn[32];
  char spin_cut_Bdn_Yup[32];
  char spin_cut_Bdn_Ydn[32];
  char pt_cut[32];
  char e12_cut[32];
  char trig_cut[32];
  char m12_cut[32];
  char z_cut[32];
  char eta_cut[32];
  //char n12_cut[32]; // done in ReduceData.C

  strcpy(spin_cut_Bup_Yup,"spin==3");
  strcpy(spin_cut_Bup_Ydn,"spin==2");
  strcpy(spin_cut_Bdn_Yup,"spin==1");
  strcpy(spin_cut_Bdn_Ydn,"spin==0");
  strcpy(trig_cut,"(TrigBits&0x200)");
  //strcpy(n12_cut,"N12==2");
  
  // CUTS
  // -- original
  /*
  sprintf(pt_cut,"Pt>1 && Pt<10");
  sprintf(e12_cut,"abs(E12-%f)<%f",70.0,25.0);
  sprintf(m12_cut,"abs(M12-%f)<%f",0.135,0.1);
  sprintf(z_cut,"Z<%f",0.8);
  sprintf(eta_cut,"abs(Eta-%f)<%f",3.7,0.2); // boundary btwn large/small is 3.3
  */
  // -- from Run11An_1.ppt
  sprintf(pt_cut,"Pt>0 && Pt<10");
  sprintf(e12_cut,"E12>45 && E12<95");
  sprintf(m12_cut,"abs(M12-%f)<%f",0.135,0.1);
  sprintf(z_cut,"Z<%f",0.8);
  sprintf(eta_cut,"Eta>3.5 && Eta<3.9");

  char cut_Bup_Yup[512];
  char cut_Bup_Ydn[512];
  char cut_Bdn_Yup[512];
  char cut_Bdn_Ydn[512];
  sprintf(cut_Bup_Yup,"%s && %s && %s && %s && %s && %s && %s",
     spin_cut_Bup_Yup,pt_cut,e12_cut,trig_cut,z_cut,m12_cut,eta_cut);
  sprintf(cut_Bup_Ydn,"%s && %s && %s && %s && %s && %s && %s",
     spin_cut_Bup_Ydn,pt_cut,e12_cut,trig_cut,z_cut,m12_cut,eta_cut);
  sprintf(cut_Bdn_Yup,"%s && %s && %s && %s && %s && %s && %s",
     spin_cut_Bdn_Yup,pt_cut,e12_cut,trig_cut,z_cut,m12_cut,eta_cut);
  sprintf(cut_Bdn_Ydn,"%s && %s && %s && %s && %s && %s && %s",
     spin_cut_Bdn_Ydn,pt_cut,e12_cut,trig_cut,z_cut,m12_cut,eta_cut);

  printf("cut_Bup_Yup=%s\n",cut_Bup_Yup);
  printf("cut_Bup_Ydn=%s\n",cut_Bup_Ydn);
  printf("cut_Bdn_Yup=%s\n",cut_Bdn_Yup);
  printf("cut_Bdn_Ydn=%s\n",cut_Bdn_Ydn);

  
  // tree projections
  printf("projecting cos_phi_Bup_Yup\n");
  tr->Project("cos_phi_Bup_Yup","cos(Phi)",cut_Bup_Yup);
  printf(" %d entries projected\n",cos_phi_Bup_Yup->GetEntries());
  printf("projecting cos_phi_Bup_Ydn\n");
  tr->Project("cos_phi_Bup_Ydn","cos(Phi)",cut_Bup_Ydn);
  printf(" %d entries projected\n",cos_phi_Bup_Ydn->GetEntries());
  printf("projecting cos_phi_Bdn_Yup\n");
  tr->Project("cos_phi_Bdn_Yup","cos(Phi)",cut_Bdn_Yup);
  printf(" %d entries projected\n",cos_phi_Bdn_Yup->GetEntries());
  printf("projecting cos_phi_Bdn_Ydn\n");
  tr->Project("cos_phi_Bdn_Ydn","cos(Phi)",cut_Bdn_Ydn);
  printf(" %d entries projected\n",cos_phi_Bdn_Ydn->GetEntries());

  // store sum of square of weights
  cos_phi_Bup_Yup->Sumw2(); // need to call Sumw2 before adding histograms
  cos_phi_Bup_Ydn->Sumw2();
  cos_phi_Bdn_Yup->Sumw2();
  cos_phi_Bdn_Ydn->Sumw2();


  // analysing power
  TH1F * all_cos_phi = new TH1F("all_cos_phi","all_cos_phi",cos_phi_bins,-1,1);
  TH1F * cos_phi_same = new TH1F("cos_phi_same","cos_phi_same",cos_phi_bins,-1,1);
  TH1F * cos_phi_diff = new TH1F("cos_phi_diff","cos_phi_diff",cos_phi_bins,-1,1);
  cos_phi_same->Add(cos_phi_Bup_Yup,cos_phi_Bdn_Ydn);
  cos_phi_diff->Add(cos_phi_Bup_Ydn,cos_phi_Bdn_Yup);
  //cos_phi_same->Sumw2(); // not necessary here since done in GetAsymmetry
  //cos_phi_diff->Sumw2(); 
  all_cos_phi = (TH1F*)cos_phi_same->GetAsymmetry(cos_phi_diff);
  all_cos_phi->SetTitle("A_{LL}");

  /*
  Float_t bcc_Bup_Yup,bcc_Bup_Ydn,bcc_Bdn_Yup,bcc_Bdn_Ydn;
  for(Int_t i=1; i<cos_phi_Bup_Yup->GetNbinsX(); i++)
  {
    bcc_Bup_Yup = cos_phi_Bup_Yup->GetBinContent(i);
    bcc_Bup_Ydn = cos_phi_Bup_Ydn->GetBinContent(i);
    bcc_Bdn_Yup = cos_phi_Bdn_Yup->GetBinContent(i);
    bcc_Bdn_Ydn = cos_phi_Bdn_Ydn->GetBinContent(i);
    if(bcc_Bup_Yup*bcc_Bup_Ydn*bcc_Bdn_Yup*bcc_Bdn_Ydn>0) 
      all_cos_phi->SetBinContent(i,
        ( bcc_Bup_Yup + bcc_Bdn_Ydn - bcc_Bup_Ydn - bcc_Bdn_Yup )/
        ( bcc_Bup_Yup + bcc_Bdn_Ydn + bcc_Bup_Ydn + bcc_Bdn_Yup ));
  };
  */

  gStyle->SetOptFit(1);
  all_cos_phi->Fit("pol1","b","",-1,1);
  c1->Close();

  
  // draw output
  TCanvas * cc = new TCanvas("cc","cc",700,500);
  all_cos_phi->Draw("e");
  
  // write output
  TFile * outfile = new TFile("single_ALL.root","RECREATE");
  all_cos_phi->Write("all_cos_phi");
  cos_phi_Bup_Yup->Write("cos_phi_Bup_Yup");
  cos_phi_Bup_Ydn->Write("cos_phi_Bup_Ydn");
  cos_phi_Bdn_Yup->Write("cos_phi_Bdn_Yup");
  cos_phi_Bdn_Ydn->Write("cos_phi_Bdn_Ydn");

}
