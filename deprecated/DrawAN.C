// computes AN for a single Pt:E bin

void DrawAN(Float_t pt=7.0,
            Float_t e12=50.0,
            const char * filename="redset/RedOutputset107ha.root")
{
  // open Output file
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("str");


  // declare cos(phi) distributions
  Int_t cos_phi_bins=10;
  TH1F * cos_phi_up = new TH1F("cos_phi_up","cos_phi_up",cos_phi_bins,-1,1);
  TH1F * cos_phi_dn = new TH1F("cos_phi_dn","cos_phi_dn",cos_phi_bins,-1,1);


  // cuts
  char spin_cut_up[32];
  char spin_cut_dn[32];
  char pt_cut[32];
  char e12_cut[32];

  strcpy(spin_cut_up,"(spin==2 || spin==3)");
  strcpy(spin_cut_dn,"(spin==0 || spin==1)");
  
  sprintf(pt_cut,"abs(Pt-%f)<%f",7.0,1.0);
  sprintf(e12_cut,"abs(E12-%f)<%f",50.0,10.0);

  char cut_up[512];
  char cut_dn[512];
  sprintf(cut_up,"%s && %s && %s",spin_cut_up, pt_cut, e12_cut);
  sprintf(cut_dn,"%s && %s && %s",spin_cut_dn, pt_cut, e12_cut);

  
  // tree projections
  printf("projecting cos_phi_up\n");
  tr->Project("cos_phi_up","cos(Phi)",cut_up);
  printf("projecting cos_phi_dn\n");
  tr->Project("cos_phi_dn","cos(Phi)",cut_dn);

  
  // analysing power
  TH1F * an_cos_phi = new TH1F((*cos_phi_up-*cos_phi_dn)/(*cos_phi_up+*cos_phi_dn));
  an_cos_phi->Sumw2();
  an_cos_phi->Fit("pol0","b","",-1,1);

  
  // draw output
  /*
  TCanvas * cc = new TCanvas("cc","cc",700,500);
  an_cos_phi->Draw("e");
  */
}
