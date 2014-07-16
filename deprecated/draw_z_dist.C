void draw_mass_dist()
{
  TChain * TwoTr = new TChain("TwoTr");
  TwoTr->Add("/home/dilks/h5/Output/Outputset*.root");

  Float_t min_bin = 0; 
  Float_t max_bin = 1;
  char md_cut[256];
  sprintf(md_cut,"N12==2 && (TrigBits&0xFFF)>0 && M12>%f && M12<%f",min_bin,max_bin);
  TH1D * md = new TH1D("md","M_{#gamma#gamma} distribution",100,min_bin,max_bin);
  printf("mass distribution cut = %s\n",md_cut);
  TwoTr->Project("md","M12",md_cut);

  Float_t xl = 0.135-0.1;
  Float_t xh = 0.135+0.1;
  Float_t max = md->GetMaximum();
  TLine * ll = new TLine(xl,0,xl,max);
  TLine * lh = new TLine(xh,0,xh,max);
  ll->SetLineColor(kRed);
  lh->SetLineColor(kRed);
  ll->SetLineWidth(2);
  lh->SetLineWidth(2);

  TCanvas * cc = new TCanvas("cc","cc",700,500);
  md->Draw();
  ll->Draw();
  lh->Draw();
};
