// draws pT-dependent A_LL plots for three energy bins

void DrawEnergySeparation(Int_t asym=3, const char * filename="spin.root")
{
  TFile * infile = new TFile(filename,"READ");
  char dir[4][16];
  strcpy(dir[1],"A_L_yellow");
  strcpy(dir[2],"A_L_blue");
  strcpy(dir[3],"A_LL");
  char plot[3][32];
  for(Int_t i=0; i<3; i++) 
  {
    sprintf(plot[i],"/%s/pt_dep_a%d_e%d",dir[asym],asym,i);
    printf("reading %s\n",plot[i]);
  };
  TGraphErrors * a[3];
  TCanvas * c[3];
  char c_n[3][8];
  for(Int_t i=0; i<3; i++)
  {
    a[i] = (TGraphErrors*) infile->Get(plot[i]);
    sprintf(c_n[i],"c%d",i);
    c[i] = new TCanvas(c_n[i],c_n[i],1200,1000);
    c[i]->SetGrid(1,1);
    a[i]->SetLineWidth(4);
    a[i]->SetMarkerSize(4);
    a[i]->Draw("APE");
    a[i]->GetXaxis()->SetLimits(0,15);
  };



  if(asym==3)
  {
    a[0]->GetHistogram()->SetMinimum(-0.01);
    a[0]->GetHistogram()->SetMaximum(0.01);
    a[1]->GetHistogram()->SetMinimum(-0.01);
    a[1]->GetHistogram()->SetMaximum(0.06);
    a[2]->GetHistogram()->SetMinimum(-0.08);
    a[2]->GetHistogram()->SetMaximum(0.08);
  }
  else if(asym==2)
  {
    a[0]->GetHistogram()->SetMinimum(-0.004);
    a[0]->GetHistogram()->SetMaximum(0.004);
    a[1]->GetHistogram()->SetMinimum(-0.01);
    a[1]->GetHistogram()->SetMaximum(0.02);
    a[2]->GetHistogram()->SetMinimum(-0.1);
    a[2]->GetHistogram()->SetMaximum(0.06);
  }

  for(Int_t i=0; i<3; i++) c[i]->Update();

};
