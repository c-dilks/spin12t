// draws kin_dep plot from spin analysis output file
// and names it
//
// -- this is currently being used to see how changing the number of 
//    pt bins affects A_LL

void draw_kin_dep(const char * filename="spinset/spin_all.root",Int_t pt_bins0=10)
{
  TFile * infile = new TFile(filename,"READ");
  TGraphErrors * g = (TGraphErrors*) infile->Get("kin_dep");
  g->Fit("pol0","","",8,10);
  gStyle->SetOptFit(1);
  TCanvas * c = new TCanvas("c","c",1000,800);
  c->SetGrid(1,1);
  g->Draw("ape");
  char printname[128];
  sprintf(printname,"kin_dep_open/pt_%d.png",pt_bins0);
  c->Print(printname,"png");
};
