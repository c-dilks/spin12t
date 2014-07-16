void compare_isolation()
{
  TFile * infile1 = new TFile("spinset/spin_all_100mr.root","READ");
  TFile * infile2 = new TFile("spinset/spin_all_35mr.root","READ");

  TGraphErrors * k1 = (TGraphErrors*) infile1->Get("kin_dep");
  TGraphErrors * k2 = (TGraphErrors*) infile2->Get("kin_dep");

  k1->SetTitle("#epsilon_{LL} vs. p_{T} -- 100mr");
  k2->SetTitle("#epsilon_{LL} vs. p_{T} -- 35mr");

  TCanvas * cc_1 = new TCanvas("cc_1","cc_1",1000,800);
  cc_1->SetGrid(1,1);
  k1->Draw("ape");
  cc_1->Print("cc_1.png","png");

  TCanvas * cc_2 = new TCanvas("cc_2","cc_2",1000,800);
  cc_2->SetGrid(1,1);
  k2->Draw("ape");
  cc_2->Print("cc_2.png","png");

  TCanvas * cc_both = new TCanvas("cc_both","cc_both",1000,800);
  cc_both->SetGrid(1,1);
  TMultiGraph * mg = new TMultiGraph();
  mg->SetTitle("#epsilon_{LL} vs. p_{T} -- red:100mr,  blue:35mr");
  k2->SetMarkerColor(kBlue);
  mg->Add(k1);
  mg->Add(k2);
  mg->Draw("ape");
  cc_both->Print("cc_both.png","png");

}
