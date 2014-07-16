// draws plots in diag.root, which is produced by Diagnostics.C

void DrawDiagnostics(const char * filename="diag.root")
{
  TFile * diag_file = new TFile("diag.root","READ");
  TH1F * mass_dist = (TH1F*) diag_file->Get("mass_dist");
  TH1F * z_dist = (TH1F*) diag_file->Get("z_dist");
  TH1F * trig_dist = (TH1F*) diag_file->Get("trig_dist");

  TH2F * sph_pt_vs_eta = (TH2F*) diag_file->Get("sph_pt_vs_eta");
  TH2F * sph_en_vs_eta = (TH2F*) diag_file->Get("sph_en_vs_eta");
  TH2F * sph_pt_vs_phi = (TH2F*) diag_file->Get("sph_pt_vs_phi");
  TH2F * sph_en_vs_phi = (TH2F*) diag_file->Get("sph_en_vs_phi");
  TH2F * sph_eta_vs_phi = (TH2F*) diag_file->Get("sph_eta_vs_phi");
  TH2F * sph_en_vs_pt = (TH2F*) diag_file->Get("sph_en_vs_pt");

  TH2F * pi0_pt_vs_eta = (TH2F*) diag_file->Get("pi0_pt_vs_eta");
  TH2F * pi0_en_vs_eta = (TH2F*) diag_file->Get("pi0_en_vs_eta");
  TH2F * pi0_pt_vs_phi = (TH2F*) diag_file->Get("pi0_pt_vs_phi");
  TH2F * pi0_en_vs_phi = (TH2F*) diag_file->Get("pi0_en_vs_phi");
  TH2F * pi0_eta_vs_phi = (TH2F*) diag_file->Get("pi0_eta_vs_phi");
  TH2F * pi0_en_vs_pt = (TH2F*) diag_file->Get("pi0_en_vs_pt");

  TH2F * thr_pt_vs_eta = (TH2F*) diag_file->Get("thr_pt_vs_eta");
  TH2F * thr_en_vs_eta = (TH2F*) diag_file->Get("thr_en_vs_eta");
  TH2F * thr_pt_vs_phi = (TH2F*) diag_file->Get("thr_pt_vs_phi");
  TH2F * thr_en_vs_phi = (TH2F*) diag_file->Get("thr_en_vs_phi");
  TH2F * thr_eta_vs_phi = (TH2F*) diag_file->Get("thr_eta_vs_phi");
  TH2F * thr_en_vs_pt = (TH2F*) diag_file->Get("thr_en_vs_pt");

  TH2F * pi0_z_vs_eta = (TH2F*) diag_file->Get("pi0_z_vs_eta");
  TH2F * pi0_z_vs_phi = (TH2F*) diag_file->Get("pi0_z_vs_phi");


  // draw output
  gStyle->SetOptStat(0);

  TCanvas * cc0 = new TCanvas("cc0","cc0",1000,1200);
  cc0->Divide(1,3);
  cc0->cd(1); mass_dist->Draw();
  cc0->cd(2); z_dist->Draw();
  cc0->cd(3); cc0->GetPad(3)->SetLogy(); trig_dist->Draw();

  TCanvas * sph_cc = new TCanvas("sph_cc","sph_cc",1000,1200);
  sph_cc->Divide(2,3);
  for(Int_t i=1; i<=6; i++) sph_cc->GetPad(i)->SetLogz();
  sph_cc->cd(1); sph_pt_vs_eta->Draw("colz");
  sph_cc->cd(2); sph_pt_vs_phi->Draw("colz");
  sph_cc->cd(3); sph_en_vs_eta->Draw("colz");
  sph_cc->cd(4); sph_en_vs_phi->Draw("colz");
  sph_cc->cd(5); sph_eta_vs_phi->Draw("colz");
  sph_cc->cd(6); sph_en_vs_pt->Draw("colz");

  TCanvas * pi0_cc = new TCanvas("pi0_cc","pi0_cc",1000,1200);
  pi0_cc->Divide(2,3);
  for(Int_t i=1; i<=6; i++) pi0_cc->GetPad(i)->SetLogz();
  pi0_cc->cd(1); pi0_pt_vs_eta->Draw("colz");
  pi0_cc->cd(2); pi0_pt_vs_phi->Draw("colz");
  pi0_cc->cd(3); pi0_en_vs_eta->Draw("colz");
  pi0_cc->cd(4); pi0_en_vs_phi->Draw("colz");
  pi0_cc->cd(5); pi0_eta_vs_phi->Draw("colz");
  pi0_cc->cd(6); pi0_en_vs_pt->Draw("colz");

  TCanvas * pi0_ccz = new TCanvas("pi0_ccz","pi0_ccz",1000,1200);
  pi0_ccz->Divide(1,2);
  for(Int_t i=1; i<=2; i++) pi0_ccz->GetPad(i)->SetLogz();
  pi0_ccz->cd(1); pi0_z_vs_eta->Draw("colz");
  pi0_ccz->cd(2); pi0_z_vs_phi->Draw("colz");

  TCanvas * thr_cc = new TCanvas("thr_cc","thr_cc",1000,1200);
  thr_cc->Divide(2,3);
  for(Int_t i=1; i<=6; i++) thr_cc->GetPad(i)->SetLogz();
  thr_cc->cd(1); thr_pt_vs_eta->Draw("colz");
  thr_cc->cd(2); thr_pt_vs_phi->Draw("colz");
  thr_cc->cd(3); thr_en_vs_eta->Draw("colz");
  thr_cc->cd(4); thr_en_vs_phi->Draw("colz");
  thr_cc->cd(5); thr_eta_vs_phi->Draw("colz");
  thr_cc->cd(6); thr_en_vs_pt->Draw("colz");


  cc0->Print("diag.pdf(","pdf");
  sph_cc->Print("diag.pdf","pdf");
  pi0_cc->Print("diag.pdf","pdf");
  pi0_ccz->Print("diag.pdf","pdf");
  thr_cc->Print("diag.pdf)","pdf");
  printf("\ndiag.pdf diagnostic kinematics file produced\n");
};
