// shows the dependence of various kinematics on geometry
// e.g., pt vs. eta
// --> outputs all plots to diag.pdf

void Diagnostics()
{
  gSystem->Load("src/RunData12.so");
  RunData12 * RD = new RunData12();
  const Int_t NBINS=100; // NUMBER OF BINS

  // open chain
  TChain * tr = new TChain("str");
  tr->Add("./redset/Red*.root");
  Float_t E12,Pt,Eta,Phi,M12,Z,b_pol,y_pol;
  Bool_t kicked,isConsistent;
  Int_t TrigBits,runnum,bx;
  Float_t N12;

  Float_t E12_min,Pt_min,Eta_min,Phi_min;
  Float_t E12_max,Pt_max,Eta_max,Phi_max;
  E12_min=Pt_min=Eta_min=Phi_min=1000;
  E12_max=Pt_max=Eta_max=Phi_max=0;

  str->SetBranchAddress("runnum",&runnum);
  str->SetBranchAddress("Bunchid7bit",&bx);
  str->SetBranchAddress("E12",&E12);
  str->SetBranchAddress("Pt",&Pt);
  str->SetBranchAddress("Eta",&Eta);
  str->SetBranchAddress("Phi",&Phi);
  str->SetBranchAddress("M12",&M12);
  str->SetBranchAddress("Z",&Z);
  str->SetBranchAddress("TrigBits",&TrigBits);
  str->SetBranchAddress("N12",&N12);

  
  // get kinematics ranges for eta,Pt,Phi
  Int_t runnum_tmp=0;
  for(Int_t x=0; x<tr->GetEntries(); x++)
  {
    if((x%100000)==0) printf("loop 1: %.2f%%\n",100*((Float_t)x)/((Float_t)tr->GetEntries()));
    tr->GetEntry(x);
    kicked = RD->Kicked(runnum,bx);
    if(runnum!=runnum_tmp)
    {
      b_pol = RD->BluePol(runnum);
      y_pol = RD->YellPol(runnum);
      isConsistent = RD->RellumConsistent(runnum);
      runnum_tmp=runnum;
    }
    // n photon cut
    //if(N12==1 && kicked==0 && isConsistent==1 && b_pol*y_pol!=0)
    // pi0 cut
    if(fabs(M12-0.135)<0.1 && Z<0.8 && (TrigBits&0x200) && N12==2 && kicked==0 && isConsistent==1 && b_pol*y_pol!=0 && Pt<15)
    {
      E12_min = (E12 < E12_min) ? E12:E12_min;
      Pt_min  = (Pt  < Pt_min)  ? Pt:Pt_min;
      Eta_min = (Eta < Eta_min) ? Eta:Eta_min;
      Phi_min = (Phi < Phi_min) ? Phi:Phi_min;

      E12_max = (E12 > E12_max) ? E12:E12_max;
      Pt_max  = (Pt  > Pt_max)  ? Pt:Pt_max;
      Eta_max = (Eta > Eta_max) ? Eta:Eta_max;
      Phi_max = (Phi > Phi_max) ? Phi:Phi_max;
    };
  };


  TH2F * sph_pt_vs_eta = new TH2F("sph_pt_vs_eta","single #gamma :: p_{T} vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,Pt_max);
  TH2F * sph_en_vs_eta = new TH2F("sph_en_vs_eta","single #gamma :: E vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,E12_max);
  TH2F * sph_pt_vs_phi = new TH2F("sph_pt_vs_phi","single #gamma :: p_{T} vs. #phi",NBINS,-3.14,3.14,NBINS,0,Pt_max);
  TH2F * sph_en_vs_phi = new TH2F("sph_en_vs_phi","single #gamma :: E vs. #phi",NBINS,-3.14,3.14,NBINS,0,E12_max);
  TH2F * sph_eta_vs_phi = new TH2F("sph_eta_vs_phi","single #gamma :: #eta vs. #phi",NBINS,-3.14,3.14,NBINS,Eta_min,Eta_max);
  TH2F * sph_en_vs_pt = new TH2F("sph_en_vs_pt","single #gamma :: E vs. p_{T}",NBINS,0,Pt_max,NBINS,0,E12_max);

  TH2F * pi0_pt_vs_eta = new TH2F("pi0_pt_vs_eta","#pi^{0} :: p_{T} vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,Pt_max);
  TH2F * pi0_en_vs_eta = new TH2F("pi0_en_vs_eta","#pi^{0} :: E vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,E12_max);
  TH2F * pi0_pt_vs_phi = new TH2F("pi0_pt_vs_phi","#pi^{0} :: p_{T} vs. #phi",NBINS,-3.14,3.14,NBINS,0,Pt_max);
  TH2F * pi0_en_vs_phi = new TH2F("pi0_en_vs_phi","#pi^{0} :: E vs. #phi",NBINS,-3.14,3.14,NBINS,0,E12_max);
  TH2F * pi0_eta_vs_phi = new TH2F("pi0_eta_vs_phi","#pi^{0} :: #eta vs. #phi",NBINS,-3.14,3.14,NBINS,Eta_min,Eta_max);
  TH2F * pi0_en_vs_pt = new TH2F("pi0_en_vs_pt","#pi^{0} :: E vs. p_{T}",NBINS,0,Pt_max,NBINS,0,E12_max);

  TH2F * thr_pt_vs_eta = new TH2F("thr_pt_vs_eta","N_{#gamma}>2 :: p_{T} vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,Pt_max);
  TH2F * thr_en_vs_eta = new TH2F("thr_en_vs_eta","N_{#gamma}>2 :: E vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,E12_max);
  TH2F * thr_pt_vs_phi = new TH2F("thr_pt_vs_phi","N_{#gamma}>2 :: p_{T} vs. #phi",NBINS,-3.14,3.14,NBINS,0,Pt_max);
  TH2F * thr_en_vs_phi = new TH2F("thr_en_vs_phi","N_{#gamma}>2 :: E vs. #phi",NBINS,-3.14,3.14,NBINS,0,E12_max);
  TH2F * thr_eta_vs_phi = new TH2F("thr_eta_vs_phi","N_{#gamma}>2 :: #eta vs. #phi",NBINS,-3.14,3.14,NBINS,Eta_min,Eta_max);
  TH2F * thr_en_vs_pt = new TH2F("thr_en_vs_pt","N_{#gamma}>2 :: E vs. p_{T}",NBINS,0,Pt_max,NBINS,0,E12_max);

  TH2F * pi0_z_vs_eta = new TH2F("pi0_z_vs_eta","#pi^{0} :: Z vs. #eta",NBINS,Eta_min,Eta_max,NBINS,0,1);
  TH2F * pi0_z_vs_phi = new TH2F("pi0_z_vs_phi","#pi^{0} :: Z vs. #phi",NBINS,-3.14,3.14,NBINS,0,1);

  TH1F * mass_dist = new TH1F("mass_dist","M_{#gamma#gamma} distribution (N12==2, jet1, M12>0, Z<0.8, kicked==0)",NBINS,0,1);
  TH1F * z_dist = new TH1F("z_dist","Z distribution (N12==2, jet1, abs(M12-0.135)<0.1, kicked==0)",NBINS,0,1);
  TH1F * trig_dist = new TH1F("trig_dist","TrigBits distribution (N12==2)",NBINS,7500,33500);

  char cut[256];
  sprintf(cut,"abs(M12-0.135)<0.1 && Z<0.8 && (TrigBits&0x200) && kicked==0 && isConsistent==1 && b_pol*y_pol!=0");
  runnum_tmp=0;
  for(Int_t x=0; x<tr->GetEntries(); x++)
  {
    if((x%100000)==0) printf("loop 2: %.2f%%\n",100*((Float_t)x)/((Float_t)tr->GetEntries()));
    tr->GetEntry(x);
    kicked = RD->Kicked(runnum,bx);
    if(runnum!=runnum_tmp)
    {
      b_pol = RD->BluePol(runnum);
      y_pol = RD->YellPol(runnum);
      isConsistent = RD->RellumConsistent(runnum);
      runnum_tmp=runnum;
    }
    // rellum / pol cut
    if( kicked==0 && isConsistent==1 && b_pol>0 && y_pol>0)
    {
      // IF YOU CHANGE THE CUTS HERE, CHANGE THEM IN THE PLOT TITLES TOO!!!!!
      if(fabs(N12-2)<0.01 && (TrigBits&0x200) && M12>0 && Z<0.8) mass_dist->Fill(M12);
      if(fabs(N12-2)<0.01 && (TrigBits&0x200) && fabs(M12-0.135)<0.1) z_dist->Fill(Z);
      if(fabs(N12-2)<0.01) trig_dist->Fill(TrigBits);

      // single photon cut
      if(fabs(N12-1)<0.01)
      {
        sph_pt_vs_eta->Fill(Eta,Pt);
        sph_en_vs_eta->Fill(Eta,E12);
        sph_pt_vs_phi->Fill(Phi,Pt);
        sph_en_vs_phi->Fill(Phi,E12);
        sph_eta_vs_phi->Fill(Phi,Eta);
        sph_en_vs_pt->Fill(Pt,E12);
      };

      // pi0 cut
      if( (TrigBits&0x200) && fabs(N12-2)<0.01 && Z<0.8 && fabs(M12-0.135)<0.1)
      {
        pi0_pt_vs_eta->Fill(Eta,Pt);
        pi0_en_vs_eta->Fill(Eta,E12);
        pi0_pt_vs_phi->Fill(Phi,Pt);
        pi0_en_vs_phi->Fill(Phi,E12);
        pi0_eta_vs_phi->Fill(Phi,Eta);
        pi0_en_vs_pt->Fill(Pt,E12);

        pi0_z_vs_eta->Fill(Eta,Z);
        pi0_z_vs_phi->Fill(Phi,Z);
      };

      // three or more photons cut
      if(N12>2.5)
      {
        thr_pt_vs_eta->Fill(Eta,Pt);
        thr_en_vs_eta->Fill(Eta,E12);
        thr_pt_vs_phi->Fill(Phi,Pt);
        thr_en_vs_phi->Fill(Phi,E12);
        thr_eta_vs_phi->Fill(Phi,Eta);
        thr_en_vs_pt->Fill(Pt,E12);
      };
    };
  };


  /*
  // draw output -- MOVED TO DrawDiagnostics.C
  gStyle->SetOptStat(0);

  TCanvas * cc0 = new TCanvas("cc0","cc0",1000,1200);
  cc0->Divide(1,3);
  cc0->cd(1); mass_dist->Draw();
  cc0->cd(2); z_dist->Draw();
  cc0->cd(3); cc0->GetPad(3)->SetLogy(); trig_dist->Draw();


  TCanvas * cc1 = new TCanvas("cc1","cc1",1000,1200);
  cc1->Divide(2,3);
  cc1->cd(1); pt_vs_eta->Draw("colz");
  cc1->cd(2); pt_vs_phi->Draw("colz");
  cc1->cd(3); en_vs_eta->Draw("colz");
  cc1->cd(4); en_vs_phi->Draw("colz");
  cc1->cd(5); z_vs_eta->Draw("colz");
  cc1->cd(6); z_vs_phi->Draw("colz");

  TCanvas * cc2 = new TCanvas("cc2","cc2",1000,1200);
  cc2->Divide(2,3);
  cc2->cd(1); eta_vs_phi->Draw("colz");
  cc2->cd(2); en_vs_pt->Draw("colz");


  cc0->Print("diag_lin.pdf(","pdf");
  cc1->Print("diag_lin.pdf","pdf");
  cc2->Print("diag_lin.pdf)","pdf");
  printf("\ndiag_lin.pdf diagnostic kinematics file produced\n");

  for(Int_t i=1; i<=3; i++) cc0->GetPad(i)->SetLogy();
  for(Int_t i=1; i<=6; i++)
  {
    cc1->GetPad(i)->SetLogz();
    cc2->GetPad(i)->SetLogz();
  };

  cc0->Print("diag_log.pdf(","pdf");
  cc1->Print("diag_log.pdf","pdf");
  cc2->Print("diag_log.pdf)","pdf");
  printf("\ndiag_log.pdf diagnostic kinematics file produced\n");

  printf("E12 range: %f -- %f\n",E12_min,E12_max);
  printf("Pt range: %f -- %f\n",Pt_min,Pt_max);
  printf("Eta range: %f -- %f\n",Eta_min,Eta_max);
  printf("Phi range: %f -- %f\n",Phi_min,Phi_max);
  */

  // write output
  TFile * outfile = new TFile("diag.root","RECREATE");

  mass_dist->Write();
  z_dist->Write();
  trig_dist->Write();

  sph_pt_vs_eta->Write();
  sph_en_vs_eta->Write();
  sph_pt_vs_phi->Write();
  sph_en_vs_phi->Write();
  sph_eta_vs_phi->Write();
  sph_en_vs_pt->Write();

  pi0_pt_vs_eta->Write();
  pi0_en_vs_eta->Write();
  pi0_pt_vs_phi->Write();
  pi0_en_vs_phi->Write();
  pi0_eta_vs_phi->Write();
  pi0_en_vs_pt->Write();

  thr_pt_vs_eta->Write();
  thr_en_vs_eta->Write();
  thr_pt_vs_phi->Write();
  thr_en_vs_phi->Write();
  thr_eta_vs_phi->Write();
  thr_en_vs_pt->Write();

  pi0_z_vs_eta->Write();
  pi0_z_vs_phi->Write();


}
  
