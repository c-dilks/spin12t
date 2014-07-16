// calculates AN for various E:Pt bins according to Eta vs. Cos(phi) distributions
// - the AN calculation is not implemented yet

void DrawAN3(const char * filename="RedOutputset107ha.root", Int_t nbins_E=15, Int_t nbins_Pt=15)
{
  const Int_t COS_BINS=10;

  // read str tree from reduced file
  char root_file[256];
  sprintf(root_file,"redset/%s",filename);
  TFile * infile = new TFile(root_file,"READ");
  
  TTree * tr = (TTree*) infile->Get("str");
  //TChain * tr = new TChain("str");
  //tr->Add("redset/*.root");
  
  Float_t M12,E12,Phi,Eta,Pt;
  Int_t spin;
  tr->SetBranchAddress("M12",&M12);
  tr->SetBranchAddress("E12",&E12);
  tr->SetBranchAddress("Phi",&Phi);
  tr->SetBranchAddress("Eta",&Eta);
  tr->SetBranchAddress("Pt",&Pt);
  tr->SetBranchAddress("spin",&spin);


  // define energy and pt bins
  const Int_t Ebins = nbins_E;
  const Int_t Ptbins = nbins_Pt;
  Float_t E_arr[Ebins];
  Float_t Pt_arr[Ptbins];
  Float_t E_div = 100./((Float_t)Ebins);
  Float_t Pt_div = 10./((Float_t)Ptbins);
  printf("E bins:\n");
  for(Int_t n=0; n<Ebins; n++) { E_arr[n] = n*E_div; printf(" %.1f",E_arr[n]); };
  printf("\nPt bins:\n");
  for(Int_t n=0; n<Ptbins; n++) { Pt_arr[n] = n*Pt_div; printf(" %.1f",Pt_arr[n]); };
  printf("\n");



  // initialise cos(phi) distributions
  printf("initialising histograms\n");
  TH2F * cos_phi_up[Ebins][Ptbins];
  TH2F * cos_phi_dn[Ebins][Ptbins];
  char cos_phi_up_n[Ebins][Ptbins][32];
  char cos_phi_dn_n[Ebins][Ptbins][32];
  char cos_phi_up_t[Ebins][Ptbins][128];
  char cos_phi_dn_t[Ebins][Ptbins][128];
  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      sprintf(cos_phi_up_n[n_E][n_Pt],"c_up_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_dn_n[n_E][n_Pt],"c_dn_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_up_t[n_E][n_Pt],
              "B#uparrow  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_dn_t[n_E][n_Pt],
              "B#downarrow  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      cos_phi_up[n_E][n_Pt] = new TH2F(cos_phi_up_n[n_E][n_Pt],cos_phi_up_t[n_E][n_Pt],COS_BINS,-1,1,COS_BINS,2.5,4.5);
      cos_phi_dn[n_E][n_Pt] = new TH2F(cos_phi_dn_n[n_E][n_Pt],cos_phi_dn_t[n_E][n_Pt],COS_BINS,-1,1,COS_BINS,2.5,4.5);
      cos_phi_up[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_dn[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_up[n_E][n_Pt]->GetYaxis()->SetTitle("#eta");
      cos_phi_dn[n_E][n_Pt]->GetYaxis()->SetTitle("#eta");
    };
  };


  // str tree loop - fill appropriate cos_phi_* TH1Fs
  printf("reading str\n");
  Int_t n_E_set;
  Int_t n_Pt_set;
  Int_t ent = tr->GetEntries();
  //ent = 1000;
  for(Int_t i=0; i<ent; i++)
  {
    if(!(i%10000)) printf("%d/%d\n",i,ent-1);
    tr->GetEntry(i);
    n_E_set=0;
    n_Pt_set=0;
    while(!(E12>=E_arr[n_E_set] && E12<E_arr[n_E_set]+E_div)) n_E_set++;
    while(!(Pt>=Pt_arr[n_Pt_set] && Pt<Pt_arr[n_Pt_set]+Pt_div)) n_Pt_set++;
    //printf("E bin: %.1f<=%.1f<%.1f\n",E_arr[n_E_set],E12,E_arr[n_E_set]+E_div);
    //printf("Pt bin: %.1f<=%.1f<%.1f\n",Pt_arr[n_Pt_set],Pt,Pt_arr[n_Pt_set]+Pt_div);
    if(spin==0 || spin==1) cos_phi_dn[n_E_set][n_Pt_set]->Fill(cos(Phi),Eta);
    else if(spin==2 || spin==3) cos_phi_up[n_E_set][n_Pt_set]->Fill(cos(Phi),Eta);
  };
  printf("end of str read loop\n");

  
  // calculate AN
  TH2F * blue_an[Ebins][Ptbins];
  TH2F * blue_an_2d_pos = 
    new TH2F("blue_an_2d_pos","Blue A_{N} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * blue_an_2d_neg = 
    new TH2F("blue_an_2d_neg","Blue A_{N} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  blue_an_2d_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  blue_an_2d_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  blue_an_2d_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  blue_an_2d_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  TH1F * an_dist = new TH1F("an_dist","Blue A_{N} distribution",20,-0.001,0.001);
  TF1 * fit_func;
  char blue_an_n[Ebins][Ptbins][32];
  char blue_an_t[Ebins][Ptbins][64];
  Float_t an;
  Bool_t full_cos_dist;
  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      full_cos_dist=1;
      for(Int_t nn=1; nn<=Ebins; nn++)
      {
        if(cos_phi_up[n_E][n_Pt]->GetBinContent(nn)==0 != 
           cos_phi_dn[n_E][n_Pt]->GetBinContent(nn)==0) full_cos_dist=0;
      };
      sprintf(blue_an_n[n_E][n_Pt],"blue_an_E%d_Pt%d",n_E,n_Pt);
      sprintf(blue_an_t[n_E][n_Pt],"Blue A_{N}  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      if(full_cos_dist && cos_phi_up[n_E][n_Pt]->GetEntries() * cos_phi_dn[n_E][n_Pt]->GetEntries())
      {
        blue_an[n_E][n_Pt] =
         new TH2F((*(cos_phi_up[n_E][n_Pt])-*(cos_phi_dn[n_E][n_Pt]))/(*(cos_phi_up[n_E][n_Pt])+*(cos_phi_dn[n_E][n_Pt])));
        blue_an[n_E][n_Pt]->SetName(blue_an_n[n_E][n_Pt]);
        blue_an[n_E][n_Pt]->SetTitle(blue_an_t[n_E][n_Pt]);
        blue_an[n_E][n_Pt]->Sumw2();
        /*
        blue_an[n_E][n_Pt]->Fit("pol0","","",-1,1);
        fit_func = (TF1*) blue_an[n_E][n_Pt]->GetFunction("pol0");
        if(fit_func!=NULL)an = fit_func->GetParameter(0);
        else an=0;
        if(an>0) blue_an_2d_pos->Fill(Pt_arr[n_Pt]+Pt_div/2,E_arr[n_E]+E_div/2,an);
        else if(an<0) blue_an_2d_neg->Fill(Pt_arr[n_Pt]+Pt_div/2,E_arr[n_E]+E_div/2,-1*an);
        if(an!=0) an_dist->Fill(an);
        */
      }
      else blue_an[n_E][n_Pt] = new TH2F(blue_an_n[n_E][n_Pt],blue_an_t[n_E][n_Pt],COS_BINS,-1,1,COS_BINS,2.5,4.5);
    };
  };
  printf("done calculating AN\n");
   

  // write to root file and draw pdfs
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  char outfilename[64];
  char pdfname[64];
  sscanf(filename,"RedOutputset%s",outfilename);
  sprintf(pdfname,"pdfset/cosphi%s.pdf",outfilename);
  sprintf(outfilename,"outset/SpinSet%s",outfilename);
  char pdfnameL[64]; sprintf(pdfnameL,"%s(",pdfname);
  char pdfnameR[64]; sprintf(pdfnameR,"%s)",pdfname);
  TFile * outfile = new TFile(outfilename,"RECREATE");
  TCanvas * cc = new TCanvas("cc","cc",1200,900);
  /*
  cc->SetLogz();
  cc->SetGrid(1,1);
  blue_an_2d_pos->Draw("colz");
  cc->Print(pdfnameL,"pdf");
  blue_an_2d_neg->Draw("colz");
  cc->Print(pdfname,"pdf");
  gStyle->SetOptStat(1);
  an_dist->Draw();
  cc->Print(pdfname,"pdf");
  gStyle->SetOptStat(0);
  */
  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      cc->Clear();
      cc->Divide(2,2);
      cc->cd(1); cos_phi_up[n_E][n_Pt]->Draw("colz");
      cc->cd(2); cos_phi_dn[n_E][n_Pt]->Draw("colz");
      cc->cd(3); blue_an[n_E][n_Pt]->Draw("colz");
      cc->Update();
      if(n_E+n_Pt==0) cc->Print(pdfnameL,"pdf");
      else if(n_E+n_Pt==Ebins+Ptbins-2) cc->Print(pdfnameR,"pdf");
      else cc->Print(pdfname,"pdf");
      //printf("cos_phi_up[%d][%d]->GetEntries()=%d\n",n_E,n_Pt,cos_phi_up[n_E][n_Pt]->GetEntries());
      //printf("cos_phi_dn[%d][%d]->GetEntries()=%d\n",n_E,n_Pt,cos_phi_dn[n_E][n_Pt]->GetEntries());
      cos_phi_up[n_E][n_Pt]->Write();
      cos_phi_dn[n_E][n_Pt]->Write();
    };
  };
  //blue_an_2d_pos->Write();
  //blue_an_2d_neg->Write();

}
