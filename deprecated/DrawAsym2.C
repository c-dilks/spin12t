// draws cos(phi) distributions for various E:Pt bins and computes A_LL
// - The A_LL vs. cos(phi) plot is computed manually rather than by making a new histogram
//   from the old ones using mathematical operations (as done in DrawAN2.C)

void DrawAsym2(const char * filename="RedOutputset107ha.root", Int_t nbins_E=5, Int_t nbins_Pt=5)
{
  const Int_t COS_BINS=10;
  Bool_t use_chain=true;

  // read str tree from reduced file
  char root_file[256];
  sprintf(root_file,"redset/%s",filename);
  TFile * infile = new TFile(root_file,"READ");
  
  // use single file "filename" or use chain constructed from all files (pick only one)
  //TTree * tr = (TTree*) infile->Get("str"); use_chain=false;
  TChain * tr = new TChain("str"); tr->Add("redset/*.root");
  
  // set reduced twotr branch addresses
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
  TH1F * cos_phi_bu[Ebins][Ptbins]; // blue up
  TH1F * cos_phi_bd[Ebins][Ptbins]; // blue dn
  TH1F * cos_phi_yu[Ebins][Ptbins]; // yellow up
  TH1F * cos_phi_yd[Ebins][Ptbins]; // yellow dn

  TH1F * cos_phi_buyu[Ebins][Ptbins]; // blue up & yellow up
  TH1F * cos_phi_buyd[Ebins][Ptbins]; // blue up & yellow dn
  TH1F * cos_phi_bdyu[Ebins][Ptbins]; // blue dn & yellow up
  TH1F * cos_phi_bdyd[Ebins][Ptbins]; // blue dn & yellow dn

  char cos_phi_bu_n[Ebins][Ptbins][32];
  char cos_phi_bd_n[Ebins][Ptbins][32];
  char cos_phi_yu_n[Ebins][Ptbins][32];
  char cos_phi_yd_n[Ebins][Ptbins][32];
  char cos_phi_buyu_n[Ebins][Ptbins][32];
  char cos_phi_buyd_n[Ebins][Ptbins][32];
  char cos_phi_bdyu_n[Ebins][Ptbins][32];
  char cos_phi_bdyd_n[Ebins][Ptbins][32];

  char cos_phi_bu_t[Ebins][Ptbins][128];
  char cos_phi_bd_t[Ebins][Ptbins][128];
  char cos_phi_yu_t[Ebins][Ptbins][128];
  char cos_phi_yd_t[Ebins][Ptbins][128];
  char cos_phi_buyu_t[Ebins][Ptbins][128];
  char cos_phi_buyd_t[Ebins][Ptbins][128];
  char cos_phi_bdyu_t[Ebins][Ptbins][128];
  char cos_phi_bdyd_t[Ebins][Ptbins][128];
  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      sprintf(cos_phi_bu_n[n_E][n_Pt],"c_bu_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_bd_n[n_E][n_Pt],"c_bd_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_yu_n[n_E][n_Pt],"c_yu_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_yd_n[n_E][n_Pt],"c_yd_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_buyu_n[n_E][n_Pt],"c_buyu_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_buyd_n[n_E][n_Pt],"c_buyd_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_bdyu_n[n_E][n_Pt],"c_bdyu_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_bdyd_n[n_E][n_Pt],"c_bdyd_E%d_Pt%d",n_E,n_Pt);

      sprintf(cos_phi_bu_t[n_E][n_Pt],
              "B#uparrow :: E_{#gamma#gamma}#in[%.1f,%.1f) :: P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_bd_t[n_E][n_Pt],
              "B#downarrow :: E_{#gamma#gamma}#in[%.1f,%.1f) :: P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_yu_t[n_E][n_Pt],
              "Y#uparrow :: E_{#gamma#gamma}#in[%.1f,%.1f) :: P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_yd_t[n_E][n_Pt],
              "Y#downarrow :: E_{#gamma#gamma}#in[%.1f,%.1f) :: P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_buyu_t[n_E][n_Pt],
              "B#uparrow Y#uparrow  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_buyd_t[n_E][n_Pt],
              "B#uparrow Y#downarrow  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_bdyu_t[n_E][n_Pt],
              "B#downarrow Y#uparrow  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(cos_phi_bdyd_t[n_E][n_Pt],
              "B#downarrow Y#downarrow  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);

      cos_phi_bu[n_E][n_Pt] = new TH1F(cos_phi_bu_n[n_E][n_Pt],cos_phi_bu_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_bd[n_E][n_Pt] = new TH1F(cos_phi_bd_n[n_E][n_Pt],cos_phi_bd_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_yu[n_E][n_Pt] = new TH1F(cos_phi_yu_n[n_E][n_Pt],cos_phi_yu_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_yd[n_E][n_Pt] = new TH1F(cos_phi_yd_n[n_E][n_Pt],cos_phi_yd_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_buyu[n_E][n_Pt] = new TH1F(cos_phi_buyu_n[n_E][n_Pt],cos_phi_buyu_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_buyd[n_E][n_Pt] = new TH1F(cos_phi_buyd_n[n_E][n_Pt],cos_phi_buyd_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_bdyu[n_E][n_Pt] = new TH1F(cos_phi_bdyu_n[n_E][n_Pt],cos_phi_bdyu_t[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_bdyd[n_E][n_Pt] = new TH1F(cos_phi_bdyd_n[n_E][n_Pt],cos_phi_bdyd_t[n_E][n_Pt],COS_BINS,-1,1);

      cos_phi_bu[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_bd[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_yu[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_yd[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_buyu[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_buyd[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_bdyu[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
      cos_phi_bdyd[n_E][n_Pt]->GetXaxis()->SetTitle("cos #phi");
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
    
    if(spin==0 || spin==1) cos_phi_bd[n_E_set][n_Pt_set]->Fill(cos(Phi));
    if(spin==2 || spin==3) cos_phi_bu[n_E_set][n_Pt_set]->Fill(cos(Phi));
    if(spin==0 || spin==2) cos_phi_yd[n_E_set][n_Pt_set]->Fill(cos(Phi));
    if(spin==1 || spin==3) cos_phi_yu[n_E_set][n_Pt_set]->Fill(cos(Phi));

    if(spin==0) cos_phi_bdyd[n_E_set][n_Pt_set]->Fill(cos(Phi));
    if(spin==1) cos_phi_bdyu[n_E_set][n_Pt_set]->Fill(cos(Phi));
    if(spin==2) cos_phi_buyd[n_E_set][n_Pt_set]->Fill(cos(Phi));
    if(spin==3) cos_phi_buyu[n_E_set][n_Pt_set]->Fill(cos(Phi));
  };
  printf("end of str read loop\n");


  // store sum of squares of weights
  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      //cos_phi_bu[n_E][n_Pt]->Sumw2();  // not needed since GetAsymmetry does it
      //cos_phi_bd[n_E][n_Pt]->Sumw2();
      //cos_phi_yu[n_E][n_Pt]->Sumw2();
      //cos_phi_yd[n_E][n_Pt]->Sumw2();
      cos_phi_buyu[n_E][n_Pt]->Sumw2(); // needed since these are added...
      cos_phi_buyd[n_E][n_Pt]->Sumw2();
      cos_phi_bdyu[n_E][n_Pt]->Sumw2();
      cos_phi_bdyd[n_E][n_Pt]->Sumw2();
    };
  };

  
  // calculate A_N and A_LL
  TH1F * blue_an[Ebins][Ptbins];    // distribution of blue A_N vs. cos(Phi)
  TH1F * yell_an[Ebins][Ptbins];    // distribution of yellow A_N vs. cos(Phi)
  TH1F * double_all[Ebins][Ptbins]; // distribution of A_LL vs. cos(Phi)

  // TH2Fs of asym for display (shows dependence on En and Pt)
  // -- one is for positive asym and the other for negative, so that 
  //    lack of asym in a particular En:Pt bin shows up as white in 
  //    both plots -- there might be a better way to do this.....
  TH2F * blue_an_2d_pos =
    new TH2F("blue_an_2d_pos","Blue A_{N} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * blue_an_2d_neg =
    new TH2F("blue_an_2d_neg","Blue A_{N} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * yell_an_2d_pos =
    new TH2F("yell_an_2d_pos","Yellow A_{N} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * yell_an_2d_neg =
    new TH2F("yell_an_2d_neg","Yellow A_{N} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * double_all_2d_pos = 
    new TH2F("double_all_2d_pos","A_{LL} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * double_all_2d_neg = 
    new TH2F("double_all_2d_neg","A_{LL} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);

  // error bars on AN and ALL
  TH2F * blue_an_err_pos =
    new TH2F("blue_an_err_pos","Blue #sigma_{A_{N}} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * blue_an_err_neg =
    new TH2F("blue_an_err_neg","Blue #sigma_{A_{N}} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * yell_an_err_pos =
    new TH2F("yell_an_err_pos","Yellow #sigma_{A_{N}} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * yell_an_err_neg =
    new TH2F("yell_an_err_neg","Yellow #sigma_{A_{N}} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * double_all_err_pos = 
    new TH2F("double_all_err_pos","#sigma_{A_{LL}} (positive)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);
  TH2F * double_all_err_neg = 
    new TH2F("double_all_err_neg","#sigma_{A_{LL}} (negative)",Ptbins,Pt_arr[0],Pt_arr[Ptbins-1]+Pt_div,
             Ebins,E_arr[0],E_arr[Ebins-1]+E_div);

  blue_an_2d_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  blue_an_2d_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  blue_an_2d_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  blue_an_2d_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  yell_an_2d_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  yell_an_2d_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  yell_an_2d_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  yell_an_2d_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  double_all_2d_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  double_all_2d_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  double_all_2d_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  double_all_2d_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  blue_an_err_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  blue_an_err_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  blue_an_err_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  blue_an_err_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  yell_an_err_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  yell_an_err_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  yell_an_err_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  yell_an_err_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  double_all_err_pos->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  double_all_err_pos->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");
  double_all_err_neg->GetXaxis()->SetTitle("P_{T} (Gev/c)");
  double_all_err_neg->GetYaxis()->SetTitle("E_{#gamma#gamma} (Gev/c^{2})");

  // histogram of A_N and A_LL
  TH1F * blue_an_dist = new TH1F("blue_an_dist","Blue A_{N} distribution",20,-0.001,0.001);
  TH1F * yell_an_dist = new TH1F("yell_an_dist","Yellow A_{N} distribution",20,-0.001,0.001);
  TH1F * double_all_dist = new TH1F("double_all_dist","A_{LL} distribution",20,-0.001,0.001);


  TF1 * blue_an_fit_func;
  TF1 * yell_an_fit_func;
  TF1 * double_all_fit_func;


  char blue_an_n[Ebins][Ptbins][64];
  char blue_an_t[Ebins][Ptbins][128];
  char yell_an_n[Ebins][Ptbins][64];
  char yell_an_t[Ebins][Ptbins][128];
  char double_all_n[Ebins][Ptbins][64];
  char double_all_t[Ebins][Ptbins][128];

  Float_t a,l,a_err,l_err;
  Float_t bc_bu,bc_bd,bc_yu,bc_yd;
  Float_t bc_buyu,bc_buyd,bc_bdyu,bc_bdyd;
  Float_t blue_asym,yell_asym,double_asym;

  TH1F * cos_phi_same[Ebins][Ptbins];
  TH1F * cos_phi_diff[Ebins][Ptbins];
  char cos_phi_same_n[Ebins][Ptbins][64];
  char cos_phi_diff_n[Ebins][Ptbins][64];

  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      sprintf(blue_an_n[n_E][n_Pt],"blue_an_E%d_Pt%d",n_E,n_Pt);
      sprintf(yell_an_n[n_E][n_Pt],"yell_an_E%d_Pt%d",n_E,n_Pt);
      sprintf(double_all_n[n_E][n_Pt],"double_all_E%d_Pt%d",n_E,n_Pt);

      sprintf(blue_an_t[n_E][n_Pt],"Blue A_{N}  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(yell_an_t[n_E][n_Pt],"Yellow A_{N}  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      sprintf(double_all_t[n_E][n_Pt],"A_{LL}  ::  E_{#gamma#gamma}#in[%.1f,%.1f)  ::  P_{T}#in[%.1f,%.1f)",
              E_arr[n_E],E_arr[n_E]+E_div,Pt_arr[n_Pt],Pt_arr[n_Pt]+Pt_div);
      
      // blue A_N
      blue_an[n_E][n_Pt] = new TH1F(blue_an_n[n_E][n_Pt],blue_an_t[n_E][n_Pt],COS_BINS,-1,1);
      if(cos_phi_bu[n_E][n_Pt]->GetEntries() * cos_phi_bd[n_E][n_Pt]->GetEntries())
      {
        blue_an[n_E][n_Pt] = (TH1F*)cos_phi_bu[n_E][n_Pt]->GetAsymmetry(cos_phi_bd[n_E][n_Pt]); // here
        /*
        for(Int_t bn=1; bn<=COS_BINS; bn++)
        {
          bc_bu = cos_phi_bu[n_E][n_Pt]->GetBinContent(bn);
          bc_bd = cos_phi_bd[n_E][n_Pt]->GetBinContent(bn);
          if(bc_bu>0 && bc_bd>0)
          {
            blue_asym = (bc_bu - bc_bd) / (bc_bu + bc_bd);
            blue_an[n_E][n_Pt]->SetBinContent(bn,blue_asym);
          };
        };
        */
        blue_an[n_E][n_Pt]->SetName(blue_an_n[n_E][n_Pt]);
        blue_an[n_E][n_Pt]->SetTitle(blue_an_t[n_E][n_Pt]);
        blue_an[n_E][n_Pt]->Fit("pol1","","",-1,1);
        blue_an_fit_func = (TF1*) blue_an[n_E][n_Pt]->GetFunction("pol1");
        if(blue_an_fit_func!=NULL) 
        {
          l = blue_an_fit_func->GetParameter(0);
          l_err = blue_an_fit_func->GetParError(0);
          a = blue_an_fit_func->GetParameter(1);
          a_err = blue_an_fit_func->GetParError(1);
        }
        else { a=0; a_err=0; };


        if(a>0)
        {
          blue_an_2d_pos->SetBinContent(n_Pt+1,n_E+1,a);
          blue_an_2d_pos->SetBinError(n_Pt+1,n_E+1,a_err);
          blue_an_err_pos->SetBinContent(n_Pt+1,n_E+1,a_err);
        }
        else if(a<0) 
        {
          blue_an_2d_neg->SetBinContent(n_Pt+1,n_E+1,-1*a);
          blue_an_2d_neg->SetBinError(n_Pt+1,n_E+1,a_err);
          blue_an_err_neg->SetBinContent(n_Pt+1,n_E+1,a_err);
        };

        if(a!=0) blue_an_dist->Fill(a);
      };

      // yellow A_N
      yell_an[n_E][n_Pt] = new TH1F(yell_an_n[n_E][n_Pt],yell_an_t[n_E][n_Pt],COS_BINS,-1,1);
      if(cos_phi_yu[n_E][n_Pt]->GetEntries() * cos_phi_yd[n_E][n_Pt]->GetEntries())
      {
        yell_an[n_E][n_Pt] = (TH1F*)cos_phi_yu[n_E][n_Pt]->GetAsymmetry(cos_phi_yd[n_E][n_Pt]);
        /*
        for(Int_t bn=1; bn<=COS_BINS; bn++)
        {
          bc_yu = cos_phi_yu[n_E][n_Pt]->GetBinContent(bn);
          bc_yd = cos_phi_yd[n_E][n_Pt]->GetBinContent(bn);
          if(bc_yu>0 && bc_yd>0)
          {
            yell_asym = (bc_yu - bc_yd) / (bc_yu + bc_yd);
            yell_an[n_E][n_Pt]->SetBinContent(bn,yell_asym);
          };
        };
        */
        yell_an[n_E][n_Pt]->SetName(yell_an_n[n_E][n_Pt]);
        yell_an[n_E][n_Pt]->SetTitle(yell_an_t[n_E][n_Pt]);
        yell_an[n_E][n_Pt]->Fit("pol1","","",-1,1);
        yell_an_fit_func = (TF1*) yell_an[n_E][n_Pt]->GetFunction("pol1");
        if(yell_an_fit_func!=NULL) 
        {
          l = yell_an_fit_func->GetParameter(0);
          l_err = yell_an_fit_func->GetParError(0);
          a = yell_an_fit_func->GetParameter(1);
          a_err = yell_an_fit_func->GetParError(1);
        }
        else { a=0; a_err=0; };

        if(a>0) 
        {
          yell_an_2d_pos->SetBinContent(n_Pt+1,n_E+1,a);
          yell_an_2d_pos->SetBinError(n_Pt+1,n_E+1,a_err);
          yell_an_err_pos->SetBinContent(n_Pt+1,n_E+1,a_err);
        }
        else if(a<0) 
        {
          yell_an_2d_neg->SetBinContent(n_Pt+1,n_E+1,-1*a);
          yell_an_2d_neg->SetBinError(n_Pt+1,n_E+1,a_err);
          yell_an_err_neg->SetBinContent(n_Pt+1,n_E+1,a_err);
        };

        if(a!=0) yell_an_dist->Fill(a);
      };
      

      // A_LL
      double_all[n_E][n_Pt] = new TH1F(double_all_n[n_E][n_Pt],double_all_t[n_E][n_Pt],COS_BINS,-1,1);
      sprintf(cos_phi_same_n[n_E][n_Pt],"cos_phi_same_E%d_Pt%d",n_E,n_Pt);
      sprintf(cos_phi_diff_n[n_E][n_Pt],"cos_phi_diff_E%d_Pt%d",n_E,n_Pt);
      cos_phi_same[n_E][n_Pt] = new TH1F(cos_phi_same_n[n_E][n_Pt],cos_phi_same_n[n_E][n_Pt],COS_BINS,-1,1);
      cos_phi_diff[n_E][n_Pt] = new TH1F(cos_phi_diff_n[n_E][n_Pt],cos_phi_diff_n[n_E][n_Pt],COS_BINS,-1,1);
      if(cos_phi_buyu[n_E][n_Pt]->GetEntries() * 
         cos_phi_buyd[n_E][n_Pt]->GetEntries() * 
         cos_phi_bdyu[n_E][n_Pt]->GetEntries() * 
         cos_phi_bdyd[n_E][n_Pt]->GetEntries())
      {
        cos_phi_same[n_E][n_Pt]->Add(cos_phi_buyu[n_E][n_Pt],cos_phi_bdyd[n_E][n_Pt]);
        cos_phi_diff[n_E][n_Pt]->Add(cos_phi_buyd[n_E][n_Pt],cos_phi_bdyu[n_E][n_Pt]);
        double_all[n_E][n_Pt] = (TH1F*)cos_phi_same[n_E][n_Pt]->GetAsymmetry(cos_phi_diff[n_E][n_Pt]);
        /*
        for(Int_t bn=1; bn<=COS_BINS; bn++)
        {
          bc_buyu = cos_phi_buyu[n_E][n_Pt]->GetBinContent(bn);
          bc_buyd = cos_phi_buyd[n_E][n_Pt]->GetBinContent(bn);
          bc_bdyu = cos_phi_bdyu[n_E][n_Pt]->GetBinContent(bn);
          bc_bdyd = cos_phi_bdyd[n_E][n_Pt]->GetBinContent(bn);
          // require at least one entry in bin for both up and down dists
          if(bc_buyu>0 && bc_buyd>0 && bc_bdyu>0 && bc_bdyd>0)
          {
            double_asym = (bc_buyu + bc_bdyd - bc_buyd - bc_bdyu) / (bc_buyu + bc_bdyd + bc_buyd + bc_bdyu);
            double_all[n_E][n_Pt]->SetBinContent(bn,double_asym);
          };
        };
        */

        double_all[n_E][n_Pt]->SetName(double_all_n[n_E][n_Pt]);
        double_all[n_E][n_Pt]->SetTitle(double_all_t[n_E][n_Pt]);
        double_all[n_E][n_Pt]->Fit("pol1","","",-1,1);
        double_all_fit_func = (TF1*) double_all[n_E][n_Pt]->GetFunction("pol1");
        if(double_all_fit_func!=NULL)
        {
          l = double_all_fit_func->GetParameter(0);
          l_err = double_all_fit_func->GetParError(0);
          a = double_all_fit_func->GetParameter(1);
          a_err = double_all_fit_func->GetParError(1);
        }
        else { a=0; a_err=0; };

        if(a>0) 
        {
          double_all_2d_pos->SetBinContent(n_Pt+1,n_E+1,a);
          double_all_2d_pos->SetBinError(n_Pt+1,n_E+1,a_err);
          double_all_err_pos->SetBinContent(n_Pt+1,n_E+1,a_err);
        }
        else if(a<0) 
        {
          double_all_2d_neg->SetBinContent(n_Pt+1,n_E+1,-1*a);
          double_all_2d_neg->SetBinError(n_Pt+1,n_E+1,a_err);
          double_all_err_neg->SetBinContent(n_Pt+1,n_E+1,a_err);
        };

        if(a!=0) double_all_dist->Fill(a);
      };
    };
  };
  printf("done calculating A_N and A_LL\n");

  // set scales
  Float_t maxim;
  maxim = (blue_an_2d_pos->GetMaximum() > blue_an_err_pos->GetMaximum()) ?
           blue_an_2d_pos->GetMaximum() : blue_an_err_pos->GetMaximum();
        blue_an_2d_pos->SetMaximum(maxim); blue_an_err_pos->SetMaximum(maxim);
  maxim = (blue_an_2d_neg->GetMaximum() > blue_an_err_neg->GetMaximum()) ?
           blue_an_2d_neg->GetMaximum() : blue_an_err_neg->GetMaximum();
        blue_an_2d_neg->SetMaximum(maxim); blue_an_err_neg->SetMaximum(maxim);
  maxim = (yell_an_2d_pos->GetMaximum() > yell_an_err_pos->GetMaximum()) ?
           yell_an_2d_pos->GetMaximum() : yell_an_err_pos->GetMaximum();
        yell_an_2d_pos->SetMaximum(maxim); yell_an_err_pos->SetMaximum(maxim);
  maxim = (yell_an_2d_neg->GetMaximum() > yell_an_err_neg->GetMaximum()) ?
           yell_an_2d_neg->GetMaximum() : yell_an_err_neg->GetMaximum();
        yell_an_2d_neg->SetMaximum(maxim); yell_an_err_neg->SetMaximum(maxim);
  maxim = (double_all_2d_pos->GetMaximum() > double_all_err_pos->GetMaximum()) ?
           double_all_2d_pos->GetMaximum() : double_all_err_pos->GetMaximum();
        double_all_2d_pos->SetMaximum(maxim); double_all_err_pos->SetMaximum(maxim);
  maxim = (double_all_2d_neg->GetMaximum() > double_all_err_neg->GetMaximum()) ?
           double_all_2d_neg->GetMaximum() : double_all_err_neg->GetMaximum();
        double_all_2d_neg->SetMaximum(maxim); double_all_err_neg->SetMaximum(maxim);
   

  // write to root file and draw pdfs
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  char outfilename[64];
  sscanf(filename,"RedOutputset%s",outfilename);

  char blue_pdfname[64];
  char yell_pdfname[64];
  char double_pdfname[64];
  if(!use_chain) 
  {
    sprintf(blue_pdfname,"pdfset/blue_an_%s.pdf",outfilename);
    sprintf(yell_pdfname,"pdfset/yell_an_%s.pdf",outfilename);
    sprintf(double_pdfname,"pdfset/double_all_%s.pdf",outfilename);
    sprintf(outfilename,"outset/SpinSet%s",outfilename);
  }
  else 
  {
    sprintf(blue_pdfname,"pdfset/blue_an_chain.pdf");
    sprintf(yell_pdfname,"pdfset/yell_an_chain.pdf");
    sprintf(double_pdfname,"pdfset/double_all_chain.pdf");
    sprintf(outfilename,"outset/chain.root");
  };

  char blue_pdfnameL[128]; sprintf(blue_pdfnameL,"%s(",blue_pdfname);
  char blue_pdfnameR[128]; sprintf(blue_pdfnameR,"%s)",blue_pdfname);
  char yell_pdfnameL[128]; sprintf(yell_pdfnameL,"%s(",yell_pdfname);
  char yell_pdfnameR[128]; sprintf(yell_pdfnameR,"%s)",yell_pdfname);
  char double_pdfnameL[128]; sprintf(double_pdfnameL,"%s(",double_pdfname);
  char double_pdfnameR[128]; sprintf(double_pdfnameR,"%s)",double_pdfname);

  TFile * outfile = new TFile(outfilename,"RECREATE");
  TCanvas * cc = new TCanvas("cc","cc",1200,900);

  cc->Divide(2,2);
  for(Int_t i=1; i<=4; i++)
  {
    cc->GetPad(i)->SetGrid(1,1);
    if(i>2) cc->GetPad(i)->SetLogz();
  };

  gStyle->SetOptStat(0);

  cc->cd(1); blue_an_2d_pos->Draw("colz");
  cc->cd(2); blue_an_err_pos->Draw("colz");
  cc->cd(3); blue_an_2d_pos->Draw("colz");
  cc->cd(4); blue_an_err_pos->Draw("colz");
  cc->Print(blue_pdfnameL,"pdf");
  cc->cd(1); blue_an_2d_neg->Draw("colz");
  cc->cd(2); blue_an_err_neg->Draw("colz");
  cc->cd(3); blue_an_2d_neg->Draw("colz");
  cc->cd(4); blue_an_err_neg->Draw("colz");
  cc->Print(blue_pdfname,"pdf");

  cc->cd(1); yell_an_2d_pos->Draw("colz");
  cc->cd(2); yell_an_err_pos->Draw("colz");
  cc->cd(3); yell_an_2d_pos->Draw("colz");
  cc->cd(4); yell_an_err_pos->Draw("colz");
  cc->Print(yell_pdfnameL,"pdf");
  cc->cd(1); yell_an_2d_neg->Draw("colz");
  cc->cd(2); yell_an_err_neg->Draw("colz");
  cc->cd(3); yell_an_2d_neg->Draw("colz");
  cc->cd(4); yell_an_err_neg->Draw("colz");
  cc->Print(yell_pdfname,"pdf");

  cc->cd(1); double_all_2d_pos->Draw("colz");
  cc->cd(2); double_all_err_pos->Draw("colz");
  cc->cd(3); double_all_2d_pos->Draw("colz");
  cc->cd(4); double_all_err_pos->Draw("colz");
  cc->Print(double_pdfnameL,"pdf");
  cc->cd(1); double_all_2d_neg->Draw("colz");
  cc->cd(2); double_all_err_neg->Draw("colz");
  cc->cd(3); double_all_2d_neg->Draw("colz");
  cc->cd(4); double_all_err_neg->Draw("colz");
  cc->Print(double_pdfname,"pdf");

  cc->Clear();
  gStyle->SetOptStat(1);
  blue_an_dist->Draw();
  cc->Print(blue_pdfname,"pdf");
  yell_an_dist->Draw();
  cc->Print(yell_pdfname,"pdf");
  double_all_dist->Draw();
  cc->Print(double_pdfname,"pdf");

  //gStyle->SetOptStat(0);
  for(Int_t n_E=0; n_E<Ebins; n_E++)
  {
    for(Int_t n_Pt=0; n_Pt<Ptbins; n_Pt++)
    {
      if(blue_an[n_E][n_Pt]->GetEntries())
      {
        cc->Clear();
        cc->Divide(2,2);
        cc->cd(1); cos_phi_bu[n_E][n_Pt]->Draw("e");
        cc->cd(2); cos_phi_bd[n_E][n_Pt]->Draw("e");
        cc->cd(3); blue_an[n_E][n_Pt]->Draw("e");
        cc->Update();
        cc->Print(blue_pdfname,"pdf");
      };
      if(yell_an[n_E][n_Pt]->GetEntries())
      {
        cc->Clear();
        cc->Divide(2,2);
        cc->cd(1); cos_phi_yu[n_E][n_Pt]->Draw("e");
        cc->cd(2); cos_phi_yd[n_E][n_Pt]->Draw("e");
        cc->cd(3); yell_an[n_E][n_Pt]->Draw("e");
        cc->Update();
        cc->Print(yell_pdfname,"pdf");
      };
      if(double_all[n_E][n_Pt]->GetEntries())
      {
        cc->Clear();
        cc->Divide(3,2);
        cc->cd(1); cos_phi_buyu[n_E][n_Pt]->Draw("e");
        cc->cd(4); cos_phi_bdyd[n_E][n_Pt]->Draw("e");
        cc->cd(2); cos_phi_buyd[n_E][n_Pt]->Draw("e");
        cc->cd(5); cos_phi_bdyu[n_E][n_Pt]->Draw("e");
        cc->cd(3); double_all[n_E][n_Pt]->Draw("e");
        cc->Update();
        cc->Print(double_pdfname,"pdf");
        //if(n_E+n_Pt==0) cc->Print(pdfnameL,"pdf");
        //if(n_E+n_Pt==Ebins+Ptbins-2) cc->Print(pdfnameR,"pdf");
        //else cc->Print(pdfname,"pdf");
        //printf("cos_phi_up[%d][%d]->GetEntries()=%d\n",n_E,n_Pt,cos_phi_up[n_E][n_Pt]->GetEntries());
        //printf("cos_phi_dn[%d][%d]->GetEntries()=%d\n",n_E,n_Pt,cos_phi_dn[n_E][n_Pt]->GetEntries());
      };

      cos_phi_bu[n_E][n_Pt]->Write();
      cos_phi_bd[n_E][n_Pt]->Write();
      cos_phi_yu[n_E][n_Pt]->Write();
      cos_phi_yd[n_E][n_Pt]->Write();
      cos_phi_buyu[n_E][n_Pt]->Write();
      cos_phi_buyd[n_E][n_Pt]->Write();
      cos_phi_bdyu[n_E][n_Pt]->Write();
      cos_phi_bdyd[n_E][n_Pt]->Write();
    };
  };
  cc->Clear();
  cc->Print(blue_pdfnameR,"pdf");
  cc->Print(yell_pdfnameR,"pdf");
  cc->Print(double_pdfnameR,"pdf");
  blue_an_2d_pos->Write();
  blue_an_2d_neg->Write();
  yell_an_2d_pos->Write();
  yell_an_2d_neg->Write();
  double_all_2d_pos->Write();
  double_all_2d_neg->Write();
}
