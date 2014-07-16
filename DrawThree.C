// draws A_LL, blue A_L, yellow A_L on one plot

void DrawThree(const char * jtype="pi0", const char * filetype="png", const char * asym_file="spin.root")
{
  // open root file
  TFile * asym_tfile = new TFile(asym_file,"READ");

  // set binning
  Int_t phi_bins0, eta_bins0, pt_bins0, en_bins0;
  if(gSystem->Getenv("PHI_BINS")==NULL){fprintf(stderr,"ERROR: source env vars\n"); return;};
  sscanf(gSystem->Getenv("PHI_BINS"),"%d",&phi_bins0);
  sscanf(gSystem->Getenv("ETA_BINS"),"%d",&eta_bins0);
  sscanf(gSystem->Getenv("PT_BINS"),"%d",&pt_bins0);
  sscanf(gSystem->Getenv("EN_BINS"),"%d",&en_bins0);

  // set jet type
  char jtype_str[32];
  if(!strcmp(jtype,"sph")) strcpy(jtype_str,"single #gamma :: ");
  else if(!strcmp(jtype,"pi0")) strcpy(jtype_str,"#pi^{0} :: ");
  else if(!strcmp(jtype,"thr")) strcpy(jtype_str,"N_{#gamma}>2 :: ");
  else strcpy(jtype_str,"");

  // define asymmetry kinematic dependence plots
  const Int_t asym_bins=4;
  char plotname[asym_bins][128];
  char plottitle[256];
  char xaxistitle[64];
  if(pt_bins0==1 && en_bins0!=1) 
  {
    strcpy(plotname[1],"/A_L_yellow/en_dep_a1_p0");
    strcpy(plotname[2],"/A_L_blue/en_dep_a2_p0");
    strcpy(plotname[3],"/A_LL/en_dep_a3_p0");
    strcpy(plottitle,"Helicity Asymmetries vs. E_{#gamma#gamma}");
    strcpy(xaxistitle,"E_{#gamma#gamma} (GeV)");
  }
  else if(pt_bins0!=1 && en_bins0==1) 
  {
    strcpy(plotname[1],"/A_L_yellow/pt_dep_a1_e0");
    strcpy(plotname[2],"/A_L_blue/pt_dep_a2_e0");
    strcpy(plotname[3],"/A_LL/pt_dep_a3_e0");
    strcpy(plottitle,"Helicity Asymmetries vs. p_{#perp}");
    strcpy(xaxistitle,"P_{#perp}  (GeV/c)");
  }
  else if(pt_bins0==1 && en_bins0==1)
  {
    strcpy(plotname[1],"/A_L_yellow/pt_dep_a1_e0");
    strcpy(plotname[2],"/A_L_blue/pt_dep_a2_e0");
    strcpy(plotname[3],"/A_LL/pt_dep_a3_e0");
    strcpy(plottitle,"Helicity Asymmetries -- single bin");
    strcpy(xaxistitle,"single bin");
  }
  else
  {
    printf("\n<><><><><><><><><><>\n\n");
    printf("pt_bins>1 && en_bins>1 ----------> three.png NOT DRAWN\n");
    printf("spin.root produced\n");
    return;
  };
  sprintf(plottitle,"%s%s",jtype_str,plottitle);


  TGraphErrors * A_LL_gr = (TGraphErrors*) asym_tfile->Get(plotname[3]);
  TGraphErrors * A_Lb_gr = (TGraphErrors*) asym_tfile->Get(plotname[2]);
  TGraphErrors * A_Ly_gr = (TGraphErrors*) asym_tfile->Get(plotname[1]);


  // fits
  Double_t min_kin, max_kin,tmp;
  A_LL_gr->GetPoint(0,min_kin,tmp);
  A_LL_gr->GetPoint(A_LL_gr->GetN() - 1,max_kin,tmp);
  min_kin-=1;
  max_kin+=1;
  printf("min_kin=%f max_kin=%f\n",min_kin,max_kin);
  TF1 * A_LL_cons_fit = new TF1("A_LL_cons_fit","pol0",min_kin,max_kin);
  TF1 * A_Lb_cons_fit = new TF1("A_Lb_cons_fit","pol0",min_kin,max_kin);
  TF1 * A_Ly_cons_fit = new TF1("A_Ly_cons_fit","pol0",min_kin,max_kin);
  TF1 * A_LL_zero_fit = new TF1("A_LL_zero_fit","pol0",min_kin,max_kin);
  TF1 * A_Lb_zero_fit = new TF1("A_Lb_zero_fit","pol0",min_kin,max_kin);
  TF1 * A_Ly_zero_fit = new TF1("A_Ly_zero_fit","pol0",min_kin,max_kin);
  A_LL_zero_fit->FixParameter(0,0);
  A_Lb_zero_fit->FixParameter(0,0);
  A_Ly_zero_fit->FixParameter(0,0);
  A_LL_gr->Fit("A_LL_cons_fit","0","",min_kin,max_kin);
  A_Lb_gr->Fit("A_Lb_cons_fit","0","",min_kin,max_kin);
  A_Ly_gr->Fit("A_Ly_cons_fit","0","",min_kin,max_kin);
  A_LL_gr->Fit("A_LL_zero_fit","0","",min_kin,max_kin);
  A_Lb_gr->Fit("A_Lb_zero_fit","0","",min_kin,max_kin);
  A_Ly_gr->Fit("A_Ly_zero_fit","0","",min_kin,max_kin);
  gStyle->SetOptFit(0);


  // define multigraph and legend
  TMultiGraph * mg = new TMultiGraph();
  A_LL_gr->SetMarkerColor(kRed);
  A_Lb_gr->SetMarkerColor(kBlue);
  A_Ly_gr->SetMarkerColor(kOrange-3);
  A_LL_gr->SetLineColor(kRed);
  A_Lb_gr->SetLineColor(kBlue);
  A_Ly_gr->SetLineColor(kOrange-3);
  A_LL_gr->SetMarkerStyle(33);
  A_Lb_gr->SetMarkerStyle(33);
  A_Ly_gr->SetMarkerStyle(33);

  mg->SetTitle(plottitle);

  TLegend * leg = new TLegend(0.1,0.9,0.3,0.7);
  leg->AddEntry(A_LL_gr,"A_{LL}","EP");
  leg->AddEntry(A_Lb_gr,"A_{L} (blue)","EP");
  leg->AddEntry(A_Ly_gr,"A_{L} (yellow)","EP");


  // get minimum and maximum range to view
  Double_t min_range, max_range;
  min_range=10000;
  max_range=-10000;
  Double_t min_tmp,max_tmp;
  Double_t pp;
  for(Int_t p=0; p<A_LL_gr->GetN(); p++)
  {
    A_LL_gr->GetPoint(p,pp,tmp);
    min_tmp = tmp - A_LL_gr->GetErrorY(p);
    max_tmp = tmp + A_LL_gr->GetErrorY(p);
    min_range = (min_tmp<min_range)?min_tmp:min_range;
    max_range = (max_tmp>max_range)?max_tmp:max_range;
  };
  for(Int_t p=0; p<A_Lb_gr->GetN(); p++)
  {
    A_Lb_gr->GetPoint(p,pp,tmp);
    min_tmp = tmp - A_Lb_gr->GetErrorY(p);
    max_tmp = tmp + A_Lb_gr->GetErrorY(p);
    min_range = (min_tmp<min_range)?min_tmp:min_range;
    max_range = (max_tmp>max_range)?max_tmp:max_range;
  };
  for(Int_t p=0; p<A_Ly_gr->GetN(); p++)
  {
    A_Ly_gr->GetPoint(p,pp,tmp);
    min_tmp = tmp - A_Ly_gr->GetErrorY(p);
    max_tmp = tmp + A_Ly_gr->GetErrorY(p);
    min_range = (min_tmp<min_range)?min_tmp:min_range;
    max_range = (max_tmp>max_range)?max_tmp:max_range;
  };
  printf("min_range = %f\n",min_range);
  printf("max_range = %f\n",max_range);

  
  // draw output
  Int_t sf=1;
  TCanvas * cc = new TCanvas("cc","cc",sf*1200,sf*1000);
  cc->SetGrid(1,1);
  A_Lb_gr->SetLineWidth(1);
  A_Ly_gr->SetLineWidth(1);
  A_LL_gr->SetLineWidth(1);
  A_Lb_gr->SetMarkerSize(2);
  A_Ly_gr->SetMarkerSize(2);
  A_LL_gr->SetMarkerSize(2);
  A_LL_gr->Print();
  mg->Add(A_Lb_gr);
  mg->Add(A_Ly_gr);
  mg->Add(A_LL_gr);
  mg->Draw("E1AP");
  //mg->GetXaxis()->SetRangeUser(0,6);
  mg->GetYaxis()->SetRangeUser(min_range-0.006,max_range+0.006); 
  //mg->GetYaxis()->SetRangeUser(-0.05,0.15); // for bx<60
  //mg->GetYaxis()->SetRangeUser(-0.07,0.07); // for bx>=60
  if(en_bins0==1 && pt_bins0==1) mg->GetYaxis()->SetRangeUser(-0.004,0.004);
  //leg->Draw();
  mg->GetXaxis()->SetTitle(xaxistitle);
  //mg->GetYaxis()->SetTitle("Asymmetry");
  cc->Update();

  char print_file[32];
  sprintf(print_file,"three.%s",filetype);
  cc->Print(print_file,filetype);
  printf("\n<><><><><><><><><><>\n\n");
  printf("three.png produced\n");
  printf("spin.root produced\n");


  // print output
  Double_t asym_LL,asym_LL_err;
  Double_t asym_Lb,asym_Lb_err;
  Double_t asym_Ly,asym_Ly_err;
  Double_t kin_pt,kin_pt_err;
  Double_t chi2_zero_LL,ndf_zero_LL,prob_zero_LL;
  Double_t chi2_cons_LL,ndf_cons_LL,prob_cons_LL,cons_LL,cons_err_LL;
  Double_t chi2_zero_Lb,ndf_zero_Lb,prob_zero_Lb;
  Double_t chi2_cons_Lb,ndf_cons_Lb,prob_cons_Lb,cons_Lb,cons_err_Lb;
  Double_t chi2_zero_Ly,ndf_zero_Ly,prob_zero_Ly;
  Double_t chi2_cons_Ly,ndf_cons_Ly,prob_cons_Ly,cons_Ly,cons_err_Ly;
  for(Int_t n=0; n<A_LL_gr->GetN(); n++)
  {
    A_LL_gr->GetPoint(n,kin_pt,asym_LL);
    A_Lb_gr->GetPoint(n,kin_pt,asym_Lb);
    A_Ly_gr->GetPoint(n,kin_pt,asym_Ly);
    asym_LL_err = A_LL_gr->GetErrorY(n);
    asym_Lb_err = A_Lb_gr->GetErrorY(n);
    asym_Ly_err = A_Ly_gr->GetErrorY(n);
    kin_pt_err = A_LL_gr->GetErrorX(n);

    chi2_zero_LL = A_LL_zero_fit->GetChisquare();
    ndf_zero_LL  = A_LL_zero_fit->GetNDF();
    prob_zero_LL = TMath::Prob(chi2_zero_LL,ndf_zero_LL);

    chi2_zero_Lb = A_Lb_zero_fit->GetChisquare();
    ndf_zero_Lb  = A_Lb_zero_fit->GetNDF();
    prob_zero_Lb = TMath::Prob(chi2_zero_Lb,ndf_zero_Lb);
    
    chi2_zero_Ly = A_Ly_zero_fit->GetChisquare();
    ndf_zero_Ly  = A_Ly_zero_fit->GetNDF();
    prob_zero_Ly = TMath::Prob(chi2_zero_Ly,ndf_zero_Ly);

    chi2_cons_LL = A_LL_cons_fit->GetChisquare();
    ndf_cons_LL  = A_LL_cons_fit->GetNDF();
    prob_cons_LL = TMath::Prob(chi2_cons_LL,ndf_cons_LL);
    cons_LL = A_LL_cons_fit->GetParameter(0);
    cons_err_LL = A_LL_cons_fit->GetParError(0);

    chi2_cons_Lb = A_Lb_cons_fit->GetChisquare();
    ndf_cons_Lb  = A_Lb_cons_fit->GetNDF();
    prob_cons_Lb = TMath::Prob(chi2_cons_Lb,ndf_cons_Lb);
    cons_Lb = A_Lb_cons_fit->GetParameter(0);
    cons_err_Lb = A_Lb_cons_fit->GetParError(0);
    
    chi2_cons_Ly = A_Ly_cons_fit->GetChisquare();
    ndf_cons_Ly  = A_Ly_cons_fit->GetNDF();
    prob_cons_Ly = TMath::Prob(chi2_cons_Ly,ndf_cons_Ly);
    cons_Ly = A_Ly_cons_fit->GetParameter(0);
    cons_err_Ly = A_Ly_cons_fit->GetParError(0);

    if(n==0) gSystem->RedirectOutput("printout.dat","w");
    else gSystem->RedirectOutput("printout.dat","a");
    printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
      kin_pt, kin_pt_err,
      asym_LL, asym_LL_err,
      asym_Lb, asym_Lb_err,
      asym_Ly, asym_Ly_err,
      chi2_zero_LL, ndf_zero_LL, prob_zero_LL,
      chi2_zero_Lb, ndf_zero_Lb, prob_zero_Lb,
      chi2_zero_Ly, ndf_zero_Ly, prob_zero_Ly,
      chi2_cons_LL, ndf_cons_LL, prob_cons_LL, cons_LL, cons_err_LL,
      chi2_cons_Lb, ndf_cons_Lb, prob_cons_Lb, cons_Lb, cons_err_Lb,
      chi2_cons_Ly, ndf_cons_Ly, prob_cons_Ly, cons_Ly, cons_err_Ly
    ); 
    gSystem->RedirectOutput(0);
  };
}
