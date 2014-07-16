// prints pdfs to "pdfset" from asymmetry calculation

void PrintPDFs(Int_t eta_bin, Int_t pt_bin, Int_t en_bin)
{
  TFile * infile = new TFile("spin.root","READ");
  
  // read TObjArrays
  TObjArray * phi_dist_arr[4];
  TObjArray * phi_dist_arr_dP[4];
  TObjArray * asym_arr;
  char phi_dist_arr_n[4][64];
  char phi_dist_arr_dP_n[4][64];
  char asym_arr_n[64];

  sprintf(asym_arr_n,"asym_g%d_p%d_e%d",eta_bin,pt_bin,en_bin);
  asym_arr = (TObjArray*) infile->Get(asym_arr_n);

  for(Int_t s=0; s<4; s++)
  {
    sprintf(phi_dist_arr_n[s],"phi_dist_s%d_g%d_p%d_e%d",s,eta_bin,pt_bin,en_bin);
    sprintf(phi_dist_arr_dP_n[s],"phi_dist_dP_s%d_g%d_p%d_e%d",s,eta_bin,pt_bin,en_bin);
    phi_dist_arr[s] = (TObjArray*) infile->Get(phi_dist_arr_n[s]);
    phi_dist_arr_dP[s] = (TObjArray*) infile->Get(phi_dist_arr_dP_n[s]);
    printf("phi_dist_arr[%d] @ %p\n",s,(void*)phi_dist_arr[s]);
  };


  // print pdfs
  char pdfname[256];
  char pdfname_L[256];
  char pdfname_R[256];
  TCanvas * canv = new TCanvas("canv","canv",2000,1500);
  gStyle->SetOptFit(1);
  sprintf(pdfname,"pdfset/asym_g%d_p%d_e%d.pdf",eta_bin,pt_bin,en_bin);
  sprintf(pdfname_L,"pdfset/asym_g%d_p%d_e%d.pdf(",eta_bin,pt_bin,en_bin);
  sprintf(pdfname_R,"pdfset/asym_g%d_p%d_e%d.pdf)",eta_bin,pt_bin,en_bin);
  Int_t NRUNS_tmp = phi_dist_arr[0]->GetEntries();
  const Int_t NRUNS = NRUNS_tmp;
  for(Int_t r=0; r<NRUNS; r++)
  {
    canv->Clear();
    canv->Divide(3,2);
    for(Int_t cc=1; cc<=6; cc++) canv->GetPad(cc)->SetGrid(1,1);
    canv->cd(1); phi_dist_arr[0]->At(r)->Draw();
    canv->cd(2); phi_dist_arr[3]->At(r)->Draw();
    canv->cd(4); phi_dist_arr[1]->At(r)->Draw();
    canv->cd(5); phi_dist_arr[2]->At(r)->Draw();
    canv->cd(6); asym_arr->At(r)->Draw("ape");
    if(r==0) canv->Print(pdfname_L,"pdf");
    else if(r+1==NRUNS) canv->Print(pdfname_R,"pdf");
    else canv->Print(pdfname,"pdf");
  };
};
