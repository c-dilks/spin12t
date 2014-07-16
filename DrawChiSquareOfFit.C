// draws chisquare / ndf vs. run (or fill) set number;

void DrawChiSquareOfFit(const char * filename="chisq")
{
  TTree * tr = new TTree();
  tr->ReadFile(filename,"num/I:chisq/F:ndf/F");
  Int_t num;
  Float_t chi,ndf;
  tr->SetBranchAddress("num",&num);
  tr->SetBranchAddress("chisq",&chi);
  tr->SetBranchAddress("ndf",&ndf);
  Int_t max_num = tr->GetMaximum("num");
  Int_t min_num = tr->GetMinimum("num");

  // get max chi2/ndf and add point to chi2/ndf vs. region plot
  Float_t max_cn = 0;
  Int_t cn_reg_pt=0;
  TGraph * cn_reg = new TGraph();
  cn_reg->SetTitle("#chi^{2}/ndf vs. region no.");
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    if(ndf>0)
    {
      max_cn = (chi/ndf > max_cn) ? chi/ndf:max_cn;
      cn_reg->SetPoint(cn_reg_pt,num,chi/ndf);
      cn_reg_pt++;
    };
  };

  // fit
  gStyle->SetOptFit(1);
  cn_reg->Fit("pol0","","",min_num,max_num);

  // chi2/ndf distribution
  TH1F * cn_dist = new TH1F("cn_dist","#chi^{2}/ndf vs. region no.",20,0,max_cn);
  tr->Project("cn_dist","chisq/ndf","ndf>0");
  
  // draw
  TCanvas * canv = new TCanvas("canv","canv",700,1000);
  canv->Divide(1,2);
  canv->cd(1); cn_dist->Draw();
  cn_reg->SetMarkerStyle(kFullCircle);
  canv->cd(2); cn_reg->Draw("APE");
}


// code below fits the chi2/ndf for all points with chi2/ndf<10

/*
void fffff(const char * filename="chisq")
{
  TTree tr;
  Float_t chisq;
  Int_t num;
  Float_t ndf;
  tr.ReadFile(filename,"num/I:chisq/F:ndf/F");
  tr.SetBranchAddress("chisq",&chisq);
  tr.SetBranchAddress("ndf",&ndf);
  tr.SetBranchAddress("num",&num);
  TGraph gr;
  Int_t pt=0;
  for(Int_t i=0; i<tr.GetEntries(); i++)
  {
    tr.GetEntry(i);
    if(ndf>0)
    {
      if(chisq/ndf < 10)
      {
        gr.SetPoint(pt,num,chisq/ndf);
        printf("%d\n",pt);
        pt++;
      };
    };
  };
  gr.Print();
  gr.SetMarkerStyle(kFullCircle);
  gr.SetMarkerColor(kRed);
  gr.Draw("AP");
  gr.Fit("pol0","","",1,21);
};
*/
