// draws pL vs. E to show how similar they are

void CompareEtoPL()
{
  TChain * tr = new TChain("str");
  tr->Add("./redset/Red*.root");
  
  TH2D * f = new TH2D("f","E_{#gamma#gamma}-p_{#parallel} vs. p_{T}",
      100,0,37,100,0,2.5);
  str->Project("f","(E12-Pt*sinh(Eta)):Pt",
    "N12==2 && (TrigBits&0x200) && Z<0.8 && abs(M12-0.135)<0.1 && kicked==0");


  TCanvas * canv = new TCanvas("canv","canv",1200,1000);
  canv->SetLogz();
  f->Draw("colz");
};
