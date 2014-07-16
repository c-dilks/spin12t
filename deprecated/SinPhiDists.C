// -- create sin phi distributions with various kinematic cuts (eta,pt,E) for each spinbit
//    (kinematic cuts are set as const Double_t's below)
// -- phi distributions for "+ -" and "- +" are weighted by relative luminosity
// -- phi distributions are written to phiset/ directory with similar name; 
//    they are named: phi_s[spinbit]_g[eta bin]_p[pt bin]_e[en bin]

void SinPhiDists(const char * filename="RedOutputset120ha.root",Int_t phi_bins0=16,
              Int_t eta_bins0=1, Int_t pt_bins0=3, Int_t en_bins0=3)
{
  // set output file name and read input tree
  char setname[32];
  sscanf(filename,"RedOutputset%s",setname);
  sprintf(filename,"redset/%s",filename);
  TFile * infile = new TFile(filename,"READ");
  TTree * tree = (TTree*) infile->Get("str");
  char outname[256];
  sprintf(outname,"phiset/phi%s",setname);


  // set ranges (set # bins = 1 to average over all values from low to high)
  const Double_t pi=3.1415;
  const Double_t phi_bins=phi_bins0; 
    const Double_t phi_low=-1; const Double_t phi_high=1;
  const Double_t eta_bins=eta_bins0;
    const Double_t eta_low=2.6; const Double_t eta_high=4.2;
  const Double_t pt_bins=pt_bins0; 
    const Double_t pt_low=0; const Double_t pt_high=10;
  const Double_t en_bins=en_bins0; 
    const Double_t en_low=0; const Double_t en_high=100;


  // define spinbit
  char spinbit_t[4][4];
  char spinbit_n[4][4];
  strcpy(spinbit_t[0],"--"); strcpy(spinbit_n[0],"nn");
  strcpy(spinbit_t[1],"-+"); strcpy(spinbit_n[1],"np");
  strcpy(spinbit_t[2],"+-"); strcpy(spinbit_n[2],"pn");
  strcpy(spinbit_t[3],"++"); strcpy(spinbit_n[3],"pp");


  // set binning
  // -- *_div = lower limit of each bin; last one is the upper limit
  // -- *_width = bin width
  Double_t eta_div[eta_bins+1];
  Double_t eta_width = (eta_high - eta_low)/eta_bins;
  for(Int_t i=0; i<eta_bins; i++) eta_div[i] = eta_low + i * eta_width;
  Double_t pt_div[pt_bins+1];
  Double_t pt_width = (pt_high - pt_low)/pt_bins;
  for(Int_t i=0; i<pt_bins; i++) pt_div[i] = pt_low + i * pt_width;
  Double_t en_div[en_bins+1];
  Double_t en_width = (en_high - en_low)/en_bins;
  for(Int_t i=0; i<en_bins; i++) en_div[i] = en_low + i * en_width;

  eta_div[eta_bins] = eta_high;
  pt_div[pt_bins] = pt_high;
  en_div[en_bins] = en_high;

  // print ranges
  /*
  for(Int_t g=0; g<eta_bins+1; g++) printf("eta_div[%d] = %.2f\n",g,eta_div[g]);
  for(Int_t p=0; p<pt_bins+1; p++) printf("pt_div[%d] = %.2f\n",p,pt_div[p]);
  for(Int_t e=0; e<en_bins+1; e++) printf("en_div[%d] = %.2f\n",e,en_div[e]);
  */


  // define phi distributions 
  TH1D * phi_dist[4][eta_bins][pt_bins][en_bins]; // [spinbit] [eta] [pt] [en]
  char phi_dist_n[4][eta_bins][pt_bins][en_bins][64];
  char phi_dist_t[4][eta_bins][pt_bins][en_bins][256];
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(phi_dist_n[s][g][p][e],"phi_s%d_g%d_p%d_e%d",s,g,p,e);
          sprintf(phi_dist_t[s][g][p][e],
           "#phi distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f)",
           spinbit_t[s],eta_div[g],eta_div[g+1],pt_div[p],pt_div[p+1],en_div[e],en_div[e+1]);
          phi_dist[s][g][p][e] = new TH1D(phi_dist_n[s][g][p][e],phi_dist_t[s][g][p][e],
            phi_bins,phi_low,phi_high);

        };
      };
    };
  };


  // define cuts
  char kinematic_cut[4][eta_bins][pt_bins][en_bins][512];
  char spin_cut[4][16];
  char eta_cut[eta_bins][32];
  char pt_cut[pt_bins][32];
  char en_cut[en_bins][32];
  char trig_cut[32]; sprintf(trig_cut,"(TrigBits&0x200)"); // other cuts
  char z_cut[16]; sprintf(z_cut,"Z<0.8");
  char mass_cut[32]; sprintf(mass_cut,"abs(M12-0.135)<0.1");
  //gROOT->ProcessLine(".! touch kinematic_cuts; rm kinematic_cuts; touch kinematic_cuts");
  for(Int_t s=0; s<4; s++)
  {
    sprintf(spin_cut[s],"spin==%d",s);
    for(Int_t g=0; g<eta_bins; g++)
    {
      sprintf(eta_cut[g],"Eta>=%.2f && Eta<%.2f",eta_div[g],eta_div[g+1]);
      for(Int_t p=0; p<pt_bins; p++)
      {
        sprintf(pt_cut[p],"Pt>=%.2f && Pt<%.2f",pt_div[p],pt_div[p+1]);
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(en_cut[e],"E12>=%.2f && E12<%.2f",en_div[e],en_div[e+1]);
          sprintf(kinematic_cut[s][g][p][e],"%s && %s && %s && %s && %s && %s && %s && isConsistent",
            spin_cut[s],eta_cut[g],pt_cut[p],en_cut[e],trig_cut,z_cut,mass_cut);
          // weight + - and - + distributions by rellum
          if(s==1||s==2) sprintf(kinematic_cut[s][g][p][e],"R3_zdc*(%s)",kinematic_cut[s][g][p][e]);
          /*
          gSystem->RedirectOutput("kinematic_cuts","a");
          printf("%s\n",kinematic_cut[s][g][p][e]);
          gSystem->RedirectOutput(0);
          */
        };
      };
    };
  };


  // projections
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          phi_dist[s][g][p][e]->Sumw2();
          tree->Project(phi_dist_n[s][g][p][e],"sin(Phi)",kinematic_cut[s][g][p][e]);
          printf("phi_dist[%d][%d][%d][%d] entries projected: %d\n",s,g,p,e,
            phi_dist[s][g][p][e]->GetEntries());
        };
      };
    };
  };


  // write output
  TFile * outfile = new TFile(outname,"RECREATE");
  for(Int_t e=0; e<en_bins; e++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t s=0; s<4; s++)
        {
          phi_dist[s][g][p][e]->Write();
        };
      };
    };
  };
};
