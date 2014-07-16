// generates list of random errors e_i[] according to selected
// RNG procedure (gaussian, poissonian, et al) and 
// propagates the uncertainty through average, to
// see how the propagated uncertainty depends on the number 
// of items in the list

void average_error_test()
{
  const Int_t TOTAL = 10000;
  Float_t e_i[TOTAL];
  Float_t e[TOTAL];
  Float_t sum_of_squares=0;
  Float_t i_arr[TOTAL];

  TRandom * rng = new TRandom();

  TH1F * e_i_hist = new TH1F("e_i_hist","uncertainties distribution",100,0,10);

  for(Int_t i=0; i<TOTAL; i++)
  {
    e_i[i] = rng->Gaus(5,5);
    //e_i[i] = rng->Poisson(5);
    sum_of_squares += pow(e_i[i],2);

    i_arr[i] = (Float_t)(i+1);
    e[i] = 1/i_arr[i]*sqrt(sum_of_squares);

    e_i_hist->Fill(e_i[i]);

    printf("%d %f %f \n",i+1,sum_of_squares,e[i]);
  };

  TGraph * gr = new TGraph(TOTAL,i_arr,e);
  gr->SetTitle("propagated uncertainty vs. N");
  gr->SetMarkerStyle(kFullCircle);

  TCanvas * c1 = new TCanvas("c1","c1",700,500);
  c1->SetGrid(1,1);
  e_i_hist->Draw();

  TCanvas * c2 = new TCanvas("c2","c2",700,500);
  c2->SetGrid(1,1);
  gr->Draw("ape");
}
