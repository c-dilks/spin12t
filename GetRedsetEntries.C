// prints the number of entries in the redset trees; appends a line to txt file "corrupted"
// with the filename if the number of entries is zero
//
// it's best to loop this over all root files in the directory redset, i.e.,
//
// touch corrupted; rm corrupted; 
// for file in $(ls redset/Red*.root); do root -b -q -l 'GetRedsetEntries.C('\"$file\"'); done

void GetRedsetEntries(const char * filename)
{
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("str");
  printf("%d events in %s\n",tr->GetEntries(),filename);
  if(tr->GetEntries() == 0)
  {
    gSystem->RedirectOutput("corrupted","a");
    printf("%s\n",filename);
    gSystem->RedirectOutput(0);
  };
};
