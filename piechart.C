void piechart()
{
  Float_t vals[] = {24271,7625,4291,1569,87,85};
  Int_t colors[] = {kYellow-9,kGreen-9,kRed-9,kBlue-9,kMagenta-9,kCyan-9};
  //Int_t colors[] = {500,516,732,700,716,532};
  Int_t nvals = sizeof(vals)/sizeof(vals[0]);
 
  TCanvas *cpie = new TCanvas("cpie","TPie test",1500,1000);
  //cpie->Divide(2,2);
 
  TPie *pie4 = new TPie("pie4", "2l1d event selection",nvals,vals,colors);

  TString filenames[] = {"QCD","TTBar","WJets","DY","WGamma","ZGamma"};

  for (Int_t i = 0; i < nvals; ++i)
    {
      // Set custom labels for each slice
      pie4->SetEntryLabel(i, filenames[i]);
    }
 
  pie4->SetLabelsOffset(-0.1);
  pie4->SetRadius(.4);
  //pie4->GetSlice(1)->SetRadiusOffset(0.1);
  //pie4->GetSlice(4)->SetRadiusOffset(0.05);   
  //pie4->SetLabelFormat("#splitline{%val (%perc)}{%txt}");
  pie4->SetLabelFormat("%txt");
  //pie4->Draw("3d");
  pie4->Draw("3d");
}

