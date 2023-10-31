#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;


// Declare all the other functions needed. See the function implementations for more
// comments about what the functions do.
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill);
void decorate(TLegend *g, float textSize, TString legendheader);
float get_nevents(TH1F *hst, float bin_lo, float bin_hi);
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi);

// This is the primary function. It opens two histogram files, extracts two histograms from
// them, and plots them on top of each other (i.e. overlays them)
void overlay()
{

  //First declare the file names (just strings)
  
  TString dyfile = "cluster_hst_output/DYJetsToLL_oct19.root";
  TString datafile = "cluster_hst_output/2016Data_oct19.root";

  //Declare other constants, strings that you might need here.

  //Declare the name of the plot that you want to overlay
  //(you can open the histogram file to see other names)
  
  //TString plotname1 = "mass_mumu_2";
  //TString plotname2 = "mass_mumu_1"; //(here we are picking the same plot from other file)
  
  TString plotname1 = "zcr_invmass";
  TString plotname2 = "zcr_met"; 


  //Also give fancy name for the axis titles
  //TString xtitle = "M_{#mu#mu} [GeV]"; // Can use latex-type commands
                                               // in strings. For example
                                               // "#mu^{1} p_{T}"
                                               
  TString xtitle = "M_{ll} [GeV]";                                         
  TString ytitle = "Events"; // Or "Entries"

  //Now let us open the files
  TFile *filedy = new TFile(dyfile);
  TFile *filedata = new TFile(datafile);

  //Now open the respective histograms from the files
  TH1F *hst1 = (TH1F*)filedy->Get(plotname1);
  TH1F *hst2 = (TH1F*)filedata->Get(plotname1);


  //Decorate the histograms using function decorate 
  // See function definition below for syntax
  // See https://root.cern.ch/root/html/TColor.html for color names.
  decorate(hst1,xtitle,ytitle,"",kBlue,2,kBlue-9,20,1);
  decorate(hst2,xtitle,ytitle,"",kRed,2,kBlack,21,0);


  //Now let us set the last bin as the overflow bin
  
    int nbins = hst1->GetNbinsX();
    hst1->SetBinContent(nbins,hst1->GetBinContent(nbins+1)+hst1->GetBinContent(nbins));
    hst2->SetBinContent(nbins,hst2->GetBinContent(nbins+1)+hst2->GetBinContent(nbins));
  

  // Now rebin the histograms if needed
  // Group nrebins bins together
  /*
    int nrebins = 2;
    hdy->Rebin(nrebins); htt->Rebin(nrebins);
  */

  //Now we scale the histograms.
  /* Scaling can be done either by some outside measure (such as size of samples)
     or by normalizing the histograms such that they have the same integral.
     Conventionally, when comparing shapes, we normalize to 1, so let us 
     try that here.
  */

  //hst1->Scale(1.0/hst1->Integral());
  //hst2->Scale(1.0/hst2->Integral());

  
  float dtlumi=36.3*1000; //pb
  float dylumi=70114140/699.0; //pb

  float dysf = dtlumi/dylumi;
  
  hst1->Scale(dysf);

  float sf = hst2->Integral() / hst1->Integral();
  hst1->Scale(sf);
  cout<<"DY Scale Factor: "<<sf<<endl;
  
  //Now let us declare a canvas 
  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  // Now we define a legend
  TLegend *lg1 = new TLegend(0.55,0.50,0.85,0.75,NULL,"NDC");
  decorate(lg1,0.05,""); // Decorated the legend using function below.
  lg1->AddEntry(hst1,"DY","l"); // Added the two entries for the two histograms
  lg1->AddEntry(hst2,"Data","ep");           // we shall be drawing.

  //c1->SetLogy();
  
  //We set the stat box for the first one to be invisible first
  hst1->SetStats(0);

  // Now to draw the histograms and the legend
  hst1->Draw("hist");            
  hst2->Draw("ep same");
  lg1->Draw();

/*  
To draw with stat boxes for each histogram
     -- Dont use SetStats(0)
     -- Then draw them, first one with option Draw(), next ones with option Draw("sames")
     -- The s at the end is for stats box
 */

}

/*
This function decorates the histogram, by setting the titles of the X and Y axes,
by setting the line and marker colors, the marker style, and by deciding whether
to fill in the histogram with a color or not.
 */


void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);

  h->SetMarkerSize(1.0);
  h->SetTitle(title);
}
 

/*
This function decorates the legend, by setting the textsize, and filling in the
legend with a pure white color.
 */
void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(10);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //but if you want one, uncomment the next line.
  //g->SetHeader(legendheader);
}


// Here are a couple of other utility functions

// For a given histogram hst, return the number of entries between bin_lo and bin_hi
float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents += hst->GetBinContent(i);

  return nevents;
}
// Partner function for above, returning the error for the above nevents
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents_err = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents_err += pow(hst->GetBinError(i),2);
  nevents_err = sqrt(nevents_err);

  return nevents_err;
}
