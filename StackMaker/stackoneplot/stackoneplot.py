import ROOT
from ROOT import *


def SetOverflowBin(histo):
    nbins = histo.GetNbinsX()
    histo.SetBinContent(nbins, histo.GetBinContent(nbins) + histo.GetBinContent(nbins+1)); ## Overflow
    histo.SetBinContent(1, histo.GetBinContent(1)+ histo.GetBinContent(0));                ## Underflow

def DrawText(X,Y,txt):
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextFont(42)
    text.SetTextSize(0.04)
    text.DrawLatex(X,Y,txt)

    return text
    
def decorate(h,color):
    h.SetLineColor(color)
    h.SetFillColor(color)
    SetOverflowBin(h) ## overflow bin is must

def PadStyling(pad,rpad):
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.20)
    pad.SetTopMargin(0.09)
    pad.SetBottomMargin(0.01)
    pad.SetTickx(1)
    pad.SetTicky(1)

    rpad.SetLeftMargin(0.15)
    rpad.SetRightMargin(0.20)
    rpad.SetTopMargin(0.02)
    rpad.SetBottomMargin(0.40)
    rpad.SetTickx(1)
    rpad.SetTicky(1)
    rpad.SetGrid(1)

def SetLegendStyle(legend):
    legend.SetTextFont(62)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.024)    
    
####################################
def main():
    
    # Read the file from the disk    
    inputDir = "../finalhaddOutput/nov15/"
    
    file_DY     =  TFile.Open(inputDir + "DY.root",     "READ")
    file_ttbar  =  TFile.Open(inputDir + "TTBar.root",  "READ")
    file_QCD    =  TFile.Open(inputDir + "QCD.root",    "READ")
    file_WJets  =  TFile.Open(inputDir + "WJets.root",  "READ")
    file_ZGamma =  TFile.Open(inputDir + "ZGamma.root", "READ")
    file_WGamma =  TFile.Open(inputDir + "WGamma.root", "READ")
    file_Data   =  TFile.Open(inputDir + "Data.root",   "READ")
    
    print(" \n file opened in ROOT successfully..")
    
    
    ## Get the histograms

    plotname = "1l2d_l0_reliso03"
    h_dy     = file_DY.Get(plotname);      decorate(h_dy,     kBlue-9);
    h_tt     = file_ttbar.Get(plotname);   decorate(h_tt,     kGreen-9); 
    h_qcd    = file_QCD.Get(plotname);     decorate(h_qcd,    kYellow-9);
    h_wjets  = file_WJets.Get(plotname);   decorate(h_wjets,  kRed-9);
    h_zgamma = file_ZGamma.Get(plotname);  decorate(h_zgamma, kCyan-9);
    h_wgamma = file_WGamma.Get(plotname);  decorate(h_wgamma, kMagenta-9);
    
    h_data = file_Data.Get(plotname);
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(0.6)
    h_data.SetLineColor(kBlack)

    SetOverflowBin(h_dy)
    SetOverflowBin(h_tt)
    SetOverflowBin(h_qcd)
    SetOverflowBin(h_wjets)
    SetOverflowBin(h_zgamma)
    SetOverflowBin(h_wgamma)
    SetOverflowBin(h_data)
   
    print("Histograms are ready...")

    '''
    ##Scale
    dtlumi = 59.8*1000;
    dylumi = 99717900/5765.0; 
    ttlumi = 28701360/88.29;    
    wzlumi = 10527550/5.052;
    
    h_dy.Scale(dtlumi/dylumi)
    h_tt.Scale(dtlumi/ttlumi)
    h_wz.Scale(dtlumi/wzlumi)
    '''

    
    ##Rebining
    rebin = 1
    h_dy.Rebin(rebin)
    h_tt.Rebin(rebin)
    h_qcd.Rebin(rebin)
    h_wjets.Rebin(rebin)
    h_zgamma.Rebin(rebin)
    h_wgamma.Rebin(rebin)   
    
    h_data.Rebin(rebin)

    # Create a list of tuples with histograms and their integrals
    
    histograms = [(h_dy,     h_dy.Integral(),     "DY"),
                  (h_tt,     h_tt.Integral(),     "TTbar"),
                  (h_qcd,    h_qcd.Integral(),    "QCD"),
                  (h_wjets,  h_wjets.Integral(),  "WJets"),
                  (h_zgamma, h_zgamma.Integral(), "ZGamma"),
                  (h_wgamma, h_wgamma.Integral(), "WGamma")]
    
    # Sort the list in descending order based on integrals
    histograms.sort(key=lambda x: x[1], reverse=True)

    #Scaling

    print("number of bins", h_data.GetNbinsX())
    bin_lo = 5
    bin_hi = h_data.GetNbinsX()
    
    h_otherbkg = h_dy.Clone()
    for hist, integral, name in histograms:
        if name!="QCD" and name!="DY":
            h_otherbkg.Add(hist)
            
            
    data_integral = h_data.Integral(bin_lo,bin_hi)
    qcd_integral = h_qcd.Integral(bin_lo,bin_hi)
    otherbkg_integral = h_otherbkg.Integral(bin_lo,bin_hi)
    print("data", data_integral)
    print("qcd", qcd_integral)
    print("other bkg", otherbkg_integral)

    qcd_sf = (data_integral-otherbkg_integral)/qcd_integral

    print("qcd before scaling", h_qcd.Integral(1,bin_hi))

    #h_qcd.Scale(qcd_sf)
   
    print("qcd after scaling", h_qcd.Integral(1,bin_hi))
    
    # Create THStack 
    h_stack = THStack()
    h_stack.Add(h_dy)
    h_stack.Add(h_tt)
    h_stack.Add(h_wjets)
    h_stack.Add(h_qcd)
    h_stack.Add(h_zgamma)
    h_stack.Add(h_wgamma)
    #for hist, integral, name in histograms:
     #   h_stack.Add(hist)
        
    h_bkg = histograms[0][0].Clone()
    for hist, integral, name in histograms[1:]:
        h_bkg.Add(hist)

    h_ratio = h_data.Clone()
    h_ratio.Divide(h_bkg)
    #h_ratio.GetXaxis().SetRangeUser(0,1)

    print("Stacking done....time for plotting!")
            
    # Legend
    legend = TLegend(0.95, 0.50, 0.80, 0.86)
    ratioleg = TLegend(0.90, 0.90, 0.81, 0.87)
    ratioleg.SetHeader(f"obs/exp={h_data.Integral()/h_bkg.Integral():.2f}   exp: {h_bkg.Integral():.0f}")

    legend.AddEntry(h_data, f"Data[{h_data.Integral():.0f}]", 'ep')
    for hist, integral, name in histograms:
        legend.AddEntry(hist, f"{name}[{integral:.0f}]", 'lf')
                
    SetLegendStyle(ratioleg)
    SetLegendStyle(legend)
                

    ######################################################################################
    ##      PLOTTING START
    ######################################################################################
    canvas = TCanvas("c","canvas",650,600)
    gStyle.SetOptStat(0)
    
    ratioPadSize = 0.3
    mainPad  = TPad("pad","pad",0,ratioPadSize,1,1)
    ratioPad = TPad("pad2","pad2",0,0,1.0,ratioPadSize)
    PadStyling(mainPad,ratioPad)
    mainPad.Draw()
    ratioPad.Draw()
    
    mainPad.cd()
    mainPad.SetLogy(1)
    h_stack.SetMinimum(0.001)
    h_stack.SetMaximum(1e8)
    
   
   
   # h_stack.GetYaxis().SetRangeUser(0,100)
    
    #h_data.GetYaxis().SetRangeUser(0.001,1e8)
    
    

    #h_data.Draw('ep')
    
    h_stack.Draw("HIST")
    h_data.Draw('ep same')

    h_stack.GetXaxis().SetRangeUser(0,1)
    h_data.GetXaxis().SetRangeUser(0,1)

    #mainPad.Modified()
    #mainPad.Update()
    
    print("h_stack address:", h_stack)  # Check if it's not None
    
    ratioleg.Draw()
    legend.Draw()
    
    h_stack.GetYaxis().SetTitle('Events')
    #h_stack.GetXaxis().SetTitle("Invariant Mass(M_{12})")
    #h_stack.GetXaxis().CenterTitle()
    h_stack.GetYaxis().CenterTitle()
    
    #beautification
    h_stack.GetXaxis().SetTitleFont(43)
    h_stack.GetXaxis().SetTitleSize(20)
    h_stack.GetXaxis().SetTitleOffset(0.8)
    h_stack.GetXaxis().SetLabelFont(43)
    h_stack.GetXaxis().SetNdivisions(513)
    h_stack.GetYaxis().SetTitleFont(43)
    h_stack.GetYaxis().SetTitleSize(20)
    h_stack.GetYaxis().SetTitleOffset(1.2)
    h_stack.GetYaxis().SetLabelFont(43)
    h_stack.GetYaxis().SetLabelSize(12)
    h_stack.GetYaxis().SetNdivisions(513)
    mainPad.SetTickx(1)
    
    ## A few text
    #text1 = DrawText(0.15,0.92,"IISER Pune Analysis")
    #text1 = DrawText(0.66,0.92,"13 TeV(2016)")

    ##plot ratiopad
    ratioPad.cd()
    h_ratio.GetXaxis().SetTitle("l0 reliso03")
    h_ratio.GetYaxis().SetTitle("obs/exp")
    h_ratio.GetXaxis().CenterTitle()
    h_ratio.GetYaxis().CenterTitle()
    #beautification
    h_ratio.GetXaxis().SetTitleFont(43)
    h_ratio.GetXaxis().SetTitleSize(20)
    h_ratio.GetXaxis().SetTitleOffset(1.2)
    h_ratio.GetXaxis().SetLabelFont(43)
    h_ratio.GetXaxis().SetLabelSize(12)
    h_ratio.GetXaxis().SetNdivisions(513)
    h_ratio.GetYaxis().SetTitleFont(43)
    h_ratio.GetYaxis().SetTitleSize(20)
    h_ratio.GetYaxis().SetTitleOffset(1.2)
    h_ratio.GetYaxis().SetLabelFont(43)
    h_ratio.GetYaxis().SetLabelSize(12)
    h_ratio.GetYaxis().SetNdivisions(503)
    h_ratio.GetYaxis().SetRangeUser(0,2)
    h_ratio.GetXaxis().SetRangeUser(0,1)
    
    h_ratio.SetTitle('')
    h_ratio.Draw("ep")
    #ratioPad.Update()


    canvas.Draw()

    canvas.SaveAs('stackoneoutput/stackplot.pdf')
    #canvas.SaveAs('qcdscaled_plots/qcdscaled_stackplot.pdf')
    

## Execute the main function
if __name__ == '__main__':
    main()
