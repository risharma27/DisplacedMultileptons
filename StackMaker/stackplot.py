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
    
######################################################################################################################################################################

def main():
    
    # Read the file from the disk    
    inputDir = "../cluster_hst_output/"

    file_names = [
        
        "ttbar3LCR2018_DY_MC_May13_v5.root",
        "ttbar3LCR2018_ttbar_MC_May13_v5.root",
        "ttbar3LCR2018_ZZ_MC_May13_v5.root",
        "ttbar3LCR2018_WZ_MC_May13_v5.root",
        "ttbar3LCR2018_TTZ_MC_May13_v5.root",
        "ttbar3LCR2018_TTW_MC_May13_v5.root",
        "ttbar3LCR2018_Data_May13_v5.root",
        
    ]  # this is a list of names of all root files
    

    # open all the ROOT files

    files = []
    for file_name in file_names:
     myfile = TFile.Open(inputDir + file_name, "READ")
     files.append(myfile)
        
    print(" \n files opened in ROOT successfully..")

    # get a list of histogram names from one of the files (all the root files have the same set of histograms)
    
    file0 = files[0]
    list_of_histograms = file0.GetListOfKeys()

    hist_names = []
    for hist in list_of_histograms:
     histname = hist.GetName()
     hist_names.append(histname)

    hist_colors = [kRed, kGreen, kBlue, kMagenta, kOrange, kCyan, kViolet, kYellow, kGray]
    
    #########################
    #  Stacking Histograms 
    #########################

    # Create a THStack to stack histograms
    h_stack = THStack()

    # Create a TLegend to hold legend entries
    legend = TLegend(0.75, 0.7, 0.9, 0.9)
    
    # Loop through each histogram name
    for i, plotname in enumerate(hist_names):
        # Initialize a list to store histograms with the same name from all files
        histograms = []

        # Loop through all the files and retrieve histograms with the same name
        for j, file_handle in enumerate(files):
            histo = file_handle.Get(plotname)
            if histo:
                color = hist_colors[j % len(colors)]  # get a color from the list
                decorate(histo, color)  # customize the color
                SetOverflowBin(histo)   # handle overflow and underflow bins
                histograms.append(histo)

        # Add histograms with the same name to the THStack
        for histogram in histograms:
            h_stack.Add(histogram)
            legend.AddEntry(histogram, f"{file_names[i]}", "lf")

    # create a TCanvas to display the stacked histogram
    canvas = TCanvas("c", "canvas", 650, 600)
    gStyle.SetOptStat(0)

    # draw the stacked histogram
    h_stack.Draw("HIST")

    # add legends, axis labels, and other decorations as needed
    legend.Draw()
            
    # save the stacked histogram to a file
    output_filename = "stackoutput/stacked_histograms.pdf"
    canvas.SaveAs(output_filename)

    # close all ROOT files
    for file_handle in file_handles:
        file_handle.Close()

    ## Execute the main function
if __name__ == '__main__':
    main()

######################################################################################################################################################################

    '''
    file_DY    = TFile.Open(inputDir + "ttbar3LCR2018_DY_MC_May13_v5.root","READ")
    file_ttbar = TFile.Open(inputDir + "ttbar3LCR2018_ttbar_MC_May13_v5.root","READ")
    file_ZZ    = TFile.Open(inputDir + "ttbar3LCR2018_ZZ_MC_May13_v5.root","READ")
    file_WZ    = TFile.Open(inputDir + "ttbar3LCR2018_WZ_MC_May13_v5.root","READ")
    file_TTZ   = TFile.Open(inputDir + "ttbar3LCR2018_TTZ_MC_May13_v5.root","READ")
    file_TTW   = TFile.Open(inputDir + "ttbar3LCR2018_TTW_MC_May13_v5.root","READ")
    file_Data  = TFile.Open(inputDir + "ttbar3LCR2018_Data_May13_v5.root","READ")
    
    print(" \n file opened in ROOT successfully..")
    
    
    ## Get the histograms

    plotname = "mass12_all"
    h_dy   = file_DY.Get(plotname);    decorate(h_dy, kRed-9);
    h_tt   = file_ttbar.Get(plotname); decorate(h_tt, kGreen-9);
    h_zz   = file_ZZ.Get(plotname);    decorate(h_zz, kMagenta-9);
    h_wz   = file_WZ.Get(plotname);    decorate(h_wz, kOrange);
    h_ttz  = file_TTZ.Get(plotname);   decorate(h_ttz,kCyan-10);
    h_ttw  = file_TTW.Get(plotname);   decorate(h_ttw,kBlue-9);

    h_data = file_Data.Get(plotname);
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(0.6)
    h_data.SetLineColor(kBlack)
    SetOverflowBin(h_data)

    
    print("Histograms are ready...")

    '''

    '''
    ##Scale
    dtlumi = 59.8*1000;
    ttlumi = 28701360/88.29;
    ttzlumi = 13280000/0.2432;
    ttwlumi = 4911941/0.2149;
    dylumi = 99717900/5765.0;  
    wzlumi =10527550/5.052;
    zzlumi =98613000/1.325;
    
    h_dy.Scale(dtlumi/dylumi)
    h_tt.Scale(dtlumi/ttlumi)
    h_wz.Scale(dtlumi/wzlumi)
    h_zz.Scale(dtlumi/zzlumi)
    h_ttw.Scale(dtlumi/ttwlumi)
    h_ttz.Scale(dtlumi/ttzlumi)

    '''
    
    '''
    
    ##Rebining
    rebin = 50
    h_dy.Rebin(rebin)
    h_tt.Rebin(rebin)
    h_wz.Rebin(rebin)
    h_zz.Rebin(rebin)
    h_ttw.Rebin(rebin)
    h_ttz.Rebin(rebin)
    h_data.Rebin(rebin)

    '''

    '''
    ##stack
    h_stack = THStack()
    h_stack.Add(h_ttw)
    h_stack.Add(h_ttz)
    h_stack.Add(h_tt)
    h_stack.Add(h_zz)
    h_stack.Add(h_wz)
    h_stack.Add(h_dy)

    ## RatioHisto (data/allbkg)
    h_bkg = h_dy.Clone()
    h_bkg.Add(h_tt)
    h_bkg.Add(h_wz)
    h_bkg.Add(h_zz)
    h_bkg.Add(h_ttz)
    h_bkg.Add(h_ttw)
    
    h_ratio = h_data.Clone()
    h_ratio.Divide(h_bkg)
    
    print("Scaling, stacking done....time for plotting!")

    ## Legend
    legend = TLegend(0.95,0.50,0.80,0.86)
    ratioleg = TLegend(0.90, 0.90, 0.81,0.87)
    ratioleg.SetHeader(f"obs/exp={h_data.Integral()/h_bkg.Integral():.2f}   exp: {h_bkg.Integral():.0f}")
    legend.AddEntry(h_data, f"Data[{h_data.Integral():.0f}]",'ep')
    legend.AddEntry(h_dy  , f"DY[{h_dy.Integral():.0f}]"    ,'lf')
    legend.AddEntry(h_wz  , f"WZ[{h_wz.Integral():.0f}]"    ,'lf')
    legend.AddEntry(h_tt  , f"TTbar[{h_tt.Integral():.0f}]" ,'lf')
    legend.AddEntry(h_zz  , f"ZZ[{h_zz.Integral():.0f}]"    ,'lf')
    legend.AddEntry(h_ttz , f"TTZ[{h_ttz.Integral():.0f}]"  ,'lf')
    legend.AddEntry(h_ttw , f"TTW[{h_ttw.Integral():.0f}]"  ,'lf')

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
    h_stack.SetMaximum(1e6)
    
    h_stack.Draw("HIST")
    h_data.Draw('ep same')

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
    text1 = DrawText(0.15,0.92,"IISER Pune Analysis")
    text1 = DrawText(0.66,0.92,"13 TeV(2018)")

    ##plot ratiopad
    ratioPad.cd()
    h_ratio.GetXaxis().SetTitle("Invariant Mass(M_{12})")
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
    
    h_ratio.SetTitle('')
    h_ratio.Draw("ep")
    #ratioPad.Update()


    canvas.Draw()
    
    canvas.SaveAs('../output/exampleStackPlot.pdf')
'''

