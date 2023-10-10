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
    inputDir = "cluster_hst_output/"    
    file_Data1  = TFile.Open(inputDir + "Data2016_SingleElectron.root","READ")
    file_Data2  = TFile.Open(inputDir + "Data2016_SingleMuon.root","READ")
    print(" \n file opened in ROOT successfully..")
    
    
    ## Get the histograms

    plotname = "2l1d_imass_3l"
    
    '''
    h_dy   = file_DY.Get(plotname);    decorate(h_dy, kRed-9);
    h_tt   = file_ttbar.Get(plotname); decorate(h_tt, kGreen-9);
    h_zz   = file_ZZ.Get(plotname);    decorate(h_zz, kMagenta-9);
    h_wz   = file_WZ.Get(plotname);    decorate(h_wz, kOrange);
    h_ttz  = file_TTZ.Get(plotname);   decorate(h_ttz,kCyan-10);
    h_ttw  = file_TTW.Get(plotname);   decorate(h_ttw,kBlue-9);
    '''

    h_data1 = file_Data1.Get(plotname); decorate(h_data1, kRed-9);
    h_data2 = file_Data2.Get(plotname); decorate(h_data2, kGreen+9);
    
    print("Histograms are ready...")

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
    
    ##stack
    h_stack = THStack()
    h_stack.Add(h_data1)
    h_stack.Add(h_data2)
    
    print("Scaling, stacking done....time for plotting!")

    ## Legend
    legend = TLegend(0.95,0.50,0.80,0.86)
    legend.AddEntry(h_data1, f"Data[{h_data1.Integral():.0f}]",'lf')
    legend.AddEntry(h_data2, f"Data[{h_data2.Integral():.0f}]",'lf')    
    
    SetLegendStyle(legend)

    ######################################################################################
    ##      PLOTTING START
    ######################################################################################
    canvas = TCanvas("c","canvas",650,600)
    gStyle.SetOptStat(0)

    mainPad  = TPad("pad","pad",0,0.3,1,1)
    mainPad.Draw()
    mainPad.cd()
    mainPad.SetLogy(1)
    h_stack.SetMinimum(0.001)
    h_stack.SetMaximum(1e6)
    
    h_stack.Draw("HIST")
    #h_data.Draw('ep same')
    
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
    text1 = DrawText(0.66,0.92,"13 TeV(2016)")


    canvas.Draw()
    
    canvas.SaveAs('stackoutput/exampleStackPlot.pdf')


## Execute the main function
if __name__ == '__main__':
    main()
