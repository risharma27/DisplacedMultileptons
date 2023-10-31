import ROOT
from ROOT import *
import sys
import os
import warnings
warnings.filterwarnings('ignore')

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

def decorate_hstack(h_stack):
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


def decorate_hratio(h_ratio):
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

    # Set ROOT's verbosity level to suppress info and warnings
    gROOT.SetBatch(True)
    gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")  # suppress info messages
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")  # suppress warning messages    

    inputDir = "../cluster_hst_output/"

    MCfile_names = [
        
        "DYJetsToLL_oct20.root",
        "TTTo2L2Nu_oct20.root",
        "WZTo3LNu_oct20.root",
        
    ]  # this is a list of names of all MC samples root files

    # open all the ROOT files
    files_mc = [TFile.Open(inputDir + file_name, "READ") for file_name in MCfile_names] #storing the MC root files in a list
    file_data = TFile.Open(inputDir + "2016Data_oct20.root", "READ")

    print("\nFiles opened in ROOT successfully..")
    
    # get a list of histogram names from one of the files (since all the root files have the same set of histograms)
    '''
    mc_file = files_mc[0]
    list_of_histograms = mc_file.GetListOfKeys()
    hist_names = [hist.GetName() for hist in list_of_histograms] #storing all the histogram names (plotnames) in a list

    '''

    hist_names = [ "flavor", "met", "imass_3l", "delR_l0l1", "delPhi_l0l1", "delPhi_l0met", "imass_l0l1", "mt0", "delR_l1l2", "delPhi_l1l2", "delPhi_l1met", "imass_l1l2", "mt1", "delR_l2l0", "delPhi_l2l0", "delPhi_l2met", "imass_l2l0", "mt2" ]
    hist_prefix = ["2l1d_", "1l2d_", "3d_"]

    hists = []

    for prefix in hist_prefix:
        for hist in hist_names:
            hist_name = prefix + hist
            hists.append(hist_name)
    
    hist_colors = [kBlue-9, kGreen-9, kRed-9]
    
    #######################
    # Luminosity Scaling #
    #######################
    
    dtlumi = 36.3*1000;
    dylumi = 616760700/5765.0
    ttlumi = 243094700/88.29
    wzlumi = 60224600/5.052

    lumi_sf = [dtlumi/dylumi, dtlumi/ttlumi, dtlumi/wzlumi]
   
    #########################
     # Stacking Histograms # 
    #########################
    
    # Loop through each histogram name
    for plotname in hists:
        hst_data = file_data.Get(plotname)
        hst_data.SetMarkerStyle(20)
        hst_data.SetMarkerSize(0.6)
        hst_data.SetLineColor(kBlack)
        #if plotname != "nEvents" and plotname != "nEvSel":
        if "flavor" in plotname:
            rebin = 1
        else:
            rebin = 5
        hst_data.Rebin(rebin)
        SetOverflowBin(hst_data)
        
        # initialize a list to store the same histograms from all files
        histograms = []
        hstMC_integral = []
        
        # loop through all the files and retrieve histograms with the same name
        for i, file_mc in enumerate(files_mc):
            hst_MC = file_mc.Get(plotname)
            if hst_MC:
                hst_MC.Scale(lumi_sf[i])
                hstMC_integral.append(hst_MC.Integral())
                color = hist_colors[i % len(hist_colors)] #filling each histogram with a different color
                decorate(hst_MC, color)  
                hst_MC.Rebin(rebin)
                SetOverflowBin(hst_MC)
                histograms.append(hst_MC)
                
        if histograms:
            #setting the xtitle 
            hst_title = plotname.split('_',1)
            if len(hst_title)>1:
             xtitle = hst_title[1]
            else:
             xtitle = plotname
            
            # Create a THStack to stack the histograms
            hst_stack = THStack()
            hst_stack.SetTitle(plotname)
          
            # Add histograms with the same name to the THStack
            for histogram in histograms:  
                hst_stack.Add(histogram)
                
            ## RatioHisto (data/allbkg)
            hst_bkg = histograms[0].Clone()
            for hbkg in range(1,len(histograms)):
                hst_bkg.Add(histograms[hbkg])
                
            hst_ratio = hst_data.Clone()
            hst_ratio.Divide(hst_bkg)

            #legend
            legend   = TLegend(0.95,0.50,0.80,0.86)
            ratioleg = TLegend(0.90,0.90,0.81,0.87)
            ratioleg.SetHeader(f"obs/exp={hst_data.Integral()/hst_bkg.Integral():.2f}   exp: {hst_bkg.Integral():.2f}")
            legend.AddEntry(hst_data,f"Data[{hst_data.Integral():.0f}]",'ep')
            for i, file_name in enumerate(MCfile_names):
                legend.AddEntry(histograms[i], f"{file_name}[{hstMC_integral[i]:.0f}]", "lf")

            SetLegendStyle(ratioleg)
            SetLegendStyle(legend)

            ###########################################################

            #                 PLOTTING START                          
            
            ###########################################################
            
            # Create a TCanvas to display the stacked histogram
            canvas = TCanvas("c", "canvas", 750, 600)
            gStyle.SetOptStat(0)

            ratioPadSize = 0.3
            mainPad  = TPad("pad","pad",0,ratioPadSize,1,1)
            ratioPad = TPad("pad2","pad2",0,0,1.0,ratioPadSize)
            PadStyling(mainPad,ratioPad)
            mainPad.Draw()
            ratioPad.Draw()

            mainPad.cd()    #ensures that the subsequent drawing
                            #happens in the mainPad
            mainPad.SetLogy(1)

            # Draw the stacked histogram
            hst_stack.SetMinimum(0.001)
            hst_stack.SetMaximum(1e6)
            hst_stack.Draw("HIST")
            hst_stack.GetYaxis().SetTitle('Events')
            hst_stack.GetYaxis().CenterTitle()
            #h_stack.GetXaxis().SetTitle(xtitle)
            #h_stack.GetXaxis().CenterTitle()
            hst_data.Draw('ep same')

            decorate_hstack(hst_stack)

            mainPad.SetTickx(1)
            
            legend.Draw()
            ratioleg.Draw()

            #Plotting ratioPad

            decorate_hratio(hst_ratio)
            ratioPad.cd()  #subsequent plotting will happen in the ratioPad
            hst_ratio.GetXaxis().SetTitle(xtitle)
            hst_ratio.GetYaxis().SetTitle("obs/exp")
            hst_ratio.GetXaxis().CenterTitle()
            hst_ratio.GetYaxis().CenterTitle()

            hst_ratio.SetTitle('')
            hst_ratio.Draw("ep")

            # Display the histogram
            canvas.Draw()
           
            # Saving the stacked histograms
            if plotname.startswith("2l1d_"):
                output_filename = f"stackoutput/evsel_2l1d/{plotname}.png"
                canvas.SaveAs(output_filename)
            elif plotname.startswith("1l2d_"):
                output_filename = f"stackoutput/evsel_1l2d/{plotname}.png"
                canvas.SaveAs(output_filename)
            elif plotname.startswith("3d_"):
                output_filename = f"stackoutput/evsel_3d/{plotname}.png"
                canvas.SaveAs(output_filename)
            else:
                output_filename = f"stackoutput/{plotname}.png"
                canvas.SaveAs(output_filename)
            

    # Close all ROOT files
    file_data.Close()
    for file_mc in files_mc:
        file_mc.Close()

if __name__ == '__main__':
    main()

