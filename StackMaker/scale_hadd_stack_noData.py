import os
import sys
import ROOT
from ROOT import *
import subprocess
import warnings
warnings.filterwarnings('ignore')

from bkgsamples import bkg_samples

####################################################################################
#                    User Defined functions                                        #
####################################################################################

#function to scale the histograms to data luminosity
def scale_histograms(input_file, output_file, scale_factor):
   
    in_file = TFile.Open(input_file, "READ")
    if not in_file or in_file.IsZombie():
        print("Error: Unable to open input file:", input_file)
        return

    output_directory = os.path.dirname(output_file)
    os.makedirs(output_directory, exist_ok=True)

    out_file = TFile(output_file, "RECREATE") # output root file which will store the scaled histograms

    # loop over all histograms in the input file
    for key in in_file.GetListOfKeys():
        obj = key.ReadObj()
        if obj.IsA().InheritsFrom("TH1"):  #check if the object is a histogram
            hist = obj.Clone()  # clone to avoid modifying the original histogram
            hist.Scale(scale_factor)
            hist.Write()  # write the scaled histogram to the output file

    in_file.Close()
    out_file.Close()


#function to hadd the all files of data and every bkg
#the output files with the scaled histograms are hadded
def hadd_files(output_dir, setname, input_files): 
    os.makedirs(output_dir, exist_ok=True) 
    output_file = os.path.join(output_dir, f"{setname}.root")
    hadd_command = ["hadd", "-f",  output_file] + input_files

    subprocess.call(hadd_command)

    
#this functions uses the above two functions
#each bkg sample has subdivisions for different pT/HT/mass bins with different xsec, this function scales all of them to their respective value
#the info is taken from the bkg_samples dictionary
#and final output is one file for every bkg and data
def process_samples(samples, input_dir, output_dir, datalumi_pre, datalumi_post):
    for sample_group, sample_info in samples.items():
        scaled_files = []  
        for sub_sample_name, sub_sample_info in sample_info.items():
            input_filename = os.path.join(input_dir, sub_sample_info["filename"])
            output_filename = os.path.join(output_dir, f"scaled_{sub_sample_info['filename']}")
            
            if sub_sample_info["data"]==0:
                if "preVFP" in sub_sample_info['filename']:
                    print(sub_sample_info['filename'])
                    scale_factor = datalumi_pre / (sub_sample_info["nevents"] / sub_sample_info["xsec"])  # calculating the scale factor
                elif "postVFP" in sub_sample_info['filename']:
                    print(sub_sample_info['filename'])
                    scale_factor = datalumi_post / (sub_sample_info["nevents"] / sub_sample_info["xsec"])  
            elif sub_sample_info["data"]==1:
                scale_factor = 1
                    
            scale_histograms(input_filename, output_filename, scale_factor)          
            scaled_files.append(output_filename)

        # merge the scaled files into one final file for each sample
        hadd_files(output_dir, sample_group, scaled_files)

        
#following functions are for decoration
def SetOverflowBin(histo):
    nbins = histo.GetNbinsX()
    histo.SetBinContent(nbins, histo.GetBinContent(nbins) + histo.GetBinContent(nbins+1)); ## Overflow
    histo.SetBinContent(1, histo.GetBinContent(1)+ histo.GetBinContent(0));                ## Underflow

def DrawText(X,Y,txt,style):
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextFont(style)
    text.SetTextSize(0.04)
    text.DrawLatex(X,Y,txt)

    return text
    
def decorate(h,color):
    h.SetLineColor(kGray+3)
    h.SetFillColor(color)
    SetOverflowBin(h) ## overflow bin is must

def decorate_hstack(h_stack):
    h_stack.GetXaxis().SetTitleFont(43)
    h_stack.GetXaxis().SetTitleSize(20)
    h_stack.GetXaxis().SetTitleOffset(0.8)
    h_stack.GetXaxis().SetLabelFont(43)
    h_stack.GetXaxis().SetLabelSize(12)
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

'''
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
'''

def PadStyling(pad):
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.20)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.1)
    pad.SetTickx(1)
    pad.SetTicky(1)


def SetLegendStyle(legend):
    legend.SetTextFont(62)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.027)    
           
            
####################################################################################                                                                             
####################################################################################



###########################
# main function begins here
###########################

def main():
    # Set ROOT's verbosity level to suppress info and warnings
    gROOT.SetBatch(True)
    gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")  # suppress info messages
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")  # suppress warning messages

    inputDir = "../cluster_hst_output/mar15/"
    outputDir = "finalhaddOutput/mar15/"
    stackDir = "stackoutput/mar15/"
    
    preVFP_lumi  = 19.3 * 1000
    postVFP_lumi = 17 * 1000

    # process all samples, scale histograms, and hadd files
    process_samples(bkg_samples, inputDir, outputDir, preVFP_lumi, postVFP_lumi)

    scaled_files_to_delete = [f for f in os.listdir(outputDir) if f.startswith("scaled")]
    for file_to_delete in scaled_files_to_delete:
        os.remove(os.path.join(outputDir, file_to_delete))
    
        
    ###########################
     # Stacking begins
    ###########################

    
    MC_files = ["DY.root","TTBar.root","WJets.root","QCD.root","WGamma.root","ZGamma.root","ZZ.root","WZ.root","WW.root","SingleTop.root"]
    #MC_files = ["DY.root"]
    #hist_colors = [kBlue-7,kTeal+9,kOrange+1,kYellow-7,kMagenta-10,kViolet+1,kAzure,kMagenta-7,kPink+1,kCyan-3]
    #hist_colors = [kBlue-7,kTeal+9,kOrange+1,kYellow-7,kViolet+6,kAzure-4,kAzure,kMagenta-7,kPink+1,kCyan-3]
    hist_colors = [kBlue-7,kGreen-2,kOrange+1,kYellow-7,kRed-10,kAzure-9,kAzure-3,kViolet-9,kRed-9,kTeal-8]
    files_mc = [TFile.Open(outputDir + file_name, "READ") for file_name in MC_files] #storing the MC root files in a list
    file_data = TFile.Open(outputDir + "Data.root", "READ")

    print("\nFiles opened in ROOT successfully..") 
    
    hist_names = ["flavor", "met", "l0_dxy", "l0_ip3d", "l0_sip3d", "l1_dxy", "l1_ip3d", "l1_sip3d", "l2_dxy", "l2_ip3d", "l2_sip3d", "imass_3l", "delR_l0l1", "delPhi_l0l1", "delPhi_l0met", "imass_l0l1", "mt0", "delR_l1l2", "delPhi_l1l2", "delPhi_l1met", "imass_l1l2", "mt1", "delR_l2l0", "delPhi_l2l0", "delPhi_l2met", "imass_l2l0", "mt2", "njet", "bjets", "l0_reliso03", "l1_reliso03", "l2_reliso03", "l2_dxy", "l2_dz", "l2_ip3d", "l2_sip3d", "l0_|dxy|", "l1_|dxy|", "l2_|dxy|"]

    #hist_2LSS = ["2LSS_flavor", "2LSS_invmass_ll", "2LSS_met", "mumu_invmass_ll", "ee_invmass_ll", "emu_or_mue_invmass_ll", "2LSS_ht"]

    #hist_3L = ["3L_invmass_3l", "3L_invmass_l0l1", "3L_met", "3L_lt", "3L_njet", "3L_ht"]

    hist_cr_ttbar = ["cr_ttbar_met", "cr_ttbar_pt0", "cr_ttbar_pt1", "cr_ttbar_lt", "cr_ttbar_njet", "cr_ttbar_bjet", "cr_ttbar_ht", "cr_ttbar_mt0", "cr_ttbar_mt1","cr_ttbar_l0iso", "cr_ttbar_l1iso", "cr_ttbar_dphil0l1", "cr_ttbar_invmassl0l1", "cr_ttbar_invmassj0j1"]

    hist_cr_ttbar_2l1d = ["cr_ttbar_2l1d_met", "cr_ttbar_2l1d_pt0", "cr_ttbar_2l1d_pt1", "cr_ttbar_2l1d_lt", "cr_ttbar_2l1d_njet", "cr_ttbar_2l1d_bjet", "cr_ttbar_2l1d_ht", "cr_ttbar_2l1d_mt0", "cr_ttbar_2l1d_mt1","cr_ttbar_2l1d_l0iso", "cr_ttbar_2l1d_l1iso", "cr_ttbar_2l1d_dphil0l1", "cr_ttbar_2l1d_invmassl0l1", "cr_ttbar_2l1d_invmassj0j1"]

    hist_vr_ttbar = ["vr_ttbar_met", "vr_ttbar_pt0", "vr_ttbar_pt1", "vr_ttbar_lt", "vr_ttbar_njet", "vr_ttbar_bjet", "vr_ttbar_ht", "vr_ttbar_mt0", "vr_ttbar_mt1","vr_ttbar_l0iso", "vr_ttbar_l1iso", "vr_ttbar_dphil0l1", "vr_ttbar_invmassl0l1", "vr_ttbar_invmassj0j1"]

    hist_vr_ttbar_2l1d = ["vr_ttbar_2l1d_met", "vr_ttbar_2l1d_pt0", "vr_ttbar_2l1d_pt1", "vr_ttbar_2l1d_lt", "vr_ttbar_2l1d_njet", "vr_ttbar_2l1d_bjet", "vr_ttbar_2l1d_ht", "vr_ttbar_2l1d_mt0", "vr_ttbar_2l1d_mt1","vr_ttbar_2l1d_l0iso", "vr_ttbar_2l1d_l1iso", "vr_ttbar_2l1d_dphil0l1", "vr_ttbar_2l1d_invmassl0l1", "vr_ttbar_2l1d_invmassj0j1"]

    hist_cr_wjets = ["cr_wjets_met", "cr_wjets_pt0", "cr_wjets_lt", "cr_wjets_njets", "cr_wjets_ht", "cr_wjets_mt0", "cr_wjets_l0iso", "cr_wjets_dphil0j0", "cr_wjets_j0pt"]

    hist_cr_wjets_1l2d = ["cr_wjets_1l2d_met", "cr_wjets_1l2d_pt0", "cr_wjets_1l2d_lt", "cr_wjets_1l2d_njets", "cr_wjets_1l2d_ht", "cr_wjets_1l2d_mt0", "cr_wjets_1l2d_l0iso", "cr_wjets_1l2d_dphil0j0", "cr_wjets_1l2d_j0pt"]

    hist_prefix = ["2l1d_", "1l2d_", "3d_"]    
    hists = []   
    for prefix in hist_prefix:
        for hist in hist_names:
            hist_name = prefix + hist
            hists.append(hist_name)

    '''
    for hist in hist_2LSS:
        hists.append(hist)

    for hist in hist_3L:
        hists.append(hist)

    for hist in hist_cr_ttbar:
        hists.append(hist)

    for hist in hist_vr_ttbar:
        hists.append(hist)

    for hist in hist_cr_wjets:
        hists.append(hist)   

    for hist in hist_cr_ttbar_2l1d:
        hists.append(hist)

    for hist in hist_vr_ttbar_2l1d:
        hists.append(hist)  

    for hist in hist_cr_wjets_1l2d:
        hists.append(hist)
    '''

    for plotname in hists:
        hst_data = file_data.Get(plotname)
        hst_data.SetMarkerStyle(20)
        hst_data.SetMarkerSize(0.6)
        hst_data.SetLineColor(kBlack)
        #if plotname != "nEvents" and plotname != "nEvSel":
        if "flavor" in plotname or "njet" in plotname:
            rebin = 1
        elif "dxy" in plotname or "dz" in plotname:
            rebin = 10
        elif "ip3d" in plotname:
            rebin = 5
        else:
            rebin = 25
        hst_data.Rebin(rebin)
        SetOverflowBin(hst_data)
          
        # initialize a list to store the same histograms from all files
        histograms = []
        hstMC_integral = []

        for i, file_mc in enumerate(files_mc):
            hst_MC = file_mc.Get(plotname)
            if hst_MC:
               hstMC_integral.append(hst_MC.Integral())
               color = hist_colors[i] 
               decorate(hst_MC, color)  
               hst_MC.Rebin(rebin)
               SetOverflowBin(hst_MC)
               histograms.append(hst_MC)

        sorted_histograms = sorted(histograms, key=lambda x: x.Integral())
                
        if histograms:
            #setting the xtitle
            
            hst_title = plotname.split('_',1)
            if len(hst_title)>1:
               xtitle = hst_title[1]
            else:
               xtitle = plotname
            
            # Create a THStack to stack the histograms
            hst_stack = THStack()
            #hst_stack.SetTitle(plotname)
            hst_stack.SetTitle(" ")
          
            # Add histograms with the same name to the THStack
            for histogram in sorted_histograms:  
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
            #print(hst_bkg , hst_bkg.Integral())
            ratioleg.SetHeader(f"obs/exp={hst_data.Integral()/hst_bkg.Integral():.5f} | exp: {hst_bkg.Integral():.0f}")

            # create pairs of samples and corresponding integrals
            file_integral_pairs = zip(MC_files, hstMC_integral)
            
            # sort based on the integral values (using the second element of each pair)
            sorted_file_integral_pairs = sorted(file_integral_pairs, key=lambda x: x[1], reverse=True)
            
            #legend.AddEntry(hst_data,f"Data[{hst_data.Integral():.0f}]",'ep')
            
            for filename, integral in sorted_file_integral_pairs:
                fileleg = filename.split('.',1)[0]
                legend.AddEntry(histograms[MC_files.index(filename)], f"{fileleg}[{integral:.0f}]", "f")
                

            SetLegendStyle(ratioleg)
            SetLegendStyle(legend)

            
            
            ###########################################################

            #                 PLOTTING START                          
            
            ###########################################################
            
            # create a TCanvas to display the stacked plots
            canvas = TCanvas("c", "canvas", 800, 600)
            gStyle.SetOptStat(0)

            ratioPadSize = 0
            mainPad  = TPad("pad","pad",0,ratioPadSize,1,1)
            #ratioPad = TPad("pad2","pad2",0,0,1.0,ratioPadSize)
            PadStyling(mainPad)
            mainPad.Draw()
            #ratioPad.Draw()

           

            # draw the stacked histogram
            mainPad.cd()    #ensures that the subsequent drawing
                            #happens in the mainPad
            mainPad.SetLogy(1)
            
            hst_stack.Draw("HIST")
              
            decorate_hstack(hst_stack)
      
            #hst_stack.GetXaxis().SetTitle(xtitle)
            #hst_stack.GetXaxis().CenterTitle()          
            hst_stack.GetYaxis().SetTitle('Events')
            hst_stack.GetYaxis().CenterTitle()
            hst_stack.GetXaxis().SetTitle(xtitle)
            hst_stack.GetXaxis().CenterTitle()
            #hst_stack.GetXaxis().SetRangeUser(0.05, 10)
            hst_stack.SetMinimum(0.001)
            hst_stack.SetMaximum(1e5)
            #hst_data.Draw('ep same')        

            mainPad.SetTickx(1)
            
            legend.Draw()
            #ratioleg.Draw()

            ## A few text
            #text1 = DrawText(0.15,0.92,"IISER Pune Analysis")
            text2 = DrawText(0.64,0.92,"13 TeV(2016)",52)

            DrawText(0.15,0.92,"CMS",62)
            DrawText(0.22,0.92,"preliminary",52)

            if "2l1d_" in plotname:
                DrawText(0.45,0.92,"2L1D",62)

            if "1l2d_" in plotname:
                DrawText(0.45,0.92,"1L2D",62)

            if "3d_" in plotname:
                DrawText(0.45,0.92,"0L3D",62)

            #plotting ratioPad

            '''
            decorate_hratio(hst_ratio)
           
            ratioPad.cd()  #subsequent plotting will happen in the ratioPad
            hst_ratio.GetXaxis().SetTitle(xtitle)
            hst_ratio.GetYaxis().SetTitle("obs/exp")
            hst_ratio.GetXaxis().CenterTitle()
            hst_ratio.GetYaxis().CenterTitle()

            hst_ratio.SetTitle('')
            hst_ratio.Draw("ep")
            '''

            # display the plot
            canvas.Draw()

            
            # saving the stacked plots
            if plotname.startswith("2l1d_"):
                output_filename = os.path.join(stackDir, "evsel_2l1d", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            elif plotname.startswith("1l2d_"):
                output_filename = os.path.join(stackDir, "evsel_1l2d", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            elif plotname.startswith("3d_"):
                output_filename = os.path.join(stackDir, "evsel_3d", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            elif plotname.startswith("cr_ttbar_"):
                output_filename = os.path.join(stackDir, "cr_ttbar", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            elif plotname.startswith("vr_ttbar_"):
                output_filename = os.path.join(stackDir, "vr_ttbar", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            elif plotname.startswith("cr_wjets_"):
                output_filename = os.path.join(stackDir, "cr_wjets", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            elif plotname.startswith("3L_"):
                output_filename = os.path.join(stackDir, "3L", f"{plotname}.png")
                canvas.SaveAs(output_filename)
            else:
                output_filename = os.path.join(stackDir, f"{plotname}.png")
                canvas.SaveAs(output_filename)

                
    print("Stacking Complete")
    
    # now close all root files
    file_data.Close()
    for file_mc in files_mc:
        file_mc.Close()
        

        
#########################           
# main function ends here
#########################



    
if __name__ == "__main__":
    main()

        
