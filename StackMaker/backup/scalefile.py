# main_script.py
import os
import ROOT
from bkgsamples import bkg_samples

def scale_histograms(input_file, output_file, scale_factor):
    # Open the input root file
    in_file = ROOT.TFile.Open(input_file, "READ")

    # Check if the file is open successfully
    if not in_file or in_file.IsZombie():
        print("Error: Unable to open input file:", input_file)
        return

    # Create the output directory if it doesn't exist
    output_directory = os.path.dirname(output_file)
    os.makedirs(output_directory, exist_ok=True)

    # Create an output root file
    out_file = ROOT.TFile(output_file, "RECREATE")

    # Loop over all histograms in the input file
    for key in in_file.GetListOfKeys():
        obj = key.ReadObj()

        # Check if the object is a histogram
        if obj.IsA().InheritsFrom("TH1"):
            hist = obj.Clone()  # Clone to avoid modifying the original histogram
            hist.Scale(scale_factor)

            # Write the scaled histogram to the output file
            hist.Write()

    # Close the input and output files
    in_file.Close()
    out_file.Close()

def process_samples(samples, input_dir, output_dir, scale_factor):
    for sample_group, sample_info in samples.items():
        for sub_sample_name, sub_sample_info in sample_info.items():
            #print(f"Processing: {sample_group} - {sub_sample_name}")
            input_filename = os.path.join(input_dir, sub_sample_info["filename"])
            output_filename = os.path.join(output_dir, f"scaled_{sub_sample_info['filename']}")
            #print(f"Input File: {input_filename}")
            #print(f"Output File: {output_filename}")
            scale_histograms(input_filename, output_filename, scale_factor)

if __name__ == "__main__":
    scale_factor = 2.0  # Change this to your desired scaling factor
    inputDir = "../cluster_hst_output/nov09/"
    outputDir = "finalhaddOutput/"

    # Process all samples
    process_samples(bkg_samples, inputDir, outputDir, scale_factor)
