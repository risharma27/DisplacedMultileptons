import os
import ROOT
import subprocess
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

def hadd_files(output_dir, set_name, input_files):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Construct the hadd command
    output_file = os.path.join(output_dir, f"{set_name}_merged.root")
    hadd_command = ["hadd", output_file] + input_files

    # Print the hadd command for debugging
    print("Running hadd command:", " ".join(hadd_command))

    # Run the hadd command
    subprocess.run(hadd_command)

def process_samples(samples, input_dir, output_dir, scale_factor):
    for sample_group, sample_info in samples.items():
        for sub_sample_name, sub_sample_info in sample_info.items():
            input_filename = os.path.join(input_dir, sub_sample_info["filename"])
            output_filename = os.path.join(output_dir, f"scaled_{sub_sample_info['filename']}")
            scale_histograms(input_filename, output_filename, scale_factor)

def merge_samples(samples, input_dir, output_dir):
    for sample_group, sample_info in samples.items():
        input_files = [
            os.path.join(input_dir, f"scaled_{sub_sample_info['filename']}")
            for sub_sample_info in sample_info.values()
        ]
        hadd_files(output_dir, sample_group, input_files)

if __name__ == "__main__":
    scale_factor = 2.0  # Change this to your desired scaling factor
    inputDir = "../cluster_hst_output/nov09/"
    outputDir = "finalhaddOutput/"

    # Process all samples and scale histograms
    process_samples(bkg_samples, inputDir, outputDir, scale_factor)

    # Merge the scaled files for each sample
    merge_samples(bkg_samples, outputDir, "finalMergedOutput/")
