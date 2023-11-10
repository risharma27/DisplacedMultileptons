import os
import ROOT
import subprocess
from bkgsamples import bkg_samples

####################################################################################

def scale_histograms(input_file, output_file, scale_factor):
    # Open the input root file
    in_file = ROOT.TFile.Open(input_file, "READ")
    if not in_file or in_file.IsZombie():
        print("Error: Unable to open input file:", input_file)
        return

    # Create the output directory if it doesn't exist
    output_directory = os.path.dirname(output_file)
    os.makedirs(output_directory, exist_ok=True)

    # Create an output root file
    out_file = ROOT.TFile(output_file, "RECREATE")

    # loop over all histograms in the input file
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
    output_file = os.path.join(output_dir, f"{set_name}.root")
    hadd_command = ["hadd", "-f",  output_file] + input_files

    subprocess.call(hadd_command)

def process_samples(samples, input_dir, output_dir, data_lumi):
    for sample_group, sample_info in samples.items():
        scaled_files = []  # List to store paths of scaled files for hadding

        for sub_sample_name, sub_sample_info in sample_info.items():
            input_filename = os.path.join(input_dir, sub_sample_info["filename"])
            output_filename = os.path.join(output_dir, f"scaled_{sub_sample_info['filename']}")
            
            # Calculate scale factor for each sub-sample
            scale_factor = data_lumi / (sub_sample_info["nevents"] / sub_sample_info["xsec"])

            # Scale histograms
            scale_histograms(input_filename, output_filename, scale_factor)

            # Append the scaled file path to the list
            scaled_files.append(output_filename)

        # Merge the scaled files into one final file for each sample
        hadd_files(output_dir, sample_group, scaled_files)

       
            
####################################################################################


###########################
# main function begins here
###########################
        
if __name__ == "__main__":
    inputDir = "../cluster_hst_output/nov09/"
    outputDir = "finalhaddOutput/"
    data_lumi = 36.3*1000

    # Process all samples, scale histograms, and hadd files
    process_samples(bkg_samples, inputDir, outputDir, data_lumi)

    scaled_files_to_delete = [f for f in os.listdir(outputDir) if f.startswith("scaled")]
    for file_to_delete in scaled_files_to_delete:
        os.remove(os.path.join(outputDir, file_to_delete))
