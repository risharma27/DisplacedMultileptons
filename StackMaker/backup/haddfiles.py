import os
import ROOT
import subprocess
from bkgsamples import bkg_samples

def hadd_files(output_dir, set_name, input_files):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Construct the hadd command
    output_file = os.path.join(output_dir, f"{set_name}_merged.root")
    hadd_command = ["hadd", output_file] + input_files

    try:
        # Run the hadd command and capture the output and errors
        result = subprocess.run(hadd_command, capture_output=True, text=True, check=True)
        
        # Print the output for debugging
        print("hadd command output:", result.stdout)
    except subprocess.CalledProcessError as e:
        # Print the error if the command failed
        print("Error running hadd command:", e.stderr)

def process_samples(samples, input_dir, output_dir, data_lumi):
    for sample_group, sample_info in samples.items():
        scaled_files = []  # List to store paths of scaled files for hadding

        for sub_sample_name, sub_sample_info in sample_info.items():
            input_filename = os.path.join(input_dir, sub_sample_info["filename"])
            
            # Calculate scale factor for each sub-sample
            scale_factor = data_lumi / (sub_sample_info["nevents"] * sub_sample_info["xsec"])

            # Create a temporary scaled histogram
            temp_hist = ROOT.TH1F("temp_hist", "", 1, 0, 1)

            # Scale the histogram and add it to the temporary histogram
            temp_hist.Add(ROOT.TFile.Open(input_filename).Get("your_histogram_name"), scale_factor)

            # Append the temporary histogram to the list
            scaled_files.append(temp_hist)

        # Merge the scaled histograms into one final histogram for each sample
        hadd_files(output_dir, sample_group, scaled_files)

if __name__ == "__main__":
    inputDir = "../cluster_hst_output/nov09/"
    outputDir = "finalhaddOutput/"
    data_lumi = 1000.0  # Replace with your actual data luminosity

    # Process all samples and hadd files
    process_samples(bkg_samples, inputDir, outputDir, data_lumi)
