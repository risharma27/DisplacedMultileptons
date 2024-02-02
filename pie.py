import ROOT

def create_pie_chart(file_names):
    # Create a TPie object
    pie = ROOT.TPie("invmass_pie", "Invariant Mass Pie Chart", len(file_names))

    # Create a legend
    #legend = ROOT.TLegend(0.7, 0.2, 0.9, 0.5)
    #legend.SetHeader("File Entries")

    # Assign colors for each file
    file_colors = { "DY.root": 2, "TTBar.root": 3, "WJets.root": 4 }

    # Loop through each file
    for idx, file_name in enumerate(file_names):
        # Construct the full path to the file
        full_path = "StackMaker/finalhaddOutput/jan18/" + file_name

        # Open the root file
        file = ROOT.TFile(full_path, "READ")

        # Check if the file is open
        if not file or file.IsZombie():
            print(f"Error opening file: {full_path}")
            continue

        # Get the histogram from the file
        invmass_hist = file.Get("2l1d_imass_3l")

        # Check if the histogram is found
        if not invmass_hist:
            print(f"Error: Could not find 'invmass' histogram in file: {full_path}")
            file.Close()
            continue

        # Add the entry to the pie chart
        pie.SetEntryLabel(idx, invmass_hist.GetTitle())
        pie.SetEntryVal(idx, invmass_hist.GetEntries())
        pie.SetEntryFillColor(idx, file_colors.get(file_name, 1))  # Default color is 1 (black)

        # Add entry to the legend
        #legend.AddEntry(invmass_hist.GetTitle(), f"{invmass_hist.GetTitle()}: {invmass_hist.GetEntries()}")

        # Close the file
        file.Close()

    # Create a canvas to draw the pie chart
    canvas = ROOT.TCanvas("canvas", "Pie Chart Canvas", 800, 600)

    # Draw the pie chart
    pie.Draw("nol")

    # Draw the legend
    #legend.Draw()

    # Display the canvas
    canvas.Draw()

    canvas.SaveAs("output_pie_chart.pdf")

# Example usage with three files
file_names = ["DY.root", "TTBar.root", "WJets.root"]
create_pie_chart(file_names)
