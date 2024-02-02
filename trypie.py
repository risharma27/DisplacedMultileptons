import ROOT

def create_pie_chart(file_names):
    # Create a TPie object
    pie = ROOT.TPie("invmass_pie", "Invariant Mass Pie Chart", len(file_names))

    # Create a legend
    legend = ROOT.TLegend(0.7, 0.2, 0.9, 0.5)
    legend.SetHeader("Background Percentage")

    # Assign colors for each file
    file_colors = {"DY.root": 2, "TTBar.root": 3, "WJets.root": 4}

    total_entries = 0

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
        entry_label = file_name.replace(".root", "")  # Use file name without extension as the label
        entry_value = invmass_hist.GetEntries()
        pie.SetEntryLabel(idx, entry_label)
        pie.SetEntryVal(idx, entry_value)
        pie.SetEntryFillColor(idx, file_colors.get(file_name, 1))  # Default color is 1 (black)

        # Add entry to the legend
        legend.AddEntry(entry_label, f"{entry_label}: {entry_value:.2f}")

        total_entries += entry_value

        # Close the file
        file.Close()

    # Calculate total percentage using the sum of all slice values
    total_percentage = sum(pie.GetSum())

    # Create a canvas to draw the pie chart
    canvas = ROOT.TCanvas("canvas", "Pie Chart Canvas", 800, 600)

    # Draw the pie chart
    pie.Draw("nol")

    # Draw the legend
    legend.Draw()

    # Set the radius, labels offset, and label format for the last slice
    pie.SetRadius(.2)
    pie.SetLabelsOffset(.01)
    pie.SetLabelFormat("#splitline{%val (%perc)}{%txt}")

    # Display the canvas
    canvas.Draw()

# Example usage with three files
file_names = ["DY.root", "TTBar.root", "WJets.root"]
create_pie_chart(file_names)
