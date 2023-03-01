import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class Steady_state_analysis:

    # Constructor to initialize the instance variables
    def __init__(self, bin_file_id):

        # Store the file ID of the CSV file containing the data
        self.bin_file_id = bin_file_id

        # Read the CSV file as a pandas DataFrame
        df = pd.read_csv(self.bin_file_id,
                         delimiter=',',
                         index_col=0,
                         header=None,
                         skiprows=1,
                         names=['Molecule_name', 'File_ID'])

        # Convert the DataFrame into a dictionary where keys are molecule names and values are file IDs
        self.file_id = df.squeeze('columns').to_dict()

    # Method to plot the photoluminescence spectra for all molecules
    def plot_PL(self):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, file_id in self.file_id.items():
            # Load the PL data from the file using numpy
            PLdata = np.loadtxt(file_id, skiprows=1)

            # Extract the wavelength and intensity values from the PL data
            wavelength = PLdata[:, 0]
            Intensity = PLdata[:, 1]

            # Normalize the intensity values by dividing them by the maximum intensity
            NormIntensity = Intensity / max(Intensity)

            # Plot the normalized PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, NormIntensity, label=molecule_name)

        # Add a legend to the plot with font size of 12
        ax.legend(fontsize=12)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')

        # Show the plot
        plt.show()

    # Method to plot the smoothed photoluminescence spectra for all molecules
    def smooth_plot_PL(self):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, file_id in self.file_id.items():
            # Load the PL data from the file using numpy
            PLdata = np.loadtxt(file_id, skiprows=1)

            # Extract the wavelength and intensity values from the PL data
            wavelength = PLdata[:, 0]
            Intensity = PLdata[:, 1]

            # Normalize the intensity values by dividing them by the maximum intensity
            NormIntensity = Intensity / max(Intensity)

            # Smooth the intensity values using a rolling mean with a window size of 20
            smoothed_Intensity = pd.DataFrame(NormIntensity).rolling(20, center=True).mean().values.flatten()

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, smoothed_Intensity, label=molecule_name)

        # Add a legend to the plot with font size of 12
        ax.legend(fontsize=12)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')

        # Show the plot
        plt.show()
