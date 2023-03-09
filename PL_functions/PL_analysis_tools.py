import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import os


class Photoluminescence:
    def __init__(self, bin_file_id):
        self.bin_file_id = bin_file_id

        # Create an empty dictionary to store the PL data for all molecules
        self.PLdata = {}

        # Read the CSV file as a pandas DataFrame
        df = pd.read_csv(self.bin_file_id,
                         delimiter=',',
                         index_col=0,
                         header=None,
                         skiprows=1,
                         names=['Molecule_name', 'File_ID'])

        # Iterate through all the molecules and their file IDs in the DataFrame
        for molecule_name, file_id in df.squeeze('columns').to_dict().items():

            # Load the PL data from the file using numpy and store it in the PLdata dictionary
            PLdata = np.loadtxt(file_id, skiprows=1)
            self.PLdata[molecule_name] = {
                'wavelength': PLdata[:, 0],
                'intensity': PLdata[:, 1]
            }

    def unnorm_plot_PL(self):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules in the PLdata dictionary
        for molecule_name, data in self.PLdata.items():

            # Extract the wavelength and intensity values from the PL data
            wavelength = data['wavelength']
            intensity = data['intensity']

            # Plot the normalized PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, intensity, label=molecule_name)

        # Add a legend to the plot with font size of 12
        ax.legend(fontsize=8)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Un-normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')

        # Show the plot
        plt.show(dpi=500)


        # Method to plot the photoluminescence spectra for all molecules
    def unnorm_smooth_plot_PL(self):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules in the PLdata dictionary
        for molecule_name, data in self.PLdata.items():
            # Extract the wavelength and intensity values from the PL data
            wavelength = data['wavelength']
            intensity = data['intensity']

            # Smooth the intensity values using a rolling mean with a window size of 20
            smoothed_Intensity = pd.DataFrame(intensity).rolling(20, center=True).mean().values.flatten()



            # Plot the normalized PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, smoothed_Intensity, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Un-normalized smoothed PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Intensity')

        # Show the plot
        plt.show(dpi=500)

    # Method to plot the smoothed photoluminescence spectra for all molecules

    def norm_plot_PL(self):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLdata.items():
            # Extract the wavelength and intensity values from the PL data
            wavelength = data['wavelength']
            intensity = data['intensity']

            # Normalize the intensity values by dividing them by the maximum intensity
            normIntensity = intensity / max(intensity)

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, normIntensity, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')

        # Show the plot
        plt.show(dpi=500)



    def norm_smooth_plot_PL(self):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLdata.items():

            # Extract the wavelength and intensity values from the PL data
            wavelength = data['wavelength']
            intensity = data['intensity']

            # Normalize the intensity values by dividing them by the maximum intensity
            NormIntensity = intensity / max(intensity)

            # Smooth the intensity values using a rolling mean with a window size of 20
            smoothed_Intensity = pd.DataFrame(NormIntensity).rolling(20, center=True).mean().values.flatten()

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, smoothed_Intensity, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Smoothed normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')

        # Show the plot
        plt.show(dpi=500)

    def norm_AUC_plot_PL(self):
        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLdata.items():
            # Extract the wavelength and intensity values from the PL data
            wavelength = data['wavelength']
            intensity = data['intensity']

            # Normalize the intensity values by dividing them by the AUC (Area Under the Curve)
            x_auc = auc(wavelength, intensity)
            if x_auc == 0:
                normIntensity = intensity
            else:
                normIntensity = intensity / x_auc

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, normIntensity, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')

        # Get the path of the bins.csv file and create the output directory
        bins_file_path = os.path.dirname(self.bin_file_id)
        output_dir = os.path.join(bins_file_path, 'plots')
        os.makedirs(output_dir, exist_ok=True)

        # Save the plot with DPI = 500 in the output directory
        output_path = os.path.join(output_dir, 'normalized_PL_plot.png')
        plt.savefig(output_path, dpi=500)

        # Show the plot
        plt.show()




