import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import os


class Photoluminescence:
    def __init__(self, bin_file_id):
        self.bin_file_id = bin_file_id

        # Create an empty dictionary to store the PL data for all molecules
        self.PLData = {}

        # Read the CSV file as a pandas DataFrame
        df = pd.read_csv(self.bin_file_id,
                         delimiter=',',
                         index_col=0,
                         header=None,
                         skiprows=1,
                         names=['Molecule_name', 'File_ID'])

        # Iterate through all the molecules and their file IDs in the DataFrame
        for molecule_name, file_id in df.squeeze('columns').to_dict().items():
            PLData = np.loadtxt(file_id, skiprows=1)
            pl_intensity = np.trapz(PLData[:, 1], PLData[:, 0])
            self.PLData[molecule_name] = {
                'wavelength': PLData[:, 0],
                'counts': PLData[:, 1],
                'PL_intensity': pl_intensity
            }
    def calculate_PL_Quench(self, molecule1, molecule2):
        # Get the PL intensity of the molecule before and after treatment
        intensity1 = self.PLData[molecule1]['PL_intensity']
        intensity2 = self.PLData[molecule2]['PL_intensity']

        if intensity2 < intensity1:
            quenched_molecule = molecule2
            unquenched_molecule = molecule1
            quenching_ratio = 1 - intensity2 / intensity1
        else:
            quenched_molecule = molecule1
            unquenched_molecule = molecule2
            quenching_ratio = 1 - intensity1 / intensity2

        print(f'{quenched_molecule} is quenched compared to {unquenched_molecule} by {quenching_ratio * 100:.2f}%')

    def unnorm_plot_PL(self, savefig=False):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules in the PLData dictionary
        for molecule_name, data in self.PLData.items():

            # Extract the wavelength and counts values from the PL data
            wavelength = data['wavelength']
            counts = data['counts']
            #counts = counts - np.min(counts)

            # Plot the normalized PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, counts, label=molecule_name)

        # Add a legend to the plot with font size of 12
        ax.legend(fontsize=9, frameon=False, loc='upper right')
        # Add horizontal line at y=0
        #ax.axhline(0, color='black', linestyle='-', linewidth=1)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Un-normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Counts')

        ax.tick_params(labelsize=10, direction='in', axis='both', which='major', length=8, width=0.5)
        ax.tick_params(labelsize=10, direction='in', axis='both', which='minor', length=4, width=0.5)
        ax.minorticks_on()


        if savefig != False:
            # Get the path of the bins.csv file and create the output directory
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Save the plot with DPI = 200 in the output directory
            output_path = os.path.join(output_dir, 'Un-normalized_PL_plot.png')
            plt.savefig(output_path, dpi=200)

        # Show the plot
        plt.show()


        # Method to plot the photoluminescence spectra for all molecules
    def unnorm_smooth_plot_PL(self, smooth_level=20, savefig=False):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules in the PLData dictionary
        for molecule_name, data in self.PLData.items():
            # Extract the wavelength and counts values from the PL data
            wavelength = data['wavelength']
            counts = data['counts']
            #counts = counts - np.min(counts)

            # Smooth the counts values using a rolling mean with a window size of 20
            smoothed_Counts = pd.DataFrame(counts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()



            # Plot the normalized PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, smoothed_Counts, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=9, frameon=False, loc='upper right')
        # Add horizontal line at y=0
        #ax.axhline(0, color='black', linestyle='-', linewidth=1)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Un-normalized smoothed PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Counts')

        ax.tick_params(labelsize=10, direction='in', axis='both', which='major', length=8, width=0.5)
        ax.tick_params(labelsize=10, direction='in', axis='both', which='minor', length=4, width=0.5)
        ax.minorticks_on()

        if savefig != False:
            # Get the path of the bins.csv file and create the output directory
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Save the plot with DPI = 200 in the output directory
            output_path = os.path.join(output_dir, 'Un-normalized_smoothed_PL_plot.png')
            plt.savefig(output_path, dpi=200)

        # Show the plot
        plt.show()


    # Method to plot the smoothed photoluminescence spectra for all molecules

    def norm_plot_PL(self, savefig=False):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLData.items():
            # Extract the wavelength and counts values from the PL data
            wavelength = data['wavelength']
            counts = data['counts']

            # Normalize the counts values by dividing them by the maximum counts
            normCounts = counts / max(counts)

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, normCounts, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)
        # Add horizontal line at y=0
        #ax.axhline(0, color='black', linestyle='-', linewidth=1)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Counts')

        if savefig != False:
            # Get the path of the bins.csv file and create the output directory
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Save the plot with DPI = 200 in the output directory
            output_path = os.path.join(output_dir, 'Normalized_Un-smoothed_PL_plot.png')
            plt.savefig(output_path, dpi=200)

        # Show the plot
        plt.show()





    def norm_smooth_plot_PL(self, smooth_level=20, savefig=False):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLData.items():

            # Extract the wavelength and counts values from the PL data
            wavelength = data['wavelength']
            counts = data['counts']

            # Normalize the counts values by dividing them by the maximum counts
            NormCounts = counts / max(counts)

            # Smooth the counts values using a rolling mean with a window size of 20
            smoothed_Counts = pd.DataFrame(NormCounts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, smoothed_Counts, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)
        # Add horizontal line at y=0
        #ax.axhline(0, color='black', linestyle='-', linewidth=1)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Smoothed normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Counts')

        if savefig != False:
            # Get the path of the bins.csv file and create the output directory
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Save the plot with DPI = 200 in the output directory
            output_path = os.path.join(output_dir, 'Smoothed_normalized_PL_plot.png')
            plt.savefig(output_path, dpi=200)

        # Show the plot
        plt.show()

    def smooth_norm_plot_PL(self, smooth_level=20, savefig=False):

        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLData.items():

            # Extract the wavelength and counts values from the PL data
            wavelength = data['wavelength']
            counts = data['counts']

            # Smooth the counts values using a rolling mean with a window size of 20
            smoothed_Counts = pd.DataFrame(counts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()


            # Normalize the counts values by dividing them by the maximum counts
            Norm_smooth_Counts = smoothed_Counts / max(smoothed_Counts)

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, Norm_smooth_Counts, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)
        # Add horizontal line at y=0
        #ax.axhline(0, color='black', linestyle='-', linewidth=1)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized smoothed PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Counts')

        if savefig != False:
            # Get the path of the bins.csv file and create the output directory
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Save the plot with DPI = 200 in the output directory
            output_path = os.path.join(output_dir, 'Normalized_smoothed_PL_plot.png')
            plt.savefig(output_path, dpi=200)

        # Show the plot
        plt.show()

    def norm_AUC_plot_PL(self, smooth_level=20, savefig=False):
        # Create a figure and an axis object using matplotlib
        fig, ax = plt.subplots()

        # Iterate through all the molecules and their file IDs in the dictionary
        for molecule_name, data in self.PLData.items():
            # Extract the wavelength and counts values from the PL data
            wavelength = data['wavelength']
            counts = data['counts']

            # Normalize the counts values by dividing them by the AUC (Area Under the Curve)
            x_auc = auc(wavelength, counts)
            if x_auc == 0:
                normCounts = counts
            else:
                normCounts = counts / x_auc

                # Smooth the counts values using a rolling mean with a window size of 20
            smoothedCounts = pd.DataFrame(normCounts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()

                # Plot the smoothed PL spectrum for the molecule with its name as the label

            # Plot the smoothed PL spectrum for the molecule with its name as the label
            ax.plot(wavelength, smoothedCounts, label=molecule_name)

        # Add a legend to the plot with font size of 8
        ax.legend(fontsize=8)
        # Add horizontal line at y=0
        #ax.axhline(0, color='black', linestyle='-', linewidth=1)

        # Set the title, xlabel, and ylabel of the plot
        ax.set_title('Normalized (AUC) PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Counts')

        if savefig != False:

            # Get the path of the bins.csv file and create the output directory
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Save the plot with DPI = 500 in the output directory
            output_path = os.path.join(output_dir, 'normalized_AUC_smoothed_PL_plot.png')
            plt.savefig(output_path, dpi=200)

        # Show the plot
        plt.show()






