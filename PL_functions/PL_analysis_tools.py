import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from scipy.signal import find_peaks
import os

class Photoluminescence:
    def __init__(self, bin_file_id):
        self.bin_file_id = bin_file_id
        self.PLData = {}  # Dictionary to store PL data

        # Read the CSV file as a pandas DataFrame
        df = pd.read_csv(self.bin_file_id,
                         delimiter=',',
                         index_col=0,
                         header=None,
                         skiprows=1,
                         names=['Molecule_name', 'File_ID'])

        # Iterate through molecules in the DataFrame and load PL data
        for molecule_name, file_id in df.squeeze('columns').to_dict().items():
            PLData = np.loadtxt(file_id, skiprows=1)
            pl_intensity = np.trapz(PLData[:, 1], PLData[:, 0])
            self.PLData[molecule_name] = {
                'wavelength': PLData[:, 0],
                'counts': PLData[:, 1],
                'PL_intensity': pl_intensity
            }

    def convert_to_eV(self, wavelength):
        h = 4.1357 * 10 ** (-15)  # Planck constant eV
        c = 299792458  # speed of light in m/s
        energy = (h * c) / (wavelength * 10 ** (-9))
        return energy

    def plot_PL(self, x_values, y_values, x_label, y_label, title, legend_labels, plot_energy=False, plot_peaks=False, savefig=False):
        fig, ax = plt.subplots()

        for idx, y_data in enumerate(y_values):
            if plot_energy:
                x_data = self.convert_to_eV(x_values[idx])
                x_label = 'Energy (eV)'
            else:
                x_data = x_values[idx]
                x_label = 'Wavelength (nm)'

            ax.plot(x_data, y_data, label=legend_labels[idx])

            if plot_peaks:
                max_intensity = np.max(y_data)
                prominence_fraction = 0.3  # You can adjust this value as needed
                prominence = prominence_fraction * max_intensity
                peaks, _ = find_peaks(y_data, prominence, distance=15)
                peak_wavelengths = x_data[peaks]
                peak_heights = y_data[peaks]
                if plot_energy:
                    for peak_wavelength, peak_height in zip(peak_wavelengths, peak_heights):
                        print(legend_labels[idx], "Peak energy (eV):", np.round(peak_wavelength, decimals=3), "Peak height:",
                          np.round(peak_height, decimals=2))
                    ax.plot(peak_wavelengths, peak_heights, "x", label="Peaks")
                else:
                    for peak_wavelength, peak_height in zip(peak_wavelengths, peak_heights):
                        print(legend_labels[idx], "Peak wavelength (nm):", np.round(peak_wavelength, decimals=2), "Peak height:",
                          np.round(peak_height, decimals=2))
                    ax.plot(peak_wavelengths, peak_heights, "x", label="Peaks")

        ax.legend(fontsize=9, frameon=False, loc='best')
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.tick_params(labelsize=10, direction='in', axis='both', which='major', length=8, width=0.5)
        ax.tick_params(labelsize=10, direction='in', axis='both', which='minor', length=4, width=0.5)
        ax.minorticks_on()

        if savefig:
            bins_file_path = os.path.dirname(self.bin_file_id)
            output_dir = os.path.join(bins_file_path, 'plots')
            os.makedirs(output_dir, exist_ok=True)

            # Get current timestamp
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

            # Create a unique filename with timestamp
            filename = f"{title.replace(' ', '_')}_{timestamp}.png"
            output_path = os.path.join(output_dir,filename)
            plt.savefig(output_path, dpi=200)

        plt.show()

    def plot_unnorm_PL(self, plot_energy=False, plot_peaks=False, savefig=False):
        x_values = []
        y_values = []
        legend_labels = []

        for molecule_name, data in self.PLData.items():
            x_values.append(data['wavelength'])
            y_values.append(data['counts'])
            legend_labels.append(molecule_name)

        self.plot_PL(x_values, y_values, 'Wavelength (nm)', 'Counts', 'Un-normalized PL spectra', legend_labels, plot_energy, plot_peaks, savefig)

    def plot_unnorm_smooth_PL(self, smooth_level=20, plot_energy=False, plot_peaks=False, savefig=False):
        x_values = []
        y_values = []
        legend_labels = []

        for molecule_name, data in self.PLData.items():
            x_values.append(data['wavelength'])
            counts = data['counts']
            smoothed_counts = pd.DataFrame(counts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()
            y_values.append(smoothed_counts)
            legend_labels.append(molecule_name)

        self.plot_PL(x_values, y_values, 'Wavelength (nm)', 'Smoothed Counts', 'Un-normalized smoothed PL spectra', legend_labels, plot_energy, plot_peaks, savefig)

    def plot_norm_PL(self, plot_energy=False, plot_peaks=False, savefig=False):
        x_values = []
        y_values = []
        legend_labels = []

        for molecule_name, data in self.PLData.items():
            x_values.append(data['wavelength'])
            counts = data['counts']
            norm_counts = counts / max(counts)
            y_values.append(norm_counts)
            legend_labels.append(molecule_name)

        self.plot_PL(x_values, y_values, 'Wavelength (nm)', 'Normalized Counts', 'Normalized PL spectra', legend_labels, plot_energy, plot_peaks, savefig)

    def plot_norm_smooth_PL(self, smooth_level=20, plot_energy=False, plot_peaks=False, savefig=False):
        x_values = []
        y_values = []
        legend_labels = []

        for molecule_name, data in self.PLData.items():
            x_values.append(data['wavelength'])
            counts = data['counts']
            norm_counts = counts / max(counts)
            smoothed_counts = pd.DataFrame(norm_counts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()
            y_values.append(smoothed_counts)
            legend_labels.append(molecule_name)

        self.plot_PL(x_values, y_values, 'Wavelength (nm)', 'Smoothed Normalized Counts', 'Smoothed normalized PL spectra', legend_labels, plot_energy, plot_peaks, savefig)
    def plot_smooth_norm_PL(self, smooth_level=20, plot_energy=False, plot_peaks=False, savefig=False):
        x_values = []
        y_values = []
        legend_labels = []

        for molecule_name, data in self.PLData.items():
            x_values.append(data['wavelength'])
            counts = data['counts']

            smoothed_counts = pd.DataFrame(counts).rolling(smooth_level, center=True,
                                                           min_periods=1).mean().values.flatten()
            smoothed_norm_counts = smoothed_counts / max(smoothed_counts)
            y_values.append(smoothed_norm_counts)
            legend_labels.append(molecule_name)

        self.plot_PL(x_values, y_values, 'Wavelength (nm)', 'Smoothed Normalized Counts', 'Smoothed normalized PL spectra', legend_labels, plot_energy, plot_peaks, savefig)

    def plot_norm_AUC_PL(self, smooth_level=20, plot_energy=False, plot_peaks=False, savefig=False):
        x_values = []
        y_values = []
        legend_labels = []

        for molecule_name, data in self.PLData.items():
            x_values.append(data['wavelength'])
            counts = data['counts']
            x_auc = auc(data['wavelength'], counts)
            if x_auc == 0:
                norm_counts = counts
            else:
                norm_counts = counts / x_auc
            smoothed_counts = pd.DataFrame(norm_counts).rolling(smooth_level, center=True, min_periods=1).mean().values.flatten()
            y_values.append(smoothed_counts)
            legend_labels.append(molecule_name)

        self.plot_PL(x_values, y_values, 'Wavelength (nm)', 'Normalized then Smoothed Counts (AUC)', 'Normalized (AUC) Smoothed PL spectra', legend_labels, plot_energy, plot_peaks, savefig)



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





