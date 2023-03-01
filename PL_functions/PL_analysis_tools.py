import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d




class Steady_state_analysis:
    def __init__(self, bin_file_id):
        self.bin_file_id = bin_file_id
        self.file_id = pd.read_csv(self.bin_file_id,
                              delimiter=',',
                              index_col=0,
                              header=None,
                              skiprows=1,
                              names=['Molecule_name', 'File_ID']).squeeze('columns').to_dict()

    def plot_PL(self):
        fig, ax = plt.subplots()
        for molecule_name, file_id in self.file_id.items():

            PLdata = np.loadtxt(file_id, skiprows=1)
            wavelength = PLdata[:, 0]
            Intensity = PLdata[:, 1]

            NormIntensity = Intensity / max(Intensity)
            ax.plot(wavelength, NormIntensity, label=molecule_name)
        ax.legend(fontsize=12)
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')
        plt.show()

    def smooth_plot_PL(self):
        fig, ax = plt.subplots()
        for molecule_name, file_id in self.file_id.items():

            PLdata = np.loadtxt(file_id, skiprows=1)
            wavelength = PLdata[:, 0]
            Intensity = PLdata[:, 1]

            NormIntensity = Intensity / max(Intensity)
            smoothed_Intensity = pd.DataFrame(NormIntensity).rolling(20, center=True).mean().values.flatten()
            ax.plot(wavelength, smoothed_Intensity, label=molecule_name)
        ax.legend(fontsize=12)
        ax.set_title('Normalized PL spectra')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Normalized Intensity')
        plt.show()





