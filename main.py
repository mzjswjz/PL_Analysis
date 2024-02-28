from PL_functions import PL_analysis_tools as PATools
import os


def main():
    # Asking user to enter file path for the bins file
    q1 = 'Where the bins file is located? (including .csv)'

    # Checking if the entered file path is valid
    while True:
        bins_file_id = input(q1)
        if os.path.isfile(bins_file_id):
            break
        print('Invalid file path. Please enter a valid file path.')

    try:
        # Creating an instance of the Steady_state_analysis class
        PLA = PATools.Photoluminescence(bins_file_id)

        # Plotting the unnormalized PL spectra
        #PLA.plot_unnorm_PL(plot_energy=True, plot_peaks=True, savefig=False)

        # Plotting the unnormalized smoothed PL spectra
        #PLA.plot_unnorm_smooth_PL(smooth_level=20, plot_energy=True, plot_peaks=True, savefig=False)

        # Plotting the normalized PL spectra
        #PLA.plot_norm_PL(plot_energy=True, plot_peaks=True, savefig=False)

        # Plotting smoothed then normalized PL spectra
        PLA.plot_smooth_norm_PL(smooth_level=25, plot_energy=False, plot_peaks=True, savefig=True)

        # Plotting the normalized then smoothed PL spectra
        #PLA.plot_norm_smooth_PL(smooth_level=20, plot_energy=True, plot_peaks=True, savefig=False)

        # Plotting the AUC normalized then smoothed PL spectra
        #PLA.plot_norm_AUC_PL(smooth_level=20, plot_energy=True, plot_peaks=True, savefig=False)

        #PLA.calculate_PL_Quench('DTDCPB:C60 110C', 'DTDCPB 110C')



    except Exception as e:
        # If any error occurs during the execution of the code, it will be printed
        print(f'Error: {e}')


if __name__ == '__main__':
    main()


