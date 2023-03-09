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
        PLA.unnorm_plot_PL()

        # Plotting the unnormalized smoothed PL spectra
        PLA.unnorm_smooth_plot_PL()

        # Plotting the normalized PL spectra
        PLA.norm_plot_PL()

        # Plotting the smoothed normalized PL spectra
        PLA.norm_smooth_plot_PL()

        PLA.norm_AUC_plot_PL()



    except Exception as e:
        # If any error occurs during the execution of the code, it will be printed
        print(f'Error: {e}')


if __name__ == '__main__':
    main()


