from PL_functions import PL_analysis_tools as PATools
import os


def main():

    q1 = 'Where the bins file is located? (including .csv)'

    while True:
        bins_file_id = input(q1)
        if os.path.isfile(bins_file_id):
            break
        print('Invalid file path. Please enter a valid file path.')

    try:
        PLA = PATools.Steady_state_analysis(bins_file_id)
        PLA.plot_PL()
        PLA.smooth_plot_PL()

    except Exception as e:
        print(f'Error: {e}')




if __name__ == '__main__':
    main()


