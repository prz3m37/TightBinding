import datetime
import warnings 
import config as cfg
warnings.filterwarnings("ignore")
from helpers import TightBindingHelpers
from execute_tb_calculations import ExecuteTightBindingCalculations


def main()->None:

    """
    Function where user calls all required methods and classes,
    Returns: None

    """

    print('\n[INFO]: Tight Binding calculations have started \n')
    execute = ExecuteTightBindingCalculations()
    helpers = TightBindingHelpers(execute.parametrization)

    start = datetime.datetime.now()
    lattice_df_format, dimension = None, None

    try: 
        energies, wave_function = execute.call_tight_binding_calculation(dimension, lattice_df_format)
    except RuntimeError:
        print("[ERROR]: Factor is exactly singular \n")
        print("[INFO]: Calculations have stopped \n")
        return
    except TypeError:
        print("[ERROR]: No data to calculate. Please check your configuration or input \n")
        print("[INFO]: Calculations have stopped \n")
        return

    density_of_states = execute.calculate_DOS(energies)
    end = datetime.datetime.now()

    print('[INFO]: Saving results \n')

    helpers.create_saving_folder()
    file_name = helpers.get_file_name(dimension)
    helpers.save_numerical_results(file_name, energies)
    helpers.plot_DOS(file_name, dimension, density_of_states)

    print('[INFO]: Calculation time for ' + str(dimension) + ' atoms: ', end - start, '\n')
    print('[INFO]: Calculations have sucesfully finished \n')

    return


if __name__ == '__main__':
    exit(main())
