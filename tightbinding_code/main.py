import sys
import datetime
import warnings
import config as cfg
warnings.filterwarnings("ignore")
from helpers import TightBindingHelpers
from execute_tb_calculations import ExecuteTightBindingCalculations

def main() -> None:

    """
    Function where user calls all required methods and classes,
    Returns: None

    """

    sys.path.insert(1, '..')
    settings = cfg.settings
    parametrization = sys.argv[1]
    helpers = TightBindingHelpers()
    execute = ExecuteTightBindingCalculations(parametrization, helpers, settings)

    helpers.create_logfile()
    helpers.save_log('\n[INFO]: Tight Binding calculations have started \n')
    lattice_df_format, dimension = None, None
    helpers.save_log('\n[INFO]: Input DataFrame loaded into memory \n')

    # TODO:  Execute...
    start = datetime.datetime.now()
    try:
        energies, wave_function = execute.call_tight_binding_calculation(dimension, lattice_df_format)
        density_of_states = execute.calculate_DOS(energies)
    except RuntimeError:
        helpers.save_log("[ERROR]: Factor is exactly singular \n")
        helpers.save_log("[INFO]: Calculations have stopped \n")
        helpers.close_logfile()
        return
    except TypeError:
        helpers.save_log("[ERROR]: No data to calculate. Please check your configuration or input \n")
        helpers.save_log("[INFO]: Calculations have stopped \n")
        helpers.close_logfile()
        return

    end = datetime.datetime.now()
    # TODO:  data managment.
    helpers.save_log('[INFO]: Saving results \n')

    helpers.create_saving_folder()
    helpers.save_numerical_results('1', energies)
    # helpers.plot_DOS(file_name, dimension, density_of_states)

    helpers.save_log('[INFO]: Calculation time for ' + str(dimension) + ' atoms: ' + str(end - start) + '\n')
    helpers.save_log('[INFO]: Calculations have successfully finished \n')
    helpers.close_logfile()

    return


if __name__ == '__main__':
    exit(main())
