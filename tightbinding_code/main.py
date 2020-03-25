import sys
import datetime
import warnings
import config as cfg
from data_manager import DataManager
from helpers import TightBindingHelpers
from execute_tb_calculations import ExecuteTightBindingCalculations
warnings.filterwarnings("ignore")


def main() -> None:
    """
    Function where user calls all required methods and classes,
    Returns: None
    """

    sys.path.insert(1, '..')
    settings = cfg.settings
    config_title = sys.argv[1]
    configuration = cfg.configuration[config_title]
    helpers = TightBindingHelpers()
    data_manager = DataManager(helpers, settings)
    execute = ExecuteTightBindingCalculations(helpers, settings, configuration)

    helpers.create_logfile()
    helpers.save_log('\n[INFO]: Tight Binding calculations have started \n')
    input_data = settings["input data"]
    dimension = input_data.shape[0]
    helpers.save_log('\n[INFO]: Input DataFrame loaded into memory \n')

    start = datetime.datetime.now()
    energies, wave_functions, interaction_matrix, dos, p_dos = execute.execute_tb_calculations(
        input_data)
    end = datetime.datetime.now()

    helpers.save_log('[INFO]: Saving results \n')
    helpers.create_saving_folder()
    data_manager.save_data(energies, wave_functions,
                           interaction_matrix, dos, p_dos, configuration)
    helpers.save_log('[INFO]: Results saved \n')

    if settings["load file"] is not None:
        helpers.save_log('[INFO]: Loading results \n')
        data = data_manager.load_data()
        # Please, process you data according to your needs.
        helpers.save_log('[INFO]: Results visualized \n')
    else:
        pass

    helpers.save_log('[INFO]: Calculation time for ' +
                     str(dimension) + ' atoms: ' + str(end - start) + '\n')
    helpers.save_log('[INFO]: Calculations have successfully finished \n')
    helpers.close_logfile()

    return


if __name__ == '__main__':
    exit(main())
