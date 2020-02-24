import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

class Plotter(object):

    def __init__(self, data_menager):
        self.__helpers = None
        self.__data_menager = data_menager

    def data_parser(self, data):
        with open(data, "r") as results:
            results = results.read().split('\n')
            results = [float(res) for res in results[:-1]]
        eigenenergies = results
        return eigenenergies



    def plot_dos(self, , gauss_width, start, stop, step):

        with open(results_file, "r") as results:
            results = results.read().split('\n')
            results = [float(res) for res in results[:-1]]

        eigenenergies = results
        gauss_width = gauss_width

        D_E = 0
        E = np.arange(start, stop, step)
        for eigenenergy in eigenenergies:
            D_E = D_E + np.exp(-(E - eigenenergy)**2 / (2 * gauss_width**2)) / (np.pi * gauss_width * np.sqrt(2))

        font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 16}
        plt.figure(figsize=(13.66, 7.68))
        plt.plot(E, D_E)
        plt.xlabel('\nEnergy [a.u.]', fontsize=15,fontdict=font)
        section = np.arange(-1, 1, 1/20.)
        plt.fill_between(E, D_E, color='blue', alpha=0.3)
        plt.ylabel('DOS\n', fontsize=15,fontdict=font)
        plt.title('Density of states\n', fontsize=15,fontdict=font)
        plt.xlim(start, stop)
        plt.ylim(bottom=0)
        plt.subplots_adjust(left=0.15)
        plt.xticks(fontsize=11)
        plt.yticks(fontsize=11)
        #plt.gca().spines['right'].set_position(('data',0))
        #plt.gca().spines['top'].set_position(('data',0))
        plt.savefig(results_file + '.png', dpi=400)
        plt.grid(False)
        plt.close()
        return

def main():


    sns.set()
    start = [-7,-6,-1.1,-6]#-7,-5.5,-5,-7,-0.1,-7,-5.,-6.6,-7,-0.5,-6.5,-7,-5,-7,-6,-7,-7,-7,0.1,0.5,-6,-0.5,-7,-7,-0.6,-7,-5.5,-6,-7,-7,-7,-7,-7,-6.5,-7,-7,-7

    stop = [7,6,10.1,6] #7,14.5,5,7,14.5,7,13.5,6.5,7,15.5,15,7,14.,7,6,7,7,7,14.5,14.5,6,10,7,7,15.5,7,13.7,6,7,7,7,7,7,6.5,7,7,7

    step = 0.01
    gauss_width = 0.06
    path = "/home/przemek/Documents/Modeling/tight_binding/results_diploma"
    results = []
    print(len(start), len(stop))
    os.chdir(path)
    for file in glob.glob("*.txt"):
        input_file = path + '/' + file
        ready_input_file = open(input_file, 'r')
        num_list = [float(num) for num in ready_input_file.read().split()]
        max_val = max(num_list)
        min_val = min(num_list)
        results.append([max_val, min_val, file])

    for num, result in enumerate(results):
        print(result[2])
        print(start[num], stop[num])
        quick_plot(path + '/' + result[2], gauss_width, start[num], stop[num], step)

    return

if __name__ == '__main__':
    exit(main())

