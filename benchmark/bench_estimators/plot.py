from pylab import *


class benchmark_integrate_nmc:

    def __init__(self, filename, label):
        self.label = label
        self.data = {}

        with open(filename) as bmfile:
            for line in bmfile:

                lsplit = line.split()

                if len(lsplit) != 9:
                    continue

                if lsplit[0][0:9] == 't/element':
                    self.data[(lsplit[2][5:-1] + "-dim", lsplit[3])] = (float(lsplit[5]), float(lsplit[7]))


def plot_compare_nmc(benchmark_list, **kwargs):
    nbm = len(benchmark_list)
    print(nbm)
    print(benchmark_list)
    print(benchmark_list[0].data)

    fig = figure()
    fig.suptitle('Estimators benchmark, comparing averaging techniques\n ', fontsize=12)
    for i in range(4):
        ax = fig.add_subplot(2, 2, i + 1)
        xlabels = [tup[1] for tup in list(benchmark_list[0].data.keys())[
                                     i * 3:(i + 1) * 3]]  # get the xlabels from first entry in data dict

        for benchmark in benchmark_list:
            values = [benchmark.data[key][0] for key in list(benchmark.data.keys())[i * 3:(i + 1) * 3]]
            errors = [benchmark.data[key][1] for key in list(benchmark.data.keys())[i * 3:(i + 1) * 3]]
            ax.errorbar(xlabels, values, xerr=None, yerr=errors, **kwargs)
        ax.set_title(list(benchmark_list[0].data.keys())[i * 3][0])
        if (i % 2 == 0):
            ax.set_ylabel('Time per element [ns]')
        ax.set_ylim(0.1, 100.)
        ax.set_yscale('log')
        ax.legend([bench.label for bench in benchmark_list])

    return fig


# Script

benchmark_list = []
for benchmark_file in sys.argv[1:]:
    try:
        benchmark = benchmark_integrate_nmc(benchmark_file, benchmark_file.split('_')[1].split('.')[0])
        benchmark_list.append(benchmark)
    except(OSError):
        print("Warning: Couldn't load benchmark file " + benchmark_file + "!")

if len(benchmark_list) < 1:
    print("Error: Not even one benchmark loaded!")
else:
    fig1 = plot_compare_nmc(benchmark_list, fmt='o--')

tight_layout()
show()
