from pylab import *

class benchmark_integrate_nmc:

    def __init__(self, filename, label):
        self.label = label
        self.data = {}

        with open(filename) as bmfile:
            for line in bmfile:

                lsplit = line.split()

                if len(lsplit) != 7:
                    continue

                if lsplit[0][0:6] == 't/step':
                    self.data[lsplit[1][1:]] = (float(lsplit[3]), float(lsplit[5]))


def plot_compare_nmc(benchmark_list, **kwargs):
    nbm = len(benchmark_list)
    print(nbm)
    print(benchmark_list)
    print(benchmark_list[0].data)
    xlabels = list(benchmark_list[0].data.keys()) # get the xlabels from first entry in data dict

    fig = figure()
    fig.suptitle('MCIntegrate benchmark, comparing different Nmc',fontsize=14)
    ax = fig.add_subplot(1, 1, 1)

    for benchmark in benchmark_list:
        values = [ benchmark.data[key][0] for key in benchmark.data.keys() ]
        errors = [ benchmark.data[key][1] for key in benchmark.data.keys() ]
        ax.errorbar(xlabels, values, xerr=None, yerr=errors, **kwargs)

    ax.set_ylabel('Time per sample [$\mu s$]')
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

if len(benchmark_list)<1:
    print("Error: Not even one benchmark loaded!")
else:
    fig1 = plot_compare_nmc(benchmark_list, fmt='o--')

show()
