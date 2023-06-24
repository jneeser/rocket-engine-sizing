import numpy as np
from matplotlib import pyplot as plt
from itertools import cycle

from CEA import BipropCEA, MonopropCEA
import PropLibrary as pl


def isp_plots(fuels, oxidisers, Pc, eps, MR_range):

	isp_arr = np.ndarray((len(fuels), len(MR_range)))
	
	# matplotlib linestyles
	lines = ["-","--","-.",":","-"]
	linecycler = cycle(lines)
	
	for i in range(len(fuels)):
		cea = BipropCEA(fuels[i], oxidisers[i], Pc)
		
		for j in range(len(MR_range)):
			cea.metric_cea_output('chamber', MR_range[j], eps)
			# calculate vacuum isp
			isp_arr[i][j] = cea.isp

		max_idx = np.where(isp_arr[i] == max(isp_arr[i]))[0][0]
		print(str(fuels[i]), ' + ' , str(oxidisers[i]), ' maximum Isp vac: ', np.round(max(isp_arr[i]),3), ' s at MR: ', np.round(MR_range[max_idx], 3))
	
	# plot isp data
	for i in range(len(fuels)):
		plt.plot(MR_range, isp_arr[i], next(linecycler), label=str(fuels[i]) + ' + ' + str(oxidisers[i]))
	
	plt.xlabel('MR')
	plt.ylabel('vac Isp [s]')
	plt.legend(loc='best')
	plt.title('Pc = ' + str(round(Pc/1e5,2)) + ' bar,' + ' eps = ' + str(eps))
	plt.show()


if __name__ == "__main__":
	fuels     = ['RP-1', pl.TMPDA, 'Isopropanol', 'C2H5OH']
	oxidisers = [pl.peroxide98, pl.peroxide98, pl.peroxide90, pl.peroxide90]
	Pc        = 10e5                    # chamber pressure [Pa]
	eps       = 200                      # expansion ratio
	MR_range  = np.linspace(3.5, 10, 99)  # mixture ratio range

	isp_plots(fuels, oxidisers, Pc, eps, MR_range)