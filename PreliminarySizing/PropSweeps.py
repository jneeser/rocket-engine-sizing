import numpy as np
from matplotlib import pyplot as plt
from itertools import cycle

from CEA import BipropCEA, MonopropCEA
import PropLibrary as pl


def isp_plots(fuels, oxidisers, Pc, eps, MR_range, P_amb):

	isp_arr = np.ndarray((len(fuels), len(MR_range)))
	
	# matplotlib linestyles
	lines = ["-","--","-.",":","-"]
	linecycler = cycle(lines)
	
	for i in range(len(fuels)):
		cea = BipropCEA(fuels[i], oxidisers[i], Pc, P_amb)
		
		for j in range(len(MR_range)):
			cea.metric_cea_output('chamber', MR_range[j], eps)
			# calculate vacuum isp
			isp_arr[i][j] = cea.ispAmb[0]

		max_idx = np.where(isp_arr[i] == max(isp_arr[i]))[0][0]
		print(str(fuels[i]), ' + ' , str(oxidisers[i]), ' maximum Isp amb: ', np.round(max(isp_arr[i]),3), ' s at MR: ', np.round(MR_range[max_idx], 3))
	
	# plot isp data
	for i in range(len(fuels)):
		plt.plot(MR_range, isp_arr[i], next(linecycler), label=str(fuels[i]) + ' + ' + str(oxidisers[i]))
	
	plt.xlabel('MR')
	plt.ylabel('Isp [s]')
	plt.legend(loc='best')
	plt.title('Pc = ' + str(round(Pc/1e5,2)) + ' bar,' + ' eps = ' + str(eps) + ', P_amb = ' + str(round(P_amb/1e5,2)) + ' bar')
	plt.savefig('Isp.png', dpi=600)
	plt.show()
	


def temperature_plots(fuels, oxidisers, Pc, eps, MR_range, P_amb):

	T_arr = np.ndarray((len(fuels), len(MR_range)))
	
	# matplotlib linestyles
	lines = ["-","--","-.",":","-"]
	linecycler = cycle(lines)
	
	for i in range(len(fuels)):
		cea = BipropCEA(fuels[i], oxidisers[i], Pc, P_amb)
		
		for j in range(len(MR_range)):
			cea.metric_cea_output('chamber', MR_range[j], eps)
			# calculate vacuum isp
			T_arr[i][j] = cea.Tc

	# plot isp data
	for i in range(len(fuels)):
		plt.plot(MR_range, T_arr[i], next(linecycler), color='r')
	
	plt.vlines(42.17, 1000, 3200, linestyles='--', color='k', label='lean extinction limit')
	plt.vlines(1.27, 1000, 3200, linestyles=':', color='k', label='rich extinction limit')
	plt.vlines(3.945, 1000, 3200, linestyles='-.', color='b', label='stoichiometric')
	plt.xlabel('O/F')
	plt.ylabel('Adiabatic Flame Temperature (K)')
	plt.legend(loc='best')
	#plt.title('Pc = ' + str(round(Pc/1e5,2)) + ' bar,' + ' eps = ' + str(eps) + ', P_amb = ' + str(round(P_amb/1e5,2)) + ' bar')
	plt.savefig('Temperature.png', dpi=600)
	plt.show()



if __name__ == "__main__":
	fuels     = ['MMH', pl.TMPDA, 'Propylene', 'Propylene', 'C2H6']
	oxidisers = ['N2O4', pl.peroxide98, 'N2O', pl.peroxide98, pl.peroxide98]
	Pc        = 10e5                    # chamber pressure [Pa]
	P_amb     = 10
	eps       = 220                         # expansion ratio
	MR_range  = np.linspace(1.5, 12, 99)  # mixture ratio range

	isp_plots(fuels, oxidisers, Pc, eps, MR_range, P_amb)
	temperature_plots(['CH4'], ['GOX'], 3e5, 1, np.linspace(1,45,500), 1e5)