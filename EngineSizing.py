import numpy as np
import scipy.optimize 

from CEA import CEA 
import PropLibrary as pl


class Isentropic():
	def __init__(self, total_pressure, total_temperature, gamma):
		self.p_t = total_pressure
		self.T_t = total_temperature
		self.gamma = gamma

	def mach(self, area_ratio):
		# for diverging section only 
		mach = lambda M:  1/(M*M) * (2/(self.gamma+1) * (1 + (self.gamma-1)/2*M*M))**((self.gamma+1)/(self.gamma-1)) - (area_ratio)**2
		mach_number = scipy.optimize.fsolve(mach, 10)
		return mach_number

	def pressure(self,mach):
		return self.p_t/((1 + (self.gamma-1)/2 * mach**2)**(self.gamma/(self.gamma-1)))

	def temperature(self,mach):
		return self.T_t/(1 + (self.gamma-1)/2 * mach**2)


class LPRE():
	def __init__(self, fuel, oxidiser, chamber_pressure, expansion_ratio, contraction_ratio, divergence_angle, L_star, thrust, mixture_ratio):
		self.chamber_pressure  = chamber_pressure
		self.expansion_ratio   = expansion_ratio
		self.fuel              = fuel
		self.oxidiser          = oxidiser
		self.divergence_angle  = divergence_angle
		self.L_star            = L_star
		self.contraction_ratio = contraction_ratio
		self.mixture_ratio     = mixture_ratio
		self.thrust            = thrust
		self.mass_flow         = 0.03
		self.R   			   = 8314.5							# universal gas constant J/g


	def get_gas_properties(self):
		# get hot gas properties from CEA
		self.cea = CEA(self.fuel, self.oxidiser, self.chamber_pressure)
		self.cea.metric_cea_output('throat', self.mixture_ratio, self.expansion_ratio)
		
		self.V = np.sqrt(self.cea.gamma * ((1+self.cea.gamma)/2) ** ((1+self.cea.gamma)/(1-self.cea.gamma)))

		# isentropic gas object
		self.static = Isentropic(self.chamber_pressure, self.cea.Tc, self.cea.gamma)

	def get_efficiency(self):
		# basic divergence efficiency for rocket nozzles
		eta_nozzle = np.sin(self.divergence_angle) / self.divergence_angle
		# assumed combustion efficency of 1
		self.eta = eta_nozzle


	def get_throat_area(self):
		# basic chamber geometry
		self.throat_area = self.mass_flow * np.sqrt(self.R / self.cea.MW * self.cea.Tc) / self.V / self.chamber_pressure
		self.throat_diameter = 2 * np.sqrt(self.throat_area/np.pi)
		self.exit_area = self.expansion_ratio * self.throat_area
		self.exit_diameter = 2 * np.sqrt(self.exit_area/np.pi)
		
		V_c = self.L_star * self.throat_area
		# assuming converging section is 15% of chamber volume 
		self.chamber_area = self.contraction_ratio * self.throat_area
		self.chamber_diameter = 2 * np.sqrt(self.chamber_area/np.pi)
		self.L_cyl = 0.85 * V_c / self.chamber_area


	def claculate(self, ambient_pressure):
		self.get_gas_properties()
		self.get_efficiency()
		
		# exist gas properties
		m_exit = self.static.mach(self.expansion_ratio)
		p_exit = self.static.pressure(m_exit)
		T_exit = self.static.temperature(m_exit)

		# exsit velocity and Isp
		self.v_exit = m_exit * np.sqrt(self.cea.gamma * 8314.5 / self.cea.MW * T_exit)
		self.isp = self.eta * self.v_exit / 9.80665

		def func(mass_flow):
			# optimisation function for estimating mass flow for a given thrust
			throat_area      = mass_flow * np.sqrt(self.R / self.cea.MW * self.cea.Tc) / self.V / self.chamber_pressure
			exit_area        = self.expansion_ratio * throat_area
			estimated_thrust = self.eta * (self.mass_flow * self.v_exit + (p_exit - ambient_pressure) * exit_area)
			
			return abs(estimated_thrust - self.thrust)
		
		
		self.mass_flow = scipy.optimize.fsolve(func, 0.1)
		self.get_throat_area()



if __name__ == '__main__':

	target_thrust = 200				# N
	fuel = pl.HIP11
	ox = pl.peroxide98
	Pc = 10e5						# Pa
	MR = 3.8						
	ER = 30
	CR = 14
	phi_div = 15 * np.pi/180
	P_amb = 10 					    # Pa
	L_star = 1.2
	

	engine = LPRE(fuel, ox, Pc, ER, CR, phi_div, L_star, target_thrust, MR)
	engine.claculate(P_amb)

	
	print('Thrust	     ', engine.thrust, ' N')
	print('mass flow     ', engine.mass_flow, ' kg/s')
	print('c* efficiency ', engine.eta)
	print('Dt            ', engine.throat_diameter*1000, ' mm')
	print('Dc            ', engine.chamber_diameter*1000, ' mm')
	print('De            ', engine.exit_diameter*1000, ' mm')
	print('Lcyl          ', engine.L_cyl*1000, ' mm')
	 
