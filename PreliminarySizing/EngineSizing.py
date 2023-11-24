import numpy as np
import scipy.optimize 
import thermo

from CEA import BipropCEA, MonopropCEA
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
	def __init__(self, fuel, oxidiser, chamber_pressure, expansion_ratio, contraction_ratio, divergence_angle, L_star, thrust, mixture_ratio, Pamb):
		self.chamber_pressure  = chamber_pressure
		self.expansion_ratio   = expansion_ratio
		self.fuel              = fuel
		self.oxidiser          = oxidiser
		self.divergence_angle  = divergence_angle
		self.L_star            = L_star
		self.contraction_ratio = contraction_ratio
		self.mixture_ratio     = mixture_ratio
		self.thrust            = thrust
		self.ambient_pressure  = Pamb
		self.mass_flow         = 0.622
		self.R   			   = 8314.5							# universal gas constant J/mol

		# get hot gas properties from CEA
		self.cea = BipropCEA(self.fuel, self.oxidiser, self.chamber_pressure, self.ambient_pressure)
		self.cea.metric_cea_output('throat', self.mixture_ratio, self.expansion_ratio)

		# isentropic gas object
		self.static = Isentropic(self.chamber_pressure, self.cea.Tc, self.cea.gamma)
		self.V = np.sqrt(self.cea.gamma * ((1+self.cea.gamma)/2) ** ((1+self.cea.gamma)/(1-self.cea.gamma)))


	def get_efficiency(self):
		# basic divergence efficiency for rocket nozzles
		eta_nozzle = np.sin(self.divergence_angle) / self.divergence_angle
		# assumed combustion efficency of 1
		self.eta = eta_nozzle


	def get_chamber_dimensions(self):
		self.throat_area = self.mass_flow * np.sqrt(self.R / self.cea.MW * self.cea.Tc) / self.V / self.chamber_pressure
		self.throat_diameter = 2 * np.sqrt(self.throat_area/np.pi)
		self.exit_area = self.expansion_ratio * self.throat_area
		self.exit_diameter = 2 * np.sqrt(self.exit_area/np.pi)
		
		V_c = self.L_star * self.throat_area
		# assuming converging section is 15% of chamber volume 
		self.chamber_area = self.contraction_ratio * self.throat_area
		self.chamber_diameter = 2 * np.sqrt(self.chamber_area/np.pi)
		self.L_cyl = 0.85 * V_c / self.chamber_area


	def claculate_m_dot(self):
		self.get_efficiency()
		
		# exist gas properties
		m_exit = self.static.mach(self.expansion_ratio)
		p_exit = self.static.pressure(m_exit)
		T_exit = self.static.temperature(m_exit)

		# exsit velocity and Isp
		self.v_exit = m_exit * np.sqrt(self.cea.gamma * self.R / self.cea.MW * T_exit)
		self.ispAmb = self.cea.ispAmb[0] * self.eta
		
		def func(mass_flow):
			# optimisation function for estimating mass flow for a given thrust
			throat_area      = mass_flow * np.sqrt(self.R / self.cea.MW * self.cea.Tc) / self.V / self.chamber_pressure
			exit_area        = self.expansion_ratio * throat_area
			estimated_thrust = self.eta * (mass_flow * self.v_exit + (p_exit - self.ambient_pressure) * exit_area)
			
			return estimated_thrust - self.thrust
		
		
		self.mass_flow = scipy.optimize.fsolve(func, 0.001)
		self.get_chamber_dimensions()	
		self.c_star = self.chamber_pressure * self.throat_area / self.mass_flow


if __name__ == '__main__':

	target_thrust = 40						# N
	fuel    = pl.TMPDA				# CEA fuel or custom from PropLibrary.py
	ox      = pl.peroxide95					# CEA oxidiser or custom from PropLibrary.py
	Pc      = 15e5							# Chamber pressure in Pa
	MR      = 6							# Mixture ratio (O/F)
	ER      = 3								# expansion ratio
	CR      = 9.5								# contraction ratio
	phi_div = np.radians(15)				# divergence angle of nozzle
	Pamb    = 101325 					    # Pa
	L_star  = 1
	

	engine = LPRE(fuel, ox, Pc, ER, CR, phi_div, L_star, target_thrust, MR, Pamb)
	engine.claculate_m_dot()

	
	print('Thrust         ', engine.thrust, ' N')
	print('Isp            ', engine.cea.ispAmb, 'NOTE Psep or Pe given in Psi')
	print('mass flow      ', np.round(engine.mass_flow, 4), ' kg/s')
	print('fuel mass flow   ', np.round(engine.mass_flow / (MR + 1), 5), ' kg/s')
	print('ox mass flow ', np.round(- engine.mass_flow / (MR + 1) + engine.mass_flow, 4), ' kg/s')
	print('c*             ', np.round(engine.c_star, 3), ' m/s')
	print('Dt             ', np.round(engine.throat_diameter*1000, 3), ' mm')
	print('Dc             ', np.round(engine.chamber_diameter*1000, 3), ' mm')
	print('De             ', np.round(engine.exit_diameter*1000, 3), ' mm')
	print('Lcyl           ', np.round(engine.L_cyl*1000, 3), ' mm')

	m_dot_f = engine.mass_flow / (MR + 1)
	m_dot_ox = engine.mass_flow - m_dot_f

	v_dot_f = m_dot_f / 780
	v_dot_ox = m_dot_ox / 1420

	v_dot = v_dot_f + v_dot_ox
	m_dot_n2 = v_dot * 23.11
	print(v_dot)
	print(m_dot_n2)