import numpy as np
import scipy.optimize 

from CEAClass import CEA 
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
    def __init__(self, fuel, oxidiser, chamber_pressure, mixture_ratio, mass_flow, expansion_ratio, contraction_ratio, divergence_angle, L_star):
        self.chamber_pressure = chamber_pressure
        self.mixture_ratio = mixture_ratio
        self.mass_flow = mass_flow
        self.expansion_ratio = expansion_ratio
        self.fuel_mass_flow = mass_flow / (1+mixture_ratio)
        self.ox_mass_flow = mass_flow - self.fuel_mass_flow
        self.fuel = fuel
        self.ox = oxidiser
        self.divergence_angle = divergence_angle
        self.L_star = L_star
        self.contraction_ratio = contraction_ratio

        # get hot gas properties from CEA
        self.cea = CEA(fuel, oxidiser, self.chamber_pressure)
        self.cea.metric_cea_output('throat', self.mixture_ratio, self.expansion_ratio)


    def get_efficiency(self):
        eta_nozzle = np.sin(self.divergence_angle) / self.divergence_angle
        eta_combustion = 0.97
        self.eta = eta_combustion * eta_nozzle


    def get_throat_area(self):
        g = self.cea.gamma
        V = np.sqrt(g * ((1+g)/2) ** ((1+g)/(1-g)))
        R_a = 8314.5

        self.throat_area = self.mass_flow * np.sqrt(R_a / self.cea.MW * self.cea.Tc) / V / self.chamber_pressure
        self.throat_diameter = 2 * np.sqrt(self.throat_area/np.pi)
        self.exit_area = self.expansion_ratio * self.throat_area
        self.exit_diameter = 2 * np.sqrt(self.exit_area/np.pi)
		
        V_c = self.L_star * self.throat_area
        # assuming converging section is 15% of chamber volume 
        self.chamber_area = self.contraction_ratio * self.throat_area
        self.chamber_diameter = 2 * np.sqrt(self.chamber_area/np.pi)
        self.L_cyl = 0.85 * V_c / self.chamber_area
	
    def get_thrust(self, ambient_pressure):
        local = Isentropic(self.chamber_pressure, self.cea.Tc, self.cea.gamma)
        self.get_throat_area()
        self.get_efficiency()

        m_exit = local.mach(self.expansion_ratio)
        p_exit = local.pressure(m_exit)
        T_exit = local.temperature(m_exit)
        v_exit = m_exit * np.sqrt(self.cea.gamma * 8314.5 / self.cea.MW * T_exit)
        self.isp = self.eta * v_exit / 9.80665
        self.thrust = self.eta * (self.mass_flow * v_exit + (p_exit - ambient_pressure) * self.exit_area)


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
	

	def func(m_dot):
		engine = LPRE(fuel, ox, Pc, MR, m_dot, ER, CR, phi_div, L_star)
		engine.get_thrust(P_amb)
		
		return (engine.thrust - target_thrust)

	m_dot = scipy.optimize.fsolve(func, 18e-3)
	

	engine = LPRE(fuel, ox, Pc, MR, m_dot, ER, CR, phi_div, L_star)
	engine.get_thrust(P_amb)
	print('Thrust	     ', engine.thrust, ' N')
	print('mass flow     ', m_dot, ' kg/s')
	print('c* efficiency ', engine.eta)
	print('Dt            ', engine.throat_diameter*1000, ' mm')
	print('Dc            ', engine.chamber_diameter*1000, ' mm')
	print('De            ', engine.exit_diameter*1000, ' mm')
	print('Lcyl          ', engine.L_cyl*1000, ' mm')
     
