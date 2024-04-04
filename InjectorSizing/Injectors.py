import numpy as np
import thermo
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

class LiquidInjector():
    def __init__(self, fluid, mixture, temperature, pressure, length, massflow, pressuredrop, inletangle):
        self.fluid = thermo.Mixture(fluid, ws=mixture,  T=temperature, P=pressure)
        self.length = length 
        self.massflow = massflow
        self.pressuredrop = pressuredrop
        self.inletangle = inletangle

    def xi1c(self, Re):
        # curfve fit from "Liquid Rocket Trhust Chambers" page 48 realting inlet efficiency to Reynolds number
        Re = np.log10(Re)
        return 3.55378*np.exp(-Re*0.647016) - 0.103358

    def xiinlet(self):
        return 0.5 + 2.5/np.pi*self.inletangle         # min 0.5 for coaxial flow before injection, max 0.9 for flow at pi/6 rad realtive to faceplate. orignially 0.5 + 1.2/np.pi*self.inletangle 
    
    def injector(self, maxiter=100, tol=1e-6):
        it = 0 
        difference = 1 
        xiinlet = self.xiinlet()
        xi = xiinlet
        while difference > tol:
            mu = 1/np.sqrt(1 + xi)
            diameter = 0.95*self.massflow**0.5 * mu**(-0.5) * (self.fluid.rho*self.pressuredrop)**(-0.25)
            vel = 1.273*self.massflow*self.fluid.rho**(-1)*diameter**(-2)
            Re = self.fluid.rho*vel*diameter/self.fluid.mu
            lam = 0.3164*Re**(-0.25)
            friction = lam*self.length/diameter
            xi1c = self.xi1c(Re)

            xi = xiinlet + xi1c + friction 
            newmu = 1/np.sqrt(1 + xi)

            it += 1
            if it > maxiter:
                raise ValueError("Not converged after ", maxiter, " iterations")
            difference = abs(mu - newmu)

        self.mu = mu
        self.diameter = diameter
        self.velocity = vel

class GasInjector():
    def __init__(self, gas, temperature, pressure, length, massflow, pressuredrop, upstreamdiameter, inletangle):
        self.gas = thermo.Chemical(gas, T=temperature, P=pressure)
        self.length = length 
        self.massflow = massflow
        self.pressuredrop = pressuredrop
        self.chamberpressure = pressure - pressuredrop
        self.upstreamdiameter = upstreamdiameter
        self.inletangle = inletangle

    def xiinlet(self):
        return 0.5 + 1.2/np.pi*self.inletangle
    
    # from "Liquid Rocket Trhust Chambers" page 52
    def injector(self, maxiter=100, tol=1e-10):
        if self.gas.phase == 'l':
            raise ValueError("Fluid is not gaseous before injection")
        it = 0 
        difference = 1 
        mu = 0.9                        # initial guess
        gamma = self.gas.Cpg/self.gas.Cvg
        R = self.gas.Cpg-self.gas.Cvg

        while difference > tol:
            lam2 = np.sqrt((gamma+1)/(gamma-1) * (1 - (self.chamberpressure/self.gas.P)**((gamma-1)/gamma)))
            q = ((gamma+1)/2)**(1/(gamma-1)) * lam2*(1 - (gamma-1)/(gamma+1)*lam2*lam2)**(1/(gamma-1))
            c = np.sqrt(gamma*R*self.gas.T) / (gamma*np.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1))))
            An = self.massflow*c / (mu*self.gas.P*q)
            diameter = 1.128*np.sqrt(self.massflow*c / (mu*self.gas.P*q))
            vel = lam2 * np.sqrt(2*gamma/(gamma-1) * R*self.gas.T * (1 - (self.chamberpressure/self.gas.P)**((gamma-1)/gamma)))
            Re = self.gas.rho*vel*diameter/self.gas.mug
            lam = 0.3164*Re**(-0.25)
            xi = lam*self.length/diameter + self.xiinlet()*(1-diameter**2/self.upstreamdiameter**2)
            newmu = 1/np.sqrt(1 + xi)

            a = np.sqrt(gamma*R*self.gas.T)
            if vel > a:
                self.pressuredrop -= 0.05e5
                print("flow velocity exceeding speed of sound, reducing design pressure drop to ", self.pressuredrop/1e6, " MPa")

            it += 1
            if it > maxiter:
                raise ValueError("Not converged after ", maxiter, " iterations")
            difference = abs(mu - newmu)
            mu = newmu

        self.mu = mu
        self.diameter = diameter
        self.velocity = vel


class AnnularOrifice():
    """
    annulus volume flow based on Poiseuille Flow
    in practice a little fucky so dont trust blindly 
    """
    def __init__(self, fluid, mixture, temperature, pressure, length, pressuredrop, massflow, r_inner):
        self.fluid = thermo.Mixture(fluid, ws=mixture, T=temperature, P=pressure)
        self.length = length
        self.pressuredrop = pressuredrop
        self.r_inner = r_inner
        self.massflow = massflow
   
    def injector(self):
        Q = self.massflow / self.fluid.rho
        G = self.pressuredrop / self.length

        func = lambda r_outer: G*np.pi / (8*self.fluid.mu) * (r_outer**4 - self.r_inner**4 - (r_outer**2 - self.r_inner**2)**2/np.log(r_outer/self.r_inner)) - Q
        r_outer = fsolve(func, 1)
        #diameter = (Q*6*self.fluid.mu*self.length / (np.pi*radius*self.pressuredrop*5/2))**(1/3)
        self.diameter = r_outer - self.r_inner

        A = np.pi*(r_outer**2 - self.r_inner**2)
        self.mu = self.massflow / (A * np.sqrt(2*self.fluid.rho*self.pressuredrop))
        self.velocity = self.massflow / self.fluid.rho / A


class AnnulusInjector():
    def __init__(self, fluid, mixture, temperature, pressure, length, annulusdiameter, massflow, pressuredrop):
        """UNVERIFIED model of annular injectors. Uses numerical solution to Darcey Weisbach eqaution to determine discharge coefficent with cryo test data 

        :param fluid: fluid flowing through annulus
        :type fluid: string
        :param length: length of annulus flow passage (coaxial distance)
        :param annulusdiameter: mean diameter of annulus
        :param pressuredrop: design pressure drop over injector 
        """        
        self.fluid = thermo.Mixture(fluid, ws=mixture, T=temperature, P=pressure)
        self.annulusdiameter = annulusdiameter 
        self.length = length
        self.massflow = massflow
        self.pressuredrop = pressuredrop

    def friction(self, Re, di):
        # solving Darcey Weisbach equation
        surface_roughness = 0
        fd = fsolve(lambda f: 1/(np.sqrt(f)) + 2*np.log10(surface_roughness/(3.7*di) + 2.51/(Re*np.sqrt(f))), 0.0000001)
        return fd[0]

    def xiinlet(self):
        #return 0.5 + 1.2/np.pi*self.inletangle             # min 0.5 for coaxial flow before injection, max 0.9 for flow at pi/6 rad realtive to faceplate 
        return 1.662                                        # fit from cryo test data

    def injector(self, maxiter=100, tol=1e-6):
        it = 0 
        difference = 1 
        mu = 1/np.sqrt(1+self.xiinlet())                               
        di = 1e-3
        while difference > tol:
            D = self.annulusdiameter + di

            vel = mu*np.sqrt(2*self.pressuredrop/self.fluid.rho)
            Re = self.fluid.rho*vel*di/self.fluid.mu

            di = self.massflow/(self.fluid.rho*np.pi*D*vel)

            xifriction = self.friction(Re, di)*self.length/di
            xiinlet = self.xiinlet()
            newmu = 1/np.sqrt(1+xifriction+xiinlet)

            it += 1
            if it > maxiter:
                raise ValueError("Not converged after ", maxiter, " iterations")
            difference = abs(mu - newmu)
            mu = newmu  

        #print(np.pi*D*di*newmu*np.sqrt(self.fluid.rho*2*self.pressuredrop))
        
        self.D = D
        self.mu = mu
        self.diameter = di
        self.velocity = vel


def annulus_verification(r_inner, outer_radius, fluid, pressure_range, temperature, discharge_coefficient):
    volumeflow = []
    D = 2*((outer_radius - r_inner)/2 + r_inner)
    di = outer_radius - r_inner
    for p in pressure_range:
        liquid = thermo.Chemical(fluid, T=temperature, P=p)
        vel = discharge_coefficient*np.sqrt(2*p/liquid.rho)

        q = 2*np.pi*D*vel*(di/2)
        q *= 60/0.001                               # conversion to l/min

        volumeflow.append(q)
    volumeflow = np.array(volumeflow)
        

    def experimental_volumeflow(pressuredrop):
        deltap = pressuredrop/1e5
        q_exp = (deltap/0.0012)**(1/2.1171)         # cryo test data in l/min
        return q_exp

    exp_volumeflow = experimental_volumeflow(pressure_range)
    error = np.dot((exp_volumeflow-volumeflow), (exp_volumeflow-volumeflow))/np.dot(exp_volumeflow,exp_volumeflow)
    print('error: ', abs(error))

    idx = np.where(pressure_range == 12e5)
    e = exp_volumeflow[idx] - volumeflow[idx]
    print(e)
    
    plt.plot(experimental_volumeflow(pressure_range), pressure_range/1e5, color = 'red', label = 'experimental volumeflow')
    plt.plot(volumeflow, pressure_range/1e5, color = 'blue', label = 'calculated volumeflow')
    #plt.loglog()
    plt.grid()
    plt.xlabel('volume flow [l/min]')
    plt.ylabel('pressure drop [bar]')
    plt.legend(loc='best')
    plt.show()


def ohnesorge_number(injector, n_holes, massflow, fluid_mass_fractions, fluid_mole_fractions):
    ohnesorge_nr = []
    massflows = massflow/n_holes
    for m in massflows:
        injector.massflow = m
        injector.injector()                 # calculate new injector properties based on new pressure drop 
        surface_tension = injector.fluid.SurfaceTensionMixture(ws=fluid_mass_fractions, zs=fluid_mole_fractions, T=injector.fluid.T, P=injector.fluid.P)
        We = injector.fluid.rho*injector.velocity**2*injector.diameter / surface_tension
        Re = injector.fluid.rho*injector.velocity*injector.diameter / injector.fluid.mu
        ohnesorge_nr.append(np.sqrt(We)/Re)

    plt.plot(n_holes, ohnesorge_nr)
    plt.grid()
    plt.xlabel('number of holes [-]')
    plt.ylabel('Weber Number [-]')
    plt.show()


def gas_liquid_tangential_injection(injector, rho_gas, m_dot_gas, chamber_diameter):
    '''
    based on https://doi.org/10.1016/j.ijmultiphaseflow.2022.104265 
    determines penetration length of liquid jet perpendicular to core gas flow 
    Assumes chamber centerline as penettration height and calculates penetration length
    for that y coordinate
    INPUTS 
        injector [injector object]: 
        x_coordiate [float]:             downstream length
        rho_gas [float]:                 combustion gas density 
        m_dot_gas [float]:               gas mass flow
        chamber_diameter [float]:        combustion chamber diameter         
    '''
    # momentum flux ratio J
    A_chamber = np.pi * chamber_diameter**2 / 4
    v_gas = m_dot_gas / (rho_gas * A_chamber)
    J = injector.fluid.rho * injector.velocity**2 / (rho_gas * v_gas**2)

    # penetration depth
    x = chamber_diameter * (0.5/3.14 * J**(-0.41))**(1/0.27)

    return x



if __name__ == '__main__':

    n_holes_fuel = 1
    dp = 0.1e5
    m_dot_fuel = 0.3e-3
    injector_length = 1e-3
    upstream_pressure = 9e5
    upstream_temperature = 293 
    injection_angle = 0                            

    liq_inj = LiquidInjector(['h2o2'], [1], upstream_temperature, upstream_pressure, injector_length, m_dot_fuel/n_holes_fuel, dp, injection_angle)
    liq_inj.injector()
   # print(liq_inj.fluid.rho)
   # print(liq_inj.diameter*1000)
    #print(liq_inj.velocity)
    #print(liq_inj.mu)

    # core gas flow properties
    rho_gas = 2.337             # kg/m^3
    chamber_diameter = 20e-3    # m
    m_dot_gas = 21.73e-3           # kg/s

    injection_length = gas_liquid_tangential_injection(liq_inj, rho_gas, m_dot_gas, chamber_diameter)
    #print(injection_length * 1000)


    inj = GasInjector('o2', 293.15, 18e5, 3e-3, 14.136e-3, 15e5, 10, 0)
    inj.injector()
    print(liq_inj.diameter*1000)