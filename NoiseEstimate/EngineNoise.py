import numpy as np
from matplotlib import pyplot as plt
import os
from PIL import Image

from scipy.optimize import bisect


class Images():
    def __init__(self, noise_source, safety_area, image_scaling, location, folder, save_fig=True, name="Noise"):
        self.source   = noise_source                          # object with properties noise, radius and angle from the source
        self.area     = safety_area                           # filename of image of your field
        self.scaling  = image_scaling                         # ratio of pixels to m in image of safety area 
        self.loc      = location                              # location of noise source on safety area in pixels [x,y]
        self.folder   = folder                                # directory for saving output images 
        self.save_fig = save_fig
        self.name     = name

        try:
            os.mkdir(folder) 
        except:
            pass

    def contour_plot(self):   
        r, theta = np.meshgrid(self.source.radii, self.source.angles)   
        
        fig, ax  = plt.subplots(subplot_kw = dict(projection='polar'))
        cp       = ax.contourf(theta, r, self.source.sound_levels)
        fig.colorbar(cp, label='noise [dB]')
        
        if self.save_fig:
            fig.savefig(self.folder + '/' + self.name + '.jpg', dpi = 600)
        else:
            plt.show()
    
    def overlay(self):
        noise = Image.open(self.folder + '/' + self.name + '.jpg')
        area  = Image.open(self.area)
        #noise = noise.rotate(90)

        scaling = 18.4/37 
        noise = noise.resize([int(scaling * s) for s in noise.size])
        x_n, y_n = noise.size
        x_a,y_a = 978, 1629

        area.paste(noise, (int(x_a - x_n/2),int(y_a - y_n/2)))
        area.show()
                
            
class PressureVessel():
    def __init__(self, burst_pressure, tank_volume, atmospheric_pressure, n_grid=100, max_distance=30):
        self.P_burst      = burst_pressure                                   # predicted burst presssure at minimum safety factor [Pa]
        self.V_tank       = tank_volume                                      # volume of the pressure vessel [m^3]
        self.P_atm        = atmospheric_pressure                             # local atmospheric pressure [Pa]
        self.n_grid       = n_grid                                           # grid points
        self.max_distance = max_distance                                     # maximum resolved distance from the pressure vessel [m]
        self.gamma        = 1.4                                              # ratio of specific heats of air
        self.E_tnt        = 4.148e6                                          # Energy of TNT [J/kg]

    def TNT_equivalent(self):
        # TNT equivalence of a pressure vessel assuming isentropic expansion
        # Energy release
        self.E_burst = self.P_burst * self.V_tank / (self.gamma - 1) * (1 - (self.P_atm / self.P_burst)**((self.gamma - 1) / self.gamma))

        # TNT equivalent mass
        self.M_tnt = self.E_burst / self.E_tnt

    def noise_distribution(self): 
        self.angles = np.linspace(0, 2*np.pi, self.n_grid)
        self.radii  = np.linspace(0, self.max_distance, self.n_grid)
        self.sound_levels = np.ndarray((self.n_grid, self.n_grid))

        
        # pressure rise and noise power
        for i in range(round(self.n_grid / 2)):
            for j in range(self.n_grid): 
                
                if self.radii[j] < 1:
                    # set anything below 1m distance to 1m value
                    dP = 0.95 * self.M_tnt**(1/3) + 3.9 * self.M_tnt**(2/3) + 13 * self.M_tnt
                else:
                    dP = 0.95 * self.M_tnt**(1/3)/self.radii[j] + 3.9 * self.M_tnt**(2/3)/self.radii[j]**2 + 13 * self.M_tnt/self.radii[j]**3
                
                self.sound_levels[i][j] = 20 * np.log(dP/20e-6)
                self.sound_levels[self.n_grid-1-i][j] = self.sound_levels[i][j]
        print(20 * np.log((0.95 * self.M_tnt**(1/3) + 3.9 * self.M_tnt**(2/3) + 13 * self.M_tnt)/20e-6))


class ExplosiveCharge():
    def __init__(self, NEM, TNT_equivalence=0.43, n_grid=100, max_distance=30):
        self.TNT_equivalence = TNT_equivalence                                  # TNT equivelnce factor for pressure. Black powder = 0.43
        self.NEM             = NEM                                              # net explosive mass [kg]
        self.n_grid           = n_grid                          # grid points
        self.max_distance     = max_distance                    # maximum resolved distance from the motor [m]
        self.gamma           = 1.4                                              # ratio of specific heats of air
        self.E_tnt           = 4.148e6                                          # Energy of TNT [J/kg]
        
    def noise_distribution(self): 
        self.angles = np.linspace(0, 2*np.pi, self.n_grid)
        self.radii  = np.linspace(0, self.max_distance, self.n_grid)
        self.sound_levels = np.ndarray((self.n_grid, self.n_grid))

        # NEM in TNT equivalent
        self.M_tnt = self.NEM * self.TNT_equivalence
        
        # pressure rise and noise power
        for i in range(round(self.n_grid / 2)):
            for j in range(self.n_grid): 
                
                if self.radii[j] < 1:
                    # set anything below 1m distance to 1m value
                    dP = 0.95 * self.M_tnt**(1/3) + 3.9 * self.M_tnt**(2/3) + 13 * self.M_tnt
                else:
                    dP = 0.95 * self.M_tnt**(1/3)/self.radii[j] + 3.9 * self.M_tnt**(2/3)/self.radii[j]**2 + 13 * self.M_tnt/self.radii[j]**3
                
                self.sound_levels[i][j] = 20 * np.log(dP/20e-6)
                self.sound_levels[self.n_grid-1-i][j] = self.sound_levels[i][j]


class RocketMotor():
    def __init__(self, thrust, exhaust_velocity, scaling_factor=0.0012, n_grid=100, max_distance=30):
        self.thrust           = thrust
        self.exhaust_velocity = exhaust_velocity
        self.scaling_factor   = scaling_factor                  # scaling factor for noise power in relation to engine power. 0.0012, determined experimentally
        self.n_grid           = n_grid                          # grid points
        self.max_distance     = max_distance                    # maximum resolved distance from the motor [m]

    def noise_power(self):
        self.P_noise       = 0.5 * self.scaling_factor * self.exhaust_velocity * self.thrust

        print(0.5 * self.exhaust_velocity * self.thrust * 0.00135962)
        self.max_noise_lvl = 10 * np.log10(self.P_noise * 1e12)

    def noise_distance(self, distance):
        return 10 * np.log10(self.P_noise * 1e12 / (2 * np.pi * distance**2))

    def noise_direction(self, theta, radius):
        # NASA noise model
        angular_dist =  -2.20411*theta**4 + 19.43914*theta**3 - 59.41912*theta**2 + 64.91368*theta - 18.73118 
        noise_dist = angular_dist + 10 * np.log(self.P_noise * 2.5e11 / (np.pi * radius**2)) / np.log(10)
        return noise_dist

    def noise_distribution(self):
        self.angles = np.linspace(0, 2*np.pi, self.n_grid)
        self.radii  = np.linspace(0, self.max_distance, self.n_grid)
        self.sound_levels = np.ndarray((self.n_grid, self.n_grid))

        # Noise model undefined for a distance of 0 m 
        self.sound_levels[:][0] = np.ones(self.n_grid) * self.max_noise_lvl
        self.sound_levels[:,0]  = np.ones(self.n_grid) * self.max_noise_lvl

        for i in range(round(self.n_grid / 2)):
            for j in range(1,self.n_grid): 
                
                self.sound_levels[i][j] = self.noise_direction(self.angles[i], self.radii[j])
                self.sound_levels[self.n_grid-1-i][j] = self.sound_levels[i][j]

        #print(self.noise_direction(0, 1))
        #print(self.noise_direction(np.pi/2, 2))

def validation():
    # tuning of the model constant using test data of three solid rocket motors
    # test data
    thrust           = np.array([398, 524, 683])                        # N
    exhaust_velocity = 1295                                             # m/s
    r_probe          = 20                                               # m
    theta_probe      = np.pi/2  
    measured_noise   = np.array([105.6, 107.4, 109.7])                  # dB

    constant         = np.ndarray(len(thrust))

    for i in range(len(thrust)):
        def func(const):
            sound = RocketMotor(thrust[i], exhaust_velocity, scaling_factor=const)
            sound.noise_power()

            return measured_noise[i] - sound.noise_direction(theta_probe, r_probe)
        
        constant[i] = bisect(func, 1e-4, 1e-2)

    print('mean sclaing factor from test data: ', np.mean(constant))
            

if __name__ == "__main__":
    thrust = 10000                      # [N]
    exhaust_velocity = 3000             # [m/s]

    sound = RocketMotor(thrust, exhaust_velocity, scaling_factor=0.0012)
    sound.noise_power()
    sound.noise_distribution()
    
    #validation()

    #sound = PressureVessel(255e5, 0.000516, 1e5)
    #sound.TNT_equivalent()
    #sound.noise_distribution()


    #sound = ExplosiveCharge(0.01, 0.43)
    #sound.noise_distribution()

    im = Images(sound, 'FellowshipField.jpg', 1, 1, 'Figures', name='SparrowNoise')
    im.contour_plot()

    