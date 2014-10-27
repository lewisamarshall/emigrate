import numpy as np
import ionize
import scipy.integrate as integrate

class Electrophoresis(object):
    # from constants import *
    V = 1
    E = 1
    dz = None
    z = None
    x = None


    def __init__(self, domain, ions, concentrations):
        self.x = np.array(domain)
        self.z = np.array(self.x[:])
        self.set_dz()
        self.ions = ions
        self.set_ion_properties()
        self.concentrations = np.array(concentrations)
        self.t = 0

    def first_derivative(self, input, method='dissipative'):
        if method == None:
            derivative = input

        elif method == 'dissipative':
            derivative = []

            derivative = np.pad(np.diff(input, n=1, axis=1),
                ((0,0),(0,1)), 'reflect')/\
                    np.tile(self.dz, (len(self.ions),1))/2 -\
                    np.pad(np.diff(input, n=2, axis=1),
                            ((0,0),(1,1)), 'reflect')/2/\
                    np.tile(self.dz, (len(self.ions),1))

        return derivative

    def second_derivative(self, input):
        derivative = self.first_derivative(self.first_derivative(input))
        return derivative

    def set_ion_properties(self):
        pH = 7
        self.diffusivity = np.array([[ion.absolute_mobility[0]]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(pH)]
                                    for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(pH)]
                                    for ion in self.ions])
        #temporary

    def set_dz(self):
        self.dz = np.pad(np.diff(self.z), (0,1), 'reflect')

    def set_current(self, concentrations):
        self.j = self.V/sum(self.dz/self.conductivity(concentrations))

    def conductivity(self, concentrations):
        conductivity = np.sum(
                              np.tile(self.molar_conductivity,
                                      (1, len(self.dz)))
                              * concentrations
                              , 0)
        return conductivity

    def set_E(self, concentrations):
        self.set_current(concentrations)
        self.E = self.j/self.conductivity(concentrations)

    def flux(self, concentrations):
        self.set_E(concentrations)
        diffusion = \
            self.second_derivative(
                                   np.tile(self.diffusivity,
                                           (1, len(self.z)))
                                           * concentrations
                                  )
        advection = \
            -self.first_derivative(np.tile(self.mobility,
                                          (1, len(self.z)))
                                  * concentrations*
                                  self.E
                                 )

        total_flux = diffusion + advection

        return total_flux

    def reshaped_flux(self, concentrations, t):
        concentrations = concentrations.reshape(self.concentrations.shape)
        flux  = self.flux(concentrations).flatten()
        return flux

    def solve(self, t):
        self.solution, self.solver_info =\
         integrate.odeint(
                          self.reshaped_flux,
                          self.concentrations.flatten(),
                          t,
                          full_output=True)
        self.solution = [sol.reshape(self.concentrations.shape) for sol in self.solution]

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    ions = ionize.Solution(['tris', 'hydrochloric acid', 'caproic acid'],
                           [0, 0, 0]
                           ).ions
    domain = np.linspace(0, 10)
    concentrations = np.array([1+erf(domain-5), 0.5-0.5*erf(domain-5), 1-erf(domain-5)])
    my_elec = Electrophoresis(domain,ions, concentrations)
    # print my_elec.concentrations.shape
    # print my_elec.flux(my_elec.concentrations).shape
    # print my_elec.reshaped_flux(my_elec.concentrations.flatten(), 1).shape
    # print my_elec.concentrations.flatten().shape
    # print my_elec.molar_conductivity
    # print my_elec.set_current(my_elec.concentrations)
    # print my_elec.j
    my_elec.solve(np.array(np.linspace(0, 1e7, 10)))
    for sol in my_elec.solution:
        for sub_sol in sol:
            plot.plot(my_elec.z, sub_sol)
    plot.show()
