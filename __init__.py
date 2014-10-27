import numpy as np
import ionize
import scipy.integrate as integrate
import scipy.sparse as sparse


class Electrophoresis(object):
    import constants
    V = 100
    E = 1
    dz = None
    z = None
    x = None
    A = 0
    B = 0
    diffusivity = None
    mobility = None
    molar_conductivity = None
    j = 0
    solution = None
    solver_info = None

    def __init__(self, domain, ions, concentrations):
        self.x = np.array(domain)
        self.z = np.array(self.x[:])
        self.set_dz()
        self.ions = ions
        self.set_ion_properties()
        # self.set_derivative_matrices()
        self.concentrations = np.array(concentrations)
        self.t = 0

    def first_derivative(self, x_input, method='dissipative'):
        if method is None:
            derivative = x_input

        elif method == 'dissipative':
            derivative = []

            derivative = np.pad(np.diff(x_input, n=1, axis=1),
                                ((0, 0), (0, 1)), 'reflect') / \
                np.tile(self.dz, (len(self.ions), 1))/2 -\
                np.pad(np.diff(x_input, n=2, axis=1),
                       ((0, 0), (1, 1)), 'reflect')/2 / \
                np.tile(self.dz, (len(self.ions), 1))

        elif method == '6th-Order':
            pass

        return derivative

    def second_derivative(self, x_input, method='dissipative'):
        if method is None:
            derivative = input
        elif method == 'dissipative':
            derivative = self.first_derivative(self.first_derivative(x_input))
        elif method == '6th-Order':
            pass
        return derivative

    def set_ion_properties(self):
        pH = 7
        self.diffusivity = np.array([[ion.diffusivity(pH)]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(pH)]
                                  for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(pH)]
                                            for ion in self.ions])

    def set_dz(self):
        self.dz = np.pad(np.diff(self.z), (0, 1), 'reflect')

    def set_current(self, concentrations):
        self.j = self.V/sum(self.dz/self.conductivity(concentrations))

    def conductivity(self, concentrations):
        conductivity = np.sum(np.tile(self.molar_conductivity,
                                      (1, len(self.dz)))
                              * concentrations, 0)
        return conductivity

    def set_E(self, concentrations):
        self.set_current(concentrations)
        self.E = self.j/self.conductivity(concentrations)

    def flux(self, concentrations):
        self.set_E(concentrations)
        diffusion = \
            self.second_derivative(np.tile(self.diffusivity,
                                           (1, len(self.z)))
                                   * concentrations
                                   )
        advection = \
            -self.first_derivative(np.tile(self.mobility,
                                           (1, len(self.z)))
                                   * concentrations *
                                   self.E
                                   )

        total_flux = diffusion + advection

        return total_flux

    def reshaped_flux(self, concentrations, t):
        if not t == self.t:
            self.calc_equilibrium()
        concentrations = concentrations.reshape(self.concentrations.shape)
        flux = self.flux(concentrations).flatten()
        return flux

    def solve(self, t):
        self.solution, self.solver_info =\
            integrate.odeint(self.reshaped_flux,
                             self.concentrations.flatten(),
                             t,
                             full_output=True)
        self.solution = [sol.reshape(self.concentrations.shape)
                         for sol in self.solution]

    def calc_equilibrium(self):
        pass

    def set_derivative_matrices(self):
        h = self.dz[0]
        N = len(self.z)
        # A = sparse(N, N)
        # B = sparse(N, N)

        # # Forward derivative for first(/last) grid point - 5th order
        # A_First = [1, 4]
        # B_First = (1./h) * [-37/12    2/3    3   -2/3    1/12]
        #
        # # Skewed derivative for second/(N-1) grid point - 5th order
        # A_Second = [1/6 1 1/2]
        # B_Second = (1/h)*[-10/18   -1/2    1    1/18]
        #
        # # Central derivative for internal grid points
        # alphaI = 1/3
        # A_Internal = [alphaI 1 alphaI]
        # B_Internal = [-1/3*(4*alphaI-1)*(1/(4*h))  -2/3*(alphaI+2)*(1/(2*h))
        #  0  2/3*(alphaI+2)*(1/(2*h))  1/3*(4*alphaI-1)*(1/(4*h))]
        #
        # A(1,1:length(A_First)) = A_First
        # B(1,1:length(B_First)) = B_First
        # A(2,1:length(A_Second)) = A_Second
        # B(2,1:length(B_Second)) = B_Second;
        #
        # A(end,end:-1:end-length(A_First)+1) = A_First
        # B(end,end:-1:end-length(B_First)+1) = B_First
        # A(end-1,end:-1:end-length(A_Second)+1) = A_Second
        # B(end-1,end:-1:end-length(B_Second)+1) = - B_Second
        #
        # for ij=3:N-2  # run on all rows
        #     A(ij,ij-1:ij+1)=A_Internal
        #     B(ij,ij-2:ij+2)=B_Internal

        A_vector = 0

        A_constructor = [A_vector, np.ones(self.z.size), np.fliplr(A_vector)]
        B_constructor = 1/h*np.array(0)
        self.A = sparse.spdiags(A_constructor, range(-1, 1+1), N, N)
        self.B = sparse.spdiags(B_constructor, range(-4, 4+1), N, N)

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    my_ions = ionize.Solution(['tris', 'hydrochloric acid', 'caproic acid'],
                              [0, 0, 0]
                              ).ions

    domain_length = 0.1
    interface_length = 0.01
    nodes = 100
    my_domain = np.linspace(-domain_length/2., domain_length/2., nodes)
    my_concentrations = np.array([np.ones(my_domain.shape)*.1,
                                 0.05-0.05*erf(my_domain/interface_length),
                                 0.05*erf(my_domain/interface_length)+0.05])
    my_elec = Electrophoresis(my_domain, my_ions, my_concentrations)
    my_elec.solve(np.array(np.linspace(0, 1.5e2, 10)))
    for my_sol in my_elec.solution:
        for sub_sol in my_sol:
            plot.plot(my_elec.z, sub_sol)
    plot.show()
