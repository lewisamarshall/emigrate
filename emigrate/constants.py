"""Constants for use in the Emigrate solver."""

# Physical Constants
k_water = 1E-14
lpm3 = 1000
faraday = 96485.34          # Faraday's const.[C/mol]
boltzmann = 8.617e-6        # EV/K
temperature = 25
temperature_K = temperature + 273.15
h_mobility = 362E-9/faraday   # Mobility of Hydroxide   % [m^2/s*V]/F --> [mol*s/Kg]
oh_mobility = 205E-9/faraday  # Mobility of Hydronium   % [m^2/s*V]/F --> [mol*s/Kg]
h_diffusivity = h_mobility / 1 * boltzmann * (temperature_K)
oh_diffusivity = oh_mobility / -1 * boltzmann * (temperature_K)
viscosity = 1E-3            # Dynamic viscosity of water [Pa s]

# Rmu = 8.31       # Universal gas const. [J/mol*K]
