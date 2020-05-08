import numpy as np

def get_continuum_regions(data):
    rho = data.equil_data['rho0']
    B02 = data.equil_data['B02']
    B03 = data.equil_data['B03']
    B0 = np.sqrt(B02**2 + B03**2)
    v02 = data.equil_data['v02']
    v03 = data.equil_data['v03']
    T = data.equil_data['T0']
    p = rho * T

    k2 = data.params['k2']
    k3 = data.params['k3']
    gamma = data.gamma

    # determine scale factor
    if data.geometry == 'cylindrical':
        eps = data.grid_gauss
    else:
        eps = np.ones_like(data.grid_gauss)

    # Alfven and slow continuum (squared)
    wA2 = (1 / rho) * ((k2 * B02 / eps) + k3 * B03)**2
    wS2 = (gamma * p / (gamma * p + B0**2)) * wA2
    # doppler shift equals dot product of k and v
    doppler = k2 * v02 + k3 * v03

    # retrieve thermal continuum
    wTH = _thermal_continuum(data, eps)

    # additional sanity check: if the slow continuum vanishes (all wS2 values are zero) there
    # is an analytical solution for the thermal continuum. This one should be equal to the solution
    # obtained through solving the polynomial equation.
    if (wS2 == 0).all() and (T != 0).all():
        LT = data.equil_data['dLdT']
        Lrho = data.equil_data['dLdrho']
        T = data.equil_data['T0']
        vA2 = B0 ** 2 / rho
        wTH_analytical = 1j * (gamma - 1) * (rho * Lrho / T - LT * (T + vA2)) / (vA2 + gamma * T)
        wTH_difference = abs(np.sort(wTH) - np.sort(np.imag(wTH_analytical)))
        # if any value in the difference is larger than 1e-12, simply print a warning
        if (wTH_difference > 1e-12).any():
            print("WARNING: slow continuum vanishes, but difference in analytical/numerical thermal continuum!")

    # calculate slow and Alfvén continuum, take doppler shift into account
    wS_pos = doppler + np.sqrt(wS2)
    wS_neg = doppler - np.sqrt(wS2)
    wA_pos = doppler + np.sqrt(wA2)
    wA_neg = doppler - np.sqrt(wA2)
    return wS_pos, wS_neg, wA_pos, wA_neg, wTH

def _thermal_continuum(data, eps):
    rho = data.equil_data['rho0']
    B02 = data.equil_data['B02']
    B03 = data.equil_data['B03']
    B0 = np.sqrt(B02**2 + B03**2)
    T = data.equil_data['T0']
    p = rho * T
    LT = data.equil_data['dLdT']
    Lrho = data.equil_data['dLdrho']
    kappa_para = data.equil_data['kappa_para']
    kappa_perp = data.equil_data['kappa_perp']

    # if there is no conduction/cooling, there is no thermal continuum. Set to zero and return
    if (LT == 0).all() and (Lrho == 0).all() and (kappa_para == 0).all() and (kappa_perp == 0).all():
        return np.zeros_like(data.grid_gauss)
    # if temperature is zero (no pressure), set to zero and return
    if (T == 0).all():
        return np.zeros_like(data.grid_gauss)

    k2 = data.params['k2']
    k3 = data.params['k3']
    # k parallel and perpendicular to B, uses vector projection and scale factor eps
    kpara = (k2 * B02 / eps + k3 * B03) / B0
    kperp = (k2 * B03 / eps - k2 * B02) / B0
    gamma = data.gamma

    cs2 = gamma * p / rho       # sound speed
    vA2 = B0**2 / rho           # Alfvén speed
    ci2 = p / rho               # isothermal sound speed
    sigma_A2 = kpara * vA2                      # Alfvén frequency
    sigma_c2 = cs2 * sigma_A2 / (vA2 + cs2)     # cusp frequency
    sigma_i2 = ci2 * sigma_A2 / (vA2 + ci2)     # isothermal cusp frequency

    # thermal and slow continuum are coupled in a third degree polynomial in omega,
    # coeffi means the coefficient corresponding to the term omega^i
    coeff3 = rho * (cs2 + vA2) * 1j / (gamma - 1)
    coeff2 = -(kappa_para * kpara**2 + rho * LT) * (ci2 + vA2) + rho**2 * Lrho
    coeff1 = -rho * (cs2 + vA2) * sigma_c2 * 1j / (gamma - 1)
    coeff0 = (kappa_para * kpara**2 + rho * LT) * (ci2 + vA2) * sigma_i2 - rho**2 * Lrho * sigma_A2

    # we have to solve this equation "gauss_gridpts" times.
    # thermal continuum corresponds to purely imaginary solution, slow continuum are two real solutions
    wTH = []
    for idx in range(len(data.grid_gauss)):
        solutions = np.roots([coeff3[idx], coeff2[idx], coeff1[idx], coeff0[idx]])
        # retrieve the imaginary parts, only those with a real part smaller than 1e-12 (numerical reasons)
        # to fetch the purely imaginary solutions
        imag_sol = solutions.imag[abs(solutions.real) < 1e-12]
        if (imag_sol == 0).all():
            # if all solutions are zero (no thermal continuum), then append zero
            w = 0
        else:
            # else there should be ONE value that is nonzero, this is the thermal continuum value.
            # Filter out non-zero solution and try to unpack. If there is more than one non-zero
            # value unpacking fails, this is a sanity check so we raise an error.
            try:
                w, = imag_sol[imag_sol != 0]        # this is an array, so unpack with ,
            except ValueError:
                raise ValueError("Found more than one solution for the thermal continuum!")
        wTH.append(w)
    return np.asarray(wTH)
