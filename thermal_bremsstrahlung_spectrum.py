def spectrum_integrate(volume_list, temperature_list, frequency_list, constants):

    import numpy

    C = constants['bremsstrahlung prefactor']
    kelvin = constants['kelvin']
    h = constants['planck constant']
    kb = constants['boltzmann constant']

    res = numpy.zeros(frequency_list.size)
    for V, T in zip(volume_list, temperature_list):
        res += V*C*(T/kelvin)**(-1/2)*numpy.exp(-h*frequency_list/kb/T)
    return res

def sun_spectrum():

    import numpy
    import pylab

    constants = {'speed of light':3e10,
                 'planck constant':6.626e-27,
                 'boltzmann constant':1.38e-16,
                 'solar radius':6.96e10,
                 'bremsstrahlung prefactor':6.8e-38,
                 'kelvin':1,
                 'electron volt':1.6e-12}
    frequency_list = numpy.logspace(1,20,num=100)
    spectrum = spectrum_integrate([constants['solar radius']**3],
                                  [5000],
                                  frequency_list,
                                  constants)
    pylab.loglog(frequency_list*constants['planck constant']/constants['electron volt'],
                 spectrum)
    pylab.xlabel('Photon Energy [eV]')
    pylab.ylabel(r'$F_{\nu}$ [erg/s/Hz]')
    pylab.show()

def mid_array(a):

    import numpy

    b = numpy.array(a)
    return 0.5*(b[1:]+b[:-1])

def snr_knot():

    import numpy
    import pylab

    constants = {'speed of light':3e10,
                 'planck constant':6.626e-27,
                 'boltzmann constant':1.38e-16,
                 'solar radius':6.96e10,
                 'bremsstrahlung prefactor':6.8e-38,
                 'kelvin':1,
                 'electron volt':1.6e-12,
                 'wind velocity':1000*1000*1e2,
                 'adiabatic index':5./3.,
                 'proton mass':1.67e-24}
    constants['clump radius'] = 1e-2*constants['solar radius']
    def temperature_profile(r):
        g = constants['adiabatic index']
        m = constants['proton mass']
        v = constants['wind velocity']
        k = constants['boltzmann constant']
        R = constants['clump radius']
        return 2*(g-1)*m*v**2/(k*(1+(r/R)**2)*(g+1)**2)
    frequency_list = numpy.logspace(1,20,num=100)
    edge_radius_list = constants['clump radius']*numpy.logspace(-2,3,num=100)
    mid_radius_list = mid_array(edge_radius_list)
    dr_list = numpy.diff(edge_radius_list)
    volume_list = numpy.pi*mid_radius_list**2*numpy.sqrt(1+(mid_radius_list/constants['clump radius'])**2)*dr_list
    temperature_list = [temperature_profile(r) for r in mid_radius_list]
    spectrum = spectrum_integrate(volume_list,
                                  temperature_list,
                                  frequency_list,
                                  constants)
    pylab.loglog(frequency_list*constants['planck constant']/constants['electron volt'],
                 spectrum)
    pylab.xlabel('Photon Energy [eV]')
    pylab.ylabel(r'$F_{\nu}$ [erg/s/Hz]')
    pylab.show()

if __name__ == '__main__':

    snr_knot()
        
