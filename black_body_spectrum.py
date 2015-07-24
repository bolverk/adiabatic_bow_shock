def spectrum_integrate(area_list, temperature_list, frequency_list, constants):

    import numpy

    c = constants['speed of light']
    h = constants['planck constant']
    k = constants['boltzmann constant']

    res = numpy.zeros(frequency_list.size)
    for A, T in zip(area_list, temperature_list):
        res += A*((2*h*frequency_list**3)/c**2)*1./(numpy.exp(h*frequency_list/(k*T))-1)
    return res

def sun_spectrum():

    import pylab
    import numpy

    constants = {'speed of light':3e10,
                 'planck constant':6.626e-27,
                 'boltzmann constant':1.38e-16,
                 'solar radius':6.96e10}
    frequency_list = numpy.logspace(1,20,num=100)
    spectrum = spectrum_integrate([4*numpy.pi*constants['solar radius']**2],
                                  [5000],
                                  frequency_list,
                                  constants)
    pylab.loglog(frequency_list*constants['planck constant'], spectrum)
    pylab.xlabel('Energy [erg]')
    pylab.ylabel(r'$F_\nu$ [erg/s/Hz]')
    pylab.show()

def mid_array(a):

    import numpy

    b = numpy.array(a)
    return 0.5*(b[1:]+b[:-1])

def snr_knot():

    import pylab
    import numpy

    constants = {'speed of light':3e10,
                 'planck constant':6.626e-27,
                 'boltzmann constant':1.38e-16,
                 'solar radius':6.96e10,
                 'proton mass':1.67e-24,
                 'adiabatic index':5./3.,
                 'wind velocity':1000*1e3*1e2,
                 'electron volt':1.6e-12}
    constants['clump radius'] = 1e-2*constants['solar radius']
    def temperature_profile(r):
        g = constants['adiabatic index']
        m = constants['proton mass']
        v = constants['wind velocity']
        k = constants['boltzmann constant']
        R = constants['clump radius']
        return 2*(g-1)*m*v**2/(k*(1+(r/R)**2)*(g+1)**2)
    edge_radius_list = constants['clump radius']*numpy.logspace(-2,3,num=100)
    mid_radius_list = mid_array(edge_radius_list)
    dr_list = numpy.diff(edge_radius_list)
    area_list = [2*numpy.pi*r*dr*numpy.sqrt(1+(r/constants['clump radius'])**2)
                 for r,dr in zip(mid_radius_list,edge_radius_list)]
    temperature_list = [temperature_profile(r) for r in mid_radius_list]
    frequency_list = numpy.logspace(1,20,num=100)
    spectrum = spectrum_integrate(area_list,
                                  temperature_list,
                                  frequency_list,
                                  constants)
    pylab.loglog(frequency_list*constants['planck constant']/constants['electron volt'], spectrum)
    pylab.xlabel('Photon Energy [eV]')
    pylab.ylabel(r'$F_\nu$ [erg/s/Hz]')
    pylab.show()

if __name__ == '__main__':

    snr_knot()
