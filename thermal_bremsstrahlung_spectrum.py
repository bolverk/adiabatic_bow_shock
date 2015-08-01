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

def calc_max_temperature(g=None,m=None,v=None,k=None):
        
    return 2*(g-1)*m*v**2/k/(g+1)**2

def calc_high_end_spectrum(R = None, v = None, g = None, f = None, h = None,
                           k = None, me = None, c = None, mp=None, n=None,
                           q = None):
    import math
    T0 = calc_max_temperature(g=g, m=mp, v=v, k=k)
    return (64.*math.sqrt(2.)*math.exp(-f*(g+1)**2*h/(2*(g-1)*mp*v**2))*
            f*(g-1)**2*mp**2*n**2*math.pi**2*q**6*R**3*v**4/
            (3*c**3*(g+1)**4*(f*h)**(9./2.)*me**(3./2.)))*h

def real_deal():

    import glob
    import h5py
    import numpy
    import pylab

    def extract_number(snapshot_file):

        import re

        res = re.match(r'snapshot_([0-9]*).h5',snapshot_file)

        return int(res.group(1))

    snapshot_list = sorted(glob.glob('snapshot_*.h5'),
                           key=extract_number)

    last_snapshot_file = snapshot_list[-1]

    with h5py.File(last_snapshot_file) as f:
        raw = {}
        for field in ['temperature','volumes']:
            raw[field] = numpy.array(f[field])

    constants = {'speed of light':3e10,
                 'planck constant':6.626e-27,
                 'boltzmann constant':1.38e-16,
                 'solar radius':6.96e10,
                 'bremsstrahlung prefactor':6.8e-38,
                 'kelvin':1,
                 'electron volt':1.6e-12,
                 'wind velocity':1000*1000*1e2,
                 'adiabatic index':5./3.,
                 'proton mass':1.67e-24,
                 'electron mass':9.1e-28}
    max_temperatue = calc_max_temperature(g = constants['adiabatic index'],
                                          m = constants['proton mass'],
                                          v = constants['wind velocity'],
                                          k = constants['boltzmann constant'])
    scaled = {}
    scaled['volumes'] = raw['volumes']*constants['solar radius']**3
    scaled['temperature'] = raw['temperature']*max_temperatue/numpy.max(raw['temperature'])
    frequency_list = numpy.logspace(1,20,num=100)
    high_end_spectrum = map(lambda x: calc_high_end_spectrum(R=0.01*constants['solar radius'],
                                                             v = constants['wind velocity'],
                                                             g = constants['adiabatic index'],
                                                             f = x,
                                                             h = constants['planck constant'],
                                                             k = constants['boltzmann constant'],
                                                             mp = constants['proton mass'],
                                                             me = constants['electron mass'],
                                                             c = constants['speed of light'],
                                                             n = 1,
                                                             q = 4.8e-10),
                            frequency_list)
    spectrum = spectrum_integrate(scaled['volumes'],
                                  scaled['temperature'],
                                  frequency_list,
                                  constants)
    pylab.loglog(frequency_list*constants['planck constant']/constants['electron volt'],
                 spectrum,
                 frequency_list*constants['planck constant']/constants['electron volt'],
                 high_end_spectrum)
    pylab.xlabel('Photon Energy [eV]')
    pylab.ylabel(r'$F_{\nu}$ [erg/s/Hz]')
    pylab.show()

if __name__ == '__main__':

    #snr_knot()
    real_deal()
        
