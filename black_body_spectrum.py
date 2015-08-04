def cache_eval(func, fname):

    import os
    import pickle

    if os.path.isfile(fname):
        return pickle.load(open(fname,'rb'))
    else:
        res = func()
        pickle.dump(res,open(fname,'wb'))
        return res

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

class RightRectrangle:

    def __init__(self,x_min, x_max, y_min, y_max):

        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max

    def __call__(self, x, y):

        return (self.x_max > x > self.x_min) and (self.y_max > y > self.y_min)

def calc_max_temperature(g=None,m=None,v=None,k=None):
        
    return 2*(g-1)*m*v**2/k/(g+1)**2

def real_deal():

    import pickle
    import numpy
    import pylab

    t_front = pickle.load(open('bb_t_front.pkl','rb'))
    a_front = pickle.load(open('bb_a_front.pkl','rb'))

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
    scaled['areas'] = numpy.array(a_front)*constants['solar radius']**2
    scaled['temperature'] = numpy.array(t_front)*max_temperatue/numpy.max(t_front)
    frequency_list = numpy.logspace(11,20,num=100)
    spectrum = spectrum_integrate(scaled['areas'],
                                  scaled['temperature'],
                                  frequency_list,
                                  constants)
    pylab.loglog(frequency_list*constants['planck constant']/constants['electron volt'],
                 spectrum)
    pylab.show()
        
if __name__ == '__main__':

    #snr_knot()
    real_deal()
