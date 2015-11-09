import matplotlib
matplotlib.use('Qt4Agg')

def get_user_input():

    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='fits a Pade approximant to the shock front')
    parser.add_argument('file_name',help='file with bow shock points')
    parser.add_argument('out_file',help='name of output file')
    parser.add_argument('--show',
                        help='toggles plotting data',
                        action='store_true')
    return parser.parse_args()

class PadeApproximant:

    def __init__(self, coef_list):
        self.coef_list = coef_list

    def __call__(self, x):

        a = self.coef_list['a']
        b = self.coef_list['b']
        c = self.coef_list['c']
        d = self.coef_list['d']
        e = self.coef_list['e']
        f = self.coef_list['f']
        x2 = x**2
        return (a+b*x2+c*x2**2+d*x2**3)/(1+e*x2+f*x2**2)

def residual(params, x_list, y_list):

    import numpy

    my_func = PadeApproximant(dict((key,params[key].value)
                               for key in ['a','b','c','d','e','f']))
    model = numpy.array([my_func(x) for x in x_list])
    return y_list - model

def main():

    import numpy
    from lmfit import minimize, Parameters

    args = get_user_input()
    data = numpy.loadtxt(args.file_name)
    
    coef_name_list = ['a','b','c','d','e','f']
    params = Parameters()
    for coef in coef_name_list:
        params.add(coef,value=1)

    out = minimize(residual,
                   params,
                   args=(data.T[0], data.T[1]))
    pade_approximant = PadeApproximant(out.values)
    serialized = []
    for itm in coef_name_list:
        serialized.append(out.values[itm])
    numpy.savetxt(args.out_file, serialized)
    if args.show:
        import pylab
        pade_approximant = PadeApproximant(out.values)
        y_reconstructed = numpy.array([pade_approximant(x)
                                       for x in data.T[0]])
        #pylab.subplot(211)
        pylab.plot(data.T[0], data.T[1], linewidth=3, label='numeric')
        pylab.plot(data.T[0], y_reconstructed, linewidth=3,
                   label='pade')
        pylab.xlim((numpy.min(data.T[0]),
                    numpy.max(data.T[0])))
        pylab.legend(loc='best')
        #pylab.gca().xaxis.set_major_formatter(pylab.NullFormatter())
        pylab.ylabel('y')
        #pylab.subplot(212)
        #difference = y_reconstructed - data.T[1]
        #pylab.scatter(data.T[0], difference)
        #for x,y in zip(data.T[0], difference):
        #    pylab.plot([x,x], [0,y], 'k')
        #pylab.xlim((numpy.min(data.T[0]),
        #            numpy.max(data.T[0])))
        #pylab.ylim((numpy.min(difference),
        #            numpy.max(difference)))
        pylab.xlabel('x')
        #pylab.ylabel('residue')
        pylab.tight_layout()
        pylab.show()

if __name__ == '__main__':

    main()
