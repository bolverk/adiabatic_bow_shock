def extract_contour(x_list, y_list, z_list, z_value):

    import pylab
    import numpy

    temp = pylab.tricontour(x_list, y_list, z_list, levels=[z_value])
    pylab.clf()
    return numpy.concatenate([itm.vertices for itm
                              in temp.collections[0].get_paths()])

def main():

    import numpy
    import h5py
    import pylab

    fname = 'snapshot_999.h5'

    data = {}
    with h5py.File(fname,'r+') as f:
        data['density'] = numpy.array(f['hydrodynamic']['density'])
        data['pressure'] = numpy.array(f['hydrodynamic']['pressure'])
        data['x'] = numpy.array(f['geometry']['x_coordinate'])
        data['y'] = numpy.array(f['geometry']['y_coordinate'])
    g = 5./3.
    data['entropy'] = numpy.log(data['pressure']) - g*numpy.log(data['density'])

    #pylab.tricontourf(data['x'],
    #                  data['y'],
    #                  data['entropy'])
    #pylab.axis('equal')
    #pylab.colorbar()
    #pylab.show()
    contour = extract_contour(data['x'],
                              data['y'],
                              data['entropy'],
                              -15)
    fit1 = numpy.polyfit(contour.T[0],
                         contour.T[1],
                         2)
    fit2 = numpy.polyfit(contour.T[0]**2,
                         contour.T[1],
                         1)
    pylab.plot(contour.T[0],
               contour.T[1],
               '.')
    pylab.plot(contour.T[0],
               pylab.polyval(fit1,contour.T[0]))
    pylab.plot(contour.T[0],
               pylab.polyval(fit2,contour.T[0]**2))
    #pylab.axis('equal')
    pylab.show()

if __name__ == '__main__':

    main()
