def extract_contour(x_list, y_list, z_list, z_value):

    import pylab
    import numpy

    temp = pylab.tricontour(x_list, y_list, z_list, levels=[z_value])
    pylab.clf()
    return numpy.concatenate([itm.vertices for itm
                              in temp.collections[0].get_paths()])

def get_input_from_user():

    from argparse import ArgumentParser

    parser = ArgumentParser(description='extracts points along the shock front from a snapshot')
    parser.add_argument('file_name', help='name of file')
    parser.add_argument('out_file', help='name of output file')
    parser.add_argument('--show', 
                        help='toggles plotting the points',
                        action='store_true')

    return parser.parse_args()

def main():

    import h5py
    import numpy
    import pickle

    args = get_input_from_user()
    if args.show:
        import pylab
    with h5py.File(args.file_name,'r+') as f:
        g = 5./3.
        obstacle_list = numpy.array(f['obstacle'])
        density_list = numpy.array(f['density'])
        pressure_list = numpy.array(f['pressure'])
        x_list = numpy.array(f['x_coordinate'])
        y_list = numpy.array(f['y_coordinate'])
    s_list = pressure_list/density_list**g
    threshold = numpy.min(s_list)*2
    contour_points = extract_contour(x_list[obstacle_list<0.5],
                                     y_list[obstacle_list<0.5],
                                     s_list[obstacle_list<0.5],
                                     threshold)
    pickle.dump(contour_points, open(args.out_file,'wb'))
    if args.show:
        pylab.plot(contour_points.T[0],
                   contour_points.T[1],'.')
        pylab.show()

if __name__ == '__main__':

    main()
