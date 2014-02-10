from __future__ import division
from twoport import *
from numpy import *

units_multipliers = {'hz':1, 'khz':1e3, 'mhz':1e6, 'ghz':1e9}
format_converters = {'ma': lambda m, a: m*exp(1j*radians(a)),
                     'ri': lambda r, i: r + 1j*i,
                     'db': lambda db, a: (10**(db/20))*exp(1j*radians(a))}

def load_snp(file_name):
    # default options per file specification
    options = ['ghz', 's', 'ma', 'r', '50']

    with open(file_name, mode='r') as fh:
        # need to skip comments until we get to the options line
        while True:
            line = fh.readline()

            if line[0] == '#':
                break
            elif line[0] != '!':
                raise Exception('unknown file format')

        # extract options into variables, using defaults when necessary
        options_raw = line.lower().split()[1:]
        options[:len(options_raw)] = options_raw
        units, type, format, null, z0 = options

        # only S-parameters supported
        if type != 's' or format not in format_converters.keys():
            raise Exception('unsupported data format')

        # use NumPy to load actual data
        data = loadtxt(fh, comments='!', unpack=True)

    # split data into frequency and magnitude/angle pairs
    if len(data) == 3: # one-port
        f, s11a, s11b = data
        s = zeros(len(f), dtype='complex')
    elif len(data) == 9: # two-port
        f, s11a, s11b, s21a, s21b, s12a, s12b, s22a, s22b = data
        s = zeros((len(f), 2, 2), dtype='complex')
    else:
        raise Exception('unsupported data format')

    # scale frequency units
    f *= units_multipliers[units]

    # load into array and create port
    if len(data) == 3: # one-port
        s[:] = format_converters[format](s11a, s11b)

        return OnePort(s=s, f=f, z0=z0)
    else:
        s[_11] = format_converters[format](s11a, s11b)
        s[_12] = format_converters[format](s12a, s12b)
        s[_21] = format_converters[format](s21a, s21b)
        s[_22] = format_converters[format](s22a, s22b)

        return TwoPort(s=s, f=f, z0=z0)

load_s2p = load_snp

def main():
    pass

if __name__ == '__main__':
    main()
