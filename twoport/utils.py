from __future__ import division
import numpy as np

def db(val):
    return 10*np.log10(abs(val))

def un_db(val):
    return 10**(val/10)

def to_eng_form(x):
    si_prefixes = {9:'P', 9:'G', 6:'M', 3:'k', 0:'',
                   -3:'m', -6:'u', -9:'n', -12:'p', -15:'f'}

    for e in sorted(si_prefixes, reverse=True):
        if x >= 10**e:
            return x/10**e, si_prefixes[e]

    return x, ''

def resolve_reactance(x, f):
    w = 2*np.pi*f

    if x >= 0:
        unit = 'H'
        lc = x/w
    else:
        unit = 'F'
        lc = 1/(w*x)

    return abs(lc), unit

def resolve_reactance_str(x, f):
    value, unit = resolve_reactance(x, f)
    value_eng, prefix = to_eng_form(value)

    return '{:>3.2f} {}{}'.format(value_eng, prefix, unit)

def group_delay(two_port):
    f = two_port.f
    s21 = two_port.s[_21]

    d = -diff(unwrap(angle(s21)))/diff(f)/(2*pi)

    return r_[d[0], d]

def main():
    print(resolve_reactance_str(125, 100e6))

if __name__ == '__main__':
    main()
