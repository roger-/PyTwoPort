from __future__ import division
from pylab import *
from twoport import *
from utils import *


class Series(TwoPort):
    def __init__(self, one_port):
        self.z0 = one_port.z0
        self.f = one_port.f
        self.description = ('Series ' + one_port.description) if one_port.description else None
        self.components = [one_port]

        if 0:
            z = one_port.z
            s = zeros((len(asarray(z).reshape(-1)), 2, 2), dtype='complex')
            s[_11] = s[_22] = z/(z + 2*self.z0)
            s[_12] = s[_21] = 2*self.z0/(z + 2*self.z0)

        s = one_port.s
        s_new = zeros((len(asarray(s).reshape(-1)), 2, 2), dtype='complex')
        s_new[_11] = s_new[_22] = (1 + s)/(3 - s)
        s_new[_12] = s_new[_21] = 2*(1 - s)/(3 - s)

        TwoPort.__init__(self, s=s_new)

class Shunt(TwoPort):
    def __init__(self, one_port):
        self.z0 = one_port.z0
        self.f = one_port.f
        self.description = ('Shunt ' + one_port.description) if one_port.description else None
        self.components = [one_port]

        if 0:
            y = one_port.y
            s = zeros((len(asarray(y).reshape(-1)), 2, 2), dtype='complex')
            s[_11] = s[_22] = -y/(y + 2/self.z0)
            s[_12] = s[_21] = (2/self.z0)/(y + 2/self.z0)

        s = one_port.s
        s_new = zeros((len(asarray(s).reshape(-1)), 2, 2), dtype='complex')
        s_new[_11] = s_new[_22] = (s - 1)/(3 + s)
        s_new[_12] = s_new[_21] = 2*(s + 1)/(s + 3)

        TwoPort.__init__(self, s=s_new)

class Capacitor(OnePort):
    def __init__(self, val, f=None):
        if f is not None:
            self.f = f

        y = 1j*2*pi*self.f*val
        self.description = '{:>3.2f} {}F capacitor'.format(*to_eng_form(val))

        OnePort.__init__(self, y=y)

class Inductor(OnePort):
    def __init__(self, val, f=None, q=inf):
        if f is not None:
            self.f = f

        z = (1j + 1/q)*2*pi*self.f*val
        self.description = '{:>3.2f} {}H inductor'.format(*to_eng_form(val))

        OnePort.__init__(self, z=z)

class Resistor(OnePort):
    def __init__(self, val, f=None):
        if f is not None:
            self.f = f

        # FIXME: work for multivalued R's
        if 0:
            self.description = '{:>3.2f} {}ohm resistor'.format(*to_eng_form(val))

        OnePort.__init__(self, r=val)

class Terminator(OnePort):
    def __init__(self):
        self.description = 'Terminator'

        OnePort.__init__(self, s=0)

class ImpedanceInverter(TwoPort):
    def __init__(self, gain=1):
        abcd = zeros((1, 2, 2), dtype='complex')
        abcd[_12] = 1j*gain
        abcd[_21] = 1j/gain

        self.description = 'Impedance inverter'

        TwoPort.__init__(self, abcd=abcd)

class Transformer(TwoPort):
    def __init__(self, n_pri=1, n_sec=1):
        s = zeros((1, 2, 2), dtype='complex')
        n = n_pri/n_sec

        s[_11] = n**2 - 1
        s[_22] = -s[_11]
        s[_21] = s[_12] = 2*n
        s /= (n**2 + 1)

        self.description = '{}:{} transformer'.format(n_pri, n_sec)

        TwoPort.__init__(self, s=s)

class PhaseShifter(TwoPort):
    def __init__(self, phase_rad, f=TwoPort.f):
        s = zeros((len(phase_rad), 2, 2), dtype='complex')
        s[_11] = s[_22] = 0
        s[_12] = s[_21] = exp(1j*phase_rad)

        TwoPort.__init__(self, s=s, f=f)

class TransmissionLine(TwoPort):
    '''
    Loss-less transmission line
    '''
    def __init__(self, length, vf=1, z0=TwoPort.z0, f=TwoPort.f):
        c =  2.99792458e8
        beta_l = length*2*pi*f/(c*vf)

        abcd = zeros((len(beta_l), 2, 2), dtype='complex')
        abcd[_11] = abcd[_22] = cos(beta_l)
        abcd[_12] = 1j*z0*sin(beta_l)
        abcd[_21] = 1j*sin(beta_l)/z0

        self.description = '{} m transmission line'.format(length)

        TwoPort.__init__(self, abcd=abcd, f=f)

class Amplifier(TwoPort):
    def __init__(self, gain_db, f=TwoPort.f):
        s = zeros((len(f), 2, 2), dtype='complex')
        s[_21] = un_db(gain_db)
        self.description = '{} dB amplifier'.format(gain_db)

        TwoPort.__init__(self, s=s, f=f)

def pi_attenuator(atten_db, z=50):
    k = 10**(atten_db/10)

    rp = z*(sqrt(k) + 1)/(sqrt(k) - 1)
    rs = z*(k - 1)/(2*sqrt(k))

    return Shunt(Resistor(rp))*Series(Resistor(rs))*Shunt(Resistor(rp))

def main1():
    f = linspace(1e6, 2e9, 2000)

    load = OnePort(z=75/2)
    source = OnePort(z=50)

    connector = TransmissionLine(1/100, 0.66, z0=75, f=f)
    cable = TransmissionLine(0.5, 0.66, z0=50, f=f)

    net = cable*connector

    plot(f/1e6, db(net.g_t(source, load)))

    net = cable

    plot(f/1e6, db(net.g_t(source, load)))

    ylabel('Gain (dB)')
    ylabel('Frequency (MHz)')

    show()

def main2():
    f = linspace(1e6, 1e9, 1000)

    tl = TransmissionLine(10, 0.66, z0=50, f=f)*\
        TransmissionLine(0.1, 0.54, z0=40, f=f)
    n = tl*Resistor(50)


    figure()
    plot(f, abs(n.inp().z))

    figure()
    plot(f/1e6, db(tl.s[_21]))

    show()

if __name__ == '__main__':
    main1()