from __future__ import division
import pylab as pyl
from twoport import *
from twoport.networks import *


def basic():
    f = linspace(100e6, 1000e6)

    r = Resistor(100)
    c = Capacitor(10e-12, f=f)
    terminator = Resistor(50)

    filter = Series(c) * Shunt(r)
    input_equiv = filter.inp()

    g = filter.transducer_gain(terminator, terminator)

    print filter[200e6]
    print 'Zin @ 200 MHz:', input_equiv[200e6].z

    pyl.figure()
    pyl.plot(f/1e6, utils.db(g))
    pyl.xlabel('Frequency (MHz)')
    pyl.ylabel('Gain (dB)')

    pyl.show()

def amplifier():
    file_name = 'BFR93A.s2p'
    analysis_freq = 400e6

    ## load a transistor model
    q = load_snp(file_name)
    f = linspace(q.f[0], q.f[-1], 1000) # new f axis with 1000 points
    q = q[f] # interpolate to new f axis

    ## plot input reflection coefficient (S11) on Smith chart
    pyl.figure()
    sc = SmithChart(show_cursor=False, labels=True)
    sc.plot_s_param(q.inp().s)

    # calculate various stability factors
    k = q.rollett_k()
    kt = q.rollett_kt()
    mu = q.mu_stability_source()

    pyl.figure()
    pyl.plot(f/1e6, kt, label='EL Tan\'s modified Rollett $k$')
    pyl.plot(f/1e6, mu, label='$\mu_i$')
    pyl.axhline(1.0, color='black') # line at stability limit
    pyl.legend(loc='lower right')

    ## plot gains

    # We can do this manually using:
    # terminator = Resistor(50)
    # gt = q.transducer_gain(terminator, terminator)
    # gmsg = q.max_stable_gain()
    # gmax = q.max_gain()
    # etc. or just:

    pyl.figure()
    plot_gains(q)

    ## plot stability circles

    q_desired = q[analysis_freq]

    pyl.figure()
    sc = SmithChart(show_cursor=False, labels=True)
    sc.two_port = q_desired # need this for cursors to work

    sc.draw_stability_circles(q_desired, plane='both')
    sc.plot_gain_circles(q_desired, surface=True, gain_cursor=True)

    # enable cursor
    # note green circle will show optimal termination at other (load) port
    sc.enable_cursor()

    pyl.show()


def main():
    basic()
    amplifier()

if __name__ == '__main__':
    main()
