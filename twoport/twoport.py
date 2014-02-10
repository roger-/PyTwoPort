from __future__ import division
import numpy as np
from numpy import asarray, ndarray, pi, array, dot, zeros, inf, identity, real, imag, sign, abs
from numpy.linalg import inv
#from pylab import *
import pylab as mpl
from scipy import interpolate
from utils import *

__all__ = ['OnePort', 'TwoPort',
           '_11', '_12', '_21', '_22', 'plot_gains']

# 2x2 identity matrix (for convenience)
I = identity(2, dtype='complex')

# few helper variables to simplify access to the two-port arrays. Now we can use
# y[_11] instead of y[:, 0, 0], for example
_11 = slice(None, None, None), 0, 0
_12 = slice(None, None, None), 0, 1
_21 = slice(None, None, None), 1, 0
_22 = slice(None, None, None), 1, 1

def pick_f(port1, port2):
    if port1.f is None:
        return port2.f
    elif port2.f is None:
        return port1.f
    elif (port1.f == port2.f).all():
        return port1.f
    else:
        raise Exception('frequency axes incompatible')

def plot_gains(two_port):
    mpl.plot(two_port.f/1e6, db(two_port.g_msg()), label='$MSG$')
    mpl.plot(two_port.f/1e6, db(two_port.g_max()), label='$G_\mathrm{max}$')
    mpl.plot(two_port.f/1e6, db(abs(two_port.s[_21])**2), label='$S_{21}$')

    mpl.legend()
    mpl.xlabel('Frequency (MHz)')
    mpl.ylabel('Gain (dB)')
    mpl.grid(which='both')

class NPort(object):
    z0 = 50
    f = None
    description = None
    components = []

    def series(self, one_port):
        raise NotImplementedError

    def parallel(self, one_port):
        raise NotImplementedError

    def input_equiv(self):
        raise NotImplementedError

    def output_equiv(self):
        raise NotImplementedError

    def interpolate(self, new_points):
        raise NotImplementedError

    def __getitem__(self, new_points):
        # if network is frequency independent, ignore indexing
        if self.f is None:
            return self

        if isinstance(new_points, slice):
            start = new_points.start
            stop = new_points.stop
            step = new_points.step

            if not start:
                start = self.f[0]
            if not stop:
                stop = self.f[-1]
            if not step:
                step = self.f[1] - self.f[0]

            new_points = np.arange(start, stop, step)
        else:
            new_points = np.asarray(new_points)
            new_points.shape = -1

        # interpolation requires monotonic frequency axis
        if (np.diff(self.f) <= 0).any():
            raise Exception('frequency axis must be monotonically increasing')

        # make sure we aren't extrapolating
        if new_points.min() < self.f.min() or new_points.max() > self.f.max():
            raise Exception('cannot extrapolate to specified frequency range')

        return self.interpolate(new_points)

class OnePort(NPort):
    def __init__(self, f=None, **kwargs):
        param_type = list(kwargs.keys())[0].lower()
        param_value = list(kwargs.values())[0]
        param_value = asarray(param_value, dtype='complex')

        if f is not None:
            self.f = f

        if param_type in ('z', 'r'):
            self._from_z(param_value)
        elif param_type == 'y':
            self._from_y(param_value)
        elif param_type == 's':
            self._from_s(param_value)
        else:
            raise Exception('unknown type: ' + param_type)

        '''elif param_type == 'l':
            self._from_z(2j*pi*self.f*param_value)
        elif param_type == 'c':
            self._from_y(2j*pi*self.f*param_value)'''

    def _from_s(self, s):
        s = asarray(s, dtype='complex').reshape(-1)
        self.s_params = s

    def _from_z(self, z):
        z = asarray(z, dtype='complex').reshape(-1)
        self.s_params = (z - self.z0)/(z + self.z0)

    def _from_y(self, y):
        y = asarray(y, dtype='complex').reshape(-1)
        self._from_z(1/y)

    def _to_z(self):
        s = self.s_params

        return self.z0*(1 + s)/(1 - s)

    def _to_s(self):
        return self.s_params

    def _to_y(self):
        return 1/self._to_z()

    def cascade(self, two_port):
        if not isinstance(two_port, TwoPort):
            raise Exception('can only cascade one-port with two-port network')

        # convert instance to two-port
        from networks import Shunt

        return Shunt(self).cascade(two_port)

    def series(self, one_port):
        z = self._to_z() + one_port._to_z()
        f = pick_f(self, one_port)

        return OnePort(z=z, f=f)

    def parallel(self, one_port):
        y = self._to_y() + one_port._to_y()
        f = pick_f(self, one_port)

        return OnePort(y=y, f=f)

    def input_equiv(self):
        return self

    def output_equiv(self):
        return self

    def s_in(self):
        return self.s_params

    def s_out(self):
        return self.s_params

    def conj(self):
        new_port = OnePort(s=self.s_params.conj(), f=self.f)

        return new_port

    def vswr(self):
        return (1 + abs(self.s_params))/(1 - abs(self.s_params))

    def interpolate(self, new_points):
        new_s = zeros(len(new_points), dtype='complex')
        s = self.s_params

        interp_real = interpolate.InterpolatedUnivariateSpline(self.f, s.real)
        interp_imag = interpolate.InterpolatedUnivariateSpline(self.f, s.imag)

        new_s.real = interp_real(new_points)
        new_s.imag = interp_imag(new_points)

        new_port = OnePort(s=new_s, f=new_points)
        new_port.description = self.description
        new_port.components = [c.interpolate(new_points) for c in self.components]

        return new_port

    def mismatch(self, load=None):
        ss = self.s_out()
        sl = 0 if load is None else load.s_in()

        return 1 - abs((sl - ss.conj())/(1 - ss*sl))**2

    def parallel_equiv(self):
        return 1/real(self.y), -1j/imag(self.y)

    def to_shunt(self):
        '''
        Convert one-port to shunt connected two-port
        '''
        from networks import Shunt

        return Shunt(self)

    def to_series(self):
        '''
        Convert one-port to series connected two-port
        '''
        from networks import Series

        return Series(self)

    def __repr__(self):
        desc = self.description if self.description else 'One-port network'

        if self.f is not None:
            freqs = ': f = {} {}Hz'.format(*to_eng_form(self.f[0]))

            if len(self.f) > 1:
                freqs += ' - {} {}Hz'.format(*to_eng_form(self.f[-1]))
        else:
            freqs = ': '

        if len(self.s) == 1:
            z = self.z[0]

            details = 'Z = {0:.2f}+{1:.2f}j'.format(float(z.real), float(z.imag))
            if abs(z.imag) > 0:
                rp, xp = self.parallel_equiv()
                details += ' = {0:.2f}//{1:.2f}j ('.format(float(rp.real), float(xp.imag))
                details += resolve_reactance_str(float(xp.imag), self.f[0]) + '), '
            else:
                details += ', '
            details += 'VSWR = {0:.2f}, '.format(float(self.vswr()))
            details += 'MM = {0:.2f} dB'.format(float(db(self.mismatch())))
        else:
            details = ''

        return '<' + desc + freqs + details + '>'

    s = property(_to_s, _from_s)
    y = property(_to_y, _from_y)
    z = property(_to_z, _from_z)

    inp = input_equiv
    out = output_equiv
    __mul__ = cascade
    __add__ = series
    __floordiv__ = parallel

class TwoPortBase(NPort):
    def __init__(self, f=None, **kwargs):
        '''for param_type in ('s', 'y', 'z', 't', 'abcd'):
            if param_type in kwargs:
                from_func = getattr(self, '_from_' + param_type)
                from_func(kwargs[param_type])
                break
        else:
            self.s_params = None'''

        if 's' in kwargs:
            self.s = kwargs['s']
        elif 'y' in kwargs:
            self.y = kwargs['y']
        elif 'z' in kwargs:
            self.z = kwargs['z']
        elif 't' in kwargs:
            self.t = kwargs['t']
        elif 'abcd' in kwargs:
            self.abcd = kwargs['abcd']
        else:
            self.s_params = None

        if f is not None:
            self.f = f

    def _from_s(self, s_params):
        s_params = asarray(s_params, dtype='complex').reshape(-1, 2, 2)
        self.s_params = s_params

    def _from_y(self, y_params):
        y_params = asarray(y_params, dtype='complex').reshape(-1, 2, 2)
        self.s_params = array([dot((I - self.z0*yi), inv(I + self.z0*yi)) for yi in y_params])

    def _from_z(self, z_params):
        z_params = asarray(z_params, dtype='complex').reshape(-1, 2, 2)
        self.s_params = array([(zi - self.z0*I).dot(inv(zi + self.z0*I)) for zi in z_params])

    def _from_abcd(self, abcd_params):
        abcd_params = asarray(abcd_params, dtype='complex').reshape(-1, 2, 2)
        self.s_params = zeros(abcd_params.shape, dtype='complex')

        a = abcd_params[_11]
        b = abcd_params[_12]
        c = abcd_params[_21]
        d = abcd_params[_22]

        self.s_params[_11] = a + b/self.z0 - c*self.z0 - d
        self.s_params[_12] = 2*(a*d - b*c)
        self.s_params[_21] = 2
        self.s_params[_22] = -a + b/self.z0 - c*self.z0 + d

        D = (a + b/self.z0 + c*self.z0 + d).reshape(-1, 1, 1)

        self.s_params /= D

    def _from_t(self, t_params):
        t_params = asarray(t_params, dtype='complex').reshape(-1, 2, 2)
        self.s_params = zeros(t_params.shape, dtype='complex')

        t11 = t_params[_11]
        t12 = t_params[_12]
        t21 = t_params[_21]
        t22 = t_params[_22]

        self.s_params[_11] = t12/t22
        self.s_params[_12] = (t11*t22 - t12*t21)/t22
        self.s_params[_21] = 1/t22
        self.s_params[_22] = -t21/t22

    def _to_abcd(self):
        abcd_params = zeros(self.s_params.shape, dtype='complex')

        s11 = self.s_params[_11]
        s12 = self.s_params[_12]
        s21 = self.s_params[_21]
        s22 = self.s_params[_22]

        abcd_params[_11] = ((1 + s11)*(1 - s22) + s12*s21)/(2*s21)
        abcd_params[_12] = ((1 + s11)*(1 + s22) - s12*s21)*self.z0/(2*s21)
        abcd_params[_21] = ((1 - s11)*(1 - s22) - s12*s21)/(2*s21*self.z0)
        abcd_params[_22] = ((1 - s11)*(1 + s22) + s12*s21)/(2*s21)

        return abcd_params

    def _to_t(self):
        t_params = zeros(self.s_params.shape, dtype='complex')

        s11 = self.s_params[_11]
        s12 = self.s_params[_12]
        s21 = self.s_params[_21]
        s22 = self.s_params[_22]
        d = self.det_s()

        t_params[_11] = -d/s21
        t_params[_12] = s11/s21
        t_params[_21] = -s22/s21
        t_params[_22] = 1/s21

        return t_params

    def _to_s(self):
        return self.s_params

    def _to_z(self):
        return array([inv(I - si).dot(I + si) for si in self.s_params.reshape(-1, 2, 2)])*self.z0

    def _to_y(self):
        return array([inv(I + si).dot(I - si) for si in self.s_params.reshape(-1, 2, 2)])/self.z0

    def det_s(self):
        s = self.s_params

        return s[_11]*s[_22] - s[_12]*s[_21]

    def s_in(self, load=None):
        if load is not None:
            sl = load.s_in()
        else:
            sl = 1

        s = self.s

        return s[_11] + s[_12]*s[_21]*sl/(1 - s[_22]*sl)

    def s_out(self, source=None):
        if source is not None:
            ss = source.s_out()
        else:
            ss = 1

        s = self.s

        return s[_22] + s[_12]*s[_21]*ss/(1 - s[_11]*ss)

    def flip_ports(self):
        s = self.s_params
        s_new = zeros(s.shape, dtype='complex')

        s_new[_11] = s[_22]
        s_new[_12] = s[_21]
        s_new[_21] = s[_12]
        s_new[_22] = s[_11]

        return TwoPort(s=s_new, f=self.f)

    def input_equiv(self):
        s_in = self.s_in()

        return OnePort(s=s_in, f=self.f)

    def output_equiv(self):
        s_out = self.s_out()

        return OnePort(s=s_out, f=self.f)

    def interpolate(self, new_points):
        new_s = zeros((len(new_points), 2, 2), dtype='complex')
        s = self.s_params

        for param in (_11, _12, _21, _22):
            interp_real = interpolate.InterpolatedUnivariateSpline(self.f, s[param].real)
            interp_imag = interpolate.InterpolatedUnivariateSpline(self.f, s[param].imag)

            new_s[param].real = interp_real(new_points)
            new_s[param].imag = interp_imag(new_points)

        new_port = TwoPort(s=new_s, f=new_points)
        new_port.description = self.description
        new_port.components = [c[new_points] for c in self.components]

        return new_port

    def cascade(self, two_port):
        if isinstance(two_port, OnePort):
            from networks import Shunt

            return self.cascade(Shunt(two_port))

        # handle case where one two port is defined at a single point
        if len(self.s) == 1:
            tp1 = self.t[0]
            t = array([dot(tp1, tp2) for tp2 in two_port.t])
        elif len(two_port.s) == 1:
            tp2 = two_port.t[0]
            t = array([dot(tp1, tp2) for tp1 in self.t])
        else:
            t = array([dot(tp1, tp2) for (tp1, tp2) in zip(self.t, two_port.t)])

        f = pick_f(self, two_port)

        return Cascade(self, two_port, t=t, f=f)

    def series(self, two_port):
        z = self._to_z() + two_port._to_z()
        f = pick_f(self, two_port)

        return TwoPort(z=z, f=f)

    def parallel(self, two_port):
        y = self._to_y() + two_port._to_y()
        f = pick_f(self, two_port)

        return TwoPort(y=y, f=f)

    '''def cascadexx(self, two_port):
        if isinstance(two_port, OnePort):
            y = self.y
            y[_22] += two_port.y

            return TwoPort(y=y, f=self.f)

        # handle case where one two port is defined at a single point
        if len(self.s) == 1:
            tp1 = self.abcd[0]
            abcd = array([dot(tp1, tp2) for tp2 in two_port.abcd])
        elif len(two_port.s) == 1:
            tp2 = two_port.abcd[0]
            abcd = array([dot(tp1, tp2) for tp1 in self.abcd])
        else:
            abcd = array([dot(tp1, tp2) for (tp1, tp2) in zip(self.abcd, two_port.abcd)])

        return TwoPort(abcd=abcd, f=self.f)'''

    def __repr__(self):
        desc = self.description if self.description else 'Two-port network'
        desc += ': '

        if self.f is not None:
            freqs = 'f = {} {}Hz'.format(*to_eng_form(self.f[0]))

            if len(self.f) > 1:
                freqs += ' - {} {}Hz'.format(*to_eng_form(self.f[-1]))
        else:
            freqs = ''

        details = ', '

        g = self.g_t()
        ind_max = np.argmax(g)
        details += '|S21| = {0:.2f} dB'.format(db(g[ind_max]))
        if self.f is not None and len(self.f) > 1:
            details += ' (@ {} {}Hz), '.format(*to_eng_form(self.f[ind_max]))
        else:
            details += ', '

        g = self.g_msg()
        ind_max = np.argmax(g)
        details += 'MSG = {0:.2f} dB'.format(db(g[ind_max]))
        if self.f is not None and len(self.f) > 1:
            details += ' (@ {} {}Hz), '.format(*to_eng_form(self.f[ind_max]))
        else:
            details += ', '

        if (self.mu_s() < 1).any():
            details += 'potentially unstable (k_min = {}), '.format(self.rollett_kt().min())
        else:
            g = self.g_max()
            ind_max = np.argmax(g)
            details += 'g_max = {0:.2f} dB'.format(db(g[ind_max]))
            if self.f is not None and len(self.f) > 1:
                details += ' (@ {} {}Hz)'.format(*to_eng_form(self.f[ind_max]))

        return '<' + desc + freqs + details + '>'

    s = property(_to_s, _from_s)
    y = property(_to_y, _from_y)
    z = property(_to_z, _from_z)
    t = property(_to_t, _from_t)
    abcd = property(_to_abcd, _from_abcd)

    inp = input_equiv
    out = output_equiv
    #__rrshift__ = cascade
    __mul__ = cascade
    __add__ = series
    __floordiv__ = parallel

class TwoPort(TwoPortBase):
    def vswr_in(self, load=None):
        s_in = self.s_in(load)

        return (1 + abs(s_in))/(1 - abs(s_in))

    def max_stable_gain(self):
        s12, s21 = self.s[_12], self.s[_21]

        return abs(s21/s12)

    def transducer_gain(self, source=None, load=None):
        ss = 0 if source is None else source.s_out()
        sl = 0 if load is None else load.s_in()

        s = self.s

        return ((1 - abs(ss)**2)/abs(1 - self.s_in(load)*ss)**2) * abs(self.s[_21])**2 * ((1 - abs(sl)**2)/abs(1 - self.s[_22]*sl)**2)

    def available_gain(self, source=None):
        ss = 0 if source is None else source.s_out()

        s = self.s

        return ((1 - abs(ss)**2)/abs(1 - s[_11]*ss)**2)*abs(s[_21]**2)/(1 - abs(self.s_out(source))**2)

    def power_gain(self, load=None):
        sl = 0 if load is None else load.s_in()

        s = self.s

        return ((1 - abs(sl)**2)/abs(1 - s[_22]*sl)**2)*abs(s[_21]**2)/(1 - abs(self.s_in(load))**2)

    def max_gain(self):
        k = self.rollett_k()
        g_msg = self.max_stable_gain()

        g_max = zeros(len(k))
        g_max[k >= 1] = g_msg[k >= 1]*(k[k >= 1] - np.sqrt(k[k >= 1]**2 - 1))
        g_max[k < 1] = g_msg[k < 1]

        return g_max

    def max_single_sided_matched_gain(self):
        k = self.rollett_k()
        g_msg = self.max_stable_gain()

        return 2*k*g_msg

    def max_double_sided_mismatched_gain(self):
        k = self.rollett_k()
        g_msg = self.max_stable_gain()

        return 2*(k + 1)*g_msg

    def matched_source(self):
        s = self.s
        d = self.det_s()

        b = 1 + abs(s[_11])**2 - abs(s[_22])**2 - abs(d)**2
        c = s[_11] - d*(s[_22].conj())

        ss = (b - sign(b)*np.sqrt(b**2 - 4*abs(c)**2))/(2*c)

        return OnePort(s=ss, f=self.f)

    def matched_load(self):
        s = self.s
        d = self.det_s()

        b = 1 + abs(s[_22])**2 - abs(s[_11])**2 - abs(d)**2
        c = s[_22] - d*(s[_11].conj())

        sl = (b - sign(b)*np.sqrt(b**2 - 4*abs(c)**2))/(2*c)

        return OnePort(s=sl, f=self.f)

    '''def mismatch(self, load):
        ss = self.s_out()
        sl = load.s_in()

        return 1 - abs((sl - ss.conj())/(1 - ss*sl))**2'''

    def mu_stability_source(self):
        s = self.s_params
        d = self.det_s()

        return (1 - abs(s[_11])**2)/(abs(s[_12]*s[_21]) + abs(s[_22] - d*(s[_11].conj())))

    def mu_stability_load(self):
        s = self.s_params
        d = self.det_s()

        return (1 - abs(s[_22])**2)/(abs(s[_12]*s[_21]) + abs(s[_11] - d*(s[_22].conj())))

    def rollett_k(self):
        s = self.s_params
        d = self.det_s()

        return (1 - abs(s[_11])**2 - abs(s[_22])**2 + abs(d)**2)/(2*abs(s[_12]*s[_21]))

    def rollett_kt(self):
        '''
        EL Tan's single parameter modified Rollett stability factor
        '''
        s = self.s_params
        d = self.det_s()

        return (1 - abs(s[_11])**2 - abs(s[_22])**2 + abs(d)**2 + 0.5*(1 - abs(d)**2 - abs(1 - abs(d)**2)))/(2*abs(s[_12]*s[_21]))

    def input_tunability(self, load):
        s = self.s_params
        d = self.det_s()
        sl = load.s_in()

        return abs(s[_12]*s[_21]*sl)/abs(1 - s[_22]*sl)/abs(s[_11] - d*sl)

    def source_stability_circle(self):
        s = self.s_params

        det = s[_11]*s[_22] - s[_12]*s[_21]

        c = s[_11] - det*(s[_22].conj())
        d = (abs(c)**2 - abs(s[_12]*s[_21])**2)/(1 - abs(s[_22])**2)

        center = c.conj()/d
        radius = abs(s[_12]*s[_21])/abs(d)
        circle_stable = d < 0

        return center, radius, circle_stable

    def load_stability_circle(self):
        s = self.s_params

        det = s[_11]*s[_22] - s[_12]*s[_21]

        c = s[_22] - det*(s[_11].conj())
        d = abs(s[_22])**2 - abs(det)**2

        center = c.conj()/d
        radius = abs(s[_12]*s[_21])/abs(d)
        circle_stable = d < 0

        return center, radius, circle_stable

    def operating_gain_circle(self, G):
        s = self.s

        det = s[_11]*s[_22] - s[_12]*s[_21]
        g = G/(abs(s[_21])**2)

        c = s[_22] - det*(s[_11].conj())
        d = abs(s[_22])**2 - abs(det)**2
        k = self.rollett_k()

        center = g*(c.conj())/(1 + g*d)
        radius = np.sqrt(g**2*abs(s[_12]*s[_21])**2 - 2*g*k*abs(s[_12]*s[_21]) + 1)/abs(1 + g*d)

        return center, radius

    def available_gain_circle(self, G):
        s = self.s_params

        det = s[_11]*s[_22] - s[_12]*s[_21]
        g = G/abs(s[_21])**2

        c = s[_11] - det*(s[_22].conj())
        d = abs(s[_11])**2 - abs(det)**2
        k = self.rollett_k()

        center = g*(c.conj())/(1 + g*d)
        radius = np.sqrt((g*abs(s[_12]*s[_21]))**2 - 2*g*k*abs(s[_12]*s[_21]) + 1)/abs(1 + g*d)

        return center, radius

    def change_ref_plane(self, delay_11, delay_22):
        pd_11 = np.exp(1j*2*pi*self.f*delay_11)
        pd_22 = np.exp(1j*2*pi*self.f*delay_22)

        correction = np.zeros(self.s.shape, dtype='complex')
        correction[_11] = pd_11**2
        correction[_22] = pd_22**2
        correction[_12] = pd_11*pd_22
        correction[_21] = pd_11*pd_22

        return TwoPort(s=self.s*correction, f=self.f)

    g_t = transducer_gain
    g_a = available_gain
    g_p = power_gain
    g_msg = max_stable_gain
    g_msm = max_single_sided_matched_gain
    g_mdm = max_double_sided_mismatched_gain
    g_max = max_gain
    mu_s = mu_stability_source
    mu_l = mu_stability_load
    k = rollett_k
    kt = rollett_kt

class Cascade(TwoPort):
    def __init__(self, left_port, right_port, *args, **kwargs):
        self.components = []
        if isinstance(left_port, Cascade) and left_port.components:
            self.components += left_port.components
        else:
            self.components += [left_port]
        if isinstance(right_port, Cascade) and right_port.components:
            self.components += right_port.components
        else:
            self.components += [right_port]

        self.description = 'Cascade of {} two-port networks'.format(len(self.components))

        TwoPort.__init__(self, *args, **kwargs)

    def merge(self):
        self.components = []

def main():
    r = OnePort(z=50+40j, f=[1e9]) // OnePort(z=50)
    print(repr(r))

    r = OnePort(z=25) + OnePort(z=25)
    print(r.z)

    z = TwoPort(s=[[0.4, 0.1], [3, 0.2]])

    rl = z.matched_load()
    rs = z.matched_source()

    print(z.rollett_k())
    print(z.g_max(), z.g_t(rs, rl))
    print(rs.z, rl.z)


if __name__ == '__main__':
    main()
