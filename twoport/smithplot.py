from __future__ import division
from numpy import asarray, array, dot, zeros, inf, identity
from numpy.linalg import inv
import numpy as np
from pylab import *
from twoport import *
from networks import *
from utils import *
from matplotlib.patches import Circle, FancyArrowPatch 	# for drawing smith chart
from matplotlib.lines import Line2D		# for drawing smith chart
from matplotlib.text import Text
from itertools import chain

# TODO:
# * Label points by freq
# * Label source/load cursors
# clip to rect

class SmithChart(object):
    z0 = 50
    current_plane = None
    two_port = None
    show_matched_termination_cursor = True
    gain_cursor = False

    def __init__(self, ax=None, show_cursor=True, admittance=True, labels=False):
        self.ax = ax if ax else gca()
        self.fig = self.ax.figure

        self.draw_smith_chart(admittance, labels)

        if show_cursor:
            self.enable_cursor()

    def _impedance_circle_coords(self, intercept_x_coords, intercept_angles):
        # find coords for constant R circles
        rs = 2/(1 - asarray(intercept_x_coords)) - 1 # convert to desired resistances
        r_circle_coords = [((r/(r+1), 0), 1/(r + 1)) for r in rs]

        # find coords for constant X circles
        intercept_angles = np.radians(intercept_angles)
        xs = np.sin(intercept_angles)/(1 - np.cos(intercept_angles)) # convert to desired reactances
        x_circle_coords = [((1, 1/x), abs(1/x)) for x in xs]

        return r_circle_coords, x_circle_coords, rs, xs

    def draw_impedance_circles(self, intercept_x_coords, intercept_angles, labels=False):
        r_circle_coords, x_circle_coords, rs, xs = self._impedance_circle_coords(intercept_x_coords, intercept_angles)

        for center, radius in chain(r_circle_coords, x_circle_coords):
            c = Circle(center, radius, **self.patch_options_dark)

            c.set_clip_path(self.smith_circle)
            c.set_clip_box(self.ax.bbox)
            self.ax.add_patch(c)

        if labels:
            for x, r in zip(intercept_x_coords, rs):
                self.ax.text(x + 0.04, 0.03, '%.0f' % round(self.z0*r), **self.label_options)

            for a, x in zip(intercept_angles, xs):
                r = (a - 90) if x > 0 else (a + 90)
                a = np.radians(a)
                d = 1.04

                self.ax.text(d*cos(a), d*sin(a), '%.0fj' % round(self.z0*x), rotation=r, **self.label_options)

    def draw_admittance_circles(self, intercept_x_coords, intercept_angles, labels=False):
        r_circle_coords, x_circle_coords, rs, xs = self._impedance_circle_coords(intercept_x_coords, intercept_angles)

        # admittance circles have same coords as impedance ones, except flipped
        # on the y-axis
        for (x, y), radius in chain(r_circle_coords, x_circle_coords):
            c = Circle((-x, -y), radius, **self.patch_options_light)

            c.set_clip_path(self.smith_circle)
            c.set_clip_box(self.ax.bbox)
            self.ax.add_patch(c)

        if labels:
            for x, r in zip(intercept_x_coords, rs):
                self.ax.text(-x, 0, '%.1f' % (1/(50*r)), **self.label_options)

            for a, x in zip(intercept_angles, xs):
                r = (a - 90) if x < 0 else (a + 90)
                a = np.radians(a)

                self.ax.text(cos(pi - a), sin(pi - a), '%.1f' % (1/(50*x)), rotation=r, **self.label_options)

    def draw_vswr_circles(self, vswr_radii, labels=False):
        for r in vswr_radii:
            c = Circle((0, 0), r, ls='dashed', **self.patch_options_light)

            c.set_clip_path(self.smith_circle)
            c.set_clip_box(self.ax.bbox)
            self.ax.add_patch(c)

        if labels:
            for r in vswr_radii:
                if r > 0:
                    vswr = (1 + r)/(1 - r)
                    self.ax.text(0, r, '%.1f' % vswr, **self.label_options)

    def draw_chart_axes(self):
        # make outer circle
        self.smith_circle = Circle((0, 0), 1, transform=self.ax.transData, fc='none',
                              **self.patch_options_axis)
        self.ax.add_patch(self.smith_circle)

        # make constant r=1 circle
        z0_circle = Circle((0.5, 0), 0.5, transform=self.ax.transData, fc='none',
                              **self.patch_options_axis)
        z0_circle.set_clip_path(self.smith_circle)
        z0_circle.set_clip_box(self.ax.bbox)
        self.ax.add_patch(z0_circle)

        # make x-axis
        line = Line2D([-1,1],[0,0], **self.patch_options_axis)
        line.set_clip_path(self.smith_circle)
        line.set_clip_box(self.ax.bbox)
        self.ax.add_line(line)

    def draw_smith_chart(self, admittance, labels):
        # plot options for constant z/y circles and axes
        self.patch_options_light = {'fc':'none', 'color':'#474959', 'alpha':0.2, 'lw':1}
        self.patch_options_dark = {'fc':'none', 'color':'#474959', 'alpha':0.5, 'lw':1}
        self.patch_options_axis = {'color':'black', 'alpha':0.8, 'lw':1.5}

        # options for z/y circle labels
        self.label_options = {'ha':'center', 'va':'center', 'size':'9', 'alpha':0.5}#,
                         #'bbox':dict(fc='white', ec='none', alpha=0.5)}
        #self.label_options = {'ha':'center', 'va':'center', 'size':'10', 'alpha':0.5}

        # x-axis coordinates where constant R circles will intersect
        intercept_x_coords = arange(-0.75, 1, 0.25)

        # angles where constant X circles will intersect (in degrees relative
        # to positive x-axis)
        intercept_angles = arange(40, 360, 40)

        # radii for vswr circles
        vswr_radii = arange(0, 1, 0.2)

        self.draw_chart_axes()
        self.draw_impedance_circles(intercept_x_coords, intercept_angles, labels)
        self.draw_admittance_circles(intercept_x_coords, intercept_angles, labels=0)
        self.draw_vswr_circles(vswr_radii, labels)

        self.ax.grid(0)
        self.ax.axis('equal')
        self.ax.axis(np.array([-1.1, 1.1, -1.1, 1.1]))

        self.save_background()

    def enable_cursor(self):
        self.cid_on_mouse_move = self.fig.canvas.mpl_connect('motion_notify_event', self._on_mouse_move)
        self.fig.canvas.mpl_connect('button_press_event', self._on_mouse_press)
        #self.cid_on_mouse_leave = self.fig.canvas.mpl_connect('axes_leave_event', self._on_mouse_leave)
        #self.cid_resize = self.fig.canvas.mpl_connect('resize_event', self._resize)
        self.cid_draw = self.fig.canvas.mpl_connect('draw_event', self._draw)

        self.cursor_text = Text(0.7, 0.1, 'n/a',
            horizontalalignment='left',
            verticalalignment='center',
            #bbox=dict(facecolor='white', alpha=1, ec='gray', pad=10),
            color='white',
            fontsize='small',
            alpha=0.8,
            backgroundcolor='gray',
            family='monospace',
            clip_box=self.ax.bbox,
            animated=True,
            transform = self.ax.transAxes)
        self.ax.add_artist(self.cursor_text)

        self.circle_matched_term = Circle((100,100), 0.02, alpha=0.5, color='green', animated=True)
        self.ax.add_patch(self.circle_matched_term)

    def plot_circle(self, center, radius, text='', text_color='black', circle_color='red',\
                filled=False, circle_alpha=0.2, linestyle='solid', hatch=None):
        text_alpha = 0.7

        text_x_offset = 0.15
        text_y_offset = 0.1

        # calculate position for text
        a, b = center

        # find closet point to edge of circle
        x = a*(1 - radius/sqrt(a**2 + b**2))
        y = b*(1 - radius/sqrt(a**2 + b**2))
        dist_to_circle = sqrt(x**2 + y**2)

        #print 'dist to sphere: ', sqrt(x**2 + y**2)

        if (x**2 + y**2) == 1:
            text_x_offset *= -1
            text_y_offset *= -1

            #textpos = ((a - x)/radius + text_x_offset, (b - y)/radius + text_y_offset)
            textpos = (x/dist_to_circle + text_x_offset, y/dist_to_circle + text_y_offset)
        else:
            textpos = (x + text_x_offset, y + text_y_offset)

        #print 'textpos: ', textpos

        # make actual circle
        c = Circle(center, radius, fc=circle_color, ec='white', alpha=circle_alpha, fill=filled, linestyle=linestyle, hatch=hatch)
        self.ax.add_patch(c)

        self.ax.annotate(text, (x, y), color=text_color,\
            fontsize='small', alpha=text_alpha, xycoords='data', xytext=textpos,\
            arrowprops=dict(arrowstyle="->", alpha=text_alpha,\
            connectionstyle="arc3,rad=0.17"))

    def draw_stability_circles(self, two_port, plane='both', label=None):
        filled = True

        if plane == 'source' or plane == 'both':
            for i, (center, radius, circle_stable) in enumerate(zip(*two_port.source_stability_circle())):
                if circle_stable:
                    hatch = None
                    color = 'green'
                else:
                    hatch = '/'
                    color = 'red'

                if label is None:
                    label_actual = '\n(%0.2f MHz)' % (two_port.f[i]/1e6) if two_port.f is not None else ''
                    label_actual = 'Source' + label_actual
                else:
                    label_actual = label

                self.plot_circle((center.real, center.imag), radius, text=label_actual, filled=filled,\
                                text_color='black', circle_color=color, hatch=hatch)

        if plane == 'load' or plane == 'both':
            for i, (center, radius, circle_stable) in enumerate(zip(*two_port.load_stability_circle())):
                if circle_stable:
                    hatch = None
                    color = 'green'
                else:
                    hatch = '/'
                    color = 'red'

                if label is None:
                    label_actual = '\n(%0.2f MHz)' % (two_port.f[i]/1e6) if two_port.f is not None else ''
                    label_actual = 'Load' + label_actual
                else:
                    label_actual = label

                self.plot_circle((center.real, center.imag), radius, text=label_actual, filled=filled,\
                                text_color='black', circle_color=color, hatch=hatch)

        self.save_background()

    def plot_gain_circles(self, two_port, gains_db=None, plane='source', \
                          gain_cursor=True, surface=False):
        self.two_port = two_port
        self.gain_cursor = gain_cursor
        self.current_plane = plane

        max_gain = self.two_port.g_max()

        if (self.two_port.kt() < 1).any():
            max_gain = self.two_port.max_double_sided_mismatched_gain()

        if gains_db == None:
            gains_db = np.trunc(10*db(max_gain)[0]*linspace(0.5, 1, 4))/10
        if surface:
            filled = False
            circle_alpha = 0.7
            linestyle = 'dashed'
        else:
            filled = True
            circle_alpha = 0.2
            linestyle = 'solid'

        for g in gains_db:
            if plane == 'source':
                center, radius = self.two_port.available_gain_circle(un_db(g))
                self.ax.set_title('Source plane')
            elif plane == 'load':
                center, radius = self.two_port.operating_gain_circle(un_db(g))
                self.ax.set_title('Load plane')

            text = str(g) + ' dB'

            self.plot_circle((center[0].real, center[0].imag), radius[0], text=text, filled=filled,\
                text_color='black', circle_color='orange', circle_alpha=circle_alpha,\
                linestyle=linestyle)

        if not surface:
            self.save_background()

            return

        num_segments = 1000

        i, r = meshgrid(linspace(-1.1, 1.1, num_segments), linspace(-1.1, 1.1, num_segments))
        gamma = i + 1j*r
        term = OnePort(s=gamma.reshape(-1))

        if plane == 'source':
            g = self.two_port.g_a(term)
        elif plane == 'load':
            g = self.two_port.g_p(term)

        g = db(g).reshape(num_segments, num_segments)
        g[abs(gamma) > 1] = nan

        im = self.ax.imshow(g, origin='lower', extent=(-1.1, 1.1, -1.1, 1.1), interpolation='bicubic', alpha=0.9)
        im.set_clim(gains_db[0], db(max_gain))

        self.save_background()

    def _on_mouse_press(self, event):
        if event.button == 3:
            print 'clearing'
            self.fig.canvas.restore_region(self.background, bbox=self.ax.bbox)

    def _on_mouse_move(self, event):
        if event.xdata is None or event.ydata is None:
            return

        gamma = event.xdata + 1j*event.ydata
        cursor_port = OnePort(s=gamma)

        self.fig.canvas.restore_region(self.background, bbox=self.ax.bbox)

        self.update_gain_cursor(cursor_port)
        self.update_info_box(cursor_port)

        #self.fig.canvas.draw()
        self.fig.canvas.blit(self.ax.bbox)

    def update_gain_cursor(self, cursor_port):
        if not self.gain_cursor or self.two_port is None:
            return

        if self.current_plane == 'source':
            term_matched = (cursor_port*self.two_port).out().conj()
            term_resulting = (self.two_port*term_matched).inp()
        else:
            term_matched = (self.two_port*cursor_port).inp().conj()
            term_resulting = (term_matched*self.two_port).out()

        self.circle_matched_term.center = real(term_matched.s), imag(term_matched.s)

        if abs(term_matched.s) > 1 or abs(term_resulting.s) > 1:
            self.circle_matched_term.set_color('r')
        else:
            self.circle_matched_term.set_color('g')

        self.ax.draw_artist(self.circle_matched_term)

    def update_info_box(self, cursor_port):
        z = cursor_port.z
        y = cursor_port.y

        data = '$\\Gamma$ = {:>7.2f} $\\angle$ {:>+8.3f}$^{{\\circ}}$'.format(float(abs(cursor_port.s)), float(degrees(angle(cursor_port.s))))
        data += '\nVSWR = {:>7.2f}'.format(float(cursor_port.vswr()))
        data += '\nMM = {:>7.2f} dB'.format(float(db(cursor_port.mismatch())))
        data += '\nZ = {:>7.2f}{:>+8.2f}j '.format(float(real(z)), float(imag(z)))

        if self.two_port is not None:
            if self.current_plane == 'source':
                term_matched = (cursor_port*self.two_port).out().conj()
                gain = self.two_port.g_t(cursor_port, term_matched)
            else:
                term_matched = (self.two_port*cursor_port).inp().conj()
                gain = self.two_port.g_t(term_matched, cursor_port)

            x = float(imag(z))
            w = 2*pi*self.two_port.f[0]

            if x > 0:
                hf = 'H'
                lc = x/w
            else:
                hf = 'F'
                lc = 1/(w*x)

            val, prefix = to_eng_form(abs(lc))
            data += '({:>3.2f} {}{}) '.format(val, prefix, hf)
        data += '$\\Omega$'

        data += '\n  = {:>7.2f} // {:>+8.2f}j '.format(float(1/real(y)), float(-1/imag(y)))
        if self.two_port is not None:
            x = float(-1/imag(y))
            w = 2*pi*self.two_port.f[0]
            if x > 0:
                hf = 'H'
                lc = x/w
            else:
                hf = 'F'
                lc = 1/(w*x)

            val, prefix = to_eng_form(abs(lc))
            data += '({:>3.2f} {}{}) '.format(val, prefix, hf)
        data += '$\\Omega$'

        if self.gain_cursor and self.two_port is not None:
            gain_str = '{:>7.2f} dB'.format(float(db(gain))) if gain > 0 else 'n/a'
            data += '\nG = ' + gain_str

        self.cursor_text.set_text(data)
        self.ax.draw_artist(self.cursor_text)

    def plot_s_param(self, s, color='blue'):
        self.plot_line(s, color)

        self.save_background()

    def plot_line(self, points, color='blue'):
        r, t = abs(points), angle(points)

        l = Line2D(r*cos(t), r*sin(t), color=color)
        self.ax.add_line(l)

    def save_background(self):
        self.fig.canvas.draw()
        self.background = self.fig.canvas.copy_from_bbox(self.ax.bbox)

    def _on_mouse_leave(self, event):
        pass

    def _resize(self, event):
        self.save_background()
        print 'resize\n'
        #self.ax.draw(None)
        #self.fig.canvas.draw()
        #self.fig.canvas.restore_region(self.background, bbox=self.ax.bbox)

    def _draw(self, event):
        self.fig.canvas.mpl_disconnect(self.cid_draw)
        #self.fig.canvas.restore_region(self.background, bbox=self.ax.bbox)
        self.save_background()
        self.cid_draw = self.fig.canvas.mpl_connect('draw_event', self._draw)


# TODO: fix cursor for s param plot

def main():
    f = linspace(1e6, 2e9, 2000)

    load = OnePort(z=75/2)
    source = OnePort(z=50)

    figure()


    connector = TransmissionLine(2/100, 0.66, z0=75, f=f)
    cable = TransmissionLine(1, 0.66, z0=50, f=f)

    net = cable

    plot(f/1e6, db(net.g_t(source, load)), label='no adapter')

    net = cable*connector

    plot(f/1e6, db(net.g_t(source, load)), label='w/ 2 cm 75 $\Omega$ adapter')



    connector = TransmissionLine(15/100, 0.66, z0=75, f=f)

    net = cable*connector

    plot(f/1e6, db(net.g_t(source, load)), label='w/ 15 cm 75 $\Omega$ pigtail')


    ylabel('Gain (dB)')
    xlabel('Frequency (MHz)')
    legend(loc=3)


    figure()
    sc = SmithChart(labels=1, show_cursor=0)
    sc.plot_s_param((net*load).inp().s)

    show()


if __name__ == '__main__':
    main()
