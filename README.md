PyTwoPort
=========

A Python library for analyzing linear two-port electrical network. Includes useful functionality for aiding in the design of small-signal RF amplifiers (e.g. Smith charts and stability analysis).

Documentation and comments are a bit sparse, but start by looking at test.py (and below).

Features
=========

* Support for one- and two-port networks (N-port functionality can be added at some point).
* Conversion between Z, Y, S, T and ABCD parameters, including spline based interpolation of new frequency points.
* Combine networks in parallel, series or cascade.
* Convert one-port networks to parallel or series two-port networks, and one-port networks to driving point one-port networks.
* Functions for calculating various gains (tranducer gain, maximum stable gain, etc.)
* Functions for computing various stability factors (Rollett's, mu, etc.) as well as the positions of stability circles.
* Ability to import S1P and S2P data files.
* Smith chart plotting, with interactive data cursors.
* Gain and stability circle plotting on Smith charts.

Dependencies
=========

* Python 2.7
* NumPy
* SciPy (used just for interpolation)
* Matplotlib (for plotting)

Usage
=========

A two-port network can be initialized using `TwoPort(f=f, <type>=<matrix>)` where `<type>` is one of the supported
network types (`z`, `y`, `s`, etc.), `<matrix>` is a numpy array and `f` is a numpy array of frequencies (in which case
`<matrix>` should be 3D with a matrix corresponding to each frequency point). The `load_snp(file_name)` can also be
used to import data from a S1P/S2P file.

If you want to build a network from scratch using lumped elements, then you can combine inductors, capacitors,
resistors, etc. in series and parallel like this:

```python
f = linspace(100e6, 1000e6) # need to define frequency points for some devices

L = Inductor(100e-9, f=f)
C = Capacitor(100e-12, f=f)
R = Resistor(50)

network = (L // C).to_series() * R.to_shunt()
```

A schematic of this network is shown [here](http://i.imgur.com/IlgTtqF.png).
The `to_series()` function converts a one-port to a series connected two-port (similarly for `to_shunt()`),
`*` combines two one-ports in cascade  and `//` two networks in parallel (`+` for series).
Additional elements/networks (transformers, transmission lines, etc.) are defined in networks.py.

The two-port parameters are now available in any format, e.g.:

```python
print network.z[_11] # _11 is a helper so we can avoid indexing with z[:, 0, 0]
print network.s[_21]
print network[500e6].y # returns matrix for f=500 MHz
```

If a network is only defined at a few frequency points, PyTwoPort can automatically interpolate new points using:

```python
network = network[400e6:450e6:5e6] # gives 5 MHz resolution
```

Note that the first and second indices can be omitted. Interpolation is done with SciPy's 1D spline function.

Two port networks can also be reduced to equivalent input or output (driving point) one-port networks:

```python
print network.input_equiv()
print network.output_equiv()
```

Here the opposite port remains open circuited.

Other things you can do with a two-port network include calculating gains, stability factors, etc.

```python
terminator = Resistor(150) # this is our source/load termination

print network.transducer_gain(terminator, terminator) # gain with this termination
print network.s[_21] # this would be the gain into 50 Ohm loads
print network.rollett_k() # stability factor (< 1.0 implies potential instability)
```

Sample plots
=========

A few sample plots (see test.py for the code):

* Instability regions for source and load termination for an active two-port network.
![image](http://i.imgur.com/4Jduvhq.png)
* Instability regions for different frequency, and a minimum gain circle.
![image](http://i.imgur.com/z4I0RfC.png)
* Gain surface and instability regions.
![image](http://i.imgur.com/NVDWINe.png)
* Gains vs. frequency
![image](http://i.imgur.com/xw9HkZ5.png)
* Reflection coefficient on Smith chart
![image](http://i.imgur.com/Oh4HyTW.png)
* Stability factors vs. frequency
![image](http://i.imgur.com/LwxYfIt.png)

License
=========
All of the code contained here is licensed by the GNU General Public License v3.

Copyright (C) 2014 by Roger https://github.com/roger-
