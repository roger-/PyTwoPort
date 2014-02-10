PyTwoPort
=========

A Python library for analyzing linear two-port electrical network. Includes useful functionality for aiding in the design of small-signal RF amplifiers (e.g. Smith charts and stability analysis).

Documentation and comments are a bit sparse, but start by looking at test.py (and below).

Features
=========

* Support for one- and two-port networks (N-port functionality can be added at some point).
* Conversion between Z, Y, S and ABC parameters, including spline based interpolation of new frequency points.
* Combine networks in parallel, series or cascade.
* Convert one-port networks to parallel or series two-port networks, and one-port networks to driving point one-port networks.
* Functions for calculating various gains (tranducer gain, maximum stable gain, etc.)
* Functions for computing various stability factors (Rollett's, mu, etc.) as well as the positions of stability circles.
* Ability to import S1P and S2P data files.
* Smith chart plotting, with interactive data cursors.
* Gain and stability circle plotting on Smith charts.

Samples
=========

A few sample plots (see test.py for the code):

* Instability regions for source and load termination for an active two-port network.
![image](http://i.imgur.com/4Jduvhq.png)
* Instability regions for different frequency, and minimum a gain circle.
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
