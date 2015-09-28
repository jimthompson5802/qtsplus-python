# umbrella module for QtPy
"""
Collection of queueing theory models.  The models are solved by a variety
of analytic and discrete event simulation techniques.

Unless otherwise specified, the implementation of the analytic models
is based on material found in "Fundamentals of Queueing Theory, 4th Ed"
by D. Gross, J. Shortle, J. Thompson, C. Harris, 2008 (GSTH).

The following Python packages must be installed for the correct operation
of QtPy:
    NumPy  (http://numpy.scipy.org)
    SciPy  (http://www.scipy.org)
    SimPy  (http://simpy.sourceforge.net)
"""

from SingleServer import *
from MultiServer import *
from Network import *
from Simulation import *
from QtPyCommon import *

