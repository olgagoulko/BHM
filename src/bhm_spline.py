#!/usr/bin/env python3
"""Module for evaluating and plotting a BHM spline

Usage::

    import numpy as np
    import matplotlib.pyplot as plt
    from bhm_spline import BHMSpline
    
    spline=BHMSpline("bhm_output.dat")
    x=np.linspace(*spline.domain())
    y=spline(x)
    plt.plot(x,y)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial as Poly


class BHMSpline:
    """The class represents a BHM spline. Can be evaluated and plotted"""

    def __read_file(self, bhm_out_name):
        """Read the file and fill relevant object variables"""

        with open(bhm_out_name,"r") as bhm_out:
            # skip initial comment lines
            while True:
                line=bhm_out.readline()
                if line[0]!='#' : break
                
            order, n_pieces = np.fromstring(line, sep=" ")
            self._order=int(order)
            n_pieces=int(n_pieces)

            line=bhm_out.readline()
            self._knots = np.fromstring(line, sep=" ")
            if self._knots.size != n_pieces+1:
                raise ValueError("Incorrect number of knots %d" % self._knots.size)

            self._spline_pieces=[]
            self._errorbar_pieces=[]
            for i in range(0, n_pieces):
                header=bhm_out.readline() # discarded

                line=bhm_out.readline()
                coeffs = np.fromstring(line, sep=" ")
                if coeffs.size != self._order+1:
                    raise ValueError("Spline piece #%d: incorrect number of coefficients %d" % (i,coeffs.size))
                self._spline_pieces.append(Poly(coeffs))

                line=bhm_out.readline()
                # self._bar_coeff[i,:] = np.fromstring(line, sep=" ")
                bar_coeffs = np.fromstring(line, sep=" ")
                if bar_coeffs.size != 2*self._order+1:
                    raise ValueError("Spline piece #%d: incorrect number of error bar coefficients %d" % (i,bar_coeffs.size))
                self._errorbar_pieces.append(Poly(bar_coeffs))

        return

    def __check_knots(self):
        """Check that the spline knots are ordered properly"""
        if not np.all(self._knots[:-1]<self._knots[1:]):
            raise ValueError("Some left knots are greater or equal right knots")
    
    def __init__(self, bhm_out_name):
        """Initialize the BHM object by reading the histogram from the given file name"""
        self.__read_file(bhm_out_name)
        # Sanity checks:
        self.__check_knots()
        return


    def domain(self):
        """Return the domain of the spline as a tuple (xmin,xmax)"""
        return tuple(self._knots[[0,-1]])


    def __call__(self,x):
        """Return the value(s) of the BHM spline at point(s) x"""
        if np.min(x)<self._knots[0] or np.max(x)>self._knots[-1]:
            raise ValueError("Argument is outside the spline domain")
        # prepare the boolean list for choosing spline pieces:
        cond=[np.logical_and(x>=r[0],x<=r[1]) for r in zip(self._knots[:-1], self._knots[1:])]
        return np.piecewise(x, cond, self._spline_pieces)


    def errorbar(self, x):
        """Return the value(s) of the BHM spline error bar at point(s) x"""
        if np.min(x)<self._knots[0] or np.max(x)>self._knots[-1]:
            raise ValueError("Argument is outside the spline domain")
        # prepare the boolean list for choosing spline pieces:
        cond=[np.logical_and(x>=r[0],x<=r[1]) for r in zip(self._knots[:-1], self._knots[1:])]
        return np.sqrt(np.piecewise(x, cond, self._errorbar_pieces))


    def plot(self, ref_fn=None):
        """Convenience method: plot the spline and the reference function ref_fn, if supplied"""
        x=np.linspace(*(self.domain()+(1000,)))
        y=self(x)
        yerr=self.errorbar(x)

        ax=plt.axes()

        ax.plot(x,y,label="spline")
        ax.fill_between(x, y-yerr,y+yerr, alpha=0.3)
        if ref_fn:
            ax.plot(x, ref_fn(x), label="reference")
        ax.legend()
        return

    def plot_difference(self, ref_fn=None):
        """Convenience method: plot the spline and the reference function ref_fn, if supplied"""
        x=np.linspace(*(self.domain()+(1000,)))
        y=self(x)
        yerr=self.errorbar(x)

        ax=plt.axes()
        if ref_fn:
            ax.plot(x,y-ref_fn(x),label="spline - reference")
            ax.fill_between(x, y-ref_fn(x)-yerr,y-ref_fn(x)+yerr, alpha=0.3)

        ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
        ax.legend()
        return


def Main(argv):
    """Plot the interpolation spline generated by BHM program.

    Usage: %s bhm_output
    """
    if len(argv)!=2:
        print("Usage: \n", Main.__doc__ % argv[0], file=sys.stderr)
        return 2

    spline=BHMSpline(sys.argv[1])
    
    spline.plot()
    plt.show()
    return 0


if __name__ == "__main__": sys.exit(Main(sys.argv))
