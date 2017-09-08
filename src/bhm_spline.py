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
            line=bhm_out.readline()
            self._xmin, self._xmax, order = np.fromstring(line, sep=" ")
            self._order=int(order)
            line=bhm_out.readline()
            n_pieces=int(line)
            self._knots_left=np.empty(n_pieces)
            self._knots_right=np.empty(n_pieces)
            self._spline_pieces=[]
            self._errorbar_pieces=[]
            for i in range(0, n_pieces):
                header=bhm_out.readline() # discarded

                line=bhm_out.readline()
                self._knots_left[i], self._knots_right[i] = np.fromstring(line, sep=" ")

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

    def __knots_ok(self):
        """Check that the spline pieces knots do not have gaps"""
        if np.any(self._knots_left[1:]!=self._knots_right[:-1]):
            print("ERROR! There is a gap between a right knot and a left knot",file=sys.stderr)
            print("Left knots:",self._knots_left,file=sys.stderr)
            print("Right knots:",self._knots_right,file=sys.stderr)
            return False

        if self._xmin!=self._knots_left[0]:
            print("ERROR! There is a gap between a leftmost knot and a minimum x",file=sys.stderr)
            print("Leftmost knot=",self._knots_left[0],file=sys.stderr)
            print("Minimum x=",self._xmin,file=sys.stderr)
            return False

        if self._xmax!=self._knots_right[-1]:
            print("ERROR! There is a gap between a rightmost knot and a maximum x",file=sys.stderr)
            print("Rightmost knot=",self._knots_right[-1],file=sys.stderr)
            print("Maximum x=",self._xmax,file=sys.stderr)
            return False

        if np.any(self._knots_left >= self._knots_right):
            print("ERROR! Some right knots are lesser than left knots",file=sys.stderr)
            return False

        return True



    def __init__(self, bhm_out_name):
        """Initialize the BHM object by reading the histogram from the given file name"""
        self.__read_file(bhm_out_name)
        # Sanity checks:
        if not self.__knots_ok():
            raise ValueError("Invalid data in the histogram file")
        return



    def domain(self):
        """Return the domain of the spline as a tuple (xmin,xmax)"""
        return (self._xmin, self._xmax)



    def __call__(self,x):
        """Return the value(s) of the BHM spline at point(s) x"""
        if np.min(x)<self._xmin or np.max(x)>self._xmax:
            raise ValueError("Argument is outside the spline domain")
        # prepare the boolean list for choosing spline pieces:
        cond=[np.logical_and(x>=r[0],x<=r[1]) for r in zip(self._knots_left, self._knots_right)]
        return np.piecewise(x, cond, self._spline_pieces)


    def errorbar(self, x):
        """Return the value(s) of the BHM spline error bar at point(s) x"""
        if np.min(x)<self._xmin or np.max(x)>self._xmax:
            raise ValueError("Argument is outside the spline domain")
        # prepare the boolean list for choosing spline pieces:
        cond=[np.logical_and(x>=r[0],x<=r[1]) for r in zip(self._knots_left, self._knots_right)]
        return np.piecewise(x, cond, self._errorbar_pieces)



    def plot(self, ref_fn=None):
        """Convenience method: plot the spline and the reference function ref_fn, if supplied"""
        x=np.linspace(self._xmin, self._xmax)
        y=self(x)
        yerr=self.errorbar(x)

        ax=plt.axes()

        ax.plot(x,y,label="spline")
        ax.fill_between(x, y-yerr,y+yerr, alpha=0.3)
        if ref_fn:
            ax.plot(x, ref_fn(x), label="reference")
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
