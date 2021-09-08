import numpy as np

from ..element_hdiv import ElementHdiv
from ...refdom import RefTri


# two Gauss points used in the definition of the basis
G_1 = .5 - np.sqrt(3) / 6.
G_2 = .5 + np.sqrt(3) / 6.


class ElementTriBDM1(ElementHdiv):
    """The lowest order Brezzi-Douglas-Marini element."""

    facet_dofs = 2
    maxdeg = 2
    dofnames = ['u^n', 'u^n']
    doflocs = np.array([[G_1, .0],
                        [G_2, .0],
                        [G_2, 1. - G_2],
                        [G_1, 1. - G_1],
                        [.0, G_2],
                        [.0, G_1]])
    refdom = RefTri

    def lbasis(self, X, i):
        x, y = X

        if i == 0:
            phi = -1. / (G_1 - G_2) * np.array([(G_2 - 1) * x, x + G_2 * y - G_2])
            dphi = -1. / (G_1 - G_2) * (2. * G_2 - 1) + 0. * x
        elif i == 1:
            phi = -1. / (G_2 - G_1) * np.array([(G_1 - 1) * x, x + G_1 * y - G_1])
            dphi = -1. / (G_2 - G_1) * (2. * G_1 - 1) + 0. * x
        elif i == 2:
            phi = np.sqrt(2) / (G_2 - G_1) * np.array([G_2 * x, (G_2 - 1) * y])
            dphi = np.sqrt(2) / (G_2 - G_1) * (2. * G_2 - 1) + 0. * x
        elif i == 3:
            phi = np.sqrt(2) / (G_1 - G_2) * np.array([G_1 * x, (G_1 - 1) * y])
            dphi = np.sqrt(2) / (G_1 - G_2) * (2. * G_1 - 1) + 0. * x
        elif i == 4:
            phi = 1. / (G_1 - G_2) * np.array([G_1 * x + y - G_1, (G_1 - 1) * y])
            dphi = 1. / (G_1 - G_2) * (2. * G_1 - 1) + 0. * x
        elif i == 5:
            phi = 1. / (G_2 - G_1) * np.array([G_2 * x + y - G_2, (G_2 - 1) * y])
            dphi = 1. / (G_2 - G_1) * (2. * G_2 - 1) + 0. * x
        else:
            self._index_error()

        return phi, dphi
