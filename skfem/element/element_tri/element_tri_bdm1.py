import numpy as np

from ..element_hdiv import ElementHdiv
from ...refdom import RefTri
from ..discrete_field import DiscreteField


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

    # def gbasis(self, mapping, X, i, tind=None):
    #     phi, dphi = self.lbasis(X, i)
    #     diff = 1 if i % 2 == 0 else -1
    #     Phi, dPhi = self.lbasis(X, i + diff)
    #     DF = mapping.DF(X, tind)
    #     detDF = mapping.detDF(X, tind)
    #     orient = self.orient(mapping, i, tind)
    #     value1 = np.einsum('ijkl,jl,kl->ikl', DF, phi,
    #                        1. / np.abs(detDF) * orient[:, None])
    #     value2 = np.einsum('ijkl,jl,kl->ikl', DF, Phi,
    #                         1. / np.abs(detDF) * orient[:, None])
    #     div1 = dphi / (np.abs(detDF) * orient[:, None])
    #     div2 = dPhi / (np.abs(detDF) * orient[:, None])
    #     ixs1=np.nonzero(orient == 1)[0]
    #     ixs2=np.nonzero(orient == -1)[0]
    #     Value=np.zeros_like(value1)
    #     for itr in range(2):
    #         Value[itr, ixs1] = value1[itr, ixs1]
    #         Value[itr, ixs2] = value2[itr, ixs2]
    #     Div=np.zeros_like(div1)
    #     Div[ixs1] = div1[ixs1]
    #     Div[ixs2] = div2[ixs2]
    #     #import pdb; pdb.set_trace()
    #     return (DiscreteField(
    #         value=Value,
    #         div=Div,
    #     ),)

    def lbasis(self, X, i):
        x, y = X

        if i == 0:
            phi = -1. / (G_2 - G_1) * np.array([(G_1 - 1) * x, x + G_1 * y - G_1])
            dphi = -1. / (G_2 - G_1) * (2. * G_1 - 1) + 0. * x
        elif i == 1:
            phi = -1. / (G_1 - G_2) * np.array([(G_2 - 1) * x, x + G_2 * y - G_2])
            dphi = -1. / (G_1 - G_2) * (2. * G_2 - 1) + 0. * x
        elif i == 2:
            phi = np.sqrt(2) / (G_1 - G_2) * np.array([G_1 * x, (G_1 - 1) * y])
            dphi = np.sqrt(2) / (G_1 - G_2) * (2. * G_1 - 1) + 0. * x
        elif i == 3:
            phi = np.sqrt(2) / (G_2 - G_1) * np.array([G_2 * x, (G_2 - 1) * y])
            dphi = np.sqrt(2) / (G_2 - G_1) * (2. * G_2 - 1) + 0. * x
        elif i == 4:
            phi = 1. / (G_1 - G_2) * np.array([G_1 * x + y - G_1, (G_1 - 1) * y])
            dphi = 1. / (G_1 - G_2) * (2. * G_1 - 1) + 0. * x
        elif i == 5:
            phi = 1. / (G_2 - G_1) * np.array([G_2 * x + y - G_2, (G_2 - 1) * y])
            dphi = 1. / (G_2 - G_1) * (2. * G_2 - 1) + 0. * x
        else:
            self._index_error()

        return phi, dphi
