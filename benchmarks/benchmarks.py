from skfem import *


class TimeMeshTri:

    def setup(self):
        self.m = MeshTri()

    def time_refine_tri(self):
        m = self.m.refined(8)


class PeakmemMeshTet:

    def peakmem_refine_tet(self):
        return MeshTet().refined(5)
