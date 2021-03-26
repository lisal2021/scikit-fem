from skfem import *


class TimeMeshTri:

    def setup(self):
        self.m = MeshTri()

    def time_refine_tri(self):
        m = self.m.refined(8)


class PeakmemMeshTri:

    def peakmem_refine_tri(self):
        return MeshTri().refined(8)
