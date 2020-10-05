"""Elmer import and export."""


from skfem.mesh import Mesh, MeshQuad


def to_file(mesh: Mesh, filename: str):
    """The mesh is written to four files.

    The names of the files are filename.{header,nodes,elements,boundary}.

    Parameters
    ----------
    mesh
        The mesh to export.
    filename
        The prefix of the filenames.

    """
    np = mesh.p.shape[1]
    nt = mesh.t.shape[1]

    # filename.header
    with open(filename + '.header', 'w') as handle:
        handle.write("{} {} {}\n".format(np,
                                         nt,
                                         len(mesh.boundary_facets())))
        handle.write("2\n")
        if isinstance(mesh, MeshQuad):
            handle.write("404 {}\n".format(nt))
            handle.write("202 {}\n".format(len(mesh.boundary_facets())))

    # filename.nodes
    with open(filename + '.nodes', 'w') as handle:
        for itr in range(np):
            handle.write("{} -1 {} {} 0.0\n".format(itr + 1,
                                                    *mesh.p[:, itr],
                                                    0))

    # filename.elements
    with open(filename + '.elements', 'w') as handle:
        for itr in range(nt):
            handle.write("{} 1 404 {} {} {} {}\n".format(itr + 1,
                                                         mesh.t[0, itr] + 1,
                                                         mesh.t[1, itr] + 1,
                                                         mesh.t[2, itr] + 1,
                                                         mesh.t[3, itr] + 1))

    # filename.boundary
    with open(filename + '.boundary', 'w') as handle:
        for itr in mesh.boundary_facets():
            handle.write("{} 1 {} {} 202 {} {}\n".format(
                itr + 1,
                mesh.f2t[0, itr] + 1,
                mesh.f2t[1, itr] + 1,
                mesh.facets[0, itr] + 1,
                mesh.facets[1, itr] + 1
            ))
