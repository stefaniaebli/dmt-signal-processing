import numpy as np
from scipy.sparse import coo_matrix, diags


def extract_simplices(simplex_tree):
    simplices = [dict() for _ in range(simplex_tree.dimension()+1)]
    for simplex, _ in simplex_tree.get_skeleton(simplex_tree.dimension()):
        k = len(simplex)
        simplices[k-1][frozenset(simplex)] = len(simplices[k-1])
    return simplices

def build_laplacians(boundaries):
    laplacians = list()
    up = coo_matrix(boundaries[0] @ boundaries[0].T)
    laplacians.append(up)
    for d in range(len(boundaries)-1):
        down = boundaries[d].T @ boundaries[d]
        up = boundaries[d+1] @ boundaries[d+1].T
        laplacians.append(coo_matrix(down + up))
    down = boundaries[-1].T @ boundaries[-1]
    laplacians.append(coo_matrix(down))
    return laplacians

def build_up_down_laplacians(boundaries):
    laplacians = list()
    ups=[]
    downs=[[]]
    up = coo_matrix(boundaries[0] @ boundaries[0].T)
    ups.append(up)
    laplacians.append(up)
    for d in range(len(boundaries)-1):
        down = boundaries[d].T @ boundaries[d]
        up = boundaries[d+1] @ boundaries[d+1].T
        laplacians.append(coo_matrix(down + up))
        ups.append(coo_matrix(up))
        downs.append(coo_matrix(down))
    down = boundaries[-1].T @ boundaries[-1]
    downs.append(down)
    laplacians.append(coo_matrix(down))
    return ups, downs,laplacians

def build_unweighted_boundaries(simplices):
    boundaries = list()
    for d in range(1, len(simplices)):
        idx_simplices, idx_faces, values = [], [], []
        for simplex, idx_simplex in simplices[d].items():
            for i, left_out in enumerate(np.sort(list(simplex))):
                idx_simplices.append(idx_simplex)
                values.append((-1)**i)
                face = simplex.difference({left_out})
                idx_faces.append(simplices[d-1][face])
        assert len(values) == (d+1) * len(simplices[d])
        boundary = coo_matrix((values, (idx_faces, idx_simplices)),
                                     dtype=np.float32,
                                     shape=(len(simplices[d-1]), len(simplices[d])))
        boundaries.append(boundary)


    return boundaries


def build_weighted_boundaries(simplices,weights):
    boundaries = list()
    for d in range(1, len(simplices)):
        idx_simplices, idx_faces, values = [], [], []
        for simplex, idx_simplex in simplices[d].items():
            for i, left_out in enumerate(np.sort(list(simplex))):
                idx_simplices.append(idx_simplex)
                values.append((-1)**i)
                face = simplex.difference({left_out})
                idx_faces.append(simplices[d-1][face])
        assert len(values) == (d+1) * len(simplices[d])
        boundary = coo_matrix((values, (idx_faces, idx_simplices)),
                                     dtype=np.float32,
                                     shape=(len(simplices[d-1]), len(simplices[d])))

        Wn=weights[d]

        w=weights[d-1].data[0]
        nz=np.nonzero(w)[0]
        inv=np.zeros(len(w))
        inv[nz]=(1/(w[nz]))
        inv_Wn=diags(inv)

        boundary=inv_Wn@boundary@Wn
        boundaries.append(boundary)


    return boundaries


def build_boundaries(simplices,weights=None):
    if weights:
        boundaries=build_weighted_boundaries(simplices, weights)

    else:
        boundaries=build_unweighted_boundaries(simplices)

    return boundaries