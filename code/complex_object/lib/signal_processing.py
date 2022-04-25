import copy
import numpy as np
from scipy.sparse import coo_matrix


def extract_xq(X, q, wq, dimq):
    q = frozenset(q)
    wq = frozenset(wq)
    x = copy.deepcopy(X)
    xq = [None]*len(x)

    for d in range(len(xq)):
        if d == dimq:
            del x[d][q]
            if len(x[d]) > 0:
                xq[d] = dict()
                for k, cell in enumerate(x[d].keys()):
                    xq[d][cell] = k

        elif d == dimq+1:
            del x[d][wq]
            if len(x[d]) > 0:
                xq[d] = dict()
                for k, cell in enumerate(x[d].keys()):
                    xq[d][cell] = k

        else:
            xq[d] = x[d]

    cells_xq = list(filter(None, xq))
    return cells_xq


def extract_total_order(X):
    tot_order = dict()
    k = 0
    for d in range(len(X)):
        for cell in X[d].keys():
            tot_order[cell] = k
            k = k+1

    return tot_order


def build_Q_boundary(cellsXq, cellsX, q, wq, dimq, kX):
    q = frozenset(q)
    wq = frozenset(wq)
    Q_boundaries = list()

    for d in range(1, len(cellsXq)):

        idx_cells, idx_faces, values = [], [], []
        if d == dimq and d != 0:
            for cell, idx_cell_q in cellsXq[d].items():

                idx_cell_x = cellsX[d][cell]
                kxT = kX[d-1].T.tolil()
                non_zero_idx = kxT.rows[idx_cell_x]

                for j in non_zero_idx:
                    value = kX[d-1].tocsr()[j, idx_cell_x]
                    values.append(value)
                    idx_cells.append(idx_cell_q)
                    idx_faces.append(j)

            boundary = coo_matrix((values, (idx_faces, idx_cells)), dtype=np.float32, shape=(
                len(cellsXq[d-1]), len(cellsXq[d])))

        if d == dimq+1:
            for cell, idx_cell_q in cellsXq[d].items():

                idx_cell_x = cellsX[d][cell]
                kxT = kX[d-1].T.tolil()
                non_zero_idx_cell = kxT.rows[idx_cell_x]
                non_zero_idx_wq = kxT.rows[cellsX[d][wq]]
                non_zero_idx = np.unique(
                    np.append(non_zero_idx_cell, non_zero_idx_wq))
                non_zero_idx = non_zero_idx[non_zero_idx != cellsX[d-1][q]]

                for j in non_zero_idx:
                    KXd = kX[d-1].tocsr()
                    value = KXd[j, idx_cell_x]-(KXd[cellsX[d-1][q], idx_cell_x]
                                                * KXd[j, cellsX[d][wq]])/(KXd[cellsX[d-1][q], cellsX[d][wq]])

                    if value != 0:
                        values.append(value)
                        idx_cells.append(idx_cell_q)
                        idx_faces.append(
                            cellsXq[d-1][list(cellsX[d-1].keys())[j]])

                boundary = coo_matrix((values, (idx_faces, idx_cells)),
                                      dtype=np.float32,
                                      shape=(len(cellsXq[d-1]), len(cellsXq[d])))

        if d == dimq+2:
            for cell, idx_cell_q in cellsXq[d].items():

                idx_cell_x = cellsX[d][cell]
                kxT = kX[d-1].T.tolil()
                non_zero_idx = kxT.rows[idx_cell_x]
                non_zero_idx = np.array(non_zero_idx)
                non_zero_idx = non_zero_idx[non_zero_idx != cellsX[d-1][wq]]

                for j in non_zero_idx:
                    value = kX[d-1].tocsr()[j, idx_cell_x]
                    if value != 0:
                        values.append(value)
                        idx_cells.append(idx_cell_q)
                        idx_faces.append(
                            cellsXq[d-1][list(cellsX[d-1].keys())[j]])

            boundary = coo_matrix((values, (idx_faces, idx_cells)), dtype=np.float32, shape=(
                len(cellsXq[d-1]), len(cellsXq[d])))

        if d != dimq and d != dimq+1 and d != dimq+2:
            boundary = kX[d-1]

        Q_boundaries.append(boundary)

    return Q_boundaries


def psi(X, Xq, q, wq, dimq, kX):
    q = frozenset(q)
    wq = frozenset(wq)
    tot_X = extract_total_order(X)
    tot_Xq = extract_total_order(Xq)

    ind_x, ind_q, values = [], [], []

    for cell in tot_X.keys():
        if cell == q:
            kxT = kX[dimq].T.tolil()
            idx_wq = X[dimq+1][wq]
            idx_q = X[dimq][q]
            non_zero_idx = kxT.rows[idx_wq]
            non_zero_idx = np.array(non_zero_idx)
            non_zero_idx = non_zero_idx[non_zero_idx != X[dimq][q]]

            for idx_face in non_zero_idx:
                KXd = kX[dimq].tocsr()
                value = -(KXd[idx_face, idx_wq])/KXd[idx_q, idx_wq]
                values.append(value)
                ind_x.append(tot_X[cell])
                face = tot_Xq[list(X[dimq].keys())[idx_face]]
                ind_q.append(face)

        if cell != q and cell != wq:

            values.append(1)
            ind_x.append(tot_X[cell])
            ind_q.append(tot_Xq[cell])

    psi = coo_matrix((values, (ind_q, ind_x)), dtype=np.float32,
                     shape=(len(tot_Xq), len(tot_X)))

    return psi


def phi(X, Xq, q, wq, dimq, kX):

    q = frozenset(q)
    wq = frozenset(wq)
    tot_X = extract_total_order(X)
    tot_Xq = extract_total_order(Xq)

    ind_x, ind_q, values = [], [], []

    kxT = kX[dimq].tolil()
    idx_wq = X[dimq+1][wq]
    idx_q = X[dimq][q]
    non_zero_idx = kxT.rows[idx_q]
    non_zero_idx = np.array(non_zero_idx)
    non_zero_idx = non_zero_idx[non_zero_idx != idx_wq]

    cofaces_q = [list(X[dimq+1].keys())[j] for j in non_zero_idx]

    for cell, idx_cell in zip(cofaces_q, non_zero_idx):
        KXd = kX[dimq].tocsr()
        value = -(KXd[idx_q, idx_cell])/KXd[idx_q, idx_wq]
        values.append(value)
        ind_x.append(tot_X[wq])
        ind_q.append(tot_Xq[cell])

    for cell in tot_Xq:
        values.append(1)
        ind_x.append(tot_X[cell])
        ind_q.append(tot_Xq[cell])

    phi = coo_matrix((values, (ind_x, ind_q)), dtype=np.float32,
                     shape=(len(tot_X), len(tot_Xq)))

    return phi


def energy(X, Xq, q, wq, dimq, kX):

    s = phi(X, Xq, q, wq, dimq, kX)
    r = psi(X, Xq, q, wq, dimq, kX)

    return np.linalg.norm(np.identity(s.shape[0]) - s @ r)


def sequence_collapses(X, kX, collapses, dim_collapses, short=False):
    all_psi, all_phi, all_xq, all_boundaries, all_complexes = [], [], [], [], []
    all_boundaries.append(kX)
    all_complexes.append(X)
    for j in range(len(collapses)):
        q = collapses[j][0]
        wq = collapses[j][1]
        x = all_complexes[-1]
        kx = all_boundaries[-1]
        xq = extract_xq(x, q, wq, dim_collapses[j])
        kq = build_Q_boundary(xq, x, q, wq, dim_collapses[j], kx)
        phi1 = phi(x, xq, q, wq, dim_collapses[j], kx)
        psi1 = psi(x, xq, q, wq, dim_collapses[j], kx)

        all_xq.append(xq)
        all_phi.append(phi1)
        all_psi.append(psi1)
        all_complexes.append(xq)
        all_boundaries.append(kq)

    if short == False:
        return all_psi, all_phi, all_boundaries, all_complexes, all_xq
    else:
        return all_psi, all_phi


def energy_sequence(all_psi, all_phi):

    psi1 = all_psi[-1]
    for j in range(2, len(all_psi)+1):
        psi1 = psi1 @ all_psi[-j]

    phi1 = all_phi[0]
    for j in range(1, len(all_phi)):
        phi1 = phi1 @ all_phi[j]

    energy = np.linalg.norm(np.identity(all_psi[0].shape[1]) - phi1@psi1)
    matrix_energy = coo_matrix(np.identity(all_psi[0].shape[1]) - phi1@psi1)
    return energy, matrix_energy


def phipsi(X, kX, collapses, dim_collapses, signal, type_collapse='up'):
    c = collapses
    d = dim_collapses

    all_psi, all_phi = sequence_collapses(
        X, kX, collapses=c, dim_collapses=d, short=True)

    psi1 = all_psi[-1]
    for j in range(2, len(all_psi)+1):
        psi1 = psi1 @ all_psi[-j]

    phi1 = all_phi[0]
    for j in range(1, len(all_phi)):
        phi1 = phi1 @ all_phi[j]

    matrix = phi1@psi1

    if type_collapse == 'up':
        a = 0
        if dim_collapses[0] > 0:
            a = len(X[dim_collapses[0]-1])
        b = a+len(X[dim_collapses[0]])
        matrix1 = matrix[a:b, a:b]
        phipsis = matrix1@np.array(signal).T

    if type_collapse == 'down':

        a = len(X[dim_collapses[0]])
        b = a+len(X[dim_collapses[0]+1])
        matrix1 = matrix[a:b, a:b]

        phipsis = matrix1.T@np.array(signal).T

    return phipsis


def loss_signal(X, kX, collapses, dim_collapses, signal, type_collapse='up'):
    c = collapses
    d = dim_collapses

    all_psi, all_phi, all_boundaries, all_complexes, all_xq = sequence_collapses(
        X, kX, collapses=c, dim_collapses=d)

    psi1 = all_psi[-1]
    for j in range(2, len(all_psi)+1):
        psi1 = psi1 @ all_psi[-j]

    phi1 = all_phi[0]
    for j in range(1, len(all_phi)):
        phi1 = phi1 @ all_phi[j]

    if type_collapse == "up":
        matrix_energy = np.identity(all_psi[0].shape[1]) - phi1@psi1
        a = 0
        if dim_collapses[0] > 0:
            a = len(X[dim_collapses[0]-1])
        b = a+len(X[dim_collapses[0]])
        matrix_energy1 = matrix_energy[a:b, a:b]

    if type_collapse == 'down':

        matrix_energy = np.identity(all_psi[0].shape[1]) - (phi1@psi1).T
        a = len(X[dim_collapses[0]])
        b = a+len(X[dim_collapses[0]+1])

        matrix_energy1 = matrix_energy[a:b, a:b]

    loss = np.linalg.norm(matrix_energy1@np.array(signal).T)**2
    return loss
