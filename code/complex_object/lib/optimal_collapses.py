import numpy as np
import random
from lib.signal_processing import sequence_collapses, energy_sequence, loss_signal, phipsi
from scipy import sparse
from lib.build_boundaries import build_up_down_laplacians


def best_up_collapse(X, kX, dimq, signal):
    s = np.array(signal)**2
    Bq = kX[dimq]
    BT = Bq.transpose()
    C = BT.tolil()
    ind = C.rows

    possible_min_q, possible_min_value = [], []
    for wq in range(Bq.shape[1]):

        boundary = ind[wq]
        sig_wq = s[boundary]
        d = np.array(C.data[wq])**2
        sig_wq = sig_wq/d
        notnan_sig = sig_wq[~np.isnan(sig_wq)]
        min_sq = np.min(notnan_sig)
        norm_dwq = np.linalg.norm(Bq.tocsr()[ind[wq], wq].todense())**2
        possible_min_value.append(min_sq*norm_dwq)

        minimum = random.choice(np.where(notnan_sig == min_sq)[0])

        idx_q_min = ind[wq][minimum]
        possible_min_q.append(idx_q_min)

    global_min = np.min(possible_min_value)

    Wq = random.choice(np.where(possible_min_value == global_min)[0])
    q = possible_min_q[Wq]

    return [[list(list(X[dimq].keys())[q]), list(list(X[dimq+1].keys())[Wq])]], global_min


def best_down_collapse(X, kX, dimq, signal):
    s = np.array(signal)**2
    Bq = kX[dimq]
    BT = Bq.copy()
    C = BT.tolil()
    ind = C.rows

    possible_min_WQ, possible_min_value = [], []

    for Q in range(Bq.shape[0]):
        boundary = ind[Q]

        sig_WQ = s[boundary]
        d = np.array(C.data[Q])**2
        sig_WQ = sig_WQ/d
        notnan_sig = sig_WQ[~np.isnan(sig_WQ)]
        min_sWQ = np.min(notnan_sig)

        norm_dQ = np.linalg.norm(Bq.tocsr()[Q, ind[Q]].todense())**2
        possible_min_value.append(min_sWQ*norm_dQ)

        minimum = random.choice(np.where(notnan_sig == min_sWQ)[0])

        idx_WQ_min = ind[Q][minimum]
        possible_min_WQ.append(idx_WQ_min)

    global_min = np.min(possible_min_value)
    q = random.choice(np.where(possible_min_value == global_min)[0])
    Wq = possible_min_WQ[q]

    return [[list(list(X[dimq].keys())[q]), list(list(X[dimq+1].keys())[Wq])]], global_min


def random_collapse(X, kX, dimq):
    nz = np.array(kX[dimq].nonzero())
    rand = np.random.choice(np.arange(len(nz[0])))
    Wq = nz[:, rand][1]
    q = nz[:, rand][0]
    return([[list(list(X[dimq].keys())[q]), list(list(X[dimq+1].keys())[Wq])]])


def sequence_optimal_up_collapses(X, kX, dimq, signal, steps, random=False):
    dX = kX.copy()
    all_X = [X]
    all_signals = [signal]
    all_collapses, all_losses, total_psi, total_phi, total_loss = [], [], [], [], []

    a1 = 0
    if dimq > 0:
        a1 = len(X[dimq-1])
    b1 = a1+len(X[dimq])

    for k in range(steps):
        if random == False:
            c = best_up_collapse(X, kX, dimq=dimq, signal=signal)[0]
        else:
            c = random_collapse(X, kX, dimq)

        all_psi, all_phi, all_boundaries, all_complexes, all_xq = sequence_collapses(
            X, kX, collapses=c, dim_collapses=[dimq])

        a = len(X[dimq-1])
        b = a+len(X[dimq])

        matrix_energy = energy_sequence(all_psi, all_phi)[
            1].toarray()[a:b, a:b]
        loss = np.linalg.norm(matrix_energy@np.array(signal).T)

        psi1 = all_psi[-1]
        for j in range(2, len(all_psi)+1):
            psi1 = psi1 @ all_psi[-j]

        total_phi.append(all_phi[-1])
        total_psi.append(all_psi[-1])
        totmatrix_energy = energy_sequence(total_psi, total_phi)[
            1].toarray()[a1:b1, a1:b1]

        totloss = np.linalg.norm(totmatrix_energy@np.array(all_signals[0]).T)
        total_loss.append(totloss)

        signal = psi1.toarray()[a:b-1, a:b]@np.array(signal).T

        kX = all_boundaries[-1]
        X = all_xq[-1]

        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)
        try:
            if len(kX) < (dimq) or np.sum(np.abs(kX[dimq])) == 0 or len(X) < (dimq+1):
                break
        except:
            pass

    phispsis = all_signals[0]-totmatrix_energy@np.array(all_signals[0]).T

    return(all_X, all_collapses, all_losses, total_loss, all_signals, phispsis)


def sequence_optimal_down_collapses(X, kX, dimq, signal, steps, random=False):
    dX = kX.copy()
    all_X = [X]
    all_signals = [signal]
    all_collapses, all_losses, total_psi, total_phi, total_loss = [], [], [], [], []

    a1 = len(X[dimq])
    b1 = a1+len(X[dimq+1])

    for k in range(steps):
        if random == False:
            c = best_down_collapse(X, kX, dimq=dimq, signal=signal)[0]
        if random == True:
            c = random_collapse(X, kX, dimq)

        all_psi, all_phi, all_boundaries, all_complexes, all_xq = sequence_collapses(
            X, kX, collapses=c, dim_collapses=[dimq])

        a = len(X[dimq])
        b = a+len(X[dimq+1])

        matrix_energy = energy_sequence(all_psi, all_phi)[
            1].toarray()[a:b, a:b]
        loss = np.linalg.norm(matrix_energy.T@np.array(signal).T)

        phi1 = all_phi[-1]
        for j in range(2, len(all_phi)+1):
            phi1 = phi1 @ all_phi[-j]

        total_phi.append(all_phi[-1])
        total_psi.append(all_psi[-1])
        totmatrix_energy = energy_sequence(total_psi, total_phi)[
            1].toarray()[a1:b1, a1:b1]

        totloss = np.linalg.norm(totmatrix_energy.T@np.array(all_signals[0]).T)
        total_loss.append(totloss)

        phiT = phi1.toarray().T

        signal = phiT[a-1:b-2, a:b]@np.array(signal).T

        kX = all_boundaries[-1]

        X = all_xq[-1]

        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)
        if np.sum(np.abs(kX[dimq])) == 0 or len(X) < dimq:
            break

    phispsis = all_signals[0]-totmatrix_energy.T@np.array(all_signals[0]).T

    return(all_X, all_collapses, all_losses, total_loss, all_signals, phispsis)


def sequence_given_up_collapses(X, kX, dimq, signal, collapses):
    dX = kX.copy()
    all_X = [X]
    all_signals = [signal]
    all_collapses = []
    all_losses = []

    for collapse in collapses:
        c = [collapse]
        all_psi, all_phi, all_boundaries, all_complexes, all_xq = sequence_collapses(
            X, kX, collapses=c, dim_collapses=[dimq])

        a = len(X[dimq-1])
        b = a+len(X[dimq])
        L = energy_sequence(all_psi, all_phi)[1].toarray()[a:b, a:b]
        loss = loss_signal(X, kX, collapses=c, dim_collapses=[
                           dimq], signal=signal)

        psi1 = all_psi[-1]
        for j in range(2, len(all_psi)+1):
            psi1 = psi1 @ all_psi[-j]

        signal = psi1.toarray()[a:b-1, a:b]@np.array(signal).T

        kX = all_boundaries[-1]
        X = all_xq[0]

        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)

    dimc = [dimq]*len(all_collapses)

    total_loss = loss_signal(
        all_X[0], dX, collapses=all_collapses, dim_collapses=dimc, signal=all_signals[0])
    phispsis = phipsi(all_X[0], dX, collapses=all_collapses,
                      dim_collapses=dimc, signal=all_signals[0])
    return(all_X, all_collapses, all_losses, total_loss, all_signals, phispsis)


def sequence_given_down_collapses(X, kX, dimq, signal, collapses):
    all_X = [X]
    dX = kX.copy()
    all_signals = [signal]
    all_collapses = []
    all_losses = []

    for collapse in collapses:
        c = [collapse]
        all_psi, all_phi, all_boundaries, all_complexes, all_xq = sequence_collapses(
            X, kX, collapses=c, dim_collapses=[dimq])

        a = len(X[dimq])
        b = a+len(X[dimq+1])

        loss = loss_signal(X, kX, collapses=c, dim_collapses=[
                           dimq], signal=signal, type_collapse='down')

        phi1 = all_phi[-1]
        for j in range(2, len(all_phi)+1):
            phi1 = phi1 @ all_phi[-j]

        phiT = phi1.toarray().T
        signal = phiT[a-1:b, a:b]@np.array(signal).T

        kX = all_boundaries[-1]

        X = all_xq[-1]

        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)

    dimc = [dimq]*len(all_collapses)
    total_loss = loss_signal(all_X[0], dX, collapses=all_collapses,
                             dim_collapses=dimc, signal=all_signals[0], type_collapse='down')
    phispsis = phipsi(all_X[0], dX, collapses=all_collapses,
                      dim_collapses=dimc, signal=all_signals[0], type_collapse='down')
    return(all_X, all_collapses, all_losses, total_loss, all_signals, phispsis)


def simulation_collapses(X, kX, dimq, signal, steps, random):
    total_losses = []
    for k in steps:
        all_X, collapses, all_losses, total_loss, all_signals, phispsis = sequence_optimal_up_collapses(
            X=X, kX=kX, dimq=dimq, signal=signal, steps=k, random=random)
        total_losses.append(total_loss)
    return(total_losses, phispsis)


def check_hodge_decomp(X, s1, kX, phispsis, trange, type_collapse='up'):
    boundaries = kX
    ups, downs, laplacians = build_up_down_laplacians(boundaries)

    down = downs[1]
    lap = laplacians[1]
    vh, vech = sparse.linalg.eigsh(lap, 30, which='SM')
    basis_h = vech[:, np.where(vh < 10**(-6))[0]]
    d = basis_h.shape

    if len(ups) > 1:
        up = ups[1]
        vup, vecup = np.linalg.eigh(up.toarray(), UPLO='L')
        if trange == None:
            uprange = len(np.where(vup > 10**(-6))[0])
        if trange != None:
            uprange = trange
        basis_up = vecup[:, np.where(vup > 10**(-6))[0]][:, :uprange]

    vdown, vecdown = np.linalg.eigh(down.toarray(), UPLO='L')
    if trange == None:

        downrange = len(np.where(vdown > 10**(-6))[0])
    if trange != None:

        downrange = trange

    basis_down = vecdown[:, np.where(vdown > 10**(-6))[0]][:, :downrange]

    hodge_basis = basis_h
    if type_collapse == "up":
        hodge_basis = np.hstack((basis_h, basis_down))
        c = hodge_basis.shape
        if len(ups) > 1:
            hodge_basis = np.hstack((hodge_basis, basis_up))

    if type_collapse == "down":
        if len(ups) > 1:
            hodge_basis = np.hstack((basis_h, basis_up))
        hodge_basis = np.hstack((hodge_basis, basis_down))

    h_dec = hodge_basis.T@s1
    h_dec_reconstruction = hodge_basis.T@(phispsis)

    return h_dec, h_dec_reconstruction, c, d
