#!/usr/bin/env python3
import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import inv
import gudhi as gd
import copy
import random
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors

def extract_simplices(simplex_tree):
    """Extract simplices from a gudhi simplex tree.

    Parameters
    ----------
    simplex_tree: gudhi simplex tree

    Returns
    -------
    simplices: List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-simplices. The dictionary's keys are sets (of size d
        + 1) of the 0-simplices that constitute the d-simplices. The
        dictionary's values are the indexes of the simplices in the boundary
        and Laplacian matrices.
    """
    simplices = [dict() for _ in range(simplex_tree.dimension()+1)]
    for simplex, _ in simplex_tree.get_skeleton(simplex_tree.dimension()):
        k = len(simplex)
        simplices[k-1][frozenset(simplex)] = len(simplices[k-1])
    return simplices

def build_laplacians(boundaries):
    """Build the Laplacian operators from the boundary operators.

    Parameters
    ----------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension.

    Returns
    -------
    laplacians: list of sparse matrices
       List of Laplacian operators, one per dimension: laplacian of degree i is in the i-th position
    """
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
    """Build the Laplacian operators from the boundary operators.

    Parameters
    ----------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension.

    Returns
    -------
    laplacians: list of sparse matrices
       List of Laplacian operators, one per dimension: laplacian of degree i is in the i-th position
    """
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
    """Build unweighted boundary operators from a list of simplices.

    Parameters
    ----------
    simplices: list of dictionaries
        List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-simplices. The dictionary's keys are sets (of size d
        + 1) of the 0-simplices that constitute the d-simplices. The
        dictionary's values are the indexes of the simplices in the boundary
        and Laplacian matrices.

    Returns
    -------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension: i-th boundary is in i-th position
    """
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
    """Build weighthd boundary operators from a list of simplices.

    Parameters
    ----------
    simplices: list of dictionaries
        List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-simplices. The dictionary's keys are sets (of size d
        + 1) of the 0-simplices that constitute the d-simplices. The
        dictionary's values are the indexes of the simplices in the boundary
        and Laplacian matrices.
    weights: list of sparse matrices
        List of sparse matrices for each dimension which diagonal contains the wights
    Returns
    -------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension: i-th boundary is in i-th position
    """
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
        
        Wn1=weights[d-1]
        w=weights[d].data[0]
        nz=np.nonzero(w)[0]
        inv=np.zeros(len(w))
        inv[nz]=(1/(w[nz]))
        inv_Wn=diags(inv) 
        boundary=Wn1@boundary@inv_Wn
        boundaries.append(boundary)
        
     
    return boundaries


def build_boundaries(simplices,weights=False):
    """Build weighted or unweighted boundary operators from a list of simplices.
        Default: unweighted
    Parameters
    ----------
    simplices: list of dictionaries
        List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-simplices. The dictionary's keys are sets (of size d
        + 1) of the 0-simplices that constitute the d-simplices. The
        dictionary's values are the indexes of the simplices in the boundary
        and Laplacian matrices.
    weights:None or list of sparse matrices, default None
        List of sparse matrices for each dimension which diagonal contains the wights
    Returns
    -------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension: i-th boundary is in i-th position
    """
    if weights==False:
        boundaries=build_unweighted_boundaries(simplices)
        
    else:
        boundaries=build_weighted_boundaries(simplices,weights=weights)
    
    return boundaries
    
def extract_xq(X,q,wq,dimq):
    """Collapse a pair Q W(Q) in a complex X

    Parameters
    ----------
    X: cell complex. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
    
    q: cell in the pair (q,wq) to be collapsed. q is a face of wq. List or string corresponding to the name of the cell.
    
    wq: cell in the pair (q,wq) to be collapsed. wq is a coface of q.
    
    dimq: integer dimension of the cell q
    
    Returns
    -------
    cells_xq: cells in XQ, the cell complex obtained after collapsing (q,wq).List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
    """
    
    q=frozenset(q)
    wq=frozenset(wq)
    x=copy.deepcopy(X)
    xq=[None]*len(x)
   
    for d in range(len(xq)):
  
        if d==dimq:
            
            del x[d][q]
            if len(x[d])>0:
                xq[d]=dict()
                for k,cell in enumerate(x[d].keys()):
                  
                    xq[d][cell]=k                
            
        if d==dimq +1:
            del x[d][wq]
            if len(x[d])>0:
                xq[d]=dict()
                for k,cell in enumerate(x[d].keys()):
                    xq[d][cell]=k
            
        if d!=dimq and d!=dimq+1:
            xq[d]=x[d]
           
    cells_xq=list(filter(None, xq))       
    return cells_xq

def extract_total_order(X):
    """Extract cells and total order of the cells

    Parameters
    ----------
    X: cell complex. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary,and chain homotopy matrices.
    -------
    tot_order: cells complex X stored as one dictionary. The dictionary's keys are sets corresponding to the name of the cells in X. The
        dictionary's values are the indexes of the cells in the chain homotopy matrices phi and psi.
    """    
    
    tot_order=dict()
    k=0
    for d in range(len(X)):
        for cell in X[d].keys():
            tot_order[cell]=k
            k=k+1   
            
    return tot_order


def build_Q_bounday(cellsXq,cellsX,q,wq,dimq,kX):
    """Build the boundary operators from a list of cellas and collapsing pair q,wq

    Parameters
    ----------
    cellsXq: cell complex given by collapsing (q,wq) in q. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
        
    cellsX: cell complex X. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
    
    q: cell in the pair (q,wq) to be collapsed. q is a face of wq. List or string corresponding to the name of the cell.
    
    wq: cell in the pair (q,wq) to be collapsed. wq is a coface of q.
    
    dimq: integer dimension of the cell q
    
    kX: List of boundary operators of X, one per dimension: i-th boundary is in i-th position. List of sparse matrices 

    Returns
    -------
    boundaries: List of boundary operators of XQ, one per dimension: i-th boundary is in i-th position. list of sparse matrices
    """
    
    q=frozenset(q)
    wq=frozenset(wq)
    
   
    
    Q_boundaries=list()
    

    for d in range(1, len(cellsXq)):
        
        idx_cells, idx_faces, values = [], [], []
        
        ##looking at the dimq-1 boundary
        if d==dimq and d!=0:
            for cell, idx_cell_q in cellsXq[d].items():
                    idx_cell_x=cellsX[d][cell]
                
                    kxT=kX[d-1].T.tolil()
                    non_zero_idx=kxT.rows[idx_cell_x]
                    
                    
                    for j in non_zero_idx:
                        value=kX[d-1].tocsr()[j,idx_cell_x]
                        values.append(value)
                        idx_cells.append(idx_cell_q)
                        idx_faces.append(j)
                        
            boundary=coo_matrix((values, (idx_faces, idx_cells)),dtype=np.float32,shape=(len(cellsXq[d-1]), len(cellsXq[d])))
            
            
                
        ### looking at the dimq boundary        
        if d==dimq+1:
            for cell, idx_cell_q in cellsXq[d].items():
                    idx_cell_x=cellsX[d][cell]
                
                    kxT=kX[d-1].T.tolil()
                    
                
                    non_zero_idx_cell=kxT.rows[idx_cell_x]
                    
                   
                    non_zero_idx_wq=kxT.rows[cellsX[d][wq]]
                    
                    
                    non_zero_idx=np.unique(np.append(non_zero_idx_cell,non_zero_idx_wq))
                    non_zero_idx=non_zero_idx[non_zero_idx!=cellsX[d-1][q]]
                    
                    for j in non_zero_idx:
                        
                        KXd=kX[d-1].tocsr()
                       
                        
                        value=KXd[j,idx_cell_x]-( KXd[cellsX[d-1][q],idx_cell_x]*KXd[j,cellsX[d][wq]])/(KXd[cellsX[d-1][q],cellsX[d][wq]])
    
                        
                        if value !=0:
                            values.append(value)
                            idx_cells.append(idx_cell_q)
                      
                            idx_faces.append( cellsXq[d-1][list(cellsX[d-1].keys())[j]] )
                    
                    boundary = coo_matrix((values, (idx_faces, idx_cells)),
                                     dtype=np.float32,
                                     shape=(len(cellsXq[d-1]), len(cellsXq[d])))
           
  

  ##looking at the dimq+1 boundary
        if d==dimq+2:
            for cell, idx_cell_q in cellsXq[d].items():
                    idx_cell_x=cellsX[d][cell]
                
                    kxT=kX[d-1].T.tolil()
                    
                    non_zero_idx=kxT.rows[idx_cell_x]
                    non_zero_idx=np.array(non_zero_idx)
                    
                    non_zero_idx=non_zero_idx[non_zero_idx!=cellsX[d-1][wq]]
                 
                    
                    
                    for j in non_zero_idx:
                        value=kX[d-1].tocsr()[j,idx_cell_x]
                        if value!=0:
                            values.append(value)
                            idx_cells.append(idx_cell_q)
                            idx_faces.append(cellsXq[d-1][list(cellsX[d-1].keys())[j]])
                        
            boundary=coo_matrix((values, (idx_faces, idx_cells)),dtype=np.float32,shape=(len(cellsXq[d-1]), len(cellsXq[d])))
            
            


            
        if d!=dimq and d!=dimq+1  and d!=dimq+2:
            
            boundary=kX[d-1]
                
                

        Q_boundaries.append(boundary)
    
    

    return Q_boundaries

def psi(X,Xq,q,wq,dimq,kX):
    """Build map psi from C(X) to C(XQ)

    Parameters
    ----------
    X: cell complex X. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
        
    Xq: cell complex given by collapsing (q,wq) in q. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
    
    q: cell in the pair (q,wq) to be collapsed. q is a face of wq. List or string corresponding to the name of the cell.
    
    wq: cell in the pair (q,wq) to be collapsed. wq is a coface of q.
    
    dimq: integer dimension of the cell q
    
    kX: List of boundary operators of X, one per dimension: i-th boundary is in i-th position. List of sparse matrices 

    Returns
    -------
    psi: sparse matrix of size X*XQ
    """   
    
    
    
    q=frozenset(q)
    wq=frozenset(wq)
    tot_X=extract_total_order(X)
    tot_Xq=extract_total_order(Xq)
    
    ind_x,ind_q,values=[],[],[]
   
    for cell in tot_X.keys():
        if cell== q:
            kxT=kX[dimq].T.tolil()
            idx_wq=X[dimq+1][wq]
            idx_q=X[dimq][q]
            non_zero_idx=kxT.rows[idx_wq]
            non_zero_idx=np.array(non_zero_idx)
            non_zero_idx=non_zero_idx[non_zero_idx!=X[dimq][q]]
            
            for idx_face in non_zero_idx:
                KXd=kX[dimq].tocsr()
                value= -(KXd[ idx_face,idx_wq ])/KXd[ idx_q,idx_wq ]
                values.append(value)
                ind_x.append(tot_X[cell])
                face= tot_Xq[list(X[dimq].keys())[idx_face]]
                ind_q.append(face)
                
        if cell!=q and cell!=wq:
            
            values.append(1)
            ind_x.append(tot_X[cell])
            ind_q.append(tot_Xq[cell])

    psi=coo_matrix((values, (ind_q, ind_x)),dtype=np.float32,shape=(len(tot_Xq), len(tot_X)))
            
            

            
    return psi

def phi(X,Xq,q,wq,dimq,kX):
    """Build map phi from C(XQ) to C(X)

    Parameters
    ----------
    X: cell complex X. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
        
    Xq: cell complex given by collapsing (q,wq) in q. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
            
    q: cell in the pair (q,wq) to be collapsed. q is a face of wq. List or string corresponding to the name of the cell.
    
    wq: cell in the pair (q,wq) to be collapsed. wq is a coface of q.
    
    dimq: integer dimension of the cell q
    
    kX: List of boundary operators of X, one per dimension: i-th boundary is in i-th position. List of sparse matrices 

    Returns
    -------
    psi: sparse matrix of size XQ*X
    """   
    
    q=frozenset(q)
    wq=frozenset(wq)
    tot_X=extract_total_order(X)
    tot_Xq=extract_total_order(Xq)
    
    ind_x,ind_q,values=[],[],[]
    
    
    kxT=kX[dimq].tolil()
    idx_wq=X[dimq+1][wq]
    idx_q=X[dimq][q]
    non_zero_idx=kxT.rows[idx_q]
    non_zero_idx=np.array(non_zero_idx)
    non_zero_idx=non_zero_idx[non_zero_idx!=idx_wq]

    cofaces_q=[list(X[dimq+1].keys())[j] for j in non_zero_idx]
   
    for cell,idx_cell in zip(cofaces_q,non_zero_idx):
        KXd=kX[dimq].tocsr()
        value= -(KXd[ idx_q,idx_cell])/KXd[ idx_q,idx_wq ]
        values.append(value)
        ind_x.append(tot_X[wq])
        ind_q.append(tot_Xq[cell])
        
        
    for cell in tot_Xq:
        values.append(1)
        ind_x.append(tot_X[cell])
        ind_q.append(tot_Xq[cell])
       
    
    phi=coo_matrix((values, (ind_x, ind_q)),dtype=np.float32,shape=(len(tot_X), len(tot_Xq)))

    return phi

def energy(X,Xq,q,wq,dimq,kX):
    """Return energy of collapsing (q,wq) in X
    """   
    
    s=phi(X,Xq,q,wq,dimq,kX)
    r=psi(X,Xq,q,wq,dimq,kX)
    
    return np.linalg.norm(np.identity(s.shape[0]) - s @ r)

def sequence_collpases(X,kX,collapses,dim_collapses):
    """Build map phi from C(XQ) to C(X)

    Parameters
    ----------
    X: cell complex X. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary.
    
    kX: List of boundary operators of X, one per dimension: i-th boundary is in i-th position. List of sparse matrices 

    collapses: sequences of collapses, list in which the d-entry is the dth-collapse and is given by a list of length 2 whith q and wq as entries
    Returns
    -------
    all_psi: list of sparse matrices, psi from C(XQj) to C(XQj+1)
    all_phi: list of sparse matrices, phi from C(XQj+1) to C(XQj)
    all_boundaries list of sparse matrices, boundaries of X,XQ1,XQ2..
    all_complexes list of sparse matrices, cell complexes of X,XQ1,XQ2..
    """   
    
    all_psi=[]
    all_phi=[]
    all_xq=[]
    all_boundaries=[]
    all_boundaries.append(kX)
    all_complexes=[]
    all_complexes.append(X)
    for j in range(len(collapses)):
        q=collapses[j][0]
        wq=collapses[j][1]
        x=all_complexes[-1]
        kx=all_boundaries[-1]
        xq=extract_xq(x,q,wq,dim_collapses[j])
        kq=build_Q_bounday(xq,x,q,wq,dim_collapses[j],kx)
        phi1=phi(x,xq,q,wq,dim_collapses[j],kx)
        psi1=psi(x,xq,q,wq,dim_collapses[j],kx)
        
        all_xq.append(xq)
        all_phi.append(phi1)
        all_psi.append(psi1)
        all_complexes.append(xq)
        all_boundaries.append(kq)
        
    return all_psi,all_phi, all_boundaries, all_complexes,all_xq
    

def energy_sequence(all_psi,all_phi):
    """Return ||1-psiphi|| and the matrix 1-psiphi
    """   
    psi1=all_psi[-1]
    for j in range(2,len(all_psi)+1):
        psi1=psi1 @ all_psi[-j]     
    
    phi1=all_phi[0]
    for j in range(1,len(all_phi)):
        phi1=phi1 @ all_phi[j]     
    
    energy= np.linalg.norm(np.identity(all_psi[0].shape[1])- phi1@psi1)
    matrix_energy=coo_matrix(np.identity(all_psi[0].shape[1])- phi1@psi1)
    return energy,matrix_energy


def phipsi(X,kX,collapses,dim_collapses,signal,type_collapse='up'):
    """Returns phipsi(s) of the signal
    """  
    c=collapses
    d=dim_collapses
    
  
    all_psi,all_phi, all_boundaries, all_complexes,all_xq=sequence_collpases(X,kX,collapses=c,dim_collapses=d)

 
    psi1=all_psi[-1]
    for j in range(2,len(all_psi)+1):
        psi1=psi1 @ all_psi[-j]     
    
    phi1=all_phi[0]
    for j in range(1,len(all_phi)):
        phi1=phi1 @ all_phi[j]     
    
    
    
    matrix= phi1@psi1
    ##The following lines can be taken away when doing collapses in different dimensions
                
    if type_collapse=='up':  
        a=0
        if dim_collapses[0]>0: 
            a=len(X[dim_collapses[0]-1])
        b=a+len(X[dim_collapses[0]])
        matrix1=matrix[a:b,a:b]
        phipsis= matrix1@np.array(signal).T
    
   
    if type_collapse=='down':        
        
        a=len(X[dim_collapses[0]])
        b=a+len(X[dim_collapses[0]+1])
        matrix1=matrix[a:b,a:b]
        
        phipsis= matrix1.T@np.array(signal).T
        
    return phipsis


def loss_signal(X,kX,collapses,dim_collapses,signal,type_collapse='up'):
    """Returns the topological error s-psiphi(s) of a signal after a sequence of collapsing  in X
    """  
    c=collapses
    d=dim_collapses
    
  
    all_psi,all_phi, all_boundaries, all_complexes,all_xq=sequence_collpases(X,kX,collapses=c,dim_collapses=d)

 
    psi1=all_psi[-1]
    for j in range(2,len(all_psi)+1):
        psi1=psi1 @ all_psi[-j]     
    
    phi1=all_phi[0]
    for j in range(1,len(all_phi)):
        phi1=phi1 @ all_phi[j]     
    
    if type_collapse=="up":
        matrix_energy=np.identity(all_psi[0].shape[1])- phi1@psi1
        ##The following lines can be taken away when doing collapses in different dimensions
        a=0
        if dim_collapses[0]>0: 
            a=len(X[dim_collapses[0]-1])
        b=a+len(X[dim_collapses[0]])
        matrix_energy1=matrix_energy[a:b,a:b]
    
    if type_collapse=='down':
       
        matrix_energy=np.identity(all_psi[0].shape[1])- (phi1@psi1).T         
        a=len(X[dim_collapses[0]])
        b=a+len(X[dim_collapses[0]+1])
     
        matrix_energy1=matrix_energy[a:b,a:b]
        
    
    loss= np.linalg.norm(matrix_energy1@np.array(signal).T)**2
    return loss

def best_up_collapse(X,kX,dimq,signal): 
    """Returns the up-collapse V which minimize the topological reconstructione error ||s-psi_Vphi_V(s)|^2, as described in Algorithm 1
    """  
    s=np.array(signal)**2
    Bq=kX[dimq]
    BT=Bq.transpose()
    C=BT.tolil()
    ind=C.rows

    possible_min_q=[]
    possible_min_value=[]
    for wq in range(Bq.shape[1]):
        
        boundary=ind[wq]
        sig_wq=s[boundary]
        d=np.array(C.data[wq])**2
        sig_wq=sig_wq/d
        notnan_sig=sig_wq[~np.isnan(sig_wq)]
        min_sq=np.min(notnan_sig)
        norm_dwq=np.linalg.norm(Bq.tocsr()[ind[wq],wq].todense())**2
        possible_min_value.append(min_sq*norm_dwq)
        
        
        minimum=random.choice(np.where(notnan_sig==min_sq)[0])
        
        idx_q_min=ind[wq][minimum]
        possible_min_q.append(idx_q_min)
    global_min=np.min(possible_min_value)
  
    Wq=random.choice(np.where(possible_min_value==global_min)[0])
    q=possible_min_q[Wq]
    
    return [[list(list(X[dimq].keys())[q]),list(list(X[dimq+1].keys())[Wq])]],global_min

def best_down_collapse(X,kX,dimq,signal): 
    """Returns the down-collapse V which minimize the topological reconstructione error ||s-phi*_Vpsi*_V(s)|^2 
    """  
    s=np.array(signal)**2
    Bq=kX[dimq]
    BT=Bq.copy()
    C=BT.tolil()
    ind=C.rows

    possible_min_WQ=[]
    possible_min_value=[]
  
    
    for Q in range(Bq.shape[0]):

        boundary=ind[Q]

        sig_WQ=s[boundary]
        d=np.array(C.data[Q])**2
        sig_WQ=sig_WQ/d
        notnan_sig=sig_WQ[~np.isnan(sig_WQ)]
        min_sWQ=np.min(notnan_sig)

        norm_dQ=np.linalg.norm(Bq.tocsr()[Q,ind[Q]].todense())**2
        possible_min_value.append(min_sWQ*norm_dQ)


        minimum=random.choice(np.where(notnan_sig==min_sWQ)[0])

        idx_WQ_min=ind[Q][minimum]
        possible_min_WQ.append(idx_WQ_min)
    global_min=np.min(possible_min_value)  
    q=random.choice(np.where(possible_min_value==global_min)[0])
    Wq=possible_min_WQ[q]

    return [[list(list(X[dimq].keys())[q]),list(list(X[dimq+1].keys())[Wq])]],global_min


def random_collapse(X,kX,dimq):
    """Returns a random collapses V =(q,mu(q)) with fixed dim(q)
    """ 
    nz=np.array(kX[dimq].nonzero())
    rand=np.random.choice(np.arange(len(nz[0])))
    Wq=nz[:,rand][1]
    q=nz[:,rand][0]
    return([[list(list(X[dimq].keys())[q]),list(list(X[dimq+1].keys())[Wq])]])

    
def sequence_optimal_up_collapses(X,kX,dimq,signal,steps,random=False):
    """Returns a sequence of up-collapse V by iterating single optimal up-collapses for a fixed number of steps. This imoplements Algorithm 2 of the paper
    """ 
    dX=kX.copy()
    all_X=[X]
    all_signals=[signal]
    all_collapses=[]
    all_losses=[]
    
    for k in range(steps):
        
        if random==False:
            c=best_up_collapse(X,kX,dimq=dimq,signal=signal)[0]
        if random==True:
            c=random_collapse(X,kX,dimq)
            
        all_psi,all_phi, all_boundaries, all_complexes,all_xq=sequence_collpases(X,kX,collapses=c,dim_collapses=[dimq])
        
        a=len(X[dimq-1])
        b=a+len(X[dimq])
        L=energy_sequence(all_psi,all_phi)[1].toarray()[a:b,a:b]
        loss=loss_signal(X,kX,collapses=c,dim_collapses=[dimq],signal=signal)
        
        psi1=all_psi[-1]
        for j in range(2,len(all_psi)+1):
                psi1=psi1 @ all_psi[-j]
                
        
        signal=psi1.toarray()[a:b-1,a:b]@np.array(signal).T
        
        
        kX=all_boundaries[-1]
        
        X=all_xq[-1]
        
        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)
        
    dimc=[dimq]*len(all_collapses)
    total_loss=loss_signal(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0])
    phispsis=phipsi(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0])

    return(all_X,all_collapses,all_losses,total_loss,all_signals,phispsis)
    
    
    
def sequence_optimal_down_collapses(X,kX,dimq,signal,steps,random=False):
    """Returns a sequence of down-collapse V by iterating single optimal down-collapses for a fixed number of steps.
    """
    dX=kX.copy()
    all_X=[X]
    all_signals=[signal]
    all_collapses=[]
    all_losses=[]
    
    for k in range(steps):
        
        if random==False:
            c=best_down_collapse(X,kX,dimq=dimq,signal=signal)[0]
        if random==True:
            c=random_collapse(X,kX,dimq)
        all_psi,all_phi, all_boundaries, all_complexes,all_xq=sequence_collpases(X,kX,collapses=c,dim_collapses=[dimq])
        
        a=len(X[dimq])
        b=a+len(X[dimq+1])
        
        loss=loss_signal(X,kX,collapses=c,dim_collapses=[dimq],signal=signal,type_collapse='down')
        
        phi1=all_phi[-1]
        for j in range(2,len(all_phi)+1):
                phi1=phi1 @ all_phi[-j]
                
        phiT=phi1.toarray().T
        signal=phiT[a+1:b,a:b]@np.array(signal).T
        
        
        
        kX=all_boundaries[-1]
        
        X=all_xq[-1]
        
        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)
        
    dimc=[dimq]*len(all_collapses)
    total_loss=loss_signal(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0],type_collapse='down')
    phispsis=phipsi(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0],type_collapse='down')

    return(all_X,all_collapses,all_losses,total_loss,all_signals,phispsis)
    

def sequence_given_up_collapses(X,kX,dimq,signal,collapses):
    """Returns the collapses complexes, the reconstruction error and the recontructed signal for a given up-matching
    """
    
    dX=kX.copy()
    all_X=[X]
    all_signals=[signal]
    all_collapses=[]
    all_losses=[]
    
    
    for collapse in collapses:
        c=[collapse]
        all_psi,all_phi, all_boundaries, all_complexes,all_xq=sequence_collpases(X,kX,collapses=c,dim_collapses=[dimq])
        
        a=len(X[dimq-1])
        b=a+len(X[dimq])
        L=energy_sequence(all_psi,all_phi)[1].toarray()[a:b,a:b]
        loss=loss_signal(X,kX,collapses=c,dim_collapses=[dimq],signal=signal)
        
        psi1=all_psi[-1]
        for j in range(2,len(all_psi)+1):
                psi1=psi1 @ all_psi[-j]
                
       
        signal=psi1.toarray()[a:b-1,a:b]@np.array(signal).T
        
        
        kX=all_boundaries[-1]
        X=all_xq[0]
        
        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)
        
    dimc=[dimq]*len(all_collapses)
    
    total_loss=loss_signal(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0])
    phispsis=phipsi(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0])
    return(all_X,all_collapses,all_losses,total_loss,all_signals,phispsis)


def sequence_given_down_collapses(X,kX,dimq,signal,collapses):
    """Returns the collapses complexes, the reconstruction error and the recontructed signal for a given down-matching
    """
    all_X=[X]
    dX=kX.copy()
    all_signals=[signal]
    all_collapses=[]
    all_losses=[]
    
    
    for collapse in collapses:
        c=[collapse]
        all_psi,all_phi, all_boundaries, all_complexes,all_xq=sequence_collpases(X,kX,collapses=c,dim_collapses=[dimq])
        
        a=len(X[dimq])
        b=a+len(X[dimq+1])
        
        loss=loss_signal(X,kX,collapses=c,dim_collapses=[dimq],signal=signal,type_collapse='down')
        
        phi1=all_phi[-1]
        for j in range(2,len(all_phi)+1):
                phi1=phi1 @ all_phi[-j]
                
        phiT=phi1.toarray().T
        signal=phiT[a-1:b,a:b]@np.array(signal).T
        
        
        
        kX=all_boundaries[-1]
        
        X=all_xq[-1]
        
        all_signals.append(signal)
        all_X.append(X)
        all_collapses.append(c[0])
        all_losses.append(loss)
        
    dimc=[dimq]*len(all_collapses)
    total_loss=loss_signal(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0],type_collapse='down')
    phispsis=phipsi(all_X[0],dX,collapses=all_collapses,dim_collapses=dimc,signal=all_signals[0],type_collapse='down')
    return(all_X,all_collapses,all_losses,total_loss,all_signals,phispsis)



def simulation_collapses(X,kX,dimq,signal,steps,random):
    """Computes multiple instantiations of optimal or random collapses
    """
    total_losses=[]
    for k in steps:
        all_X,collapses,all_losses,total_loss,all_signals,phispsis= sequence_optimal_up_collapses(X=X,kX=kX,dimq=dimq,signal=signal,steps=k,random=random)
        total_losses.append(total_loss)    
    return(total_losses)