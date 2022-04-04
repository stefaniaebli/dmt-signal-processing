# Library documentation

### This is a short description of every function's action and parameters. Read if any explanation needed.

## `build_boundaries.py`

* `extract_simplices`
    
    Extract simplices from a gudhi simplex tree.

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

* `build_laplacians`
    
    Build the Laplacian operators from the boundary operators.

    Parameters
    ----------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension.

    Returns
    -------
    laplacians: list of sparse matrices
       List of Laplacian operators, one per dimension: laplacian of degree i is in the i-th position
    
* `build_up_down_laplacians`

    Build the Laplacian operators from the boundary operators.

    Parameters
    ----------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension.

    Returns
    -------
    laplacians: list of sparse matrices
       List of Laplacian operators, one per dimension: laplacian of degree i is in the i-th position
    
* `build_unweighted_boundaries`

    Build unweighted boundary operators from a list of simplices.

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
    
* `build_weighted_boundaries`

    Build weighthd boundary operators from a list of simplices.

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

* `build_boundaries`

    Build weighted or unweighted boundary operators from a list of simplices.
        Default: unweighted

    Parameters
    ----------
    simplices: list of dictionaries
        List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-simplices. The dictionary's keys are sets (of size d
        + 1) of the 0-simplices that constitute the d-simplices. The
        dictionary's values are the indexes of the simplices in the boundary
        and Laplacian matrices.

    weights:
        List of sparse matrices for each dimension which diagonal contains the wights or None, default None.

    Returns
    -------
    boundaries: list of sparse matrices
       List of boundary operators, one per dimension: i-th boundary is in i-th position


## `signal_processing.py`

* `extract_xq`

    Collapse a pair Q W(Q) in a complex X

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
    
* `extract_total_order`

    Extract cells and total order of the cells

    Parameters
    ----------
    X: cell complex. List of dictionaries, one per dimension d. The size of the dictionary
        is the number of d-cells. The dictionary's keys are sets corresponding to the name of the cells in XQ. The
        dictionary's values are the indexes of the cells in the boundary,and chain homotopy matrices.

    Returns
    -------
    tot_order: cells complex X stored as one dictionary. The dictionary's keys are sets corresponding to the name of the cells in X. The
        dictionary's values are the indexes of the cells in the chain homotopy matrices phi and psi.
    
* `build_Q_boundary`

    Build the boundary operators from a list of cellas and collapsing pair q,wq

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
    
* `psi`

    Build map psi from C(X) to C(XQ)

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

* `phi`

    Build map phi from C(XQ) to C(X)

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
    
* `energy`

    Return energy of collapsing (q,wq) in X

* `sequence_collapses`

    Build map phi from C(XQ) to C(X)

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
    
* `energy_sequence`

    Return ||1-psiphi|| and the matrix 1-psiphi

* `phipsi`

    Returns phipsi(s) of the signal

* `loss_signal`

    Returns the topological error s-psiphi(s) of a signal after a sequence of collapsing  in X


## `optimal_collapses.py`

* `best_up_collapse`

    Returns the up-collapse V which minimize the topological reconstructione error ||s-psi_Vphi_V(s)|^2, as described in Algorithm 1

* `best_down_collapse`

    Returns the down-collapse V which minimize the topological reconstructione error ||s-phi*_Vpsi*_V(s)||^2

* `random_collapse`

    Returns a random collapses V =(q,mu(q)) with fixed dim(q)
    
* `sequence_optimal_up_collapses`

    Returns a sequence of up-collapse V by iterating single optimal up-collapses for a fixed number of steps. This imoplements Algorithm 2 of the paper
    
* `sequence_optimal_down_collapses`

    Returns a sequence of down-collapse V by iterating single optimal down-collapses for a fixed number of steps.

* `sequence_given_up_collapses`

    Returns the collapses complexes, the reconstruction error and the recontructed signal for a given up-matching

* `sequence_given_down_collapses`

    Returns the collapses complexes, the reconstruction error and the recontructed signal for a given down-matching

* `simulation_collapses`

    Computes multiple instantiations of optimal or random collapses

* `check_hodge_decomp`

    None.
    