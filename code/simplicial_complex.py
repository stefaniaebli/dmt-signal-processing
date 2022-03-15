from itertools import chain, combinations

from attr import frozen

class ChainComplex(simplicial_complex,cell_complex=False, chain_complex=False, weights=False):
    def __init__(self) -> None:
        self.__init__ = self

        if cell complex==True:
            self._ordered_cells=cell_complex
            self._chain_complex=chain_complex
        elif:
            self._ordered_cells=extract_simplices(simplicial_complex)
            self._chain_complex= build_boundaries(self,weights)


    def extract_simplices(simplicial_complex):
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
            simplex_tree=gd.SimplexTree()
            for s in simplicial_complex
                simplex_tree.insert(s)

            simplices = [dict() for _ in range(simplex_tree.dimension()+1)]
            #simplices = [dict() for _ in range(simplex_tree.dimension()+1)]
            for simplex, _ in simplex_tree.get_skeleton(simplex_tree.dimension()):
                k = len(simplex)
                simplices[k-1][frozenset(simplex)] = len(simplices[k-1])
            return simplices

    def build_unweighted_boundaries(self):
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
        for d in range(1, len(self._ordered_cells)):
            idx_simplices, idx_faces, values = [], [], []
            for simplex, idx_simplex in self._ordered_cells[d].items():
                for i, left_out in enumerate(np.sort(list(simplex))):
                    idx_simplices.append(idx_simplex)
                    values.append((-1)**i)
                    face = simplex.difference({left_out})
                    idx_faces.append(self._ordered_cells[d-1][face])
            assert len(values) == (d+1) * len(self._ordered_cells[d])
            boundary = coo_matrix((values, (idx_faces, idx_simplices)),
                                         dtype=np.float32,
                                         shape=(len(self._ordered_cells[d-1]), len(self._ordered_cells[d])))
            boundaries.append(boundary)


        return boundaries


    def build_weighted_boundaries(self,weights):
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
            for simplex, idx_simplex in self._ordered_cells[d].items():
                for i, left_out in enumerate(np.sort(list(simplex))):
                    idx_simplices.append(idx_simplex)
                    values.append((-1)**i)
                    face = simplex.difference({left_out})
                    idx_faces.append(self._ordered_cells[d-1][face])
            assert len(values) == (d+1) * len(simplices[d])
            boundary = coo_matrix((values, (idx_faces, idx_simplices)),
                                         dtype=np.float32,
                                         shape=(len(self._ordered_cells[d-1]), len(self._ordered_cells[d])))

            Wn=weights[d]

            w=weights[d-1].data[0]
            nz=np.nonzero(w)[0]
            inv=np.zeros(len(w))
            inv[nz]=(1/(w[nz]))
            inv_Wn=diags(inv)

            boundary=inv_Wn1@boundary@Wn
            boundaries.append(boundary)


        return boundaries


def build_boundaries(self,weights):
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
        boundaries=build_unweighted_boundaries(self)

    else:
        boundaries=build_weighted_boundaries(self,weights=weights)

    return boundaries



    def add_simplex(self, simplex):
        if simplex not in self[len(simplex)-1].keys():
            self[len(simplex)-1][frozenset(simplex)] = len(self[len(simplex)-1])
            for element in self.__subsimplices(simplex):
                if element not in self[len(element)-1].keys():
                    self[len(element)-1][frozenset(element)] = len(self[len(element)-1])


    def delete_simplex(self, simplex):
        try:
            del self[len(simplex)-1][simplex]
            for element in self.__subsimplices(simplex):
                try:
                    del self[len(element)-1][element]
                except KeyError:
                    pass
            #@TODO value reenumeration
        except KeyError:
             print(f'{simplex} does not exist', simplex)


    def __subsimplices(simplex):
        s = list(simplex)
        return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))
