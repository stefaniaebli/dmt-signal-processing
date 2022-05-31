from distutils.command import build
from itertools import chain, combinations
from attr import frozen
import gudhi as gd
from lib.build_boundaries import extract_simplices, build_boundaries


class ChainComplex:
    def __init__(self, list_of_simplices=[], kX=None, X=None):
        if kX and X != None:
            self.kX = kX
            self.X = X
        else:
            simplex_tree = gd.SimplexTree()
            for element in list_of_simplices:
                simplex_tree.insert(element)
            self.X = extract_simplices(simplex_tree)
            self.kX = build_boundaries(self.X)

    def add_simplex(self, simplex):
        simplex = frozenset(simplex)
        if simplex not in self.X[len(simplex)-1].keys():
            self.X[len(simplex)-1][frozenset(simplex)
                                   ] = len(self[len(simplex)-1])
            for element in self.__subsimplices(simplex):
                if element not in self.X[len(element)-1].keys():
                    self.X[len(element)-1][frozenset(element)
                                           ] = len(self[len(element)-1])
        self.kX = build_boundaries(self.X)

    def delete_simplex(self, simplex):
        simplex = frozenset(simplex)
        try:
            del self.X[len(simplex)-1][simplex]

            for element in ChainComplex.__subsimplices(simplex):
                try:
                    del self.X[len(element)-1][element]
                except KeyError:
                    pass
            for id, key in enumerate(self.X[len(simplex)-1]):
                self.X[len(simplex)-1][key] = id

            for iterator in range(len(simplex), len(self.X)):
                for key in self.X[iterator].keys():
                    for element in simplex:
                        if element in key:
                            del self.X[iterator][key]
                for id, key in enumerate(self.X[iterator]):
                    self.X[iterator][key] = id

            self.kX = build_boundaries(self.X)

        except KeyError:
            print(f'{simplex} does not exist')

    def __subsimplices(simplex):
        s = list(simplex)
        return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))
