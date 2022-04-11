from itertools import chain, combinations
from dataclasses import dataclass
import gudhi as gd
from lib.build_boundaries import extract_simplices, build_boundaries



class ChainComplex:
    def __init__(self, list_of_simplices):
        simplex_tree = gd.SimplexTree()
        for element in list_of_simplices:
            simplex_tree.insert(element)
        self.X = extract_simplices(simplex_tree)
        self.kX = build_boundaries(self.X)


    def add_simplex(self, simplex):
        if simplex not in self.X[len(simplex)-1].keys():
            self.X[len(simplex)-1][frozenset(simplex)] = len(self[len(simplex)-1])
            for element in self.__subsimplices(simplex):
                if element not in self.X[len(element)-1].keys():
                    self.X[len(element)-1][frozenset(element)] = len(self[len(element)-1])


    def delete_simplex(self, simplex):
        try:
            del self.X[len(simplex)-1][simplex]
            for element in self.__subsimplices(simplex):
                try:
                    del self.X[len(element)-1][element]
                except KeyError:
                    pass
            for id, key in enumerate(self[len(simplex)-1]):
                self.X[len(simplex)-1][key] = id
        except KeyError:
             print(f'{simplex} does not exist', simplex)


    def __subsimplices(simplex):
        s = list(simplex)
        return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))
