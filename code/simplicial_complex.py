from itertools import chain, combinations

class ChainComplex():
    def __init__(self) -> None:
        self.__init__ = self

        return self


    def add_simplex(self, simplex):
        if simplex not in self[len(simplex)-1].keys():
            self[len(simplex)-1][frozenset(simplex)] = len(self[len(simplex)-1])
            #subsimplices


    def delete_simplex(self, simplex):
        try:
            del self[len(simplex)-1][simplex]
            for element in self.__subsimplices(simplex):
                try:
                    del self[len(element)-1][element]
                except KeyError:
                    pass
            #reenumerate the values
        except KeyError:
             print(f'{simplex} does not exist', simplex)
        

    def __subsimplices(simplex):
        s = list(simplex)
        return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))
