#symplicial complex object file

class ChainComplex():
    def __init__(self) -> None:
        self.__init__ = self

        return self


    def add_simplex(self, simplex):
        if simplex not in self[len(simplex)-1].keys():
            self[len(simplex)-1][frozenset(simplex)] = len(self[len(simplex)-1])
    
        return self


    def delete_simplex(self, simplex):
        try:
            del self[len(simplex)-1][simplex]
            for element in self.__subsimplices(simplex):
                try:
                    del self[len(element)-1][element]
                except KeyError:
                    pass
        except KeyError:
             print(f'{simplex} does not exist', simplex)
        

    def __subsimplices(self, simplex):
        subsimplices = []
        for i in range(1, len(simplex)+1):
            #generate all subsimplices via popping elements
