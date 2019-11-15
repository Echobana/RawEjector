#! python3

class RateGraph():
    def __init__(self, rates):
        self.graph = dict()
        for orig, dest, rate in rates:
            self.add_conversion(orig, dest, rate)

    def add_conversion(self, orig, dest, rate):
        if orig not in self.graph:
            self.graph[orig] = {}
        self.graph[orig][dest] = rate

    def get_neighbors(self, node):
        if node not in self.graph:
            return None
        return self.graph[node].items()

    def get_nodes(self):
        return self.graph.keys()

    def get_conversion(self, orig, dest):
        return self.graph[orig][dest]


if __name__ == "__main__":
    rates = [['atma', 'atmg', (1, 1)],
             ['atma', 'MPa', (0.1, 0)],
             ['atmg', 'MPa', (0.1, 0.1)]]
    graph = RateGraph(rates)
    print(graph.get_conversion('atma', 'atmg'))
