from Instance import *
from pyvis.network import Network
import networkx as nx

class NetworkModel():
    def __init__(self, hostMutMatrix, fitnessParams, viralMutMatrix, infectMatrix):
        self.hostMutationAdjacencyMatrix = hostMutMatrix
        self.viralMutationAdjacencyMatrix = viralMutMatrix
        self.infectionMatrix = infectMatrix
        self.fitnessParameters = fitnessParams
        self.instances = []

        self.networkParameters = {"Host Mutation Matrix" : self.hostMutationAdjacencyMatrix,
                             "Viral Mutation Matrix" : self.viralMutationAdjacencyMatrix,
                             "Infection Matrix" : self.infectionMatrix}

        self.nxRepresentation = self.initializeNetworkXRepresentation()


    def initializeNetworkXRepresentation(self):
        return nx.from_numpy_array(self.hostMutationAdjacencyMatrix)
    
    def createInstance(self, initialPopulationValues, carryingCapacity, numIters):
        self.instances.append(Instance(initialPopulationValues, carryingCapacity, self.fitnessParameters, self.hostMutationAdjacencyMatrix, numIters))

    def runInstance(self, instanceIndex):
        self.instances[instanceIndex].runAlgorithm()
    
    def getNetworkParameter(self, networkParameter):
        return self.networkParameters[networkParameter]

    def displayModel(self):
        nt = Network("500px", "500px")
        nt.from_nx(self.nxRepresentation)
        nt.show("nx.html")

    def getInstance(self):
        return 0

    def saveModel(self):
        return 0

    def loadModel(self):
        return 0

    def loadInstance(self):
        return 0
        
#Test suite
mutationMatrix = np.zeros((2, 2))
mutationMatrix[0, 0] = 1
mutationMatrix[0, 1] = 0
mutationMatrix[1, 0] = 0.25
mutationMatrix[1, 1] = 0.75

netTest = NetworkModel(mutationMatrix, [10, 20], [], [])

print("\nTesting createInstance()\n")
netTest.createInstance([50, 50], 10000, 10000)

for inst in netTest.instances:
    inst.allPopsReportForDuty()

print("\nTesting runInstance()\n")
netTest.runInstance(0)

for inst in netTest.instances:
    inst.allPopsReportForDuty()

print("\nTesting initializeNetworkXRepresentation\n")
netTest.initializeNetworkXRepresentation()

print("\nTesting displayModel()\n")
netTest.displayModel()