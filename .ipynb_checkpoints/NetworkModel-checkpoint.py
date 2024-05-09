from Instance import *
from pyvis.network import Network
import networkx as nx
import numpy as np
import time
import random
from matplotlib.animation import FuncAnimation

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

        self.loadedTimeSeries = [[], [], []]
        
        self.frameScalingFactor = 2000
        self.figure, self.axis = plt.subplots()


    def initializeNetworkXRepresentation(self):
        networkGraph = nx.from_numpy_array(self.hostMutationAdjacencyMatrix)
        nodePositions = nx.spring_layout(networkGraph)
        return [networkGraph, nodePositions]
        
    
    def createInstance(self, initialHostPopulationValues, initialPhagePopulationValues, carryingCapacity, lysisQuantities, numIters):
        self.instances.append(Instance(initialHostPopulationValues, initialPhagePopulationValues, carryingCapacity, self.fitnessParameters,
                                       lysisQuantities, self.hostMutationAdjacencyMatrix, self.infectionMatrix, self.viralMutationAdjacencyMatrix,
                                       numIters))

    def runInstance(self, instanceIndex):
        self.instances[instanceIndex].runAlgorithm()
    
    def getNetworkParameter(self, networkParameter):
        return self.networkParameters[networkParameter]

    '''def displayModel(self):
        colorMap = [(0.5, 0.5, 0), (0.5, 0.5, 0), (0.5, 0.5, 0)]
        nodes = nx.draw_networkx_nodes(self.nxRepresentation[0], pos = self.nxRepresentation[1], node_color = colorMap)
        edges = nx.draw_networkx_edges(self.nxRepresentation[0], pos = self.nxRepresentation[1])
        
        return nodes,'''
        
    def displayFrame(self, frame):
        populationVals = self.loadedTimeSeries[1][:, self.frameScalingFactor * frame]
        colorMap = self.getColorsFromPops(populationVals)
        self.axis.clear()
        nx.draw(self.nxRepresentation[0], self.nxRepresentation[1], with_labels = True, node_color = colorMap, ax = self.axis)

    def getColorsFromPops(self, populationValues):
        carryingCapacity = self.loadedTimeSeries[2]
        colorMap = []
        for j in range(len(populationValues)):
            colorMap.append((0.5, 0.5, populationValues[j] / carryingCapacity))
            #colorMap.append((random.random(), random.random(), random.random()))
        #print("The Color Map is presently {}".format(colorMap))
        return colorMap
            
    def getInstance(self, instanceIndex):
        return self.instances[instanceIndex]

    def saveModel(self):
        return 0

    def loadModel(self):
        return 0

    def loadInstance(self, instanceIndex):
        self.loadedTimeSeries[0] = self.getInstance(instanceIndex).timeValues
        self.loadedTimeSeries[1] = self.getInstance(instanceIndex).hostPopulationValues
        self.loadedTimeSeries[2] = self.getInstance(instanceIndex).carryingCapacity

    def displayInstanceTimecourses(self, instanceIndex):
        self.instances[instanceIndex].displayResultsAsTimecourses()
        plt.show()

#Test suite
hostMutationMatrix = np.zeros((3, 3))
hostMutationMatrix[0, 0] = 0.98
hostMutationMatrix[0, 1] = 0.01
hostMutationMatrix[0, 2] = 0.01
hostMutationMatrix[1, 0] = 0.01
hostMutationMatrix[1, 1] = 0.98
hostMutationMatrix[1, 2] = 0.01
hostMutationMatrix[2, 0] = 0.01
hostMutationMatrix[2, 1] = 0.01
hostMutationMatrix[2, 2] = 0.98

phageHostInfectionMatrix = np.zeros((3, 3))
phageHostInfectionMatrix[0, 0] = 1
phageHostInfectionMatrix[0, 1] = 0
phageHostInfectionMatrix[0, 2] = 0
phageHostInfectionMatrix[1, 0] = 0
phageHostInfectionMatrix[1, 1] = 1
phageHostInfectionMatrix[1, 2] = 0
phageHostInfectionMatrix[2, 0] = 0
phageHostInfectionMatrix[2, 1] = 0
phageHostInfectionMatrix[2, 2] = 1

phageMutationMatrix = np.zeros((3, 3))
phageMutationMatrix[0, 0] = 0.998
phageMutationMatrix[0, 1] = 0.001
phageMutationMatrix[0, 2] = 0.001
phageMutationMatrix[1, 0] = 0.001
phageMutationMatrix[1, 1] = 0.998
phageMutationMatrix[1, 2] = 0.001
phageMutationMatrix[2, 0] = 0.001
phageMutationMatrix[2, 1] = 0.001
phageMutationMatrix[2, 2] = 0.998

burstSize = 2
lysisLength = 10
lysisQuantities = [burstSize, lysisLength]

carryingCapacity = 100000
numIters = 9600000

netTest = NetworkModel(hostMutationMatrix, [0.5, 0.5, 0.5], phageMutationMatrix, phageHostInfectionMatrix)

print("\nTesting createInstance()\n")
netTest.createInstance([1000, 0, 0], [1, 0, 0], carryingCapacity, lysisQuantities, numIters)

for inst in netTest.instances:
    inst.allPopsReportForDuty()

print("\nTesting runInstance()\n")
netTest.runInstance(0)

for inst in netTest.instances:
    inst.allPopsReportForDuty()

print("\nTesting initializeNetworkXRepresentation\n")
netTest.initializeNetworkXRepresentation()

print("\nTesting displayModel()\n")
#netTest.displayModel()

print("\nTesting loadInstance()\n")
netTest.loadInstance(0)

print(netTest.loadedTimeSeries)

print("\nTesting displayIteration()\n")
#netTest.displayIteration(10000)

#netTest.displayInstanceTimecourses(0)

fig, ax = plt.subplots()

G = netTest.nxRepresentation[0]
pos = nx.spring_layout(G)

def displayUnobjectFrame(frame):  
    populationVals = netTest.loadedTimeSeries[1][:, netTest.frameScalingFactor * frame]
    #print("Population Values in use are {}".format(populationVals))
    colorMap = netTest.getColorsFromPops(populationVals)
    #print("Color Map in use is {}".format(colorMap))
    ax.clear()
    nx.draw(G, pos, with_labels = True, node_color = colorMap, ax = ax)
    

print("\n Testing animation methodology\n")
rate = 50
numIterations = len(netTest.loadedTimeSeries[0])
numFrames = int(numIters / netTest.frameScalingFactor)
ani = FuncAnimation(fig, displayUnobjectFrame, frames = numFrames, interval = rate, blit=False)

#ani.save('CoevolutionwithUnidirectionalMutation.mp4', fps=30)

plt.show()