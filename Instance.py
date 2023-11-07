from Populations import *
import matplotlib.pyplot as plt

class Instance():
    def __init__(self, initialPopulationValues, carryingCap, fitnessParameters, hostMutMatrix, numIters):
        self.hostPopulations = []
        self.carryingCapacity = carryingCap
        self.numIterations = numIters
        
        self.initializePopulations(initialPopulationValues, fitnessParameters)

        self.timeValues = np.zeros(numIters + 1)
        self.hostPopulationValues = np.zeros((len(self.hostPopulations), numIters + 1))
        self.hostPopulationValues[:, 0] = [self.hostPopulations[k].populationLevel for k in range(len(self.hostPopulations))]

        self.hostMutationMatrix = hostMutMatrix
        self.completedRun = False

    def initializePopulations(self, initialPopVals, fitnessParams):
        for j in range(len(fitnessParams)):
            self.hostPopulations.append(hostPopulation(j, initialPopVals[j], fitnessParams[j]))

        hostPopVals = [pop.populationLevel for pop in self.hostPopulations]
        for pop in self.hostPopulations:
            pop.computeGrowthRate(hostPopVals, self.carryingCapacity)
            pop.regenerateTimer()

    def runAlgorithm(self):
        for k in range(self.numIterations):
            #Step 1
            (minTimer, virocell, minPopIdx) = self.findMinPopTimer()
    
            #Step 2
            self.timeValues[k + 1] = self.timeValues[k] + minTimer
    
            #Step 3
            self.updateHostPopulation(minPopIdx)
            self.hostPopulationValues[:, k + 1] = [self.hostPopulations[k].populationLevel for k in range(len(self.hostPopulations))]
    
            #Step 4
            self.recomputeGrowthRates()
    
            #Step 5
            self.updateHostTimers(minTimer)
            #self.updatePhageTimers()
    
            #Step 6
            self.regeneratePopTimers()

        self.completedRun = True
    
    def findMinPopTimer(self):
        hostPopTimers = [pop.getTimer() for pop in self.hostPopulations]
        minHostTimer = min(hostPopTimers)
        
        minHostPopTimerIdx = hostPopTimers.index(minHostTimer)

        virocell = 0
        
        return (minHostTimer, virocell, minHostPopTimerIdx)

    def updateHostPopulation(self, minPopIdx):
        #Get the growth rate's sign to see growth or death
        growth_rate = self.hostPopulations[minPopIdx].growthRate
        
        if growth_rate > 0:
            localMutationProbabilities = self.hostMutationMatrix[minPopIdx, :]
            possibilities = [n for n in range(len(self.hostPopulations))]
            #random.choice(what you're choosing from, how many you're choosing, are you replacing, weights)
            spawnChoiceIdx = np.random.choice(possibilities, None, True, localMutationProbabilities)
            self.hostPopulations[spawnChoiceIdx].changeHostPopulation(1)
        elif growth_rate < 0:
            self.hostPopulations[minPopIdx].changeHostPopulation(-1)
        else:
            self.hostPopulations[minPopIdx].changeHostPopulation(0)

        
    def recomputeGrowthRates(self):
        hostPopVals = [pop.populationLevel for pop in self.hostPopulations]
        for pop in self.hostPopulations:
            pop.computeGrowthRate(hostPopVals, self.carryingCapacity)

    def updateHostTimers(self, minPopTimer):
        for pop in self.hostPopulations:
            pop.updateTimer(minPopTimer)

    def updatePhageTimers(self):
        return 0

    def regeneratePopTimers(self):
        for pop in self.hostPopulations:
            pop.regenerateTimer()
    
    def allPopsReportForDuty(self):
        for j in range(len(self.hostPopulations)):
            self.hostPopulations[j].reportPopInfo()

    def displayResultsAsTimecourses(self):
        if (self.completedRun == False):
            print("\nYou need to generate some data first!\n")
        else:
            for j in range(len(self.hostPopulations)):
                plt.plot(self.timeValues, self.hostPopulationValues[j, :])