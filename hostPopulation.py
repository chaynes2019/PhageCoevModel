class hostPopulation():
    def __init__(self, hostPopNum, initPopLevel, fitnessParam):
        self.populationIndex = hostPopNum
        self.populationLevel = initPopLevel
        self.fitnessParameter = fitnessParam
        self.timerValue = 0
        self.growthRate = 0

    def computeGrowthRate(self, hostPopulations, carryingCapacity):
        populationSum = sum(hostPopulations)
        individualFitness = self.fitnessParameter * (1 - (populationSum / carryingCapacity))
        #This is the in-program representation of growthRate = fitness(p) * p, where p is the number of cells in a given population
        self.growthRate = individualFitness*self.populationLevel
    
    def setTimer(self, val):
        self.timerValue = val

    def getTimer(self):
        return self.timerValue

    def reportPopInfo(self):
        print("\nInformation for Population {}".format(self.populationIndex))
        print("-----------------------------")
        print("Population Level: {}".format(self.populationLevel))
        print("Fitness Parameter: {}".format(self.fitnessParameter))
        print("Would-be Exponential Growth Rate: {}".format(self.fitnessParameter * self.populationLevel))
        print("Logistic Growth Rate: {}".format(self.growthRate))
        print("Timer Value: {}".format(self.getTimer()))