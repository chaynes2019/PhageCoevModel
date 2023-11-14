import numpy as np
from random import randrange
class phagePopulation():
    def __init__(self, phagePopIdx, initialPopLevel, burstNum, lysisLength):
        self.populationIndex = phagePopIdx
        self.populationLevel = initialPopLevel
        self.burstSize = burstNum
        self.lysisTime = lysisLength
        self.inverseLysisRandomResolution = 1000
        

        if (initialPopLevel > 0):
            self.timerValue = self.generateLysisTime()
        else:
            self.timerValue = float('inf')

    def setTimer(self, val):
        self.timerValue = val

    def getTimer(self):
        return self.timerValue

    def updateTimer(self, minTimerValue):
        if ((minTimerValue == float('inf')) | (self.getTimer() == float('inf'))):
            self.setTimer(float('inf'))
        else:
            self.setTimer(self.getTimer() - minTimerValue)

    def regenerateTimer(self):   
        if (self.populationLevel > 0):
            if (self.getTimer() == 0):
                self.setTimer(self.generateLysisTime())
            elif (self.getTimer() == float('inf')):
                self.setTimer(self.generateLysisTime())
        #This block underneath covers the case of an extinct population. What if
        #circumstances change and it needs to start generating event times again?
        else:
            #Assigns extinct populations a time of positive infinity
            self.setTimer(float('inf'))

    def generateLysisTime(self):
        return self.lysisTime + randrange(-self.inverseLysisRandomResolution, self.inverseLysisRandomResolution, 1) / self.inverseLysisRandomResolution
    
    def changePopulation(self, changeVal):
        self.populationLevel += changeVal
        if (self.populationLevel < 0):
            self.populationLevel = 0
    
    
    def reportPopInfo(self):
        print("\nInformation for Phage Population {}".format(self.populationIndex))
        print("-----------------------------")
        print("Population Level: {}".format(self.populationLevel))
        print("Burst Size: {}".format(self.burstSize))
        print("Lysis Time: {}".format(self.lysisTime))
        print("Timer Value: {}".format(self.getTimer()))
        