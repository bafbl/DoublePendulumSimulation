import math
import random
import turtle
import time as thread
import pandas.core.frame
import datetime

def init_drawing():
    global screen
    screen = turtle.Screen()
    screen.clear()
    screen.tracer(0)
    global don
    don = turtle.Turtle()
    don.speed(0)
    don.width(3)
    don.hideturtle()
    don.radians()

def draw_pendulum(dP):
    don.dot(10)
    don.setheading(dP.angleUpper - math.pi/2)

    don.forward(dP.lenUpper * 100)
    don.dot(10)
    don.setheading(dP.angleLower - math.pi/2)
    don.forward(dP.lenLower * 100)
    don.dot(10)

def draw_screen(i, t, pendulum):
    global time_last_drawn

    if t > 0 and (t - time_last_drawn) < 1.0 / 30:
        return

    don.home()
    don.clear()

    time_last_drawn = t
    draw_pendulum(pendulum)
    # write loop count
    don.penup()
    don.goto((250, 250))
    don.pendown()
    don.write("loop #%d, t=%0.3fsecs" % (i, t), True, font=("Arial", 12, "normal"))
    # write time
    don.penup()
    don.goto((250, 225))
    don.pendown()
    don.write("Angles: %+6.2f, %+6.2f" % (pendulum.angleUpper / math.pi, pendulum.angleLower / math.pi), font=("Arial", 12, "normal"))
    screen.update()

class DoublePendulum:
    angleUpper: float
    angleLower: float
    velUpper: float
    velLower: float

    initialAngleUpper: float
    initialAngleLower: float
    initialVelUpper: float
    initialVelLower: float

    lenUpper: float
    lenLower: float
    massUpper: float
    massLower: float
    g: float

    def __init__(self, state: list):
        self.setPendulumState(state)
        self.setInitialPendulumState(state)

    def setPendulumState(self, state: list):
        self.angleUpper = state[0]
        self.angleLower = state[1]
        self.velUpper = state[2]
        self.velLower = state[3]

    def setInitialPendulumState(self, state: list):
        self.initialAngleUpper = state[0]
        self.initialAngleLower = state[1]
        self.initialVelUpper = state[2]
        self.initialVelLower = state[3]
        self.lenUpper = 1
        self.lenLower = 1
        self.massUpper = 1
        self.massLower = 1
        self.g = 9.81

    def tickAngleUpper(self, deltaT):
        return self.angleUpper + deltaT * self.velUpper

    def tickAngleLower(self, deltaT):
        return self.angleLower + deltaT * self.velLower

    def tickVelUpper(self, deltaT):
        return self.velUpper + deltaT * ((-self.g * (2 * self.massUpper + self.massLower) * math.sin(self.angleUpper) - self.massLower * self.g * math.sin(self.angleUpper - 2 * self.angleLower) - 2 * math.sin(self.angleUpper - self.angleLower) * self.massLower * (math.pow(self.velLower, 2) * self.lenLower + math.pow(self.velUpper, 2) * self.lenUpper * math.cos(self.angleUpper - self.angleLower))) / (self.lenUpper * (2 * self.massUpper + self.massLower - self.massLower * math.cos(2 * self.angleUpper - 2 * self.angleLower))))

    def tickVelLower(self, deltaT):
        return self.velLower + deltaT * ((2 * math.sin(self.angleUpper - self.angleLower) * (math.pow(self.velUpper, 2) * self.lenUpper * (self.massUpper + self.massLower) + self.g * (self.massUpper + self.massLower) * math.cos(self.angleUpper) + math.pow(self.velLower, 2) * self.lenLower * self.massLower * math.cos(self.angleUpper - self.angleLower))) / (self.lenLower * (2 * self.massUpper + self.massLower - self.massLower * math.cos(2 * self.angleUpper - 2 * self.angleLower))))

    def doTicks(self, deltaT):
        return [self.tickAngleUpper(deltaT), self.tickAngleLower(deltaT), self.tickVelUpper(deltaT), self.tickVelLower(deltaT)]

    def getState(self):
        return [self.angleUpper, self.angleLower, self.velUpper, self.velLower]

    def getInitialState(self):
        return [self.initialAngleUpper, self.initialAngleLower, self.initialVelUpper, self.initialVelLower]


def generateInitialState(range):
    return [random.random() * math.pi * 2, random.random() * math.pi * 2, random.random() * range, random.random() * range]


init_drawing()

pendulums = []

currTime = datetime.datetime.now()
prevTime: datetime

pendulumCurrentStates = []
pendulumPreviousStates = []
pendulumInitialStates = []

upperAngleDataFrame = pandas.DataFrame([])
lowerAngleDataFrame = pandas.DataFrame([])
upperVelDataFrame = pandas.DataFrame([])
lowerVelDataFrame = pandas.DataFrame([])

deltaT = 0.01

time = 0

for i in range(100):
    pendulums.append(DoublePendulum(generateInitialState(1/2 * math.pi)))
    pendulumCurrentStates.append(pendulums[i].getState())
    pendulumInitialStates.append(pendulums[i].getInitialState())

InitialStatesFrame = pandas.DataFrame(pendulumInitialStates, columns=["Upper Angle", "Lower Angle", "Upper Vel", "Lower Vel"])
InitialStatesFrame.to_excel('InitialData.xlsx', float_format="%.9f", index=False)

for n in range(20):
    for i in range(round(10 / deltaT)):
        #print(pendulumCurrentStates)
        #print(time)
        for j in pendulums:
            j.setPendulumState(j.doTicks(deltaT))
        #for j in pendulums:
        #    pendulumCurrentStates[pendulums.index(j)] = j.getState()
        #draw_screen(i, time, pendulums[0])
        time += deltaT

    for i in pendulums:
        pendulumCurrentStates[pendulums.index(i)] = i.getState()
        i.setPendulumState(i.getInitialState())

    pendulums[0] = DoublePendulum([math.pi/2, math.pi/2, 0 ,0])

    currentData = pandas.DataFrame(pendulumCurrentStates, columns=[deltaT.__str__() + " Upper Angle", deltaT.__str__() + " Lower Angle", deltaT.__str__() + " Upper Vel", deltaT.__str__() + " Lower Vel"])

    upperAngleDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Upper Angle")
    lowerAngleDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Lower Angle")
    upperVelDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Upper Vel")
    lowerVelDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Lower Vel")

    pendulumPreviousStates = pendulumCurrentStates
    print(deltaT)
    deltaT = deltaT / 2
    time = 0
    with pandas.ExcelWriter('Data.xlsx', mode='w') as writer:
        upperAngleDataFrame.to_excel(writer, sheet_name="UpperAngles", float_format="%.9f", index=False)
        lowerAngleDataFrame.to_excel(writer, sheet_name="LowerAngles", float_format="%.9f", index=False)
        upperVelDataFrame.to_excel(writer, sheet_name="UpperVelocities", float_format="%.9f", index=False)
        lowerVelDataFrame.to_excel(writer, sheet_name="LowerVelocities", float_format="%.9f", index=False)
    print("spreadsheet updated")
    prevTime = currTime
    currTime = datetime.datetime.now()
    delta = currTime - prevTime
    expected = currTime + 2 * delta
    print("Time taken: " + delta.__str__() + " Expected next: " + expected.__str__())

print(upperAngleDataFrame)
