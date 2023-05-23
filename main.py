#!/usr/local/bin/python3

import sys
import math
import random
import turtle
import time as tt
import pandas.core.frame
import datetime
from multiprocessing import Pool
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from tqdm import tqdm

from multiprocessing import Pool, TimeoutError, Array

# taken from https://www.ncl.ucar.edu/Document/Graphics/named_colors.shtml
colors = ("black","RoyalBlue","LightSkyBlue",\
             "PowderBlue","lightseagreen","PaleGreen","Wheat","Brown",\
             "Pink")

def init_drawing():
    global screen
    global point_history

    screen = turtle.Screen()
    screen.clear()
    screen.tracer(False)
    global don
    don = turtle.Turtle()
    don.speed(0)
    don.width(3)
    don.hideturtle()
    don.radians()

def draw_pendulum(x, y, dP, trace_history=False):
    #print('Drawing %s'%(dP,))
    don.penup()
    don.goto((x,y))
    don.pendown()
    don.dot(10)
    don.setheading(dP.angleUpper - math.pi/2)

    don.forward(dP.lenUpper * 100)
    don.dot(10)
    don.setheading(dP.angleLower - math.pi/2)
    don.forward(dP.lenLower * 100)
    don.dot(10)
    if trace_history:
        dP.point_history.append(don.position())

        for p in dP.point_history:
            don.penup()
            don.goto(p)
            don.pendown()
            don.dot(2, "blue")

time_last_drawn=None
def draw_screen(pendulums, force=False, trace_history=False):
    global time_last_drawn

    #draw at most 60fps
    if not force and time_last_drawn and (pendulums[0].t - time_last_drawn) < 1.0 / 60:
        return
    #All the pendulums should have the same time, so just pick the first one
    time_last_drawn = pendulums[0].t

    don.home()
    don.clear()

    for i in range(len(pendulums)):
      pendulum=pendulums[i]
      don.pencolor(colors[i % len(colors)])
      draw_pendulum(-animation_start_x+animation_spacing_x*i, 0, pendulum, trace_history)

      if show_animation_data:
        don.penup()
        don.goto((-400, 300-i*25))
        don.pendown()
        don.write("%s" % (pendulum,), True, font=("Arial", 10, "normal"))
    screen.update()

def normalize_angle(angle):
    if (angle >= 2*math.pi):
      return angle-2*math.pi
    if ( angle<0 ):
      return angle+2*math.pi
    return angle

class DoublePendulum:
    t: float
    myName: str
    deltaT: float
    angleUpper: float
    angleLower: float
    velUpper: float
    velLower: float
    tipX: float
    tipY: float

    initialAngleUpper: float
    initialAngleLower: float
    initialVelUpper: float
    initialVelLower: float

    lenUpper: float
    lenLower: float
    massUpper: float
    massLower: float
    g: float

    point_history=[]

    def __init__(self, myName: str, state: list, deltaT=0.001, creationMessage = True):
        self.myName=myName
        self.setInitialPendulumState(state)
        self.setDeltaT(deltaT)
        if creationMessage:
            print('Created pendulum %s' %(self,))

    def __str__(self):
        return 'P%s [t=%.1f tick=%0.2e sec (%.0f/sec)]: %.1fpi|%.1fpi/s - %.1fm - %.1fkg- %.1fpi|%.1fpi/s - %.1fm - %.1fkg' % (self.myName, self.t, self.deltaT, 1/self.deltaT, self.angleUpper/math.pi, self.velUpper/math.pi, self.lenUpper, self.massUpper, self.angleLower/math.pi, self.velLower/math.pi, self.lenLower, self.massLower)

    def setDeltaT(self, deltaT):
        self.deltaT = deltaT

    def setInitialPendulumState(self, state: list):
        self.t=0
        self.initialAngleUpper = self.angleUpper = state[0] % (2 * math.pi)
        self.initialAngleLower = self.angleLower = state[1] % (2 * math.pi)
        self.initialVelUpper = self.velUpper = state[2]
        self.initialVelLower = self.velLower = state[3]
        self.lenUpper = 1
        self.lenLower = 1
        self.massUpper = 1
        self.massLower = 1
        self.g = 9.81
        self.point_history = []

    def resetToInitialState(self):
        self.angleUpper = self.initialAngleUpper
        self.angleLower = self.initialAngleLower
        self.velUpper = self.initialVelUpper
        self.velLower = self.initialVelLower
        self.t=0
        self.point_history=[]

    def tickAngleUpper(self):
        return self.angleUpper + self.deltaT * self.velUpper

    def tickAngleLower(self):
        return self.angleLower + self.deltaT * self.velLower

    def tickVelUpper(self):
        return self.velUpper + self.deltaT * ((-self.g * (2 * self.massUpper + self.massLower) * math.sin(self.angleUpper) - self.massLower * self.g * math.sin(self.angleUpper - 2 * self.angleLower) - 2 * math.sin(self.angleUpper - self.angleLower) * self.massLower * (math.pow(self.velLower, 2) * self.lenLower + math.pow(self.velUpper, 2) * self.lenUpper * math.cos(self.angleUpper - self.angleLower))) / (self.lenUpper * (2 * self.massUpper + self.massLower - self.massLower * math.cos(2 * self.angleUpper - 2 * self.angleLower))))

    def tickVelLower(self):
        return self.velLower + self.deltaT * ((2 * math.sin(self.angleUpper - self.angleLower) * (math.pow(self.velUpper, 2) * self.lenUpper * (self.massUpper + self.massLower) + self.g * (self.massUpper + self.massLower) * math.cos(self.angleUpper) + math.pow(self.velLower, 2) * self.lenLower * self.massLower * math.cos(self.angleUpper - self.angleLower))) / (self.lenLower * (2 * self.massUpper + self.massLower - self.massLower * math.cos(2 * self.angleUpper - 2 * self.angleLower))))

    def doTick(self):
        #print('doTick(%.2f, %.5e)' % (self.t,self.deltaT))
        newAngleUpper = self.tickAngleUpper() % (2 * math.pi)
        newAngleLower = self.tickAngleLower() % (2 * math.pi)
        newVelUpper = self.tickVelUpper()
        newVelLower = self.tickVelLower()

        self.angleLower = newAngleLower
        self.angleUpper = newAngleUpper
        self.velUpper = newVelUpper
        self.velLower = newVelLower
        self.t += self.deltaT

    def doAllTicks(self, stopT):
        # print('Processing pendulum %s' % (self,))
        while self.t < stopT:
            self.doTick()

    def getState(self):
        return [self.angleUpper, self.angleLower, self.velUpper, self.velLower]

    def getInitialState(self):
        return [self.initialAngleUpper, self.initialAngleLower, self.initialVelUpper, self.initialVelLower]

    def updCartesian(self):
        self.tipX = math.sin(self.angleUpper) + math.sin(self.angleLower)
        self.tipY = -math.cos(self.angleUpper) - math.cos(self.angleLower)

    def beginTrace(self, iterations):
        self.pathX = np.zeros(iterations)
        self.pathY = np.zeros(iterations)
        self.velHisUp = np.zeros(iterations)
        self.velHisLow = np.zeros(iterations)

    def raceShadowP(self, scaleFactor, divergenceThreshold, maxTime):
        p2 = DoublePendulum("shadow", [parameter * scaleFactor for parameter in self.getInitialState()], self.deltaT, creationMessage = False)
        diffUnsigned = 0
        while (self.t <= maxTime) and (diffUnsigned < divergenceThreshold):
            for pendulum in [self , p2]:
                pendulum.doAllTicks(pendulum.t + 1.0 / 63)
                pendulum.updCartesian()
            diff = math.sqrt((self.tipX - p2.tipX) ** 2 + (self.tipY - p2.tipY) ** 2)
            diffUnsigned += diff
        return self


def generateInitialState(range):
    return [random.random() * math.pi * 2, random.random() * math.pi * 2, random.random() * range, random.random() * range]


def process_pendulum(p, scaleFactor, divergenceThreshold, maxTime):
    p.raceShadowP(scaleFactor, divergenceThreshold, maxTime)
    return p.t
#    global resultArray
#    resultArray[int(p.myName)] = p.t


def animate_pendulums(args):
    global deltaT_round0
    global pendulum_count
    global pendulums

    # expecting at least 2 arguments: pendulum#, then how many times to divide deltaT by 2
    #   3 0 1 2   <-- Pendulum[3], updated with deltaT=0.01  (0.01/2^0), and deltaT=0.005 and deltaT=0.0025
    pendulum_number = int(args[0])

    pendulums_to_animate = []
    pendulum = pendulums[pendulum_number]

    # The rest of the args specify the 2^n that divides deltaT for the pendulums
    for i in range(2, len(args)):
        deltaT_divider = int(args[i])
        pendulums_to_animate.append(
            DoublePendulum('%s/%d' % (pendulum_number, math.pow(2, deltaT_divider)), pendulum.getState(),
                           deltaT_round0 / math.pow(2, deltaT_divider)))

    workerPool = Pool(processes=cpus)  # start worker processes

    init_drawing()

    t = 0
    draw_screen(pendulums_to_animate, force=True)
    tt.sleep(5)
    draw_interval = 1.0 / 60

    prev_clock = tt.time()

    while t < simulation_duration:
        # make sure we're not looping faster than 50ms
        cur_clock = tt.time()
        loop_delay = cur_clock - prev_clock
        if loop_delay < 0.05:
            tt.sleep(0.05 - loop_delay)
        prev_clock = cur_clock

        t += draw_interval
        # move all the pendulums along in their simulation
        async_results = [workerPool.apply_async(process_pendulum, (p, t)) for p in pendulums_to_animate]
        for i in range(len(pendulums_to_animate)):
            pendulums_to_animate[i] = async_results[i].get()

        draw_screen(pendulums_to_animate, force=True)

    print('Sleeping for 90 seconds...')
    tt.sleep(90)


def run_pendulums_collect_data():
    global deltaT_round0
    global pendulum_count
    global pendulums

    currTime = datetime.datetime.now()
    prevTime: datetime

    pendulumCurrentStates = []
    pendulumPreviousStates = []
    pendulumInitialStates = []

    upperAngleDataFrame = pandas.DataFrame([])
    lowerAngleDataFrame = pandas.DataFrame([])
    upperVelDataFrame = pandas.DataFrame([])
    lowerVelDataFrame = pandas.DataFrame([])

    for i in range(pendulum_count):
        pendulumCurrentStates.append(pendulums[i].getState())
        pendulumInitialStates.append(pendulums[i].getInitialState())

    InitialStatesFrame = pandas.DataFrame(pendulumInitialStates,
                                          columns=["Upper Angle", "Lower Angle", "Upper Vel", "Lower Vel"])
    InitialStatesFrame.to_excel('InitialData.xlsx', float_format="%.9f", index=False)

    workerPool = Pool(processes=cpus)  # start worker processes

    for n in range(rounds):
        # print(pendulumCurrentStates)
        deltaT = deltaT_round0 / math.pow(2, n)
        print('Running pendulums: %d ticks/sec (deltaT=%.3e)' % (1 / deltaT, deltaT))

        # Set all the pendulums' deltaT to be the current deltaT value
        for p in pendulums:
            p.deltaT = deltaT

        # Start running through the pendulums in separate processes
        async_results = [workerPool.apply_async(process_pendulum, (p, simulation_duration)) for p in pendulums]

        # Collect the results from the workers
        for i in range(pendulum_count):
            pendulums[i] = async_results[i].get()

        # Save pendulums' states and reset them back to their initial states
        pendulumCurrentStates = []
        for p in pendulums:
            pendulumCurrentStates.append(p.getState())
            p.resetToInitialState()

        currentData = pandas.DataFrame(pendulumCurrentStates,
                                       columns=[deltaT.__str__() + " Upper Angle", deltaT.__str__() + " Lower Angle",
                                                deltaT.__str__() + " Upper Vel", deltaT.__str__() + " Lower Vel"])

        upperAngleDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Upper Angle")
        lowerAngleDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Lower Angle")
        upperVelDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Upper Vel")
        lowerVelDataFrame[deltaT] = currentData.get(deltaT.__str__() + " Lower Vel")

        pendulumPreviousStates = pendulumCurrentStates
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

def graph_evolution(iThetaA: float, iThetaB: float, iVelA: float, iVelB: float, framerate: int = 63, simTime: int = 20):
    """
    iThetaA: initial angle of first arm (in radians, so 1/2 corresponds to pi radians or 180 degrees)
    iThetaB: initial angle of second arm
    iVelA: initial angular velocity of first arm
    iVelB: initial angular velocity of second arm
    """

    init_drawing()

    localDT = 0.0001
    p = DoublePendulum("pendy", [iThetaA * math.pi, iThetaB * math.pi, iVelA * math.pi, iVelB * math.pi], localDT)

    iterations = framerate * simTime
    t = np.linspace(0, simTime, iterations)
    pathX = np.zeros(iterations)
    pathY = np.zeros(iterations)
    velUp = np.zeros(iterations)
    velLow = np.zeros(iterations)
    for i in range(iterations):
        p.doAllTicks(p.t + 1.0/framerate)
        #print(p)
        draw_screen([p])
        print("%.5f\t"%p.t, end="")
        print("%.8f\t%.8f\t%.8f\t%.8f" % tuple(p.getState()))
        pathX[i] = math.sin(p.angleUpper) + math.sin(p.angleLower)
        pathY[i] = -math.cos(p.angleUpper) - math.cos(p.angleLower)
        velUp[i] = p.velUpper
        velLow[i] = p.velLower

    pathPoints = np.array([pathX, pathY]).T.reshape(-1, 1, 2)
    pathSegments = np.concatenate([pathPoints[:-1], pathPoints[1:]], axis=1)
    velPoints = np.array([velUp, velLow]).T.reshape(-1, 1, 2)
    velSegments = np.concatenate([velPoints[:-1], velPoints[1:]], axis=1)

    fig, axs = plt.subplots(2, 1, sharex = False, sharey = False)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(0, simTime)
    pathlc = LineCollection(pathSegments, cmap = 'rainbow', norm = norm)
    vellc = LineCollection(velSegments, cmap = 'rainbow', norm = norm)
    # Set the values used for colormapping
    pathlc.set_array(t)
    pathlc.set_linewidth(2)

    vellc.set_array(t)
    vellc.set_linewidth(2)

    pathLine = axs[0].add_collection(pathlc)
    fig.colorbar(pathLine, ax=axs[0])

    velLine = axs[1].add_collection(vellc)
    fig.colorbar(velLine, ax=axs[1])

    axs[0].set_xlim(-2.125, 2.125)
    axs[0].set_ylim(-2.125, 2.125)
    axs[1].set_xlim(velUp.min()*1.05, velUp.max()*1.05)
    axs[1].set_ylim(velLow.min()*1.05, velLow.max()*1.05)
    plt.show()
    #input("Press Enter to stop...")


def compare_two_pendulums(iThetaA: float, iThetaB: float, iVelA: float, iVelB: float, secondPScaling: list = [1.01, 1.0, 1.0, 1.0], framerate: int = 63, simTime: int = 20):
    """
    iThetaA: initial angle of first arm (in radians, so 1/2 corresponds to pi radians or 180 degrees)
    iThetaB: initial angle of second arm
    iVelA: initial angular velocity of first arm
    iVelB: initial angular velocity of second arm
    secondPScaling: how much each of the initial values (thetaA, thetaB, velA, velB) should be scaled by for the second pendulum
    """

    init_drawing()

    localDT = 0.0001
    p1 = DoublePendulum("pendy", [iThetaA * math.pi, iThetaB * math.pi, iVelA * math.pi, iVelB * math.pi], localDT)
    p2 = DoublePendulum("offPendy", [value*scale for value,scale in zip([iThetaA * math.pi, iThetaB * math.pi, iVelA * math.pi, iVelB * math.pi], secondPScaling)], localDT)

    iterations = framerate * simTime
    t = np.linspace(0, simTime, iterations)
    p1.beginTrace(iterations)
    p2.beginTrace(iterations)
    diffUnsigned = np.zeros(iterations)
    for i in range(iterations):
        for pendulum in [p1, p2]:
            pendulum.doAllTicks(pendulum.t + 1.0 / framerate)

            pendulum.pathX[i] = math.sin(pendulum.angleUpper) + math.sin(pendulum.angleLower)
            pendulum.pathY[i] = -math.cos(pendulum.angleUpper) - math.cos(pendulum.angleLower)
            pendulum.velHisUp[i] = pendulum.velUpper
            pendulum.velHisLow[i] = pendulum.velLower
            """
            print("%.5f\t"%pendulum.t, end="")
            print("%.8f\t%.8f\t%.8f\t%.8f" % tuple([v/math.pi for v in pendulum.getState()]))
        print("")  # newline to make comparing values easier (pair up values)
        """
        diff = math.sqrt((p1.pathX[i] - p2.pathX[i])**2 + (p1.pathY[i] - p2.pathY[i])**2)
        # diffSquares[i] = diff**2 + diffSquares[i-1]  # for the first element, the [i-1] is the last element, which will be 0, so the wraparound is helpful because no edge case conditional is needed (it gives 0)
        diffUnsigned[i] = diff + diffUnsigned[i-1]
        # draw_screen([p1, p2])

    path1Points = np.array([p1.pathX, p1.pathY]).T.reshape(-1, 1, 2)
    path1Segments = np.concatenate([path1Points[:-1], path1Points[1:]], axis=1)
    path2Points = np.array([p2.pathX, p2.pathY]).T.reshape(-1, 1, 2)
    path2Segments = np.concatenate([path2Points[:-1], path2Points[1:]], axis=1)

    fig, axs = plt.subplots(2, 2, sharex = False, sharey = False)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(0, simTime)
    path1lc = LineCollection(path1Segments, cmap = 'rainbow', norm = norm)
    path2lc = LineCollection(path2Segments, cmap = 'rainbow', norm = norm)
    # Set the values used for colormapping
    path1lc.set_array(t)
    path1lc.set_linewidth(0.9)

    path2lc.set_array(t)
    path2lc.set_linewidth(0.9)

    path1Line = axs[0, 0].add_collection(path1lc)
    fig.colorbar(path1Line, ax=axs[0, 0])

    path2Line = axs[1, 0].add_collection(path2lc)
    fig.colorbar(path2Line, ax=axs[1, 0])

    axs[0, 0].set_xlim(-2.125, 2.125)
    axs[0, 0].set_ylim(-2.125, 2.125)
    axs[1, 0].set_xlim(-2.125 , 2.125)
    axs[1, 0].set_ylim(-2.125 , 2.125)

    axs[0 , 0].set_title('Pendulum 1')
    axs[1 , 0].set_title('Pendulum 2')
    axs[0 , 1].set_title('Abs Diff')
    axs[0 , 1].plot(t , diffUnsigned)
    axs[1 , 1].set_title('Log10(Abs Diff)')
    axs[1 , 1].plot(t , [math.log(v,10) for v in diffUnsigned])

    plt.show()
    #input("Press Enter to stop...")


def phase_diagram(iVelA: float, iVelB: float, resolution: int = 10, epsilon: float = 10.0, secondPScaling: list = [1.0, 1.001, 1.0, 1.0], framerate: int = 63, maxSimTime: int = 30):
    """
    iVelA: initial angular velocity of first arm
    iVelB: initial angular velocity of second arm
    resolution: bins of numpy array (will have shape of (resolution, resolution) )
    epsilon: when the pendulums are considered to have diverged
    secondPScaling: how much each of the initial values (thetaA, thetaB, velA, velB) should be scaled by for the second pendulum
    """

    # init_drawing()

    localDT = 0.0001

    timeToDiverge = np.empty(shape = (resolution,resolution))
    p1 = DoublePendulum("pendy" , [0 , 0 , iVelA * math.pi , iVelB * math.pi] , localDT)
    p2 = DoublePendulum("offPendy" , [value * scale for value , scale in zip([0 , 0 , iVelA * math.pi , iVelB * math.pi] , secondPScaling)] , localDT)
    for thetaAIndex,thetaA in enumerate(np.linspace(0, math.pi, resolution)):
        for thetaBIndex,thetaB in enumerate(np.linspace(0 , math.pi , resolution)):
            p1.setInitialPendulumState([thetaA, thetaB, iVelA * math.pi, iVelB * math.pi])
            p2.setInitialPendulumState([value * scale for value , scale in zip([thetaA , thetaB , iVelA * math.pi , iVelB * math.pi] , secondPScaling)])
            p1.resetToInitialState()
            p2.resetToInitialState()

            diffUnsigned = 0
            while (p1.t <= maxSimTime) and (diffUnsigned < epsilon):
                for pendulum in [p1, p2]:
                    pendulum.doAllTicks(pendulum.t + 1.0 / framerate)

                    pendulum.tipX = math.sin(pendulum.angleUpper) + math.sin(pendulum.angleLower)
                    pendulum.tipY = -math.cos(pendulum.angleUpper) - math.cos(pendulum.angleLower)
                    """
                    print("%.5f\t"%pendulum.t, end="")
                    print("%.8f\t%.8f\t%.8f\t%.8f" % tuple([v/math.pi for v in pendulum.getState()]))
                print("")  # newline to make comparing values easier (pair up values)
                """
                diff = math.sqrt((p1.tipX - p2.tipX)**2 + (p1.tipY - p2.tipY)**2)
                diffUnsigned += diff
            timeToDiverge[thetaAIndex, thetaBIndex] = p1.t
            print("%.9f\t"%timeToDiverge[thetaAIndex, thetaBIndex], end="")
            print([thetaAIndex,thetaBIndex])

    fig , axs = plt.subplots(1 , 1 , constrained_layout = True , squeeze = False)
    psm = plt.pcolormesh(timeToDiverge , rasterized = True , cmap = 'hot')
    fig.colorbar(psm)

    plt.show()
    #input("Press Enter to stop...")

def run_phase_diagram_collect_data():
    global resolution
    global pendulums
    global shadowScaleFactor
    global shadowDivergenceThreshold
    global shadowMaxRuntime
    global iVelA
    global iVelB
    global upperThetaMin
    global upperThetaMax
    global lowerThetaMin
    global lowerThetaMax

    for thetaAIndex,thetaA in enumerate(np.linspace(upperThetaMin * math.pi, upperThetaMax * math.pi, resolution)):
        for thetaBIndex,thetaB in enumerate(np.linspace(lowerThetaMin * math.pi , lowerThetaMax * math.pi , resolution)):
            pendulums.append(DoublePendulum(str(resolution*thetaAIndex+thetaBIndex), [thetaA,thetaB,iVelA * math.pi, iVelB * math.pi], deltaT = 0.001))

    workerPool = Pool(processes = cpus)  # start worker processes

    # Start running through the pendulums in separate processes
    async_results = [workerPool.apply_async(process_pendulum , (p , shadowScaleFactor, shadowDivergenceThreshold, shadowMaxRuntime)) for p in pendulums]

    # Collect the results from the workers
    # timeToDiverge = np.zeros((resolution , resolution))
    timeToDivergeDF = pandas.DataFrame(np.zeros((resolution , resolution)))
    excelSheetName = str("u=%.3f-%.3fl=%.3f-%.3ft=%d"%tuple([upperThetaMin , upperThetaMax , lowerThetaMin , lowerThetaMax , shadowMaxRuntime]))
        # t= : max time
        # u= : upper angle range (min-max)
        # l= : lower angle range (min-max)
        # resolution is not included because it is simply the number of rows/columns

    rowTimeToDiverge = np.zeros(resolution)
    print("Starting simulation at ", datetime.datetime.now().replace(microsecond = 0))
    for i in tqdm(range(resolution)):
        for j in range(resolution):
            rowTimeToDiverge[j] = async_results[i*resolution+j].get()
        timeToDivergeDF[i:(i+1)] = rowTimeToDiverge

    with pandas.ExcelWriter('AData.xlsx' , mode = 'w') as writer:
        timeToDivergeDF.columns = np.linspace(upperThetaMin , upperThetaMax , resolution)
        timeToDivergeDF.index = np.linspace(lowerThetaMin , lowerThetaMax , resolution)
        timeToDivergeDF.to_excel(writer , sheet_name = excelSheetName , float_format = "%.9f")


def plot_phase_diagram_from_excel(filepath: str, sheetName = 0):
    df = pandas.read_excel(filepath, sheet_name = sheetName, index_col = 0)

    fig , axs = plt.subplots(1 , 1 , constrained_layout = True , squeeze = False)
    psm = plt.pcolormesh(df , rasterized = True , cmap = 'hot')
    fig.colorbar(psm)
    plt.show()


collectData = True

spreadsheetFilepath = "AData.xlsx"

resolution = 16
shadowScaleFactor = 1.001
shadowDivergenceThreshold = 10
shadowMaxRuntime = 60

# in pi radians (so a value of 0.5 corresponds to 0.5pi radians or 90 degrees)
iVelA = 0.75
iVelB = 0.75
upperThetaMin = 0.0
upperThetaMax = 2.0
lowerThetaMin = 0.0
lowerThetaMax = 2.0


pendulums = []
cpus=8
deltaT_round0=0.01
rounds=10
pendulum_count=100
simulation_duration = 30
animation_start_x=275
animation_spacing_x=220
show_animation_data=True
seed=tt.time()
initialVars = [0.5 , 0.5 , 0 , 1.5]
scalingVars = [1.001, 1.001, 1.001, 1.001]


def main(argv):
    # compare_two_pendulums(float(argv[0]), float(argv[1]), float(argv[2]), float(argv[3]))
    # compare_two_pendulums(iThetaA = initialVars[0], iThetaB = initialVars[1], iVelA = initialVars[2], iVelB = initialVars[3], secondPScaling = scalingVars)
    # phase_diagram(iVelA = 1.0, iVelB = 1.0, resolution = 240, epsilon = 10, secondPScaling = scalingVars)
    if collectData:
        run_phase_diagram_collect_data()
    plot_phase_diagram_from_excel(spreadsheetFilepath)

def xmain(argv):
    global cpus
    global deltaT_round0
    global rounds
    global pendulum_count
    global simulation_duration
    global animation_start_x
    global animation_spacing_x
    global show_animation_data
    global seed
    global pendulums

    mode = 'run_and_save_data'
    usage='test.py [--pendulum_count <number>] [--cpus <number>] [--duration <secs>] [--deltaT_round0 <interval>] [--rounds <rounds>] [--seed <random seed>] [--mode <mode>] [--animation_start_x <x-coord>] [--animation_spacing_x <pixels>] [--hide_animation_data]'
    try:
        opts, args = getopt.getopt(argv,"h",["help", "cpus=","deltaT_round0=","rounds=","seed=","mode=","duration=","animation_start_x=","animation_spacing_x=","hide_animation_data"])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h' or opt == '--help':
            print(usage)
            sys.exit()
        elif opt == '--cpus':
            cpus=int(arg)
        elif opt == '--duration':
            simulation_duration=float(arg)
        elif opt == '--deltaT_round0':
            deltaT_round0=float(arg)
        elif opt == '--rounds':
            rounds=int(arg)
        elif opt == '--seed':
            seed=int(arg)
        elif opt == '--pendulum_count':
            pendulum_count=int(arg)
        elif opt == '--animation_start_x':
            animation_start_x=int(arg)
        elif opt == '--animation_spacing_x':
            animation_spacing_x=int(arg)
        elif opt == '--hide_animation_data':
            show_animation_data=False
        elif opt == '--mode':
            mode=arg

    random.seed(seed)

    print("Generating pendulums with seed %d" % seed)

    for i in range(pendulum_count):
        pendulums.append(DoublePendulum("Pendulum #%d" % i, generateInitialState(1 / 2 * math.pi)))


    if mode == 'run_and_save_data':
      run_pendulums_collect_data()
    elif mode == 'animate':
      animate_pendulums(args)
    else:
      print('Unknown mode')
      sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])


