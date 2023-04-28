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

from multiprocessing import Pool, TimeoutError

# taken from https://www.ncl.ucar.edu/Document/Graphics/named_colors.shtml
colors = ("black","RoyalBlue","LightSkyBlue",\
             "PowderBlue","lightseagreen","PaleGreen","Wheat","Brown",\
             "Pink")

def init_drawing():
    global screen
    screen = turtle.Screen()
    screen.clear()
    screen.tracer(False)
    global don
    don = turtle.Turtle()
    don.speed(0)
    don.width(3)
    don.hideturtle()
    don.radians()

def draw_pendulum(x, y, dP):
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

time_last_drawn=None
def draw_screen(pendulums, force=False):
    global time_last_drawn

    #draw at most 60fps
    if not force and time_last_drawn and (pendulum.t - time_last_drawn) < 1.0 / 60:
        return
    #All the pendulums should have the same time, so just pick the first one
    time_last_drawn = pendulums[0].t

    don.home()
    don.clear()

    for i in range(len(pendulums)):
      pendulum=pendulums[i]
      don.pencolor(colors[i % len(colors)])
      draw_pendulum(-animation_start_x+animation_spacing_x*i, 0, pendulum)

      if show_animation_data:
        don.penup()
        don.goto((-400, 300-i*25))
        don.pendown()
        don.write("%s" % (pendulum,), True, font=("Arial", 10, "normal"))
    screen.update()

def normalize_angle(angle):
    if (angle >= 2*math.pi):
      return angle-2*math.pi;
    if ( angle<0 ):
      return angle+2*math.pi;
    return angle;

class DoublePendulum:
    t: float
    myName: str
    deltaT: float
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

    def __init__(self, myName: str, state: list, deltaT=0.001):
        self.myName=myName
        self.setInitialPendulumState(state)
        self.setDeltaT(deltaT)
        print('Created pendulum %s' %(self,))

    def __str__(self):
        return 'P%s [t=%.1f tick=%0.2e sec (%.0f/sec)]: %.1fpi|%.1fpi/s - %.1fm - %.1fkg- %.1fpi|%.1fpi/s - %.1fm - %.1fkg' % (self.myName, self.t, self.deltaT, 1/self.deltaT, self.angleUpper/math.pi, self.velUpper/math.pi, self.lenUpper, self.massUpper, self.angleLower/math.pi, self.velLower/math.pi, self.lenLower, self.massLower)

    def setDeltaT(self, deltaT):
        self.deltaT = deltaT

    def setInitialPendulumState(self, state: list):
        self.t=0
        self.initialAngleUpper = self.angleUpper = state[0]
        self.initialAngleLower = self.angleLower = state[1]
        self.initialVelUpper = self.velUpper = state[2]
        self.initialVelLower = self.velLower = state[3]
        self.lenUpper = 1
        self.lenLower = 1
        self.massUpper = 1
        self.massLower = 1
        self.g = 9.81

    def resetToInitialState(self):
        self.angleUpper = self.initialAngleUpper
        self.angleLower = self.initialAngleLower
        self.velUpper = self.initialVelUpper
        self.velLower = self.initialVelLower
        self.t=0

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
        newAngleUpper = self.tickAngleUpper()
        newAngleLower = self.tickAngleLower()
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


def generateInitialState(range):
    return [random.random() * math.pi * 2, random.random() * math.pi * 2, random.random() * range, random.random() * range]

def process_pendulum(p, stopT):
    p.doAllTicks(stopT)
    return p

def process_pendulum(p, stopT):
    p.doAllTicks(stopT)
    return p


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
            p.resetToInitialState();

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


pendulums = []
cpus=1
deltaT_round0=0.01
rounds=10
pendulum_count=100
simulation_duration=10
animation_start_x=275
animation_spacing_x=220
show_animation_data=True
seed=tt.time()


def main(argv):
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


