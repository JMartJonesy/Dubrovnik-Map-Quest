"""
	Jesse Martinez
	Lab 1
"""

import array
from struct import unpack
from queue import PriorityQueue
from numpy import transpose, dot
from numpy.linalg import inv
from math import sqrt, pow, e, atan2, pi, floor, ceil
import xml.etree.ElementTree as XML

ROW_COL = 3601
latMeters = 111085.70174301133
lonMeters = 82014.96771009536

tree = None
tests = 0
features = 2

nodes = {}
graph = {}

wi = []
elevs = []
examples = []

#Returns XML tree
def readXML(fileName):
	return (XML.parse(fileName))

#Reads in walking data
def readWalks(fileName):
	global tests
	for line in open(fileName):
		example = line.strip().split(",")
		examples.append((example[:-4], [(float(example[-4]) * 60) + float(example[-3])] + [example[-1]]))
	tests = ceil(.6 * len(examples))

#Find example(s) with closest elevation and distances to the given path
#	if a single example is found then returns that examples time
#	else calculates the average of the two found examples
def nearestNeighbor(path):
	xs = [getElevChange(path), getDistChange(path)]
	minElevDiff = float("inf")
	minDistDiff = float("inf")
	minExp = None
	for example in examples[:tests]:
		exElevChange = getElevChange(example[0])
		exDistChange = getDistChange(example[0])
		if (abs(exElevChange - xs[0]) < minElevDiff) and (abs(exDistChange - xs[1]) < minDistDiff):
			minElevDiff = abs(exElevChange - xs[0])
			minDistDiff = abs(exDistChange - xs[1])
			minExp = example
	
	return minExp[1][0]

#Using linear algebra to find wi values (F^T F)^-1 F^T T
def linearAlgebraSolver():
	global wi, tests
	t = []
	f = []
	for example in examples[:tests]:
		t.append(example[1][0])
		f.append([getElevChange(example[0]), getDistChange(example[0])])
	
	wi = dot(dot(inv(dot(transpose(f),f)),transpose(f)), t)

#Runs a gradient descent to find the wi values that minimize the square error
def gradientDescent():
	global wi, tests
	#wi = [0 for i in range(features)]
	error = float("inf")
	prevError = 0
	stepper = .0001

	while abs(error - prevError) > 30:
		dWi = [0 for i in range(features)]
		for example in examples[:tests]:
			xs = [getElevChange(example[0]), getDistChange(example[0])]
			y = example[1][0]
			prevError = error
			error = 0
			for i in range(features):
				diff = y - linearRegression(example[0])
				error += diff * diff
				dWi[i] += -1 * diff * xs[i]
		for i in range(features):
			wi[i] -= stepper * dWi[i]
		print(error, prevError)
		input()

#Linear equation for predicting the time to walk a path
def linearRegression(path):
	xs = [getElevChange(path), getDistChange(path)]
	hx = 0
	for i in range(features):
		hx += wi[i] * xs[i]
	return hx

#Sums the elevation change of each node and neighbor node in the path
def getElevChange(pathNodes):
	dE = 0
	for i in range(len(pathNodes) - 1):
		dE += nodes[pathNodes[i+1]][-1] - nodes[pathNodes[i]][-1]
	return dE

#Sums the XY distance of each node and neighbor node in the path
def getDistChange(pathNodes):
	dD = 0
	for i in range(len(pathNodes) - 1):
		dD += sqrt(pow(nodes[pathNodes[i+1]][2] - nodes[pathNodes[i]][2], 2) + pow(nodes[pathNodes[i+1]][3] - nodes[pathNodes[i]][3], 2))
	return dD

#Uses interpolation to find the elevation of a given lat and lon location
def smoothDatElevation(lat, lon):
	lat, lon = getLatLonSeconds(lat, lon)

	tLeft = (floor(lat), ceil(lon), getElevation(floor(lat), ceil(lon)))
	tRight = (ceil(lat), ceil(lon), getElevation(ceil(lat), ceil(lon)))
	bLeft = (floor(lat), floor(lon), getElevation(floor(lat), floor(lon)))
	bRight = (ceil(lat), floor(lon), getElevation(ceil(lat), floor(lon)))

	tDiff = lat - tLeft[0]
	tMid = (tLeft[0] + tDiff, tLeft[1], ((1 - tDiff) * tLeft[2]) + (tDiff * tRight[2]))

	bDiff = lat - bLeft[0]
	bMid = (bLeft[0] + bDiff, bLeft[1], ((1 - bDiff) * bLeft[2]) + (bDiff * bRight[2]))

	mDiff = lon - bMid[1]
	return (((1 - mDiff) * bMid[2]) * (bDiff * tMid[2]))

#Convert lat and lon to seconds
def getLatLonSeconds(lat, lon):
	latSec = 3600 - (lat - int(lat)) * 3600
	lonSec = (lon - int(lon)) * 3600
	return (latSec, lonSec)

#Retrieves elevation given lat and lon in seconds
def getElevation(latSec, lonSec):
	return elevs[(latSec * ROW_COL) + lonSec]

#Reads in elevation file and generates a dictionary of node-key and tuple-value pairs
#The tuple representation (latitude, longitude, latitude in meters, longitude in meters, elevation
def generateNodes():
	root = tree.getroot()
	getElevs("N42E018.hgt")
	for node in root.findall("node"):
		lat = float(node.get("lat"))
		lon = float(node.get("lon"))
		nodes[node.get("id")] = (lat, lon, (lat * latMeters), (lon * lonMeters), smoothDatElevation(lat, lon))

#Generates the graph of ways represented as a dictionary with node-key to a list of neighbors-value
def generateGraph():
	root = tree.getroot()
	for way in root.findall("way"):
		prev = None 
		for node in way.findall("nd"):
			ref = node.get("ref")
			if ref not in graph:
				graph[ref] = []
			if prev != None:
				graph[ref].append(prev)
				graph[prev].append(ref)
			prev = ref

#Reads in the elevations file and creates a list of elevations
def getElevs(fileName):
	global elevs
	elevs = []
	with open(fileName, 'rb') as fila:
		elevs = array.array('h')
		elevs.fromfile(fila, ROW_COL * ROW_COL)
	elevs.byteswap()

#A* search from start to goal
def aStar(start, goal, linear = True):
	if start not in nodes:
		return start + " node not in map."
	elif goal not in nodes:
		return goal + "node not in map."

	global wi
	if linear:
		costFn = linearRegression
		if len(wi) == 0:
			linearAlgebraSolver()
	else:
		costFn = nearestNeighbor

	pq = PriorityQueue()
	pq.put((0, start))

	costs = {}
	cameFrom = {}

	costs[start] = 0
	cameFrom[start] = None

	while not pq.empty():
		cost, current = pq.get()

		if current == goal:
			minutes = 0
			seconds = costs[goal]
			while seconds > 60:
				seconds -= 60
				minutes += 1
			print("Walking Time " + str(minutes) + " minutes and " + str(seconds) + " seconds")
			print("PATH:")
			return constructPath(cameFrom, goal)

		for neighbor in graph[current]:
			newCost = costs[current] + costFn([current, neighbor])

			if neighbor not in costs or newCost < costs[neighbor]:
				costs[neighbor] = newCost
				priority = newCost + costFn([neighbor,goal])
				pq.put((priority, neighbor))
				cameFrom[neighbor] = current

	return "No path found between " + start + " and " + goal 

#xyz distance between the two nodes given converted into second using the average walking speed 1.4m/s
def heuristic(nodeFrom, nodeTo):
	xSquare = (nodes[nodeTo][2] - nodes[nodeFrom][2]) * (nodes[nodeTo][2] - nodes[nodeFrom][2])
	ySquare = (nodes[nodeTo][3] - nodes[nodeFrom][3]) * (nodes[nodeTo][3] - nodes[nodeFrom][3])
	zSquare = (nodes[nodeTo][4] - nodes[nodeFrom][4]) * (nodes[nodeTo][4] - nodes[nodeFrom][4])
	return (sqrt(xSquare + ySquare + zSquare) / 1.4) / 60

#Naismith rule used to find seconds to walk between the two nodes given
def naiSmith(nodeFrom, nodeTo):
	if nodeFrom == nodeTo:
		return 0
	
	xSquare = (nodes[nodeTo][2] - nodes[nodeFrom][2]) * (nodes[nodeTo][2] - nodes[nodeFrom][2])
	ySquare = (nodes[nodeTo][3] - nodes[nodeFrom][3]) * (nodes[nodeTo][3] - nodes[nodeFrom][3])
	xyMins = (sqrt(xSquare + ySquare) / 5000) * 60
	elevChange = nodes[nodeTo][4] - nodes[nodeFrom][4]
	slope = atan2(elevChange, sqrt(xSquare + ySquare)) * (180 / pi)
	if elevChange < 0 and slope < -12 and slope > -5:
		xyMins -= (elevChange / 300) * 10
	elif elevChange < 0 and slope > -12:
		xyMins += (elevChange / 300) * 10
	elif elevChange > 0:
		xyMins += (elevChange / 600) * 60
	return (xyMins * 60)

#Tobler's Hiking Function used to find seconds to walk between the two nodes given
def toblers(nodeFrom, nodeTo):
	if nodeFrom == nodeTo:
		return 0
	xSquare = (nodes[nodeTo][2] - nodes[nodeFrom][2]) * (nodes[nodeTo][2] - nodes[nodeFrom][2])
	ySquare = (nodes[nodeTo][3] - nodes[nodeFrom][3]) * (nodes[nodeTo][3] - nodes[nodeFrom][3])
	zSquare = (nodes[nodeTo][4] - nodes[nodeFrom][4]) * (nodes[nodeTo][4] - nodes[nodeFrom][4])
	dh = nodes[nodeTo][4] - nodes[nodeFrom][4]
	dx = sqrt(xSquare + ySquare + zSquare)

	w = 6 * pow(e, (-3.5) * (abs((dh/dx) + 0.05)))
	w = (.06 / w) * dx
	return (w*60)

#Construct the path taken after doing an A* search
def constructPath(cameFrom, node):
	path = []
	while node != None:
		path.append(node)
		node = cameFrom[node]
	return path[::-1]

#Initialize all them goodies
def initItUp():
	global tree
	tree = readXML("dbv.osm")
	readWalks("walks.txt")
	generateNodes()
	generateGraph()

if __name__ == "__main__":
	initItUp()
	#linearAlgebraSolver()
	#gradientDescent()
	#print(wi)
	#for example in examples[tests:]:
	#	print(linearRegression(example[0]), example[1][0])
	
	while(input("Search (Y/N):") == "Y"):
		start = input("Input Starting Node:")
		end = input("Input Destination Node:")
		print("LINEAR REGRESSION-----------------------")
		print(aStar(start, end))
		print("NEAREST NEIGHBOR------------------------")
		print(aStar(start, end, False))
