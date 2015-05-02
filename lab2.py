"""
	Jesse Martinez
	Lab 1
"""

import array
from struct import unpack
from queue import PriorityQueue
from math import sqrt, pow, e, atan2, pi, floor, ceil
import xml.etree.ElementTree as XML

ROW_COL = 3601
latMeters = 111085.70174301133
lonMeters = 82014.96771009536

nodes = {}
graph = {}

wi = []
elevs = []
features = 2

#Returns XML tree
def readXML(fileName):
	return (XML.parse(fileName))

def gradientDescent():
	wi = [0 for i in range(features)]
	prevWi = [1 for i in range(features)]
	stepper = 1.2

	while wi != prevWi:
		dWi = [0 for i in range(features)]
		for example in examples:
			pathNodes = example[:-3]
			xs = []
			xs.append(getElevChange(pathNodes))
			xs.append(getDistChange(pathNodes))
			y = (example[-3] * 60) + example[-2]
			for i in range(features):
				dWi[i] += -1 * (y - linearRegression(xs)) * xs[i]
		for i in range(features):
			wi[i] -= stepper * dWi[i]

def linearRegression(xs):
	hx = 0
	for i in range(features):
		hx += wi[i] * xs[i]
	return (hx * (1/1.6))

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
def aStar(start, goal):
	if start not in nodes:
		return start + " node not in map."
	elif goal not in nodes:
		return goal + "node not in map."

	pq = PriorityQueue()
	pq.put((0, start))

	costs = {}
	cameFrom = {}

	costs[start] = 0
	cameFrom[start] = None

	while not pq.empty():
		cost, current = pq.get()

		if current == goal:
			minutes = int(costs[goal])
			seconds = int(100 * (costs[goal] - minutes))
			while seconds > 60:
				seconds -= 60
				minutes += 1
			print(nodes[goal])
			print("Walking Time " + str(minutes) + " minutes and " + str(seconds) + " seconds")
			return constructPath(cameFrom, goal)

		for neighbor in graph[current]:
			newCost = costs[current] + toblers(current, neighbor)

			if neighbor not in costs or newCost < costs[neighbor]:
				costs[neighbor] = newCost
				priority = newCost + heuristic(neighbor,goal)
				pq.put((priority, neighbor))
				cameFrom[neighbor] = current

	return "No path found between " + start + " and " + goal 

#xyz distance between the two nodes given converted into second using the average walking speed 1.4m/s
def heuristic(nodeFrom, nodeTo):
	xSquare = (nodes[nodeTo][2] - nodes[nodeFrom][2]) * (nodes[nodeTo][2] - nodes[nodeFrom][2])
	ySquare = (nodes[nodeTo][3] - nodes[nodeFrom][3]) * (nodes[nodeTo][3] - nodes[nodeFrom][3])
	zSquare = (nodes[nodeTo][4] - nodes[nodeFrom][4]) * (nodes[nodeTo][4] - nodes[nodeFrom][4])
	return (sqrt(xSquare + ySquare + zSquare) / 1.4) / 60

#Naismith rule used to find minutes to walk between the two nodes given
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
	return xyMins

#Tobler's Hiking Function used to find minutes to walk between the two nodes given
def toblers(nodeFrom, nodeTo):
	if nodeFrom == nodeTo:
		return 0
	#km/h
		   #lat To     #lat From
	xSquare = (nodes[nodeTo][2] - nodes[nodeFrom][2]) * (nodes[nodeTo][2] - nodes[nodeFrom][2])
	ySquare = (nodes[nodeTo][3] - nodes[nodeFrom][3]) * (nodes[nodeTo][3] - nodes[nodeFrom][3])
	zSquare = (nodes[nodeTo][4] - nodes[nodeFrom][4]) * (nodes[nodeTo][4] - nodes[nodeFrom][4])
	     #elev To    #elev From
	dh = nodes[nodeTo][4] - nodes[nodeFrom][4]
	dx = sqrt(xSquare + ySquare + zSquare)

	w = 6 * pow(e, (-3.5) * (abs((dh/dx) + 0.05)))
	w = (.06 / w) * dx
	return w

#Construct the path taken after doing an A* search
def constructPath(cameFrom, node):
	path = []
	while node != None:
		path.append(node)
		node = cameFrom[node]
	return path[::-1]

if __name__ == "__main__":
	tree = readXML("dbv.osm")
	generateNodes()
	generateGraph()
	print(nodes)	
	#while(input("Search (Y/N):") == "Y"):
	#	print(aStar(input("Input Starting Node:"), input("Input Destination Node:")))
