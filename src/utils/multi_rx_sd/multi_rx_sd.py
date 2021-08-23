import sys
import utm
import time
import math
import cartopy
import numpy as np
import matplotlib.pyplot as plt

from area import area
from shapely.geometry import Point, LineString, Polygon

from socket import *

sock = socket(AF_INET, SOCK_STREAM)

'''
Steps for comparison
1. Get coordinates and velocity vectors for each receiver
2. Check if coordinates from receivers are same
3. Propagate ref coords
4. Check proximity of coordinates from each receiver with the respective propagated coordinated
5. Calculate perimeter and area of the polygon formed by the received coordinates 
'''
def get_coords_vel():
    while 1:
        try:
            sock.connect(('127.0.0.1',3333))
            print ("Connected")
            break
        except ConnectionRefusedError:
            time.sleep(5)
            print ("Retrying..")

    while 1:
        try:
            msg = '\rstatus\n'
            # print (type(msg))
            # msg = sys.stdin.readline()
            # print (type(msg))

            sock.send(msg.encode())
            data = sock.recv(1024).decode()
            data = data.split("\n")
            for row in data:
                if "Receiver Position" in row:
                    print (row.split(" ")[-3:])

        except KeyboardInterrupt:
            sock.close()
            break

    sock.close()

def propagate_coords(ref_coords, vel):
    return_coords = []

    for coords in ref_coords:
        dt = 1 # To do implement dt (use time.time()

        metersPerDegLat = 111111.0;
        metersPerRadLat = metersPerDegLat * 180 / math.pi;
        metersPerDegLon = metersPerDegLat * math.cos(coords[1]);
        metersPerRadLon = metersPerDegLon * 180 / math.pi;

        lat = coords[1] + vel[0] * dt / metersPerRadLat;
        lon = coords[0] + vel[1] * dt / metersPerRadLon;

        return_coords.append([lon, lat])

    return return_coords

def get_distance(coord1, coord2):
    R = 6378.137    
    lat1 = coord1[1]
    lon1 = coord1[0]
    lat2 = coord2[1]
    lon2 = coord2[0]

    dLat = lat2 * math.pi / 180 - lat1 * math.pi / 180
    dLon = lon2 * math.pi / 180 - lon1 * math.pi / 180
    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.cos(
        lat1 * math.pi / 180
    ) * math.cos(lat2 * math.pi / 180) * math.sin(dLon / 2) * math.sin(dLon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c
    return d * 1000

def check_proximity(coord1, coord2):
    distance = get_distance(coord1, coord2)
    return distance < 5 # proximity of 5 m (configurable)

def calculate_perimeter(coords):
    perimeter = 0
    for i in range(len(coords)): 
        perimeter = perimeter + get_distance(coords[i], coords[(i+1)%4])

    return perimeter

def main():
    ref_coords = [[-71.085777, 42.338863], [-71.083863, 42.338867], [-71.083864, 42.338050], [-71.085777, 42.338055], [-71.085777, 42.338863]]
    geo_json = {'type':'Polygon','coordinates':[]}
    geo_json_ref = geo_json
    geo_json_recv = geo_json

    geo_json_ref['coordinates'].append(ref_coords)
    ref_perimeter = calculate_perimeter(ref_coords)
    ref_area = area(geo_json_ref)

    while 1:
        # Step 1
        recv_coords, vel = get_coords_vel()
        geo_json_recv['coordinates'].append(recv_coords)

        # Step 2
        temp = []
        [temp.append(i) for i in recv_coords if i not in temp]
        if len(temp) != len(recv_coords):
            print ("Spoofing detected: 1 or more coordinates are same")

        # Step 3
        propagated_coords = propagate_coords(ref_coords[0:-1], vel)

        # Step 4
        for i in range(len(propagated_coords)):
            if not check_proximity(propagated_coords[i], recv_coords[i]):
                print ("Proximity check failed for RX# %d" %(i))
        
        # Step 5.1
        recv_perimeter = calculate_perimeter(recv_coords)
        perimeter_error_threshold = 1 # %

        if (abs(recv_perimeter - ref_perimeter) * 100) / ref_perimeter > perimeter_error_threshold:
            print ("Perimeter check failed")

        # Step 5.2
        recv_area = area(geo_json_recv)
        area_error_threshold = 1 # %
        if (abs(recv_area - ref_area) * 100) / ref_area > area_error_threshold:
            print ("Area check failed")

        ref_coords = propagated_coords

