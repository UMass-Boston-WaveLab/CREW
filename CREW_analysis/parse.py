"""
File created on 7/1/2018
@author: Eric Brown

This script parses the cero output.
"""

import json
import csv

'''
get files
parse json
pull out all tags and add them to data object
'''


class Position:
    """
    Class for storing test data.

    The structure of this class is determined by the formatting of the CERO JSON output.
    While running a test, start a timer on your phone and note the time you open and close the door to the shielded room.
    GNURadio saves a receive time tag and it will be preserved when making the object.
    When the test is complete, enter your time set into the crew_read function as a 2xN array where N is the number of sample
    positions taken, this will serve to pull only valid data out of the json. The object class also contains several analysis functions.

    The offset value in the header is about .007 seconds from the rx_time. Not sure what that tag is meant to represent.
    Hugh can you look into that for me?


    Attributes
    ----------
    data : dict of dicts
        Each packet has attached tags from GNURadio. They are saved into the
    run_num : int
        This is for bookeeping. You record the run # then can compare results across runs
    time_set : ndarray
        This is an array of door open and door close times, so that we can delineate data from one position from another and remove
        data from the set that is impacted by RF signals outside the room.
    """
    def __init__(self, data, pos_num, time_sets):
        self.data = data
        self.pos_num = pos_num
        self.time_set = time_set

    def get_data(self, data):
        """
        Returns a data cube containing the data from each packet received in specified time slice.
        :param time_slice:
        :return: data ndarray
        """





def condition_json(filename):
    """
    The CERO framework atm doesn't dump the json cleanly, so it has to be modified
    :param filename:
    :return:
    """
    with open(filename, 'r') as json_file:
        data = json_file.read()
        json_file.close()
    data = data.replace('}{', '},\n{')
    json_file = open('conditioned_data.json', 'w+')
    json_file.seek(0, 0)
    json_file.write('[\n' + data + '\n]')
    json_file.close()


def crew_read(run_num, run_date, run_notes, times_filename='time_sets.csv', json_filename='conditioned_data.json'):
    """
    While running test, create a CSV with a comma delimiter, (or code another delimiter). This method will
    create a set of position objects for the run. These objects will have attributes for date, time, notes...

    :param run_num, run_date, run_notes:
    :param times_filename:
    :param json_filename:
    :return:
    """
    with open(json_filename) as json_file:
        data = json.load(json_file)
        json_file.close()
    time_sets = []
    with open(times_filename) as time_file:
        times = csv.reader(time_file, delimiter=',')
        for row in times:
            time_sets = time_sets.append(row)
    return RunSet(data, run_num, run_date, run_notes, time_sets)
