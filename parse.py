"""
File created on 7/1/2018
@author: Eric Brown

This script parses the cero output.
"""

import scipy as sp
import json

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
    tags : dict of dicts
        Each packet has attached tags from GNURadio. They are saved into the
    run_num : int
        This is for bookeeping. You record the run # then can compare results across runs
    time_sets : ndarray
        This is an array of door open and door close times, so that we can delineate data from one position from another and remove
        data from the set that is impacted by RF signals outside the room.
    """
    def __init__(self, data, run_num, time_sets):
        self.data = data
        self.run_num = run_num
        self.time_sets = time_sets

    def get_tag(self, key):
        """
        Return data from particular key. As is, valid keys are: packet_num, rx_time, ofdm_sync_carrier_offset,
        ofdm_sync_chan_taps, packet_len

        :param key:
        :return: tag dict
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
    data =
    json_file = open('fixed_data.json')
    json_file.seek(0, 0)
    json_file.write('[\n' + data + '\n]')
    json_file.


def crew_read(filename = 'fixed_data.json'):
    '''

    :param filename:
    :return:
    '''
    with open(filename) as json_file:
        data = json.load(json_file)
        for p =

    return Position(tags, run_num, time_sets)
