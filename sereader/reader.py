"""The main reader, that creates the correpsonding
subclass of the given serial image. 
"""

from codecs import open
import logging
import numpy as np


class SeriesData(object):
    def __init__(self,series_version,data_type_id,tag_type_id):
        #Need to assert the following
        #byte-order is little endian
        # series ID = 0x0197

        #Make sure that the class argument values are the right type when calling the reader function

        assert series_version == '0x0220' or series_version == '0x0210', "Invalid series version: {0:s}".format(series_version)
        self.series_version = series_version

        assert data_type_id == '0x4120' or data_type_id == '0x4122', "Invalid data type ID: {0:s}".format(data_type_id)
        self.data_type_id = data_type_id

        assert tag_type_id == '' or tag_type_id == '', "Invalid tag type ID: {0:s}".format(tag_type_id)
        self.tag_type_id = tag_type_id



#Just for quick testing
if __name__ == '__main__':
    A = SeriesData('0x0220','0x4122',2)