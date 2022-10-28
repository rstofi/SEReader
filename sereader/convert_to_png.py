"""A script to read .ser [TIA (Emispec) file format] microscopy images

For the ITA or (.SER) format see: https://www.ntu.edu.sg/home/cbb/info/TIAformat/index.html

Author: Kristof Rozgonyi
e-mail: rstofi@gmail.com
License: MIT
Date: 2020
"""
#===============
#=== Imports ===
#===============
import sys;
import numpy as np;
import logging;
import struct;

#For plotting
import matplotlib
import matplotlib.pyplot as plt

#=================
#=== Functions ===
#=================
#=== Header ===
def check_series_ID_and_Version(byte_array):
    """Read the first 6 bytes in and
    Check the validity of the
    series ID and series Version
    
    I decided to work in a hex format and
    the input will be always the byte_array
    as I thinking of creating a SER class
    for the convesrion in the future
    
    :return: series_Version as it will be needed later
    """
    #Asser little endian encoding
    assert hex(int.from_bytes(byte_array[0:2], 'little')) == '0x4949', "Invalid byte ordering: {0:s}".format(
            hex(int.from_bytes(byte_array[0:2], 'little')));

    #Valid ID's and Versions
    valid_series_ID_list = ['0x197'];
    valid_series_Version_list = ['0x210','0x220'];

    #Chek ID and version validity
    series_ID = hex(int.from_bytes(byte_array[2:4], 'little'));
    series_Version = hex(int.from_bytes(byte_array[4:6], 'little'));

    log.info("Series ID: {0:s}".format(series_ID));
    log.info("Series Version: {0:s}".format(series_Version)); 

    assert series_ID in valid_series_ID_list, "Invalid series ID: {0:s}".format(series_ID);
    assert series_Version in valid_series_Version_list, "Invalid series ID: {0:s}".format(series_Version);

    return series_Version;

def get_data_and_tag_ID(byte_array):
    """Get the data type ID
    and the tag type ID

    return: 
    """
    #Valid Data and Tag type ID's
    valid_data_type_ID_list = ['0x1420','0x1422'];
    valid_tag_type_ID_list = ['0x4152','0x4142'];

    #Data type
    data_type_ID = hex(int.from_bytes(byte_array[6:10], 'little'));
    if data_type_ID == '0x4120':
        log.info("Tag Type ID: {0:s} : 1D array".format(data_type_ID));
    elif data_type_ID == '0x4122':
        log.info("Tag Type ID: {0:s} : 2D array".format(data_type_ID));
    else:
        assert data_type_ID in valid_data_type_ID_list, "Invalid data type ID: {0:s}".format(data_type_ID);

    #Tag type
    tag_type_ID = hex(int.from_bytes(byte_array[10:14], 'little'));
    if tag_type_ID == '0x4152':
        log.info("Tag Type ID: {0:s} : time tag only".format(tag_type_ID));
    elif tag_type_ID == '0x4142':
        log.info("Tag Type ID: {0:s} : time and 2D position tag".format(tag_type_ID));
    else:
        assert tag_type_ID in valid_tag_type_ID_list, "Invalid tag type ID: {0:s}".format(tag_type_ID);

    return data_type_ID, tag_type_ID


def get_total_and_valid_element_num(byte_array):
    """Get the number of total elements
    and valid elements

    return:
    """

    #Number of total elements
    N_total_element = int.from_bytes(byte_array[14:18], 'little');
    log.info("Total number of elements: {0:d}".format(N_total_element));

    N_valid_element = int.from_bytes(byte_array[18:22], 'little');
    log.info("Valid number of elements: {0:d}".format(N_valid_element));

    assert N_total_element >= N_valid_element, "Total number of elements ({0:d}) is lower than valid elements ({1:d})!".format(
            N_total_element,N_valid_element);

    #Check if all data is written
    if N_total_element != N_valid_element:
        log.info("Not all elements of the series were written!");

def get_Offset_Array_Offset(byte_array, SeriesVersion):
    """Get the Offset Array Offset

    Indicates the absolute offset (in bytes) in the Series
    Data File of the beginning of the Data Offset Array.

    4 bytes if series version <= 0x210
    8 bytes if series version >= 0x220

    :param: SeriesVersion the output from check_series_ID_and_Version
    """

    #Get the series version

    #Due the fact that only two valid series version exist in the code
    #and in the documentation so far, we only use these cases

    if SeriesVersion == "0x210":
        Offset_Array_Offset_length = int(4);
    elif SeriesVersion == "0x220":
        Offset_Array_Offset_length = int(8);

    log.info("The Offset Array Offset byte length is: {0:d}".format(Offset_Array_Offset_length));

    #Get the actual value of the offset
    Offset_Array_Offset = int.from_bytes(byte_array[22:22+Offset_Array_Offset_length], 'little');

    log.info("The Offset Array Offset is: {0:d}".format(Offset_Array_Offset));

    #Returen the Offset_Array_Offset byte length and its value
    return Offset_Array_Offset_length, Offset_Array_Offset;

def get_Number_of_Dimensions(byte_array,OffsetArrayOffset_length):
    """Get the Number of dimensions

    Indicates the number of dimensions of the Series Data File.
    This indicates the number of dimensions of the indices, NOT the number of dimensions of the data.
    
    :param: OffsetArrayOffset_length is the byte lenghth of the OffsetArrayOffset
    """
    byte_offset = 22 + OffsetArrayOffset_length;

    Number_of_Dimensions = int.from_bytes(byte_array[byte_offset:byte_offset+4], 'little');

    log.info("Number of Dimensions: {0:d}".format(Number_of_Dimensions));

    #Return the number of Dimension array dimension
    return Number_of_Dimensions;

#=== Dimension array ===
def get_Dimension_Size(byte_array,OffsetArrayOffset_length,dim_byte_offset=0):
    """Get the dimension size e.i. the number of elements


    :param: OffsetArrayOffset_length is the byte lenghth of the OffsetArrayOffset
    :param dim_byte_offset: The byte offset for the Nth element in the dimension array

    return:
    """
    byte_offset = 26 + OffsetArrayOffset_length + dim_byte_offset;

    N_dimension_size = int.from_bytes(byte_array[byte_offset:byte_offset + 4], 'little');
    log.info("Number of elements in this Dimension: {0:d}".format(N_dimension_size));

def get_Calibration_Element(byte_array,OffsetArrayOffset_length,dim_byte_offset=0):
    """Get the Calibration Offset, the Calibration Delta and the Calibration element


    :params: see above
    :return: 
    """
    byte_offset = 30 + OffsetArrayOffset_length + dim_byte_offset;


    Calibration_Offset = float(struct.unpack('<d',byte_array[byte_offset:byte_offset + 8])[0]);
    Calibration_Delta = float(struct.unpack('<d',byte_array[byte_offset + 8:byte_offset + 16])[0]);
    Calibration_Element = int.from_bytes(byte_array[byte_offset + 16:byte_offset + 20], 'little');

    log.info("Calibration Element index: {0:d} with Offset {1:.2f} and Delta {2:.2f}".format(Calibration_Element,Calibration_Offset,Calibration_Delta));

    return Calibration_Element,Calibration_Offset,Calibration_Delta;


def get_Description_Length(byte_array,OffsetArrayOffset_length,dim_byte_offset=0):
    """Get the length of the dimension element descripction

    :params: see above
    :return:
    """
    byte_offset = 50 + OffsetArrayOffset_length + dim_byte_offset;

    Element_Description_Length = int.from_bytes(byte_array[byte_offset:byte_offset + 4], 'little');

    log.info("The length of the element descripction string is {0:d}".format(Element_Description_Length));

    return Element_Description_Length;

def get_Element_Description(byte_array,OffsetArrayOffset_length,element_description_length,dim_byte_offset=0):
    """Get the description of the element

    :params: see above
    :param element_description_length: the output of the get_description_length function
    :return:
    """
    byte_offset = 54 + OffsetArrayOffset_length + dim_byte_offset;

    Element_Descrition = str(byte_array[byte_offset:byte_offset + element_description_length], 'utf-8');

    log.info("Element descrition:\n{0:s}".format(Element_Descrition));

    return Element_Descrition;

def get_Units_Length(byte_array,OffsetArrayOffset_length,element_description_length,dim_byte_offset=0):
    """Get the length of the Units string of the element

    :params: see above
    :return:
    """
    byte_offset = 54 + OffsetArrayOffset_length + element_description_length + dim_byte_offset;

    print(byte_offset)

    Element_Units_Length = int.from_bytes(byte_array[byte_offset:byte_offset + 4], 'little');

    log.info("The length of the element units string is {0:d}".format(Element_Units_Length));

    return Element_Units_Length;

def get_Units_Description(byte_array,OffsetArrayOffset_length,element_description_length,element_units_length,dim_byte_offset=0):
    """Get the units description string for the element

    :params: see above
    :param element_units_length: output of the get_Units_Length function
    :return:
    """
    if element_units_length == 0:
        log.info("No associated unit string for this element!");

        return None;
    else:
        byte_offset = 58 + OffsetArrayOffset_length + element_description_length + dim_byte_offset;

        Element_Units = str(byte_array[byte_offset:byte_offset + element_units_length], 'utf-8');

        log.info("Element units descrition:\n{0:s}".format(Element_Units));

        return Element_Units;

def get_Data_Offset_Array(byte_array,OffsetArrayOffset,SeriesVersion,TotalNumberOfElements):
    """Return a list of the byte offest for each individual elements in
    the TotalNumber of elements array

    :param byte_array :
    :param OffsetArrayOffset: the second output of the get_Offset_Array_Offset() function
    :param SeriesVersion: the series version
    :param TotalNumberOfElements: equals to the dimensions

    :return:
    """

    DataOffsetArray = np.zeros(TotalNumberOfElements,dtype=int)

    #Check series version
    if SeriesVersion == "0x210":
        OffsetArrayOffset_length = int(4)
    elif SeriesVersion == "0x220":
        OffsetArrayOffset_length = int(8)

    byte_offset = OffsetArrayOffset_length;

    #Loop through all elements (dimensions)
    for NElement in range(0,TotalNumberOfElements):

        DataOffsetArray[NElement] = int.from_bytes(byte_array[OffsetArrayOffset:OffsetArrayOffset+byte_offset], 'little');

        byte_offset += OffsetArrayOffset_length;

    return DataOffsetArray;

def get_Tag_Offset_Array(byte_array,OffsetArrayOffset,SeriesVersion,TotalNumberOfElements,DataOffsetArray_length):
    """Return a list of the byte offest for each individual elements in
    the TotalNumber of elements array

    :param byte_array :
    :param OffsetArrayOffset: the second output of the get_Offset_Array_Offset() function
    :param SeriesVersion: the series version
    :param TotalNumberOfElements: equals to the dimensions
    :param DataOffsetArray_length: the lenght of the array outputted from get_Data_Offset_Array()

    :return:
    """
    TagOffsetArray = np.zeros(TotalNumberOfElements,dtype=int)

    #Check series version
    if SeriesVersion == "0x210":
        OffsetArrayOffset_length = int(4)
    elif SeriesVersion == "0x220":
        OffsetArrayOffset_length = int(8)

    byte_offset_start = OffsetArrayOffset + (DataOffsetArray_length * OffsetArrayOffset_length)
    byte_offset =  OffsetArrayOffset_length

    #print(OffsetArrayOffset,byte_offset_start,byte_offset_start+byte_offset)

    #Loop through all elements (dimensions)
    for NElement in range(0,TotalNumberOfElements):

        TagOffsetArray[NElement] = int.from_bytes(byte_array[byte_offset_start:byte_offset_start+byte_offset], 'little');

        byte_offset += OffsetArrayOffset_length;

    return TagOffsetArray;

def get_TimeOnlyTag(byte_array,TagOffset,TagTypeID):
    """Return the time stamp of the data

    :param byte_array:
    :param TagOffset: the byte offset where the ElementTag starts
    :param TagTypeID: the Tag Type ID gained from the header
    
    :return:
    """

    tag_type_ID = hex(int.from_bytes(byte_array[TagOffset:TagOffset + 4], 'little'))#Two rouge byte at the end, but does nor matter

    assert tag_type_ID == TagTypeID, "Invalid Tag Type ID: {0:s} (not {1:s})!".format(tag_type_ID,TagTypeID)

    time_stamp = int.from_bytes(byte_array[TagOffset + 4:TagOffset + 8], 'little');

    return time_stamp;

def get_Time_and_PositionTag(byte_array,TagOffset,TagTypeID):
    """Return the time stamp and the position tags of the data

    :param byte_array:
    :param TagOffset: the byte offset where the ElementTag starts
    :param TagTypeID: the Tag Type ID gained from the header has 
    
    :return:
    """

    tag_type_ID = hex(int.from_bytes(byte_array[TagOffset:TagOffset + 4], 'little'))#Two rouge byte at the end, but does nor matter

    assert tag_type_ID == TagTypeID, "Invalid Tag Type ID: {0:s} (not {1:s})!".format(tag_type_ID,TagTypeID)

    time_stamp = int.from_bytes(byte_array[TagOffset + 4:TagOffset + 8], 'little');

    position_X = float(struct.unpack('<d',byte_array[TagOffset + 8:TagOffset + 16])[0]);
    position_Y = float(struct.unpack('<d',byte_array[TagOffset + 16:TagOffset + 16])[0]);

    return time_stamp,position_X,position_Y;

def get_1DdataElementHeader(byte_array,DataOffset,DataTypeID,TagTypeID):
    """Get the 1D data from the data element array for a single element
    Returns the calibration offset, the calibration delta, the calibration element, data type, the data array size
    and the data array itself 
    
    :param byte_array:
    :param DataOffset: the byte offset where the DataElement starts
    :param DataTypeID: the Data Type ID gained from the header
    :param TagTypeID: the Tag Type ID gained from the header
    
    :return:
    """

    #=== This seems to be unnecessary ===
    #Get the full offset based on the tag type that is befor the data element
    #Okay here it should be the tag array...
    if TagTypeID == "0x4152":
        #byte_offset = DataOffset + 8
        byte_offset = DataOffset + 0
    elif TagTypeID == "0x4142":
        #byte_offset = DataOffset + 24
        byte_offset = DataOffset + 0
    #====================================

    Calibration_Offset = float(struct.unpack('<d',byte_array[byte_offset:byte_offset + 8])[0]);
    Calibration_Delta = float(struct.unpack('<d',byte_array[byte_offset + 8:byte_offset + 16])[0]);
    Calibration_Element = int.from_bytes(byte_array[byte_offset + 16:byte_offset + 20], 'little');

    DataType = int(struct.unpack('<H',byte_array[byte_offset + 20:byte_offset + 22])[0]);

    Data_Array_Size = int.from_bytes(byte_array[byte_offset + 22:byte_offset + 26], 'little');

    return Calibration_Offset, Calibration_Delta, Calibration_Element, DataType, Data_Array_Size

def get_2DdataElementData(byte_array,DataOffset,DataTypeID,DataType,DataSize):
    """Get the actual data of the Element

    :param byte_array:
    :param DataOffset: the byte offset where the DataElement starts
    :param DataTypeID: the Data Type ID gained from the header
    :param DataType: theencoding type of the data as described in the documentation
    :param DataSize: the lenght of the data in the only axis
    
    :return:
    """

    byte_offset = DataOffset + 26

    #Get the data length => need a different function to generate this
    if DataType == 2:
        Data_Array_byte_Size = DataSizeX * DataSizeY * 2;

    #Read in as a single tuple (!) but I convert it to a numpy array with np.asarray()
    Image_Array = np.asarray(struct.unpack("<{0:d}H".format(DataSize),byte_array[byte_offset:byte_offset + Data_Array_byte_Size]))

    return Image_Array

def get_2DdataElementHeader(byte_array,DataOffset,DataTypeID,TagTypeID):
    """Get the 2D data from the data element array for a single element
    Returns the calibration offset, the calibration delta, the calibration element, data type, the data array size
    and the lenght of the data array in bytes but not the data array itself 
    
    :param byte_array:
    :param DataOffset: the byte offset where the DataElement starts
    :param DataTypeID: the Data Type ID gained from the header
    :param TagTypeID: the Tag Type ID gained from the header
    
    :return:
    """

    #=== This seems to be unnecessary ===
    #Get the full offset based on the tag type that is befor the data element
    #Okay here it should be the tag array...
    if TagTypeID == "0x4152":
        #byte_offset = DataOffset + 8
        byte_offset = DataOffset + 0
    elif TagTypeID == "0x4142":
        #byte_offset = DataOffset + 24
        byte_offset = DataOffset + 0
    #====================================

    Calibration_Offset_X = float(struct.unpack('<d',byte_array[byte_offset:byte_offset + 8])[0]);
    Calibration_Delta_X = float(struct.unpack('<d',byte_array[byte_offset + 8:byte_offset + 16])[0]);
    Calibration_Element_X = int.from_bytes(byte_array[byte_offset + 16:byte_offset + 20], 'little');

    Calibration_Offset_Y = float(struct.unpack('<d',byte_array[byte_offset + 20:byte_offset + 28])[0]);
    Calibration_Delta_Y = float(struct.unpack('<d',byte_array[byte_offset + 28:byte_offset + 36])[0]);
    Calibration_Element_Y = int.from_bytes(byte_array[byte_offset + 36:byte_offset + 40], 'little');

    #The data offset and delta should be the same in the xy direction
    assert Calibration_Offset_X == Calibration_Offset_Y, "The calibration offset is different ({0:e},{1:e}) for the x and y directions!".format(
            Calibration_Offset_X,Calibration_Offset_Y)

    assert Calibration_Delta_X == Calibration_Delta_Y, "The calibration delta is different ({0:e},{1:e}) for the x and y directions!".format(
            Calibration_Delta_X,Calibration_Delta_Y)

    DataType = int(struct.unpack('<H',byte_array[byte_offset + 40:byte_offset + 42])[0]);

    Data_Array_Size_X = int.from_bytes(byte_array[byte_offset + 42:byte_offset + 46], 'little');
    Data_Array_Size_Y = int.from_bytes(byte_array[byte_offset + 46:byte_offset + 50], 'little');

    return Calibration_Offset_X, Calibration_Offset_Y, \
            Calibration_Delta_X, Calibration_Delta_Y, \
            Calibration_Element_X,Calibration_Element_Y, \
            DataType, \
            Data_Array_Size_X,Data_Array_Size_Y

def get_2DdataElementData(byte_array,DataOffset,DataTypeID,DataType,DataSizeX,DataSizeY):
    """Get the actual data of the Element

    :param byte_array:
    :param DataOffset: the byte offset where the DataElement starts
    :param DataTypeID: the Data Type ID gained from the header
    :param DataType: theencoding type of the data as described in the documentation
    :param DataSizeX: the lenght of the data in the X axis (N pixel)
    :param DataSizeY: the lenght of the data in the Y axis (N pixel)
    
    :return:
    """

    byte_offset = DataOffset + 50

    #Get the data length => need a different function to generate this
    if DataType == 2:
        Data_Array_byte_Size = DataSizeX * DataSizeY * 2;

    #Read in as a single tuple (!) but I convert it to a numpy array with np.asarray()
    Image_Array = np.asarray(struct.unpack("<{0:d}H".format(DataSizeX * DataSizeY),byte_array[byte_offset:byte_offset + Data_Array_byte_Size]))

    return Image_Array

def save_2DdataElemwntAsImage(ImageArray,DataSizeX,DataSizeY,OffsetValue,DeltaValue,ElementIndexX,ElementIndexY,OutputName):
    """Save a 2D image element with some default settings

    :param ImageArray: the 1D (!) Image array that is the output of the get_2DdataElementData()
    :param DataSizeX: the lenght of the data in the X axis (N pixel)
    :param DataSizeY: the lenght of the data in the Y axis (N pixel)
    :param OffsetValue: the offset of the reference pixel
    :param DeltaValue: the offset of the incremet
    :param ElementIndexX: the x index of the reference pixel
    :param ElementIndexY: the y index of the reference pixel
    :param OutputName: the name of the output file 

    :return :
    """
    #Reshape it into an image
    ImageArray = ImageArray.reshape(DataSizeX,DataSizeY);

    #Cange the image type to float
    ImageArray = ImageArray.astype(float)

    #Sacle with the offset
    Calibration_Element_Value = ImageArray[ElementIndexX,ElementIndexY];
    #ImageArray[ElementIndexX,ElementIndexY] = OffsetValue #I assume this

    #Scale everything with the offset
    for i in range(0,DataSizeX):
        for j in range(0,DataSizeX):
            ImageArray[i,j] = OffsetValue + ((ImageArray[i,j] - Calibration_Element_Value) * DeltaValue)

    #Plot
    #plt.matshow(Image_format,cmap='gray',origin='lower') 

    #Save the array onlya s an image
    #plt.imsave(OutputName,ImageArray,cmap='viridis',origin='lower')
    plt.imsave(OutputName,ImageArray,cmap='gray',origin='lower')
    
    #This can also do eg .eps !

#===================================================
#            __  __          _____ _   _ 
#           |  \/  |   /\   |_   _| \ | |
#           | \  / |  /  \    | | |  \| |
#           | |\/| | / /\ \   | | | . ` |
#           | |  | |/ ____ \ _| |_| |\  |
#           |_|  |_/_/    \_\_____|_| \_|
#
#===================================================
if __name__ == '__main__':
    """Create combined image plots
    """
    #Setup logging
    log = logging.getLogger(__name__);
    log.setLevel(logging.INFO);
    log.addHandler(logging.StreamHandler(sys.stdout));

    #=== File to rad in ===
    #input_file = './14.41.47 Scanning Acquire_1.ser';
    input_file = './14.23.42 Scanning Acquire_1.ser'

    #=== Read in the binary file as strings ===
    with open(input_file, mode='rb') as file: # b is important -> binary
        image_byte_array = file.read()

    file.close();

    #=== Work with the header ===
    log.info("Reafing in .SER file: {0:s}".format(input_file));

    log.info("Reading the header in...");
    series_version = check_series_ID_and_Version(image_byte_array);
    data_type_ID, tag_type_ID = get_data_and_tag_ID(image_byte_array);
    get_total_and_valid_element_num(image_byte_array);
    offset_array_offset_length, offset_array_offset = get_Offset_Array_Offset(image_byte_array, series_version);
    N_Dimensions = get_Number_of_Dimensions(image_byte_array,offset_array_offset_length);
    log.info("...done.");

    #=== Work with the dimension array ===
    log.info("Reading the Dimension Array format...");
    DM_offset = 0;#The dimension array offset for the Nth element
    
    calibration_element_list = np.zeros(N_Dimensions,dtype=int);
    calibration_offset_list = np.zeros(N_Dimensions,dtype=float);
    calibration_delta_list = np.zeros(N_Dimensions,dtype=float);

    description_lenghth_list = np.zeros(N_Dimensions,dtype=int);
    element_description_list = [];

    units_length_list = np.zeros(N_Dimensions,dtype=int);
    units_description_list = [];

    for ND in range(0,N_Dimensions):
        log.info("Dimension no {0:d}...".format(ND));
        get_Dimension_Size(image_byte_array,offset_array_offset_length,DM_offset);
        calibration_element_list[ND], calibration_offset_list[ND], calibration_delta_list[ND] = \
        get_Calibration_Element(image_byte_array,offset_array_offset_length,DM_offset);
        description_lenghth_list[ND] = get_Description_Length(image_byte_array,offset_array_offset_length,DM_offset);
        element_description_list.append(get_Element_Description(image_byte_array,offset_array_offset_length,description_lenghth_list[ND],DM_offset));
        units_length_list[ND] = get_Units_Length(image_byte_array,offset_array_offset_length,description_lenghth_list[ND],DM_offset);
        units_description_list.append(get_Units_Description(image_byte_array,offset_array_offset_length,description_lenghth_list[ND],units_length_list[ND],DM_offset));

        DM_offset += description_lenghth_list[ND] + units_length_list[ND];
        log.info("...done.")

    dimensions_array_end_byte_offset = 26 + offset_array_offset_length + ((32 * N_Dimensions) + DM_offset);#26 + offset + the length of the dimensions array including all elements

    assert dimensions_array_end_byte_offset  == offset_array_offset, \
    "OffsetArrayOffset ({0:d} bytes) has not beer reached. The dimension array read till {1:d} bytes!".format(offset_array_offset,dimensions_array_end_byte_offset);

    #Log the Element lists
    log.info("For each dimensions the respective element, offset, delta, units and description lists:\n{}\n{}\n{}\n{}\n{}".format(
            calibration_element_list,calibration_offset_list,calibration_delta_list,units_description_list,element_description_list));

    log.info("...done.");

    #=== Work with the Data offset array ===
    log.info("Reading the Data Offset Array format...");

    #Create a list of the byte offset for the individual data elements
    Data_Offset_Array = get_Data_Offset_Array(image_byte_array,offset_array_offset,series_version,N_Dimensions)

    log.info("The Data Offset Array:\n{}".format(Data_Offset_Array))

    log.info("...done")

    #=== Work with the Tag offset array ===
    log.info("Reading the Tag Offset Array format...");

    #Create a list of the byte offset for the individual tage elements
    Tag_Offset_Array = get_Tag_Offset_Array(image_byte_array,offset_array_offset,series_version,N_Dimensions,np.size(Data_Offset_Array))

    log.info("The Tag Offset Array:\n{}".format(Data_Offset_Array))

    log.info("...done")

    #=== Working with the Data elements and Data Tags ===

    #Initialise variable lists nevertheless there are big redundancies here!

    #Tag
    time_stamp_list = np.zeros(N_Dimensions,dtype=int)
    if tag_type_ID == "0x4152":
        position_tag_list = None
    elif tag_type_ID == "0x4142":
        position_tag_list = np.zeros((N_Dimensions,2),dtype=float)

    #Data
    if data_type_ID == "0x4120":
        data_calibration_offset_list = np.zeros(N_Dimensions,dtype=float)
        data_calibration_delta_list = np.zeros(N_Dimensions,dtype=float)
        data_calibration_element_list = np.zeros(N_Dimensions,dtype=int)
        data_type_list = np.zeros(N_Dimensions,dtype=int)#Non redundant info and only 2 bytes in the original data set => numpy increases the size
        data_array_lenght_list = np.zeros(N_Dimensions,dtype=int)
        data_array_list = [] #store it in a list format as each element can have a different size
    elif data_type_ID == "0x4122":
        data_calibration_offset_list = np.zeros((N_Dimensions,2),dtype=float)
        data_calibration_delta_list = np.zeros((N_Dimensions,2),dtype=float)
        data_calibration_element_list = np.zeros((N_Dimensions,2),dtype=int)
        data_type_list = np.zeros(N_Dimensions,dtype=int)#Non redundant info and only 2 bytes in the original data set => numpy increases the size
        data_array_lenght_list = np.zeros((N_Dimensions,2),dtype=int)
        data_array_list = [] #store it in a list format as each element can have a different size

    log.info("Reading in every element:")
    for ND in range(0,N_Dimensions):
        log.info("\tElement no. {0:d}".format(ND))

        #Read Data Tag Array 
        if tag_type_ID == "0x4152":
            time_stamp_list[ND] = get_TimeOnlyTag(image_byte_array,Tag_Offset_Array[ND],tag_type_ID)
            log.info("\tANSI-standard (UNIX) time stamp: {0:d}".format(time_stamp_list[ND]))

        elif tag_type_ID == "0x4142":
            time_stamp_list[ND], position_tag_list[ND,0], position_tag_list[ND,1] = \
            get_Time_and_PositionTag(image_byte_array,Tag_Offset_Array[ND],tag_type_ID)
            log.info("\tANSI-standard (UNIX) time stamp: {0:d}".format(time_stamp_list[ND]))

        a = get_2DdataElementHeader(image_byte_array,Data_Offset_Array[ND],data_type_ID,tag_type_ID)

        #Read Data Element Array
        if data_type_ID == "0x4120":
            data_calibration_offset_list[ND], data_calibration_delta_list[ND], data_calibration_element_list[ND], data_type_list[ND], data_array_lenght_list[ND] = get_1DdataElementHeader(image_byte_array,Data_Offset_Array[ND],data_type_ID,tag_type_ID)

            log.info("\tElement calibration element at ({0:d},{1:d}) xy coordinates with {2:e} offest and +\- {3:e} delta".format(
                    data_calibration_element_list[ND], data_calibration_element_list[ND],data_calibration_offset_list[ND],data_calibration_delta_list[ND]))
            log.info("\tThe data type is {0:d}".format(data_type_list[ND]))
            log.info("\tThe element size is {0:d}x{1:d} pixels".format(data_array_lenght_list[ND],data_array_lenght_list[ND]))

            data_array_list.append(get_1DdataElementData(image_byte_array,Data_Offset_Array[ND],data_type_ID,data_type_list[ND],data_array_lenght_list[ND]))

        elif data_type_ID == "0x4122":
            data_calibration_offset_list[ND,0], data_calibration_offset_list[ND,1], data_calibration_delta_list[ND,0], data_calibration_delta_list[ND,1], data_calibration_element_list[ND,0], data_calibration_element_list[ND,1], data_type_list[ND], data_array_lenght_list[ND,0], data_array_lenght_list[ND,1] = get_2DdataElementHeader(image_byte_array,Data_Offset_Array[ND],data_type_ID,tag_type_ID)

            log.info("\tElement calibration element at ({0:d},{1:d}) xy coordinates with {2:e} offest and +\- {3:e} delta".format(
                    data_calibration_element_list[ND,0], data_calibration_element_list[ND,1],data_calibration_offset_list[ND,0],data_calibration_delta_list[ND,0]))
            log.info("\tThe data type is {0:d}".format(data_type_list[ND]))
            log.info("\tThe element size is {0:d}x{1:d} pixels".format(data_array_lenght_list[ND,0],data_array_lenght_list[ND,1]))

            data_array_list.append(get_2DdataElementData(image_byte_array,Data_Offset_Array[ND],data_type_ID,data_type_list[ND],data_array_lenght_list[ND,0],data_array_lenght_list[ND,1]))

            log.info("\tSave image...")
            save_2DdataElemwntAsImage(data_array_list[ND],data_array_lenght_list[ND,0],data_array_lenght_list[ND,1],data_calibration_offset_list[ND,0],data_calibration_delta_list[ND,0],data_calibration_element_list[ND,0], data_calibration_element_list[ND,1],'./test_output.png')
            log.info("\t...done")

    log.info("...done")

    #NOTES

    #Tag seems to be in front and that the online documentation is up to date
    #Data nd Tag type ID's are 4bytes instead of 2!
