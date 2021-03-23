# -*- coding: utf-8 -*-
"""
Spyder Editor

This script loads a .wav file, reads the data interned, splits it into chunks based on size, and stores it as a compressed
h5 file.
"""
from scipy.io import wavfile
import tables
import numpy

# read data from wav
fs, data = wavfile.read('C:/Users/TheRe/Downloads/RTL-NoDelay-test/SDRSharp_20210210_041402Z_437200000Hz_IQ.wav')

# Get metadata of IQ data for chunking purposes
rows = numpy.size(data,0)
nochunks = 256
szchunk = rows/nochunks

# Set export location & filename
folder = 'C:/Users/TheRe/Downloads/RTL-NoDelay-test/H5/'
name = '20210210_041402Z_437200000Hz_256.'+"h5"

# Set compression level
filters = tables.Filters(complevel=5, complib='zlib')

# Save to h5 format
datah5 = tables.open_file(folder+name, mode = "w", title = name)
datah5.create_earray('/','time_data', atom=None, title='', filters=filters, \
                         expectedrows=rows, \
                         byteorder=None, createparents=False, obj=data, track_times=True)
datah5.set_node_attr('/time_data','sample_freq', fs)
datah5.close()
