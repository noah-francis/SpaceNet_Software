import numpy
from numpy import array, cross, dot, sqrt
from datetime import datetime, timezone
from numpy.linalg import norm

# read STK file (rough probably)
def read_STK(filename):
    '''Read and parse an STK ECEF satellite position file (assumes data directory)
    ----------------------------------------------------------------------------
    Inputs:
    --------
    `filename` -- string with the file name
    
    Returns:
    --------
    `data` -- dictionary containing the STK file info
    '''

    data_dir = '/Users/rockm/Documents/CU_Boulder/ASEN_4018/Software/'
    
    # create dictionary
    data = {
        'Time': [],
        'r_vec': [],
        'r_mag': [],
    }
    
    # month dictionary
    month = {
    'Jan': 1,
    'Feb': 2,
    'Mar': 3,
    'Apr': 4,
    'May': 5,
    'Jun': 6,
    'Jul': 7,
    'Aug': 8,
    'Sep': 9,
    'Oct': 10,
    'Nov': 11,
    'Dec': 12
    }
    
    # parse data
    with open(data_dir + filename, 'r') as file:
        lines = list(file.readlines())
    for i, line in enumerate(lines):
        if (i >= 7):
            Line = line.split()

            # Time stuff
            year = int(Line[2])
            day = int(Line[0])
            hour = int(Line[3][0:2])
            minute = int(Line[3][3:5])
            # can add seconds here once we start caring about them

            # Space stuff
            x = float(Line[4]) 
            y = float(Line[5])
            z = float(Line[6])
            r_vec = array([x, y, z])
            r_mag = float(Line[7]) # magnitude of position vector

            # append Time and Space to the data dictionary's lists
            data['Time'].append(datetime(year, month[Line[1]], day, hour, minute, tzinfo=timezone.utc))
            data['r_vec'].append(r_vec)
            data['r_mag'].append(r_mag)
    data['r_vec'] = array(data['r_vec'])
    
    return data

# Gibb's Method 
def Gibby(r_vec, r_mag, mu):
    '''Use Gibb's Method to produce a velocity vector from three consecutive
    position vectors.
    ----------------------------------------------------------------------------
    Inputs:
    --------
    `r_vec` -- numpy array of three consecutive (in time) position vectors (3,3) [km]
    `r_mag` -- list of corresponding vector magnitudes [km]
    `mu`    -- Standard Gravitational Parameter of the central body [km**3 / s**2]
    
    Returns:
    --------
    `v_vec`        -- velocity vector (numpy array) at the second position vector's 
                      time (1,3) [km/s]
    `coplanar_val` -- dot product of unit vectors that should be coplanar for Gibb's
                      Method to work properly [~] (closer to zero is better)
    '''
    
    # compute unit vectors
    C_01 = cross(r_vec[0], r_vec[1]) / norm(cross(r_vec[0], r_vec[1]))
    C_12 = cross(r_vec[1], r_vec[2]) / norm(cross(r_vec[1], r_vec[2]))
    C_20 = cross(r_vec[2], r_vec[0]) / norm(cross(r_vec[2], r_vec[0]))
    
    # check how coplanar the position vectors are
    coplanar_val = dot((r_vec[0] / r_mag[0]), C_12)
    
    # calculate N, D, and S
    N = r_mag[0] * cross(r_vec[1], r_vec[2]) + r_mag[1] * cross(r_vec[2], r_vec[0]) + r_mag[2] * cross(r_vec[0], r_vec[1])
    D = cross(r_vec[0], r_vec[1]) + cross(r_vec[1], r_vec[2]) + cross(r_vec[2], r_vec[0])
    S = r_vec[0] * (r_mag[1] - r_mag[2]) + r_vec[1] * (r_mag[2] - r_mag[0]) + r_vec[2] * (r_mag[0] - r_mag[1])
    
    # calculate velocity vector at second position vector
    v_vec = sqrt(mu / (norm(N) * norm(D))) * ((cross(D, r_vec[1]) / r_mag[1]) + S)
    
    return v_vec, coplanar_val