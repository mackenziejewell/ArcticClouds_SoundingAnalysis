from datetime import datetime
import pandas as pd
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir

def check_for_WUA_soundings(date_list = [], station = 'PABR', suppress_prints = True):
    
    """Use siphon Wyoming Upper Air Data Request to check whether sounding data is 
    available from University of Wyoming at given station for given list of dates. 
    Returns pandas data frame containing dates where sounding data does/does not exist.
    
    ref:
    https://unidata.github.io/siphon/latest/examples/upperair/Wyoming_Request.html
    

INPUT: 
- date_list: list of datetime objects
- station: string code for desired station (string, default: 'PABR' Point Barrow)
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- df: pandas data frame containing dates where sounding data does/does not exist
    sounding_exist column: 1 indicates data does exist, 0 indicates no data found

DEPENDENCIES:
from datetime import datetime
import pandas as pd
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir

Latest recorded update:
02-28-2022
    """
    
    assert type(station) == str, f'station should be string, not {type(station)}'
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    
    
    sounding_exist = []

    # run through dates and check whether data exists
    for date in date_list:
        
        try:
            # request data using siphon
            sounding = WyomingUpperAir.request_data(date, station)
            sounding_exist.append(1)

        except:
            sounding_exist.append(0)

    
    # create pandas data frame containing dates where sounding data does/does not exist
    df = pd.DataFrame(list(zip(date_list, sounding_exist)),
               columns =['date', 'sounding_exist'])
    
    if suppress_prints == False:
        print(f' >>> ({np.sum(sounding_exist)}/{len(date_list)}) provided dates had sounding data')
    
    return df