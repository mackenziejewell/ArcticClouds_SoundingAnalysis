from time import sleep
from tqdm import tqdm
import sys
from datetime import datetime
import pandas as pd
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir

def check_for_WUA_soundings(date_list = [], station = 'PABR', show_progress = True, num_tries = 1):
    
    """Use siphon Wyoming Upper Air Data Request to check whether sounding data is 
    available from University of Wyoming at given station for given list of dates. 
    Returns pandas data frame containing dates where sounding data does/does not exist.
    
    ref:
    https://unidata.github.io/siphon/latest/examples/upperair/Wyoming_Request.html
    

INPUT: 
- date_list: list of datetime objects
- station: string code for desired station (string, default: 'PABR' Point Barrow)
- show_progress: bool, whether or not to show progress bar (default: True)
- num_tries: int, number of tries to check for data since some fail when data exists (default 1)

OUTPUT:
- df: pandas data frame containing dates where sounding data does/does not exist
    sounding_exist column: 1 indicates data does exist, 0 indicates no data found

DEPENDENCIES:
from datetime import datetime
import pandas as pd
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir
from time import sleep
from tqdm import tqdm
import sys
    
Latest recorded update:
03-21-2023
    """
    
    assert type(station) == str, f'station should be string, not {type(station)}'
    assert type(show_progress) == bool, f'show_progress should be bool, not {type(show_progress)}'
    
    
    sounding_exist = []

    
    # run with progress bar
    if show_progress == True:
    
        # use tqdm to show progress bar if desired
        with tqdm(total=len(date_list), file=sys.stdout) as pbar:
            
            for dd, date in enumerate(date_list):
                
                # request data using siphon
                # try request multiple times in case there is error with first
                for attempt in range(num_tries):
                    try:
                        # request data using siphon
                        sounding = WyomingUpperAir.request_data(date, station)
                        sounding_exist.append(1)
                        worked = True
                        break
                    except:
                        worked = False
                if worked == False:
                    sounding_exist.append(0)
                # show progress
                pbar.set_description('processed: %d' % (1 + dd))
                pbar.update(1)
                sleep(0.01)
                
            tqdm._instances.clear()
            

    # run without progress bar
    else:
        # run through dates and check whether data exists
        for date in date_list:
            
            # request data using siphon
            # try request multiple times in case there is error with first
            for attempt in range(num_tries):
                try:
                    # request data using siphon
                    sounding = WyomingUpperAir.request_data(date, station)
                    sounding_exist.append(1)
                    worked = True
                    break
                except:
                    worked = False
            if worked == False:
                sounding_exist.append(0)

    # create pandas data frame containing dates where sounding data does/does not exist
    df = pd.DataFrame(list(zip(date_list, sounding_exist)),
               columns =['date', 'sounding_exist'])
    
    print(f'\n>>> ({np.sum(sounding_exist)}/{len(date_list)}) provided dates had sounding data')
    
    return df