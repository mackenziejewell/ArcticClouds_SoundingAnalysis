import numpy as np
import pandas as pd
from datetime import datetime

def grab_sounding_dates(csv = './sounding_checks/checked_20000101Z00_20001231Z23_hourly6.csv', hours = [0,6,12,18]):
    
    """Open homemade csv files listing whether soundings exist at certain dates.
    
INPUT: 
- csv: name of csv file to open
- hours: list of hours (Z) to check for sounding data (default: [0,6,12,18])

OUTPUT:
- doy: julian day, with hour as decimal (numpy array)
- sou: numpy array corresponding to doy, with 1 == sounding exists and 0 == no sounding found
- year: year of file (int)

DEPENDENCIES:
import numpy as np
import pandas as pd
from datetime import datetime

Latest recorded update:
03-01-2023
    """
    
    # open csv file corresponding to year
    ds = pd.read_csv(csv)
    
    doy = np.array([])
    sou = np.array([])
    
    # pull out dates from soundings (as day of year)
    for dd, date in enumerate(ds.date):
    
        # convert to datetime
        date_asdatetime = pd.to_datetime(date)
        
        # grab year of file
        if dd == 0:
            year = int(date_asdatetime.strftime('%Y'))

        # find julian day, with hour as decimal
        DOY = int(date_asdatetime.strftime('%j')) + int(date_asdatetime.strftime('%H'))/24

        # save specified hours of day
        if (DOY%1)*24 in hours:
            doy = np.append(doy, DOY)
            sou = np.append(sou, ds.sounding_exist[dd])
            
    return doy, sou, year