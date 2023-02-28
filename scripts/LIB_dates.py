from datetime import datetime
from datetime import timedelta

def generate_date_list(date_i = datetime(2022, 1, 1, 0), 
                       date_f = datetime(2022, 1, 2, 0), 
                       hourly = 12, suppress_prints = True):
    
    """Generate list of dates between date_i and date_f with specified hourly steps.

INPUT: 
- date_i: initial date (datetime object)
- date_f: final date (datetime object)
- hourly: time step in hours to use between dates (int)
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- date_list: list of datetime objects

DEPENDENCIES:
from datetime import datetime
from datetime import timedelta

Latest recorded update:
02-28-2022
    """
    
    assert type(hourly) == int, f'hourly should be integer, not {type(hourly)}'
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    assert str(type(date_i)) == "<class 'datetime.datetime'>", f'date_i should be datetime object, not {type(date_i)}'
    assert str(type(date_f)) == "<class 'datetime.datetime'>", f'date_i should be datetime object, not {type(date_f)}'

    # start list
    date_list = [date_i]

    # determine number of time steps to use between dates from hourly
    hours_between_dates = (date_f - date_list[-1]).total_seconds()/3600
    time_steps = int(hours_between_dates/hourly)
    
    # add dates to list 
    for ii in range(time_steps):
        date_list.append(date_list[-1] + timedelta(hours=hourly))
        
    if suppress_prints == False:
        print(f' >>> list of lenth {len(date_list)} generated with {hourly}-hourly steps between {date_list[0]} and {date_list[-1]}')

        if date_list[-1] != date_f:
            print(f' >>> with provided values, final date {date_list[-1]} does not match provided date_f')
    
    
    return date_list