
import numpy as np
from metpy.units import units
import matplotlib.pyplot as plt
from metpy.calc import (mixing_ratio_from_specific_humidity,
                       specific_humidity_from_dewpoint)
    

    
def interp_across_bins(bin_centers = [], binned_var = [], suppress_prints = True):
    
    """Linearly interpolate across empty bins (containing nans) of data for given data array.

INPUT: 
- bin_centers: evenly spaced heights (pint quantities with units) corresponding to how data will be binned
- binned_var: binned variable, with nans where no data were found
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- new_binned_var: updated binned variable with nans replaced by interpolated values

DEPENDENCIES:
import numpy as np
from metpy.units import units

Latest recorded update:
03-24-2023
    """
    
    assert len(binned_var) > 0, f'binned_var should be non-empty'
    assert len(bin_centers) > 0, f'bin_centers should be non-empty'

    assert str(type(binned_var)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'binned_var should be pint quantity, not {type(binned_var)}'
    assert str(type(bin_centers)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'bin_centers should be pint quantity, not {type(bin_centers)}'
    
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    
    # copy binned var to create new version where nans are 
    # replace with interpolated values
    new_binned_var = binned_var

    # run through each bin
    for i, bin_i in enumerate(bin_centers):

        # identify empty bin
        if np.isnan(binned_var[i]):

            current_height = bin_centers[i]

            if suppress_prints == False:
                print(f'bin {i}: {bin_centers[i]}, {binned_var[i]:.1f}')

            # find nearest non-nan data
            #--------------------------
            # search for nearest lower-height bin with data
            for iii in range(i)[::-1]:

                # identify non-empty bin and grab its height/data
                if np.isnan(binned_var[iii]) == False:
                    low_height = bin_centers[iii]
                    low_value  = binned_var[iii]
                    if suppress_prints == False:
                        print(f' - lower height: {low_height}, {low_value:.1f}')
                    break

            # search for nearest higher-height bin with data
            for iii in range(i, len(bin_centers)):

                # identify non-empty bin and grab its height/data
                if np.isnan(binned_var[iii]) == False:
                    high_height = bin_centers[iii]
                    high_value  = binned_var[iii]
                    if suppress_prints == False:
                        print(f' - higher height: {high_height}, {high_value:.1f}')
                    break

            # Interpolate
            #------------
            # change in variable over non-nan height
            Dvar = high_value - low_value
            # non-nan height change
            Dh = high_height - low_height
            # height change to nan value
            dh = current_height - low_height
            # linearly inerpolate variable to new height
            interp_var = low_value + (Dvar/Dh)*(dh)

            if suppress_prints == False:
                print(f'   > dvar/dh = {Dvar/Dh:.1f} over {dh}')
                print(f'   > var ~= {interp_var:.1f}')
                print('')

            # replace nan with interpolated value
            new_binned_var[i] = interp_var

    return new_binned_var   
    
    
    
    






def interpolate_soundings(variable = [], heights = [], bin_width = 0.25*units.kilometer, 
                          min_height = 0*units.kilometer, max_height = 12*units.kilometer, 
                          method = 'mean', suppress_plots = False, suppress_prints = True):
    
    """Interpolate sounding data with consistent height steps for easier intercomparison.
    
    ref:
    https://unidata.github.io/siphon/latest/examples/upperair/Wyoming_Request.html
    

INPUT: 
- variable: sounding data variable to interpolate (metpy variable, i.e. with units) 
    either single variable or list of multiple variables
- heights: sounding data heights corresponding to data of variable to interpolate (metpy variable, i.e. with units) 
- bin_width: size of vertical bin steps to use (pint quantity with units, default: 0.25*units.kilometer)
- min_height: lowest bin edge (pint quantity with units, default: 0*units.kilometer)
- max_height: highest bin edge (pint quantity with units, default: 12*units.kilometer)
- method: method used to describe multiple data points falling within single data bin (string)
        'mean' --> take average of variable across all data points found within given bin
        'max' --> take maximum value of variable across all data points found within given bin
- suppress_plots: bool, whether or not to supress plot that shows binning data (default: True)
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- new_binned_var: binned/interpolated data at evenly spaced heights corresponding to bin_centers
    either single variable or list of multiple variables
- bin_centers: new evenly spaced heights corresponding to new binned/interpolated data


DEPENDENCIES:
import numpy as np
from metpy.units import units
import matplotlib.pyplot as plt

# homemade functions:
create_regular_h
bin_data
add_extra_bins
interp_across_bins

Latest recorded update:
03-24-2023
    """
    
    
    assert type(method) == str, f'method should be string, not {type(method)}'
    assert method in ['mean', 'max'], f"method should be one of 'mean', 'max', not {method}"
    assert type(suppress_plots) == bool, f'suppress_plots should be bool, not {type(suppress_plots)}'
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    assert str(type(variable)) in ["<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", "<class 'list'>"], f'variable should be pint quantity or list of variables that are pint quantity, not {type(variable)}'
    assert str(type(heights)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'heights should be pint quantity, not {type(heights)}'
    assert str(type(bin_width)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'bin_width should be pint quantity, not {type(bin_width)}'
    assert str(type(min_height)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'min_height should be pint quantity, not {type(min_height)}'
    assert str(type(max_height)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'max_height should be pint quantity, not {type(max_height)}'
   
    # check whether iterating over single variable or multiple variables
    if str(type(variable)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>":
        variable = [variable]
        if suppress_prints == False:
            print('\n >>> single variable to interpolate\n')

    elif type(variable) == list:
        if suppress_prints == False:
            print('\n >>> multiple variables to interpolate\n')
    
    
    # create bins
    #------------
    bin_centers, bin_edges = create_regular_h(bin_width = bin_width, min_height = min_height, max_height = max_height)

    # create empty array to store lists of binned variables
    binned_var = [None] * len(variable)
    num_data_inbin = [None] * len(variable)
    
    # for each variable, 
    # bin data, and set to nan where no data found
    #----------------------------------------------
    for vv, var in enumerate(variable):
        binned_var[vv], num_data_inbin[vv] = bin_data(var = var, heights = heights, bin_centers = bin_centers, bin_edges = bin_edges, method = method)
        
    # make sure that num_data_inbin is same across all variables
    # otherwise something strange may be going on
    if len(variable) > 1:
        for vv in range(len(variable)-1):
            assert (num_data_inbin[vv]==num_data_inbin[vv+1]).all(), f'num_data_inbin does not match across variables {vv} and {vv+1}'
        
        
    # if last bin is empty, temporarily keep adding points on until a 
    # non-empty bin is found. This higher height will be used to interpolate
    # across desired max_height, then we will crop to desired height later
    #------------------------------------------------------------
    for vv, var in enumerate(variable):
        
        # if last bin is empty, add extra bins
        if num_data_inbin[vv][-1] == 0:
            
            # only update bin_centers, bin_edges if last variable to run through
            if vv+1 == len(variable):
                out_data = add_extra_bins(var=var, heights=heights,bin_width=bin_width, bin_centers=bin_centers, bin_edges=bin_edges, binned_var=binned_var[vv], num_data_inbin=num_data_inbin[vv], output = ['binned_var', 'num_data_inbin', 'bin_centers', 'bin_edges'], method='mean', try_extra=100, suppress_prints=True)
                
                binned_var[vv], num_data_inbin[vv], bin_centers, bin_edges = out_data
            
            # don't update bin_centers, bin_edges if not last variable to run through
            else:
                out_data = add_extra_bins(var=var, heights=heights,bin_width=bin_width, bin_centers=bin_centers, bin_edges=bin_edges, binned_var=binned_var[vv], num_data_inbin=num_data_inbin[vv], output = ['binned_var', 'num_data_inbin'], method='mean', try_extra=100, suppress_prints=True)
                
                binned_var[vv], num_data_inbin[vv] = out_data

    # check that first and last bin are not nans    
    vv = 0
    assert np.isnan(binned_var[vv].magnitude[0]) == False, f'missing data in lowest bin: {bin_edges[0]} - {bin_edges[1]}'
    assert np.isnan(binned_var[vv].magnitude[-1]) == False, f'missing data in highest bin: {bin_edges[-2]} - {bin_edges[-1]}'
    
   
    # plot bins
    #----------
    if suppress_plots == False:
        
        vv = 0
        fig, ax = plt.subplots(figsize=(2, 7)) 
        ax.scatter(np.zeros_like(bin_edges), bin_edges, marker='_', c='k', label='edge')
        ax.scatter(np.zeros_like(bin_centers), bin_centers, marker='.', c='k', label='center')
        ax.set_xticks([]);
        ax.set_xlim(-0.1,2)
        ax.set_ylabel(f'height ({bin_width.units})')
        for cc, center in enumerate(bin_centers):
            if int(num_data_inbin[vv][cc])>0:
                ax.scatter(np.arange(num_data_inbin[vv][cc])*0.1, np.full(int(num_data_inbin[vv][cc]), center), 
                           marker='o',edgecolor='k', facecolor='orange', label='points')
            else:
                ax.scatter(0.2, center, marker='x', c='blue', label='no points')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(2.1, 1));
                
        plt.show()
        

    # interpolate across empty bins
    #------------------------------
    # now interpolate data in bins without data
    # linearly interpolate with height across nearest heights with data
    
    # duplicate list to start filling in nans
    new_binned_var = [None] * len(variable)

    if suppress_prints == False:
        print(f'\n-----------------------------\ninterpolate across empty bins\n-----------------------------')
        print('>>> Search for empty bins to interpolate data with height\n')

    for vv, var in enumerate(variable):
        new_binned_var[vv] = interp_across_bins(bin_centers=bin_centers, binned_var=binned_var[vv], suppress_prints = True)
          
            
    # crop bins
    #----------
    # if we had to add on some extra bins to interplate
    # re-crop to desired max_height
    if bin_edges[-1] != max_height:
        
        # create empty array to store lists of binned variables
        binned_var_2 = [None] * len(variable)
        bin_centers_2 = np.array([])
        bin_edges_2 = np.array([])
        
        # crop edges, centers, and variables if edge exceeds max_height
        for edge in bin_edges:
            if edge <= max_height:
                bin_edges_2 = np.append(bin_edges_2, edge)
        for cc, center in enumerate(bin_centers):
            if center < max_height:
                for vv, var in enumerate(variable):
                    if cc == 0:
                        binned_var_2[vv] = np.array([])
                    binned_var_2[vv] = np.append(binned_var_2[vv], new_binned_var[vv][cc])
                bin_centers_2 = np.append(bin_centers_2, center)
             
    else:
        binned_var_2 = new_binned_var
        bin_edges_2 = bin_edges
        bin_centers_2 = bin_centers
    
    return binned_var_2, bin_centers_2, bin_edges_2






def add_extra_bins(var = [], heights = [], bin_width = [], bin_centers = [], bin_edges = [], 
                   binned_var = [], num_data_inbin = [], output = ['binned_var', 'num_data_inbin', 'bin_centers', 'bin_edges'],
                   method = 'mean', try_extra = 100, suppress_prints = True):
    
    """Add extra bins above binned sounding data until non-empty bin is found. Continue to sort data stored at irregular heights into bins with regular height intervals. This is to support later interpolation across highest desired height interval when previously highest bin was empty.

INPUT: 
- var: data variable to be binned by height
- heights: height data (pint quantities with units) corresponding to var
- bin_width: size of vertical bin steps to use (pint quantity with units, default: 0.25*units.kilometer)
- bin_centers: evenly spaced heights (pint quantities with units) corresponding to how data will be binned
- bin_edges: edges of bins of evenly spaced heights (pint quantities with units)
- binned_var: binned variable, with nans where no data were found
- num_data_inbin: array listing number of data points that were found in each bin
- output: list of variables to output, in order. strings used to indicate name of each variable.
    default: ['binned_var', 'num_data_inbin', 'bin_centers', 'bin_edges']
    - 'binned_var': updated binned variable array, with nans where no data were found and non-nan highest value
    - 'num_data_inbin': updated array listing number of data points that were found in each bin
    - 'bin_centers': updated evenly spaced heights corresponding to binned_var
    - 'bin_edges': updated  edges of bins with evenly spaced heights corresponding to binned_var
    
- method: method used to describe multiple data points falling within single data bin (string)
        'mean' --> take average of variable across all data points found within given bin
        'max' --> take maximum value of variable across all data points found within given bin
- try_extra: int, maximum number of extra bins to try to find non-nan value
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- output_data, list of variables as specified in 'output' variable

DEPENDENCIES:
import numpy as np
from metpy.units import units

Latest recorded update:
03-24-2023
    """
    
    assert len(var) > 0, f'var should be non-empty'
    assert len(heights) > 0, f'heights should be non-empty'
    assert len(bin_centers) > 0, f'bin_centers should be non-empty'
    assert len(bin_edges) > 0, f'bin_edges should be non-empty'
    
    assert str(type(bin_width)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'bin_width should be pint quantity, not {type(bin_width)}'
    assert str(type(heights)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'heights should be pint quantity, not {type(heights)}'
    assert str(type(bin_centers)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'bin_centers should be pint quantity, not {type(bin_centers)}'
    assert str(type(bin_edges)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'bin_edges should be pint quantity, not {type(bin_edges)}'
    
    assert type(method) == str, f'method should be string, not {type(method)}'
    assert method in ['mean', 'max'], f"method should be one of 'mean', 'max', not {method}"
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    
    if suppress_prints == False:
        print(f'\n--------------------------\ntake {method} of data in bins\n--------------------------')
        
    # run through up to 'try_extra' bins higher, 
    # but break loop if non-empty bin is found
    for j in range(try_extra):

        # add on extra bin with desired height interval
        bin_edges = np.append(bin_edges, bin_edges[-1] + bin_width)
        bin_centers = np.append(bin_centers, bin_centers[-1] + bin_width)

        # find edges of highest bin
        h_min = bin_edges[-2]
        h_max = bin_edges[-1]
        
        # find indices of data points within bin
        h_index = (heights>h_min) & (heights<=h_max)
        
        # variable to determine when to break loop
        # break when non-nan value found in highest bin
        break_loop = False

        # find data of points in bin
        var_in_bin = var[h_index]
        
        ## within each bin...
        #-------------------
        # if no data fall into bin, set var to nan
        # and keep running operation
        if len(var_in_bin)==0:
            mean_var = np.nan
            num_data_inbin = np.append(num_data_inbin, 0)
            break_loop = False

        # if single data values falls into bin, set var to value
        # and break loop
        elif len(var_in_bin)==1:
            mean_var = var_in_bin
            num_data_inbin=  np.append(num_data_inbin, 1)
            break_loop = True

        # if multiple values found in bin, use method to 
        # determine how to encapsulate data in single value
        # and break loop
        else:
            if str(method) == 'mean':
                mean_var = np.nanmean(var_in_bin)
            elif str(method) == 'max':
                mean_var = np.nanmax(var_in_bin)
            num_data_inbin = np.append(num_data_inbin, len(var_in_bin))
            break_loop = True

        # save to binned variable array
        binned_var = np.append(binned_var, mean_var)

        # break loop if data were found in appended bin
        if break_loop == True:
            break
    
    if suppress_prints == False:
        print(f' >>> added {j+1} bins until non-nan value was found')
    
    
    assert len(bin_centers)  == len(bin_edges)-1, f'mismatch in bin_centers length ({len(bin_centers)}) and bin_edges length ({len(bin_edges)})'
    assert len(binned_var)  == len(bin_centers), f'mismatch in binned_var (length {len(binned_var)}) and bin_centers (length {len(bin_centers)})'
    assert len(num_data_inbin)  == len(bin_centers), f'mismatch in num_data_inbin (length {len(num_data_inbin)}) and bin_centers (length {len(bin_centers)})'
    
    # save variables to dictionary
    all_data = {}
    all_data['binned_var'] = binned_var
    all_data['num_data_inbin'] = num_data_inbin
    all_data['bin_centers'] = bin_centers
    all_data['bin_edges'] = bin_edges
    
    
    # save only desired variables from all_data, in order, to output
    output_data = [all_data[desired] for desired in output]
    
    return output_data




def bin_data(var = [], heights = [], bin_centers = [], bin_edges = [], 
             method = 'mean', suppress_prints = True):
    
    """Sort data stored at irregular heights into bins with regular height intervals.

INPUT: 
- var: data variable to be binned by height
- heights: height data (pint quantities with units) corresponding to var
- bin_centers: evenly spaced heights (pint quantities with units) corresponding to how data will be binned
- bin_edges: edges of bins of evenly spaced heights (pint quantities with units)
- method: method used to describe multiple data points falling within single data bin (string)
        'mean' --> take average of variable across all data points found within given bin
        'max' --> take maximum value of variable across all data points found within given bin
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- binned_var: array of binned variable, with nans where no data were found
- num_data_inbin: array listing number of data points that were found in each bin

DEPENDENCIES:
import numpy as np
from metpy.units import units

Latest recorded update:
03-24-2023
    """
    
    assert len(var) > 0, f'var should be non-empty'
    assert len(heights) > 0, f'heights should be non-empty'
    assert len(bin_centers) > 0, f'bin_centers should be non-empty'
    assert len(bin_edges) > 0, f'bin_edges should be non-empty'
    
    assert type(method) == str, f'method should be string, not {type(method)}'
    assert method in ['mean', 'max'], f"method should be one of 'mean', 'max', not {method}"
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    
    if suppress_prints == False:
        print(f'\n--------------------------\ntake {method} of data in bins\n--------------------------')
        
    # create empty array to store lists of binned variables
    binned_var = np.array([])
    num_data_inbin = np.array([])
    
    # bin/consolidate data, fill empty bins with nan
    #-----------------------------------------------
    # run through and average data within bins
    # that contain data
    for i, bin_i in enumerate(bin_centers):

        # find edges of bin
        h_min = bin_edges[i]
        h_max = bin_edges[i+1]

        # find indices of data points within bin
        h_index = (heights>h_min) & (heights<=h_max)

        # find data of points in bin
        var_in_bin = var[h_index]

        # within each bin...
        #-------------------
        # if no data fall into bin, set var to nan
        if len(var_in_bin)==0:
            mean_var = np.nan
            num_data_inbin = np.append(num_data_inbin, 0)
            
        # if single data values falls into bin, set var to value
        elif len(var_in_bin)==1:
            mean_var = var_in_bin
            num_data_inbin = np.append(num_data_inbin, 1)

        # if multiple values found in bin, use method to 
        # determine how to encapsulate data in single value
        else:
            if str(method) == 'mean':
                mean_var = np.nanmean(var_in_bin)
            elif str(method) == 'max':
                mean_var = np.nanmax(var_in_bin)
            num_data_inbin = np.append(num_data_inbin, len(var_in_bin))

        
        # save to binned variable array
        binned_var = np.append(binned_var, mean_var)

    return binned_var, num_data_inbin



def create_regular_h(bin_width = 0.25*units.kilometer, min_height = 0*units.kilometer, 
                     max_height = 12*units.kilometer, num_decimals = 3, suppress_prints = True):
    
    """Create regular h array given height range and interval width.

INPUT: 
- bin_width: size of vertical bin steps to use (pint quantity with units, default: 0.25*units.kilometer)
- min_height: lowest bin edge (pint quantity with units, default: 0*units.kilometer)
- max_height: highest bin edge (pint quantity with units, default: 12*units.kilometer)
- num_decimals: int, number of decimals to round height data to (default: 3)
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- bin_centers: new evenly spaced heights corresponding to new binned/interpolated data
- bin_edges: edges of bins

DEPENDENCIES:
import numpy as np
from metpy.units import units

Latest recorded update:
03-23-2023
    """

    assert str(type(bin_width)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'bin_width should be pint quantity, not {type(bin_width)}'
    assert str(type(min_height)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'min_height should be pint quantity, not {type(min_height)}'
    assert str(type(max_height)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'max_height should be pint quantity, not {type(max_height)}'
    assert type(suppress_prints) == bool, f'suppress_prints should be bool, not {type(suppress_prints)}'
    
    
    # create bins
    #------------
    # check all units match
    assert bin_width.units == min_height.units, f'units of bin_width ({bin_width.units}) and min_height ({min_height.units}) should match'
    assert max_height.units == min_height.units, f'units of max_height ({max_height.units}) and min_height ({min_height.units}) should match'
    
    # create bins
    bin_edges = np.arange(min_height.magnitude, 
                          max_height.magnitude + bin_width.magnitude, 
                          bin_width.magnitude)
    bin_centers = (bin_edges+bin_width.magnitude/2)[:-1]
    
    # round to num_decimals and add units
    bin_edges = np.round(bin_edges, num_decimals) * bin_width.units
    bin_centers = np.round(bin_centers, num_decimals) * bin_width.units
    
        
    if suppress_prints == False:
        print('\n------------\ncreate bins\n------------')
        print(f'bin_edges:\n{bin_edges}')
        print(f'\nbin_centers:\n{bin_centers}')

    return bin_centers, bin_edges



def relative_humidity_from_temperature_pressure(P, Td, T, method = 'improved_magnus_1996'):

    """
    Function to calculate relative humidity, relative to both liquid water and ice from temperature, pressure, and dewpoint. 
    Formulas used are from:
    Rogers and Yau text (1996, EBOOK ISBN: 9780080570945)
    Alduchov and Eskridge (1996, doi: 10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2) 
    Huang (2018, doi: 10.1175/JAMC-D-17-0334.1)
    
    INPUT:
    - T: temperature as metpy variable, i.e. with units. Will be converted to degC
    - Td: dewpoint as metpy variable, i.e. with units. Will be converted to degC
    - P: pressure as metpy variable, i.e. with units. Will be converted to hPa
    - method: equations to use to calculate saturation pressures from T 
        either: 'improved_magnus_1996' (default) (doi: 10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2)
            or: 'huang_2018', an improved improved magnus (doi: 10.1175/JAMC-D-17-0334.1)
            
    OUTPUT:
    - RH_liq: relative humidity with respect to liquid as metpy variable with units 'percent'
    - RH_ice: relative humidity with respect to ice as metpy variable with units 'percent'
    
    DEPENDENCIES:
    import numpy as np
    from metpy.units import units
    from metpy.calc import (mixing_ratio_from_specific_humidity,
                       specific_humidity_from_dewpoint)
    homemade function:
    saturation_vaporpressures_from_temperature

    Latest recorded update:
    03-22-2023
    """

    assert type(method) == str, f'method should be string, not {type(method)}'
    assert method in ['improved_magnus_1996', 'huang_2018'], f'method should be either "improved_magnus_1996" or "huang_2018", not {method}'
    assert str(type(T)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'T should be pint quantity (including units), not {type(T)}'
    assert str(type(Td)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'Td should be pint quantity (including units), not {type(Td)}'
    assert str(type(P)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'P should be pint quantity (including units), not {type(P)}'
    
    # convert T, Td to celsius if not already
    T_celsius = T.to('degC')
    Td_celsius = Td.to('degC')
    
    # convert P to hPa if not already
    P_hPa = P.to('hPa')
    
    # mixing ratio from dewpoint temperature, pressure
    #-------------------------------------------------

    # specific humidity from pressure, dewpoint, initially dimensionless
    q = specific_humidity_from_dewpoint(P_hPa, Td_celsius).to('g/kg')

    # mixing ratio from specific humidity, initially dimensionless
    w = mixing_ratio_from_specific_humidity(q).to('g/kg')

    # saturation mixing ratio from temperature, pressure
    #---------------------------------------------------

    # saturation vapor pressure with respect to water and ice from temperature
    es_liq, es_ice = saturation_vaporpressures_from_temperature(T_celsius, method=method)

    # saturation mixing ration from saturation pressures and pressure, initially dimensionless
    # Rogers and Yau 2.18
    ws_liq = (0.622 * (es_liq) / ( P - es_liq )).to('g/kg')
    ws_ice = (0.622 * (es_ice) / ( P - es_ice )).to('g/kg')

    # relative humidity w.r.t ice and liquid
    #---------------------------------------
    
    # Rogers and Yau 2.20
    RH_liq = (100*w/ws_liq)*units('percent')
    RH_ice = (100*w/ws_ice)*units('percent')

    return RH_liq, RH_ice


def saturation_vaporpressures_from_temperature(T, method = 'improved_magnus_1996'):
    
    """
    Function to calculate saturation vapor pressure over water and ice from temperature. Formulas used are from 
    Alduchov and Eskridge (1996, doi: 10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2) and Huang (2018, doi: 10.1175/JAMC-D-17-0334.1)
    
    INPUT:
    - T: temperature as metpy variable, i.e. with units. Will be converted to degC
    - method: equations to use to calculate saturation pressures from T 
        either: 'improved_magnus_1996' (default) (doi: 10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2)
            or: 'huang_2018', an improved improved magnus (doi: 10.1175/JAMC-D-17-0334.1)
            
    OUTPUT:
    - es_liq: saturation vapor pressure with respect to liquid as metpy variable with units 'hPa'
    - es_ice: saturation vapor pressure with respect to ice as metpy variable with units 'hPa'
    
    DEPENDENCIES:
    import numpy as np
    from metpy.units import units

    Latest recorded update:
    03-22-2023
    """
    
    assert type(method) == str, f'method should be string, not {type(method)}'
    assert method in ['improved_magnus_1996', 'huang_2018'], f'method should be either "improved_magnus_1996" or "huang_2018", not {method}'
    assert str(type(T)) == "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>", f'T should be pint quantity (including units), not {type(T)}'
    
    # convert T to celsius if not already
    T_celsius = T.to('degC').magnitude
    
    # Improved Magnus Form of saturation vapor pressures from Alduchov and Eskridge (1996)
    # from:
    # Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus Form Approximation of Saturation Vapor Pressure. 
    # J. Appl. Meteor. Climatol., 35, 601–609, https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    
    if str(method) == 'improved_magnus_1996':
        es_liq = (((610.94 * np.exp( 17.625 * T_celsius / (T_celsius + 243.04) ) ))* units('Pa')).to('hPa')
        es_ice = (((611.21 * np.exp( 22.587 * T_celsius / (T_celsius + 273.86) ) ))* units('Pa')).to('hPa')

    # Further improved saturation vapor pressures from Huang (2018)
    # from:
    # Huang, J., 2018: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice. 
    # J. Appl. Meteor. Climatol., 57, 1265–1272, https://doi.org/10.1175/JAMC-D-17-0334.1.
    elif str(method) == 'huang_2018':
        es_liq = ((np.exp( 34.494 - (4924.99/(T_celsius + 237.1)) ) / ((T_celsius + 105)**(1.57)))* units('Pa')).to('hPa')
        es_ice = ((np.exp( 43.494 - (6545.8/(T_celsius + 278)) ) / ((T_celsius + 868)**(2)))* units('Pa')).to('hPa')
    
    return es_liq, es_ice
    

    
