
import numpy as np
from metpy.units import units
import matplotlib.pyplot as plt
from metpy.calc import (mixing_ratio_from_specific_humidity,
                       specific_humidity_from_dewpoint)
                        
                        
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
    

    
def interpolate_soundings(variable = [], heights = [], bin_width = 0.25*units.kilometer, 
                          min_height = 0*units.kilometer, max_height = 12*units.kilometer, 
                          method = 'mean',
                          suppress_plots = False,
                          suppress_prints = True):
    
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
- bin_edges: edges of bins

DEPENDENCIES:
import numpy as np
from metpy.units import units
import matplotlib.pyplot as plt

Latest recorded update:
03-14-2023
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
    # check all units match
    assert bin_width.units == min_height.units, f'units of bin_width ({bin_width.units}) and min_height ({min_height.units}) should match'
    assert max_height.units == min_height.units, f'units of max_height ({max_height.units}) and min_height ({min_height.units}) should match'
    bin_edges = np.arange(min_height.magnitude,max_height.magnitude+bin_width.magnitude,bin_width.magnitude)*bin_width.units
    bin_centers = (bin_edges.magnitude+bin_width.magnitude/2)[:-1]*bin_width.units

    if suppress_prints == False:
        print('\n------------\ncreate bins\n------------')
        print(f'bin_edges:\n{bin_edges}')
        print(f'\nbin_centers:\n{bin_centers}')

    if suppress_prints == False:
        print(f'\n--------------------------\ntake {method} of data in bins\n--------------------------')
   
    
    num_data_inbin = []
    
    # create empty array to store lists of binned variables
    binned_var = [None] * len(variable)
    
    for vv, var in enumerate(variable):
        
        binned_var[vv] = np.array([])
        
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

            # take mean of data values 
            if len(var_in_bin)==0:
                mean_var = np.nan
                if vv == 0:
                    num_data_inbin.append(0)
            elif len(var_in_bin)==1:
                mean_var = var_in_bin
                if vv == 0:
                    num_data_inbin.append(1)
            else:
                if str(method) == 'mean':
                    mean_var = np.mean(var_in_bin)
                elif str(method) == 'max':
                    mean_var = np.max(var_in_bin)
                if vv == 0:
                    num_data_inbin.append(len(var_in_bin))

            # save to binned variable array
            binned_var[vv] = np.append(binned_var[vv], mean_var)
        
    # if last bin is empty 
    # keep adding points on until a non-empty bin is found
    # this higher height will be used to interpolate across
    # desired max_height, then we will crop to desired height later
    if num_data_inbin[-1] == 0:
        
        # run through up to 50 bins higher, 
        # but break loop if non-empty bin is found
        for j in range(50):
        
            bin_edges = np.append(bin_edges, bin_edges[-1] + bin_width)
            bin_centers = np.append(bin_centers, bin_centers[-1] + bin_width)

            # find edges of bin
            h_min = bin_edges[-2]
            h_max = bin_edges[-1]
            
            # find indices of data points within bin
            h_index = (heights>h_min) & (heights<=h_max)

            break_loop = False
            
            for vv, var in enumerate(variable):
                
                # find data of points in bin
                var_in_bin = var[h_index]
                
                # take mean of data values 
                if len(var_in_bin)==0:
                    mean_var = np.nan
                    if vv == 0:
                        num_data_inbin.append(0)
                    break_loop = False
                elif len(var_in_bin)==1:
                    mean_var = var_in_bin
                    if vv == 0:
                        num_data_inbin.append(1)
                    break_loop = True
                else:
                    if str(method) == 'mean':
                        mean_var = np.mean(var_in_bin)
                    elif str(method) == 'max':
                        mean_var = np.max(var_in_bin)
                    if vv == 0:
                        num_data_inbin.append(len(var_in_bin))
                    break_loop = True
                    
                # save to binned variable array
                binned_var[vv] = np.append(binned_var[vv], mean_var)

            if break_loop == True:
                break

    # check that first and last bin are not nans    
    assert np.isnan(binned_var[vv].magnitude[0]) == False, f'missing data in lowest bin: {bin_edges[0]} - {bin_edges[1]}'
    assert np.isnan(binned_var[vv].magnitude[-1]) == False, f'missing data in highest bin: {bin_edges[-2]} - {bin_edges[-1]}'
    
    # this second assert is becoming an issue when chosen highest bin height is empty
    # so now including a check for empty high bin, then temporarily adding some extra
    # bins if this is the case to interpolate across chosen highest bin, recropping later
    
    # plot bins
    #----------
    if suppress_plots == False:
        fig, ax = plt.subplots(figsize=(2, 7)) 
        ax.scatter(np.zeros_like(bin_edges), bin_edges, marker='_', c='k', label='edge')
        ax.scatter(np.zeros_like(bin_centers), bin_centers, marker='.', c='k', label='center')
        ax.set_xticks([]);
        ax.set_xlim(-0.1,2)
        ax.set_ylabel(f'height ({bin_width.units})')
        for cc, center in enumerate(bin_centers):
            if num_data_inbin[cc]>0:
                ax.scatter(np.arange(num_data_inbin[cc])*0.1, 
                           np.full(num_data_inbin[cc], center), 
                           marker='o',edgecolor='k', facecolor='orange', label='points')
            else:
                ax.scatter(0.2, center, 
                           marker='x', c='blue', label='no points')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(2.1, 1));
                
        plt.show()
        
    

    # interpolate across empty bins
    #------------------------------
    # duplicate list to start filling in nans
    new_binned_var = binned_var

    if suppress_prints == False:
        print(f'\n-----------------------------\ninterpolate across empty bins\n-----------------------------')
    
    # now interpolate data in bins without data
    # linearly interpolate with height across nearest heights with data
    if suppress_prints == False:
        print('>>> Search for empty bins to interpolate data with height\n')
    
    # for each variable
    for vv, var in enumerate(variable):
        
        # for each bin
        for i, bin_i in enumerate(bin_centers):

            # identify empty bin
            if np.isnan(binned_var[vv][i]):

                current_height = bin_centers[i]

                if suppress_prints == False:
                    print(f'bin {i}: {bin_centers[i]}, {binned_var[vv][i]:.1f}')

                # find nearest non-nan data
                #--------------------------
                # search for nearest lower-height bin with data
                for iii in range(i)[::-1]:

                    # identify non-empty bin and grab its height/data
                    if np.isnan(binned_var[vv][iii]) == False:
                        low_height = bin_centers[iii]
                        low_value  = binned_var[vv][iii]
                        if suppress_prints == False:
                            print(f' - lower height: {low_height}, {low_value:.1f}')
                        break

                # search for nearest higher-height bin with data
                for iii in range(i, len(bin_centers)):

                    # identify non-empty bin and grab its height/data
                    if np.isnan(binned_var[vv][iii]) == False:
                        high_height = bin_centers[iii]
                        high_value  = binned_var[vv][iii]
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
                new_binned_var[vv][i] = interp_var
            
          
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
