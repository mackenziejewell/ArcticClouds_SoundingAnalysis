
import numpy as np
from metpy.units import units

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

Latest recorded update:
02-28-2023
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
                num_data_inbin.append(0)
            elif len(var_in_bin)==1:
                mean_var = var_in_bin
                num_data_inbin.append(1)
            else:
                if str(method) == 'mean':
                    mean_var = np.mean(var_in_bin)
                elif str(method) == 'max':
                    mean_var = np.max(var_in_bin)
                num_data_inbin.append(len(var_in_bin))

            # save to binned variable array
            binned_var[vv] = np.append(binned_var[vv], mean_var)
        

        # check that first and last bin are not nans    
        assert np.isnan(binned_var[vv].magnitude[0]) == False, f'missing data in lowest bin: {bin_edges[0]} - {bin_edges[1]}'
        assert np.isnan(binned_var[vv].magnitude[-1]) == False, f'missing data in highest bin: {bin_edges[-2]} - {bin_edges[-1]}'
    
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
            
            
    return new_binned_var, bin_centers, bin_edges
