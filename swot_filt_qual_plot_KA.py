# swot_filt_qual_plots.py

import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import glob as glob
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset

import matplotlib.colors as colors
import matplotlib.cm as cmx
from itertools import product

import time

from mpl_toolkits.basemap import Basemap

def plot_swot(axn, var, lonv, latv, vmin, vmax, cmap, merv, parv, box, plotstyle, levels=0): 
    """
    Improve part of if lon, lat 1D or 2D
    Plots an individual plot of the data using scatter and Basemap.
    Input variables:
    - xx = x-axis variable, in this case, longitude
    - yy = y-axis variable, in this case, latitude
    - var = variable to be plotted
    - ax = axis
    - box = box region of the area wanted to be shown
    , where box is a 1 x 4 array: 
    [minimum_longitude maximum_longitude minimum_latitude maximum_latitude]
    - vmin = minimum value of the colorbar 
    - vmax = maximum value of the colorbar
    - merv = a 1 x 4 array to know if to label and where the meridians
    -- [0 0 0 0] = no labels
    -- [1 0 0 1]
    - parv = like merv, but for the parallels' labels
    - cmap = colormap to be used
    - plotstyle:
    -- '0' only scatter
    -- '1' scatter with contour
    - levels: levels of contour.  
    -- 0 if none (or just scatter)
    -- levels = np.arange(min, max, interval)
    Output variables:
    - c1 = plot object (for colorbar)
    - my_mapn = map obejct (to do more plots on same figure)
    """
    
    lomin = box[0]
    lomax = box[1]
    lamin = box[2]
    lamax = box[3]
    
    my_mapn = Basemap(projection='merc'
                      , lat_0=(lamin+lamax)/2
                      , lon_0=(lomin+lomax)/2
                      , resolution = 'l'
                      , llcrnrlon = lomin, llcrnrlat= lamin
                      , urcrnrlon = lomax
                      , urcrnrlat = lamax, ax = axn
                      , area_thresh = 1000)
    
    x, y = my_mapn(lonv, latv) # compute map proj coordinates.
    
    if plotstyle == '1':
        my_mapn.contour(x, y, var, levels, colors='w', linewidth=2) 
    c1 = my_mapn.scatter(x, y, c=var, s=8, linewidth=0
                         , vmin=vmin, vmax=vmax, cmap=cmap)
    
    #my_mapn.drawcoastlines() 
    #my_mapn.drawmapboundary()
    #my_mapn.fillcontinents(color='.75') 
    
    if box == box_p09:
        my_mapn.drawmeridians(np.arange(-160.5, 140.5, 1), labels=merv, size=10)
        my_mapn.drawparallels(np.arange(0., 70., 1), labels=parv, size=10)
    elif box == box_p22:
        my_mapn.drawmeridians(np.arange(-160., 140., 1), labels=merv, size=10)
        my_mapn.drawparallels(np.arange(0.5, 70.5, 1), labels=parv, size=10)
    
    return my_mapn, c1

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

import SWOTdenoise as swotd

# Need to load the SWOTdenoise module

def derivatives_calc_method(ssh, lon, lat, x_ac, x_al, order):
    """
    Without gaussian convolution
    """
    if np.ma.isMaskedArray(ssh) == False:
        ssh = np.ma.asarray(ssh)
        #print 'ssh had to be masked'
    
    # Fill nadir gap with masked fill values
    ssh_filled, lon_filled, lat_filled, x_ac_filled = swotd.fill_nadir_gap(ssh, lon, lat, x_ac, x_al, method='fill_value')  # fill the nadir gap with masked fill values
        
    ssh_filled_d = ssh_filled
    
    if order == '1':
        deriv_ssh = np.sqrt(swotd.gradx(ssh_filled_d)**2 + swotd.grady(ssh_filled_d)**2)
    
    elif order == '2':
        deriv_ssh   = swotd.laplacian(ssh_filled_d)
    
    deriv_ssh_ma = np.ma.array(deriv_ssh, mask = ssh_filled.mask, fill_value = ssh_filled.fill_value )       # generate masked array
    deriv_ssh_ma = swotd.empty_nadir_gap(deriv_ssh_ma, x_ac_filled, ssh, x_ac)                                 # Remove value in the gap

    return deriv_ssh_ma

def mask_borders(var):
    varm = var.copy()
    varm[:,0] = np.ma.masked
    varm[:,-1] = np.ma.masked
    varm[0,:] = np.ma.masked
    varm[-1,:] = np.ma.masked
    return varm

box_p09 = [3., 5., 37.5, 39.5]
box_p22 = [1.5, 3.5, 37., 39]


def swot_filt_qual_plot(filedir, filedir_den_root, lambdas, levels_ssh, savedir, ncyc, npass, dseason, vms):
        """
        ncyc_break = until which cycle (included) you want to generate plots.
        vms = array with colorbar vmins and maxs for each var [vmin_SSH, vmax_SSH, vmin_gra, vmax_gra, vmin_lap, vmax_lap]. If set to 'None', automatically set.
        """

        filename_a = filedir + 'all_noises/' + 'MED_1km_nogap_' + dseason + '_swotFastPhase_BOX_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_v2.nc' # modified for local from cal1
        print(filename_a)
        
        filename_k = filedir + 'just_karin/' + 'MED_1km_nogap_' + dseason + '_swotFastPhase_BOX_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_v2.nc'
        print(filename_k)
    
        plot_figrows = 3

        # Data:
        fileroot = filename_a.split('_c')[0].split('/')[-1]
        #fileroot = filename.split('.nc')[0].split('/')[-1]
        #MED_1km_nogap_JAS12_swotFastPhase_BOX_c01_p009_v2

        #ncyc  = fileroot.split('_c')[-1].split('_p')[0]
        #npass = filename.split('_p')[-1].split('_v')[0]
        print('Cycle ', ncyc, 'pass ', npass)

        if npass  == np.int(9):
            boxp = box_p09
        elif npass == np.int(22):
            boxp = box_p22
        else:
            print 'box error'

        # Load SSH_model and _obs data:
        sshm, lon, lat, x_ac, x_al = swotd.read_data(filename_a, 'ADT_model_box', 'lon_box', 'lat_box', 'x_ac', 'x_al')

        sshoa, _, _, _, _ = swotd.read_data(filename_a
              , 'ADT_obs_box', 'lon_box', 'lat_box', 'x_ac', 'x_al')

        sshok, _, _, _, _ = swotd.read_data(filename_k
              , 'ADT_obs_k_box', 'lon_box', 'lat_box', 'x_ac', 'x_al')
        
        # Calculate derivatives:

        am = derivatives_calc_method(sshm, lon, lat, x_ac, x_al, order='1')
        aoa = derivatives_calc_method(sshoa, lon, lat, x_ac, x_al, order='1')
        aok = derivatives_calc_method(sshok, lon, lat, x_ac, x_al, order='1')

        vm = derivatives_calc_method(sshm, lon, lat, x_ac, x_al, order='2')
        voa = derivatives_calc_method(sshoa, lon, lat, x_ac, x_al, order='2')
        vok = derivatives_calc_method(sshok, lon, lat, x_ac, x_al, order='2')
        vm = mask_borders(vm)
        voa = mask_borders(voa)
        vok = mask_borders(vok)
        
        if vms == 'None':
            vs = np.empty(6)
            vs[0] = round(sshm.min()-.05, -int(np.floor(np.log10(abs(sshm.min()-.05)))))
            vs[1] = round(sshm.max()+.05, -int(np.floor(np.log10(abs(sshm.max()+.05)))))

            vs[2] = am.min() #+ .05
            vs[3] = am.max() #- .05
            vs[4] = vm.min() #+ .05
            vs[5] = vmin*-1 #lapm.max()/2. #- .05
        else:
            vs = vms

        for ii in xrange (0, len(lambdas)):
           
            fig1 = plt.figure(figsize=(14, 10))  # (w,h)
            fig1.patch.set_facecolor('none')
        
        
        
            gs = gridspec.GridSpec(plot_figrows, 6, width_ratios=[.19, .19, .19, .19, .19, .01])

            ax01 = plt.subplot(gs[0,0])
            ax02 = plt.subplot(gs[0,1])
            ax03 = plt.subplot(gs[0,2])
            ax04 = plt.subplot(gs[0,3])
            ax05 = plt.subplot(gs[0,4])
            ax11 = plt.subplot(gs[1,0])
            ax12 = plt.subplot(gs[1,1])
            ax13 = plt.subplot(gs[1,2])
            ax14 = plt.subplot(gs[1,3])
            ax15 = plt.subplot(gs[1,4])
            ax21 = plt.subplot(gs[2,0])
            ax22 = plt.subplot(gs[2,1])
            ax23 = plt.subplot(gs[2,2])
            ax24 = plt.subplot(gs[2,3])
            ax25 = plt.subplot(gs[2,4])

	    print('lambda: ', lambdas[ii])

	    #filedir_den = filedirs_den #[ii]

	    #for mm in xrange(1, len(myfiles)):
	    #    
	    #    filename = myfiles[mm]

	    # One file only:
	    #print('filedir_den')
	    #print(filedir_den)
	    print('fileroot')
	    print(fileroot)

	    filename_den_a = filedir + "all_noises/" + fileroot + '_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_v2_denoised_p1_' + str(0).zfill(5) + '_p2_' + str(np.int(lambdas[ii])).zfill(7) + '_p3_' + str(0).zfill(7) + '.nc' 
        #'/mnt/meom/workdir/gomeznl/data/occigen_outputs/all_noises_f' + filedir_den + fileroot + '_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_v2_denoised_p1_' + str(0).zfill(5) + '_p2_' + str(np.int(lambdas[ii])).zfill(7) + '_p3_' + str(0).zfill(7) + '.nc' #str(lambdas[ii][2]).zfill(7) + '.nc'
        
	    #MED_1km_nogap_JAS13_swotFastPhase_BOX_c92_p009_v2_denoised_p1_00000_p2_0000705_p3_0000000.nc
	    print('filename_den_a')
	    print(filename_den_a)
	    
	    filename_den_k = filedir + "just_karin/" + fileroot + '_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_v2_denoised_p1_' + str(0).zfill(5) + '_p2_' + str(np.int(lambdas[ii])).zfill(7) + '_p3_' + str(0).zfill(7) + '.nc'
        #'/mnt/meom/workdir/gomeznl/data/occigen_outputs/just_karin_f' + filedir_den + fileroot + '_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_v2_denoised_p1_' + str(0).zfill(5) + '_p2_' + str(np.int(lambdas[ii])).zfill(7) + '_p3_' + str(0).zfill(7) + '.nc' #str(lambdas[ii][2]).zfill(7) + '.nc'
            
	    # SSH

	    # Load SSH_obs denoised data:
	    sshoaf, _, _, _, _ = swotd.read_data(filename_den_a, 'SSH'
						       , 'lon', 'lat', 'x_ac', 'x_al')
	    
            sshokf, _, _, _, _ = swotd.read_data(filename_den_k, 'SSH'
                                                       , 'lon', 'lat', 'x_ac', 'x_al')
	    
	    # Calculate derivatives:

	    aoaf = derivatives_calc_method(sshoaf, lon, lat, x_ac, x_al, order='1')

	    voaf = derivatives_calc_method(sshoaf, lon, lat, x_ac, x_al, order='2')
	    voaf = mask_borders(voaf)
	    
	    aokf = derivatives_calc_method(sshokf, lon, lat, x_ac, x_al, order='1')

            vokf = derivatives_calc_method(sshokf, lon, lat, x_ac, x_al, order='2')
            vokf = mask_borders(vokf)            
	    
	    am[:,50] = np.ma.masked
 	    aoa[:,50] = np.ma.masked
	    aok[:,50] = np.ma.masked
	    aoaf[:,50] = np.ma.masked
	    aokf[:,50] = np.ma.masked

	    savename = savedir + fileroot + '_c' + str(ncyc).zfill(2) + '_p' + str(npass).zfill(3) + '_' + str(filedir_den_root.split('/')[-3]) + '_' + str(filedir_den_root.split('/')[-2]) + '_lambd_' + str(int(lambdas[ii])).zfill(5) + '_KA.png'
	    print savename

	    cmap1 = 'viridis' #cmocean.cm.haline #same, viridis #getncvcm('banded', reverse=False, half=0)
	    cmap2 = 'inferno' #worse, magma #getncvcm('manga', reverse=False, half=0)
	    cmap3 = 'RdBu' #cmocean.cm.curl #pink green #

	    #lon_sat[lon_sat > 180] -= 360

	    #SSH_r = np.ma.masked_invalid(SSH_r)

	    vmin1 = vs[0]
	    vmax1 = vs[1]

	    vmin2 = vs[2]
	    vmax2 = vs[3]

	    vmin3 = vs[4]
	    vmax3 = vs[5]

	    ##################################################

	    my_map01, c1 = plot_swot(ax01, sshm, lon, lat, vmin1, vmax1, cmap1, merv=[0,0,0,0], parv=[1,0,0,1]
				     , box=boxp, plotstyle='0', levels=levels_ssh)
	    ax01.set_ylabel('$SSH$ $[m]$', labelpad=40, size = 16)
	    ax01.set_title('SSH_model', size=14)

	    plot_swot(ax02, sshok, lon, lat, vmin1, vmax1, cmap1, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
		      , plotstyle='0', levels=levels_ssh)
	    ax02.set_title('SSH_obs_K', size=14)

	    plot_swot(ax03, sshoa, lon, lat, vmin1, vmax1, cmap1, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
		      , plotstyle='0', levels=levels_ssh)
	    ax03.set_title('SSH_obs', size=14)
		
	    plot_swot(ax04, sshokf, lon, lat, vmin1, vmax1, cmap1, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
                      , plotstyle='0', levels=levels_ssh)
            ax04.set_title('SSH_obs_K_f', size=14)

            plot_swot(ax05, sshoaf, lon, lat, vmin1, vmax1, cmap1, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
                      , plotstyle='0', levels=levels_ssh)
            ax05.set_title('SSH_obs_f', size=14)

	    # vel.

	    my_map02, c2 = plot_swot(ax11, am, lon, lat, vmin2, vmax2, cmap2, merv=[0,0,0,0], parv=[1,0,0,1]
				     , box=boxp, plotstyle='0')
	    ax11.set_ylabel(r'$|\nabla$$SSH|$ $[m]$', labelpad=40, size = 16)
	    
	    plot_swot(ax12, aok, lon, lat, vmin2, vmax2, cmap2, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
		      , plotstyle='0')

	    plot_swot(ax13, aoa, lon, lat, vmin2, vmax2, cmap2, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
		      , plotstyle='0')
	    
            plot_swot(ax14, aokf, lon, lat, vmin2, vmax2, cmap2, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
                      , plotstyle='0')

            plot_swot(ax15, aoaf, lon, lat, vmin2, vmax2, cmap2, merv=[0,0,0,0], parv=[0,0,0,0], box=boxp
                      , plotstyle='0')

            # vort.
            my_map03, c3 = plot_swot(ax21, vm, lon, lat, vmin3, vmax3, cmap3, merv=[1,0,0,1], parv=[1,0,0,1]
                         , box=boxp, plotstyle='0')
            ax21.set_ylabel(r'$\Delta$$SSH$ $[m]$', labelpad=40, size = 16)

            plot_swot(ax22, vok, lon, lat, vmin3, vmax3, cmap3, merv=[1,0,0,1], parv=[0,0,0,0], box=boxp
                      , plotstyle='0')

            plot_swot(ax23, voa, lon, lat, vmin3, vmax3, cmap3, merv=[1,0,0,1], parv=[0,0,0,0], box=boxp
                  , plotstyle='0')

            plot_swot(ax24, vokf, lon, lat, vmin3, vmax3, cmap3, merv=[1,0,0,1], parv=[0,0,0,0], box=boxp
                      , plotstyle='0')

            plot_swot(ax25, voaf, lon, lat, vmin3, vmax3, cmap3, merv=[1,0,0,1], parv=[0,0,0,0], box=boxp
                      , plotstyle='0')

            ####################################################################################

            axC1 = plt.subplot(gs[0,-1])
            axC2 = plt.subplot(gs[1,-1])
            axC3 = plt.subplot(gs[2,-1])

            cbar1 = plt.colorbar(c1, cax=axC1, extend='both')
            cbar1.ax.tick_params(labelsize=10)
            #cbar1.formatter.set_powerlimits((0, 0))
            #cbar1.update_ticks()

            cbar2 = plt.colorbar(c2, cax=axC2, extend='max')
            cbar2.ax.tick_params(labelsize=10)
            cbar2.formatter.set_powerlimits((0, 0))
            cbar2.update_ticks()

            cbar3 = plt.colorbar(c3, cax=axC3, extend='both')
            cbar3.ax.tick_params(labelsize=10)
            cbar3.formatter.set_powerlimits((0, 0))
            cbar3.update_ticks()            
 
	    #savename = '/Users/laura/Documents/PhD_stuff/my_presentations/SWOTmeetingCanada_2018/presentation/figs/' + '3_vars_p015_DEF_redBlue_M' + '.jpeg'
	    plt.savefig(str(savename), bbox_inches='tight', dpi=300)

	    plt.show()
	    
	    #plot_comb(lon, lat, ssh_model, ssh_obs, ssh_obs_den, ssh_model_gra, ssh_obs_gra
	    #	  , ssh_obs_den_gra, ssh_model_lap, ssh_obs_lap, ssh_obs_den_lap
	    #	  , vs, box, savename, ax01, ax02, ax03, ax11, ax12, ax13, ax21, ax22, ax23, levels_ssh)
      
	    #ssh_diffs(lon, lat, ssh_model, ssh_obs, ssh_obs_den, suptitle, savename)
	    #break
	    print('Figure saved at: ', savename)
	    #plt.close('all')

	    #break
	#if ncyc_break == ncyc:
	#    break


