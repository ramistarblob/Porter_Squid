import gofish
import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube
import bettermoments as bm
import astropy.io.fits as pyfits
from astropy import wcs
import cmocean
import seaborn as sns
from astropy import units as u
from pprint import pprint
import math
from astropy.wcs import WCS
from astropy.visualization import quantity_support
from pvextractor import extract_pv_slice, Path
from pvextractor import PathFromCenter
from astropy.coordinates import Galactic
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.utils import data
from astropy.coordinates import SpectralCoord
from radio_beam import Beam
import astropy.constants as c
import astropy.units as u



'''creates moment maps and applies channel mask; firstchannel and lastchannel are saved to the dictionary'''

def fits_info(path):
    fits =  pyfits.open(path)
    header = fits[0].header
    pprint(header)
    fit_wcs = WCS(header)
    fits.close()
    return fit_wcs

def channel_map(isotop, path, cube_FOV, smooth, polyorder, clip=None, smooth_threshold_mask=0.0,vel_range=None, cmap=cmocean.cm.ice): #firstchannel=0, lastchannel=-1, vel_range=None):
    #if firstchannel != 0:
    #    firstchannel= isotop['firstchannel']
    #if lastchannel != -1:
    #    lastchannel = isotop['lastchannel'] 

    cube = imagecube(isotop['file'], FOV=cube_FOV, bunit='Jy/beam')
    data, velax = bm.load_cube(isotop['file'])
    #data = cube.data
    velax=velax/1000
    #print("velocity axis", "velax")
    smoothed_data = bm.smooth_data(data=data, smooth=smooth, polyorder=polyorder)
    data = smoothed_data
    #print(data)
    rms_smoothed = bm.estimate_RMS(data=data, N=5)  # this estimates rms from 1st a nd last 5 channels 
    print('RMS =',rms_smoothed)

    fits =  pyfits.open(path)
    header = fits[0].header

    fit_wcs = WCS(header)
    print(fit_wcs)
    fits.close()

    if clip is not None:
        threshold_mask = bm.get_threshold_mask(data=data,
                                                clip=clip, # this is a 2sigma clip to throw away noise artefacts
                                                smooth_threshold_mask=smooth_threshold_mask) 

    else:
        threshold_mask = bm.get_threshold_mask(data=data,
                                                clip=None, # this is a 2sigma clip to throw away noise artefacts
                                                smooth_threshold_mask=0.0) 
        
    if vel_range is None:
        user_mask_3d_fin = bm.get_user_mask(data=data, user_mask_path=None)
        #vel_range = [(-5.8847, 9.8319), (14.4358, 25.8661)]  # Default velocity ranges
    else:
    # Convert the velocity range into quantities with units of km/s
    #vel_range_qty = [
    #    (vel_range[0][0], vel_range[0][1]),
    #    (vel_range[1][0], vel_range[1][1])
    #]

        # Generate the velocity mask based on the provided velocity range
        mask = np.zeros_like(velax, dtype=bool)
        #pprint(mask)
        for i in range(len(mask)):
           for vmin, vmax in vel_range:
                if vmin <= vmax:
                    # Case for increasing velocities
                    mask[i] = (velax[i] >= vmin) and (velax[i] <= vmax) or mask[i]
                else:
                    # Case for decreasing velocities
                    mask[i] = (velax[i] <= vmin) and (velax[i] >= vmax) or mask[i]



        # Convert to integer (1 for included channels, 0 for excluded)
        user_mask = mask.astype(int)
        #pprint(user_mask)
        #pprint(len(user_mask))

        # Expand the user mask to match the shape of the 3D data cube
        user_mask_3d = np.expand_dims(user_mask, axis=(1,2)) 
        user_mask_3d_fin = np.broadcast_to(user_mask_3d, data.shape)  # Broadcast to match the shape of the data cube
        #print(user_mask_3d_fin)


    #channel_mask = bm.get_channel_mask(data=data,
                                       # firstchannel=firstchannel,
                                       #lastchannel=lastchannel)  # here we throw away all channels outside of the line width
    

    # Combine masks: apply the user mask to the threshold and channel masks
    #combined_mask = bm.get_combined_mask(user_mask=user_mask_3d_fin,
                                         #threshold_mask=threshold_mask,
                                         #channel_mask=channel_mask,
                                         #combine='and')

    masked_data=data*user_mask_3d_fin*threshold_mask
    #print(user_mask_3d)
    #print(user_mask_3d_fin)
    # Filter out channels where the data is entirely masked
    valid_channels = [ch for ch in range(masked_data.shape[0]) if np.any(masked_data[ch])]
    #print(valid_channels)

    col = 5  # Or any number of columns you'd like to set initially
    row = math.ceil(len(valid_channels) / col)  # Automatically calculate the required number of rows

    fig, ax = plt.subplots(row, col, figsize=(17, 3 * row))

        # Flattening ax array to handle any number of subplots (even if it is a non-rectangular grid)
    ax = ax.flatten()

    #print(masked_data.shape[0])
    # Loop over all channels and plot only those that fall within the velocity range
    for i in range(len(valid_channels)):
        # If the velocity index is within the specified range (i.e., we want to plot it)
        #if np.any(user_mask_3d_fin[i]):  # If any of the data is not masked
        im = ax[i].imshow(masked_data[valid_channels[i]], cmap=cmap, origin='lower',extent=cube.extent)
        ax[i].set_ylabel('Offset (arcsec)')
        ax[i].set_xlabel('Offset (arcsec)')
        cube.plot_beam(ax[i], x0=0.1, y0=0.1, color='white')
        #lon = ax[i].coords[0]
        #lat = ax[i].coords[1]
        #lon.set_axislabel("")
        #lat.set_axislabel("")
        ax[i].text(0.95, 0.95, s=f'{valid_channels[i]}', color='white', weight='bold', fontsize=12,
                    horizontalalignment='right', verticalalignment='top', transform=ax[i].transAxes)
        
        ax[i].text(0.05, 0.95, s=f'{velax[valid_channels[i]]:.2f} km/s', color='white', weight='bold', fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax[i].transAxes)

        #ax[i].axis('off')  # Turn off unused axe
        

    # Hide any unused subplots if there are fewer images than the grid
    for j in range(len(valid_channels), row * col):
        ax[j].axis('off')  # Turn off unused axes
    
    #fig.text(0.5, 0.09, 'RA', ha='center', va='center', fontsize=18)  # X-axis label
    #fig.text(0.02, 0.5, 'DEC', ha='center', va='center', rotation='vertical', fontsize=18)  # Y-axis label


    #plt.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.15, hspace=0.2, wspace=0.4)
    plt.tight_layout()
    plt.show()
     

def create_moment_files(isotop, smooth, polyorder, clip_mom0_8=None, clip_mom1_2_9=None, smooth_threshold_mask=0.0, vel_range=None):
    #if firstchannel != 0:
    #    isotop['firstchannel'] = firstchannel
    #if lastchannel != -1:
    #    isotop['lastchannel'] = lastchannel
    
    #ube = imagecube(isotop['file'], FOV=cube_FOV, bunit='Jy/beam')
    data, velax = bm.load_cube(isotop['file'])
    velax=velax/1000
    #print("velocity axis", "velax")
    smoothed_data = bm.smooth_data(data=data, smooth=smooth, polyorder=polyorder)
    data = smoothed_data
    #print(data)
    rms_smoothed = bm.estimate_RMS(data=data, N=5)  # this estimates rms from 1st and last 5 channels 
    print('RMS =',rms_smoothed)

    if vel_range is None:
        user_mask_3d_fin = bm.get_user_mask(data=data, user_mask_path=None)
        #vel_range = [(-5.8847, 9.8319), (14.4358, 25.8661)]  # Default velocity ranges
    else:
    # Convert the velocity range into quantities with units of km/s
    #vel_range_qty = [
    #    (vel_range[0][0], vel_range[0][1]),
    #    (vel_range[1][0], vel_range[1][1])
    #]

        # Generate the velocity mask based on the provided velocity range
        mask = np.zeros_like(velax, dtype=bool)
        #pprint(mask)
        for i in range(len(mask)):
           for vmin, vmax in vel_range:
                if vmin <= vmax:
                    # Case for increasing velocities
                    mask[i] = (velax[i] >= vmin) and (velax[i] <= vmax) or mask[i]
                else:
                    # Case for decreasing velocities
                    mask[i] = (velax[i] <= vmin) and (velax[i] >= vmax) or mask[i]



        # Convert to integer (1 for included channels, 0 for excluded)
        user_mask = mask.astype(int)
        #pprint(user_mask)
        #pprint(len(user_mask))

        # Expand the user mask to match the shape of the 3D data cube
        user_mask_3d = np.expand_dims(user_mask, axis=(1,2)) 
        user_mask_3d_fin = np.broadcast_to(user_mask_3d, data.shape)  # Broadcast to match the shape of the data cube
        #print(user_mask_3d_fin)

    if clip_mom0_8 is not None:
        threshold_mask = bm.get_threshold_mask(data=data,
                                                clip=clip_mom0_8, # this is a 2sigma clip to throw away noise artefacts
                                                smooth_threshold_mask=smooth_threshold_mask) 

    else:
        threshold_mask = bm.get_threshold_mask(data=data,
                                                clip=None, # this is a 2sigma clip to throw away noise artefacts
                                                smooth_threshold_mask=0.0) 
    
    if clip_mom1_2_9 is not None:
        threshold_mask_mom1 = bm.get_threshold_mask(data=data,
                                                    clip=clip_mom1_2_9,  # this is a 2sigma clip to throw away noise artefacts
                                                    smooth_threshold_mask=smooth_threshold_mask)
    
    else: 
        threshold_mask_mom1 = bm.get_threshold_mask(data=data,
                                                    clip=None,  # this is a 2sigma clip to throw away noise artefacts
                                                    smooth_threshold_mask=0.0)

    #channel_mask = bm.get_channel_mask(data=data,
                                       #firstchannel=isotop['firstchannel'],
                                       #lastchannel=isotop['lastchannel'])  # here we throw away all channels outside of the line width

    # Combine masks: apply the user mask to the threshold and channel masks
    #combined_mask = bm.get_combined_mask(user_mask=user_mask_3d,
                                         #threshold_mask=threshold_mask,
                                         #channel_mask=channel_mask,
                                         #combine='and')
    
    #mask_mom1 = bm.get_combined_mask(user_mask=user_mask_3d,
                                     #threshold_mask=threshold_mask_mom1,
                                     #channel_mask=channel_mask,
                                     #combine='and')

    # Apply the mask to the data
    masked_data=data*user_mask_3d_fin*threshold_mask  # Now we mask the data accordingly
    masked_data_mom1_2_9=data*user_mask_3d_fin*threshold_mask_mom1

    # Create moment maps
    moment0 = bm.collapse_zeroth(velax=velax, data=masked_data, rms=rms_smoothed)
    moment1 = bm.collapse_first(velax=velax, data=masked_data_mom1_2_9, rms=rms_smoothed)
    moment2 = bm.collapse_second(velax=velax, data=masked_data_mom1_2_9, rms=rms_smoothed)
    moment8 = bm.collapse_eighth(velax=velax, data=masked_data, rms=rms_smoothed)
    moment9 = bm.collapse_ninth(velax=velax, data=masked_data_mom1_2_9, rms=rms_smoothed)

    # Save moment maps to FITS files
    bm.save_to_FITS(moments=moment0, method='zeroth', path=isotop['file'])
    bm.save_to_FITS(moments=moment1, method='first', path=isotop['file'])
    bm.save_to_FITS(moments=moment2, method='second', path=isotop['file'])
    bm.save_to_FITS(moments=moment8, method='eighth', path=isotop['file'])
    bm.save_to_FITS(moments=moment9, method='ninth', path=isotop['file'])

'''outputs moment 0, 1, and 8, 9 images, as well as moment 0 snr plots and noise plots'''
def moment_map(isotop, cube_FOV, thresh_mom1=None, thresh_mom2=None, thresh_mom9=None, rms_mom0=None, mom_0_contours=None,show_mom0_snr_plot=False, show_noise_plots=False):
    cube = imagecube(isotop['file'], FOV = cube_FOV,bunit='Jy/beam')
    file_M0   = isotop['file'][:-5]+'_M0.fits'
    file_M1   = isotop['file'][:-5]+'_M1.fits'
    file_M2   = isotop['file'][:-5]+'_M2.fits'
    file_M8   = isotop['file'][:-5]+'_M8.fits'
    file_M9   = isotop['file'][:-5]+'_M9.fits'

    file_dM0 = isotop['file'][:-5]+'_dM0.fits'
    file_dM1 = isotop['file'][:-5]+'_dM1.fits'
    file_dM2   = isotop['file'][:-5]+'_dM2.fits'
    file_dM8 = isotop['file'][:-5]+'_dM8.fits'
    file_dM9 = isotop['file'][:-5]+'_dM9.fits'

    cubeMom_0 = imagecube(file_M0, FOV=cube_FOV)
    cubeMom_1 = imagecube(file_M1, FOV=cube_FOV)
    cubeMom_2 = imagecube(file_M2, FOV=cube_FOV)
    cubeMom_8 = imagecube(file_M8, FOV=cube_FOV)
    cubeMom_9= imagecube(file_M9, FOV=cube_FOV)

    cubeMom_dm0 = imagecube(file_dM0, FOV=cube_FOV)
    cubeMom_dm1 = imagecube(file_dM1, FOV=cube_FOV)
    cubeMom_dm2 = imagecube(file_dM2, FOV=cube_FOV)
    cubeMom_dm8 = imagecube(file_dM8, FOV=cube_FOV)
    cubeMom_dm9 = imagecube(file_dM9, FOV=cube_FOV)

    print(f'rms mom 0 = {np.nanstd(cubeMom_dm0.data)}, rms_mom0_try = {np.sqrt(np.mean(cubeMom_dm0.data**2))}, rms mom 1 = {np.nanstd(cubeMom_dm1.data)}, rms mom 2 = {np.nanstd(cubeMom_dm2.data)} ,rms mom 8 = {np.nanstd(cubeMom_dm8.data)}')

    # Plot moment 0 map
    mapPeak=np.amax(cubeMom_0.data) 
    print('mom 0 map peak', mapPeak)

    contour_levels = [level * mapPeak for level in mom_0_contours]

    fig, ax = plt.subplots(figsize=(10, 5))
    pltmom0 = ax.imshow(cubeMom_0.data, origin='lower', cmap=cmocean.cm.ice_r, extent=cubeMom_0.extent)
    cubeMom_0.plot_beam(ax, x0=0.1, y0=0.1, color='white')
    cbar = plt.colorbar(pltmom0)
    cbar.set_label(r'Jy $\mathrm{beam}^{-1}$ $\mathrm{km}$ $\mathrm{s}^{-1}$', fontsize=16, labelpad=20, rotation=270)
    plt.contour(cubeMom_0.data, origin='lower', extent=cubeMom_0.extent, levels=contour_levels,colors='k',linewidths=1)
    plt.title('Moment 0 Map')
    plt.xlabel('Offset (arcsec)')
    plt.ylabel('Offset (arcsec)')
    plt.show()

    SNR_mom0 = cubeMom_0.data/cubeMom_dm0.data
    if thresh_mom1 is not None:
        # Set a threshold for significant emission
        threshold = thresh_mom1 * rms_mom0

        # Mask irrelevant pixels where moment 0 is below the threshold
        #mom1_masked = np.where(SNR_mom0 > thresh_mom1, cubeMom_1.data, np.nan)
        mom1_masked = np.where(cubeMom_0.data > threshold, cubeMom_1.data, np.nan)

        # Plot masked moment 1 
        plt.figure(figsize=(10, 5))
        pltmom1_mask = plt.imshow(mom1_masked, origin='lower', cmap='coolwarm', extent=cubeMom_1.extent)
        cbar = plt.colorbar(pltmom1_mask)
        cbar.set_label(label='km/s',  rotation=270, fontsize=16, labelpad=20)
        plt.title('Moment 1 Map')
        plt.xlabel('Offset (arcsec)')
        plt.ylabel('Offset (arcsec)')
    
    else:
        # Plot moment 1 
        plt.figure(figsize=(10, 5))
        pltmom1 = plt.imshow(cubeMom_1.data, origin='lower', cmap='coolwarm', extent=cubeMom_1.extent)
        cbar = plt.colorbar(pltmom1)
        cbar.set_label(label='km/s', rotation=270, fontsize=16, labelpad=20)
        plt.title('Moment 1 Map')
        plt.xlabel('Offset (arcsec)')
        plt.ylabel('Offset (arcsec)')
    
    #Plot moment 2 map
    if thresh_mom2 is not None:
        threshold = thresh_mom2 * rms_mom0
        #mom2_masked = np.where(SNR_mom0 > thresh_mom2, cubeMom_2.data, np.nan)
        mom2_masked = np.where(cubeMom_0.data > threshold, cubeMom_2.data, np.nan)
        plt.figure(figsize=(10, 5))
        pltmom2_masked= plt.imshow(mom2_masked, origin='lower', cmap='coolwarm', extent=cubeMom_2.extent)
        cbar = plt.colorbar(pltmom2_masked)
        cbar.set_label(label=r'$\mathrm{km}^{2}$/$\mathrm{s}^{2}$', rotation=270, fontsize=16, labelpad=20)
        plt.title('Moment 2 Map')
        plt.xlabel('Offset (arcsec)')
        plt.ylabel('Offset (arcsec)')

    else:
        plt.figure(figsize=(10, 5))
        pltmom2 = plt.imshow(cubeMom_2.data, origin='lower', cmap='coolwarm', extent=cubeMom_2.extent)
        cbar = plt.colorbar(pltmom2)
        cbar.set_label(label=r'$\mathrm{km}^{2}$/$\mathrm{s}^{2}$', rotation=270, fontsize=16, labelpad=20)
        plt.title('Moment 2 Map')
        plt.xlabel('Offset (arcsec)')
        plt.ylabel('Offset (arcsec)')

    #Plot moment 8 map
    plt.figure(figsize=(10, 5))
    pltmom8 = plt.imshow(cubeMom_8.data, origin='lower', cmap=cmocean.cm.ice, extent=cubeMom_8.extent)
    cbar = plt.colorbar(pltmom8)
    cbar.set_label(label=r'Jy $\mathrm{beam}^{-1}$',rotation=270, fontsize=16, labelpad=20)
    plt.title('Moment 8 Map')
    plt.xlabel('Offset (arcsec)')
    plt.ylabel('Offset (arcsec)')
    plt.show()

    #Plot moment 9 map
    if thresh_mom9 is not None:
        threshold = thresh_mom9 * rms_mom0
        #mom2_masked = np.where(SNR_mom0 > thresh_mom2, cubeMom_2.data, np.nan)
        mom9_masked = np.where(cubeMom_0.data > threshold, cubeMom_9.data, np.nan)
        plt.figure(figsize=(10, 5))
        pltmom9_masked= plt.imshow(mom9_masked, origin='lower', cmap='coolwarm', extent=cubeMom_9.extent)
        cbar = plt.colorbar(pltmom9_masked)
        cbar.set_label(label=r'$\mathrm{km}$/$\mathrm{s}$', rotation=270, fontsize=16, labelpad=20)
        plt.title('Moment 9 Map')
        plt.xlabel('Offset (arcsec)')
        plt.ylabel('Offset (arcsec)')

    else:
        plt.figure(figsize=(10, 5))
        pltmom9 = plt.imshow(cubeMom_9.data, origin='lower', cmap='coolwarm', extent=cubeMom_9.extent)
        cbar = plt.colorbar(pltmom9)
        cbar.set_label(label=r'$\mathrm{km}$/$\mathrm{s}$', rotation=270, fontsize=16, labelpad=20)
        plt.title('Moment 9 Map')
        plt.xlabel('Offset (arcsec)')
        plt.ylabel('Offset (arcsec)')


    if show_mom0_snr_plot == True:
        # moment 0 snr plot 
        fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(12, 3))
        div= np.divide(cubeMom_0.data, cubeMom_dm0.data)
        im_div = ax3.imshow(div, origin='lower', cmap=cmocean.cm.ice_r)
        ax3.contour(div, origin='lower',levels=[3,5], colors='black',linewidths=1)
        fig3.colorbar(im_div, ax=ax3)
    else:
        pass

    
    if show_noise_plots == True:
        fig2, ax2 = plt.subplots(nrows=1, ncols=5, figsize=(15,5))
        im0_dm = ax2[0].imshow(cubeMom_dm0.data, origin='lower', cmap=cmocean.cm.ice_r, extent=cubeMom_dm0.extent)

        im1_dm = ax2[1].imshow(cubeMom_dm1.data, cmap='coolwarm', origin='lower',extent=cubeMom_dm1.extent)
        im2_dm = ax2[2].imshow(cubeMom_dm2.data, cmap='coolwarm', origin='lower',extent=cubeMom_dm2.extent)

        im8_dm = ax2[3].imshow(cubeMom_dm8.data, cmap=cmocean.cm.ice, origin='lower', extent=cubeMom_dm8.extent)
        im9_dm = ax2[3].imshow(cubeMom_dm9.data, cmap=cmocean.cm.ice, origin='lower', extent=cubeMom_dm9.extent)

        fig2.colorbar(im0_dm)
        fig2.colorbar(im1_dm)
        fig2.colorbar(im2_dm)
        fig2.colorbar(im8_dm)
        fig2.colorbar(im9_dm)

        fig2.suptitle('Noise Moment Maps', fontsize=16)
    else:
         pass

''' finds brightest pixel in moment 0 map and outputs transformed coordinates '''

def finding_x0_y0(isotop, cube_FOV):
    file_M0   = isotop['file'][:-5]+'_M0.fits'
    hdulist = pyfits.open(file_M0)
    print(hdulist[0].header)
    w = wcs.WCS(hdulist[0].header)


    cubeMom_0 = imagecube(file_M0, FOV=cube_FOV)

    y0_pix, x0_pix = np.unravel_index(np.argmax(cubeMom_0.data),cubeMom_0.data.shape)
    print(np.argmax(cubeMom_0.data))
    print("Potential x0, y0 pix coords:",x0_pix, y0_pix)

    x0_deg, y0_deg = w.wcs_pix2world(x0_pix,y0_pix,1)

    print("x0 deg, y0 deg:",x0_deg,y0_deg)

    fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(15,5))
    ax1.imshow(cubeMom_0.data, origin='lower', extent=cubeMom_0.extent, cmap=cmocean.cm.ice_r)
    [0,0]

    # Convert pixel coordinates to data coordinates
    x0_data = (cubeMom_0.extent[0] + ((x0_pix +0.5)/ cubeMom_0.data.shape[0]) * (cubeMom_0.extent[1] - cubeMom_0.extent[0]))
    y0_data = (cubeMom_0.extent[2] + ((y0_pix +0.5)/ cubeMom_0.data.shape[1]) * (cubeMom_0.extent[3] - cubeMom_0.extent[2]))

    print("Potential x0, y0 data coords (arcsec):", x0_data, y0_data)

    # Mark the (x0, y0) point
    ax1.plot(x0_data, y0_data, 'x', color='red',markersize=10)  # 'ro' means red circle marker
    ax1.set_xlabel('Offset (arcsec)',fontsize=16)
    ax1.set_ylabel('Offset (arcsec)',fontsize=16)

    #seeing if this actually works

    #fig, ax2 = plt.subplots(nrows=1,ncols=1,figsize=(15,5))
    #ax2.imshow(cubeMom_0.data, origin='lower', cmap=cmocean.cm.ice_r)

    # Mark the (x0, y0) point
    #ax2.plot(x0_pix, y0_pix, 'x', color='red',markersize=10)

    #ax2.set_xlabel('RA (arcsec)',fontsize=16)
    #ax2.set_ylabel('Dec (arcsec)',fontsize=16)'''

    plt.show()

    
   
'''extracts line profiles '''     
#need to fig out what to do with giant list of params
def extract_line_profiles(isotop,cube_FOV, xlim, xlim_fill, vel_range=None,r_min=None, r_max=None, inc= None, PA=None,x0=None,y0=None,mstar=None, dist=None, apply_kep=False, mom_0_contours=None,line_center=[1], line_color=["firebrick"], show_line_center=False, line_label=None):
    
    cube = imagecube(isotop['file'], FOV = cube_FOV)
    data, velax = bm.load_cube(isotop['file'])
    velax=velax/1000
    #data = cube.data

    file_M0   = isotop['file'][:-5]+'_M0.fits'
    file_dM0 = isotop['file'][:-5]+'_dM0.fits'
    cubeMom_0 = imagecube(file_M0, FOV=cube_FOV)
    cubeMom_dm0 = imagecube(file_dM0, FOV=cube_FOV)
    rms_mom0 = np.nanstd(cubeMom_dm0.data)
    
    if apply_kep==True:
         xin, yin, dyin = cube.integrated_spectrum(r_min=r_min, r_max=r_max, inc=inc, PA=PA, x0=x0, y0=y0, empirical_uncertainty=True, mstar=mstar, dist=dist)
    else:
        xin, yin, dyin = cube.integrated_spectrum(r_min=r_min, r_max=r_max, inc=inc, PA=PA, x0=x0, y0=y0,empirical_uncertainty=True)
    
        #print(len(xin),len(yin))
        #print(xin,yin)
    print(yin)
                                        # xin [m/s], yin [Jy], dyin [Jy] Note Jy as now integrated unit
    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(15,5))
    ax1.errorbar(xin/1000, yin, dyin, fmt=' ', capsize=1.5, capthick=1.5, color='firebrick', lw=1.0) #note only works when mstar and dist are called
    ax1.step(xin/1000, yin, where='mid', color='firebrick', lw=2.0,label=line_label)
    ax1.plot([-100,100],[0.0,0.0],color='k',ls=':')

    y_min, y_max = ax1.get_ylim()
    ax1.fill_between([xlim_fill[0],xlim_fill[1]],[-20,-20],[y_max+20,y_max+20],color='firebrick',alpha=0.1) #this spans the complete line width
    
    if show_line_center==True:
        for lc,col in zip(line_center,line_color):
            ax1.plot([lc,lc],[y_min-2,y_max+2],color=col,ls='--',alpha=0.5)
    else:
        pass
    
    #print(vel_range)
    #print(np.argmin(np.abs(velax[::-1] - vel_range[0])))
    #print(velax[::-1][23])
    index = np.where(velax[::-1] == vel_range[0])[0][0]
    index2 = np.where(velax[::-1] == vel_range[1])[0][0]
    print(index,index2)

    vel_range_index = [index,index2]
    vmax = vel_range_index[0]
    vmin = vel_range_index[1]

    xin_km = xin/1000

    xin_no_em = np.concatenate((xin[:vmin], xin[vmax:]))
    yin_no_em = np.concatenate((yin[:vmin], yin[vmax:]))

    # If you also want xin_no_em in km/s:
    #xin_no_em_km = xin_no_em / 1000
    dF = np.sqrt(len(xin[vmin:vmax]))*np.std(yin_no_em)

    #print(xin[vmin:vmax])
    #print(vmin,vmax)
    int_flux = (np.trapz(yin[vmin:vmax], xin[vmin:vmax]))
    err = (np.sqrt((0.05*int_flux)**2+(dF)**2))
    
    print(int_flux/1000)
    print(err/1000)

    ax1.set_title(r'Integrated Flux = {:.4f} $\pm$ {:.4f} Jy$\,$km$\,$/s'.format(int_flux/1000,err/1000), fontsize=17)
    #print(len(yin))
    print('y array to integrate (Jy)',yin[vmin:vmax])
    print('x array corresponding to y (km/s)',xin_km[vmin:vmax])
    ax1.set_xlabel(r'Velocity [km$\,$s$^{-1}$]', fontsize=17)
    ax1.set_ylabel(r'Integrated line flux [Jy]', fontsize=17)
    ax1.tick_params(axis = 'both', labelsize = 15)

    ax1.set_xlim(xlim[0],xlim[1]),ax1.set_ylim(y_min-0.5,y_max+0.5)

    ax1.legend(loc='upper left')
    im0 = ax2.imshow(cubeMom_0.data, origin='lower', extent=cubeMom_0.extent, cmap=cmocean.cm.ice_r)
    #ax[1].contour(cubeMom_0.data, origin='lower', extent=cubeMom_0.extent, levels=[500, 2500])
    mapPeak=np.amax(cubeMom_0.data)
    rms=5
    color=['red','white','green']
    mapPeak=np.amax(cubeMom_0.data) 
    print('mom 0 map peak', mapPeak)
    contour_levels = [level * mapPeak for level in mom_0_contours]
    ax2.contour(cubeMom_0.data, origin='lower', extent=cubeMom_0.extent, levels=contour_levels,color=color)
    #print(10*rms, 100*rms)
    cube.plot_mask(ax=ax2, r_min=r_min, r_max=r_max, PA_min=-180.0, PA_max=180.0, mask_frame='disk',inc=inc, PA=PA, x0=x0, y0=y0)
    

    # Mark the (x0, y0) point
    ax2.plot(x0,y0, 'x', color='red',markersize=10)  # 'ro' means red circle marker
    ax2.set_title('Moment 0 Map', fontsize = 17)
    ax2.set_xlabel('Offset (arcsec)', fontsize = 17)
    ax2.set_ylabel('Offset (arcsec)', fontsize = 17)
    ax2.tick_params(axis = 'both', labelsize=15)
    cube.plot_beam(ax=ax2)
    cbar = fig.colorbar(im0)
    cbar.set_label('Jy$\,$beam$^{-1}\,$km$\,$s$^{-1}$',rotation=270,fontsize=17, labelpad=20)
    cbar.ax.tick_params(labelsize=15) 

'''pv diagram'''
def plot_pv(isotop, cube, smooth, polyorder,thresh,rms_mom0, mom_title, path_coords, width, fit_wcs,cmap,cbar_label, shrink,moment_xlim=None, moment_ylim=None, pv_xlim=None,pv_ylim=None):
    
    #path1 = PathFromCenter(center=Galactic(180.0289276656*u.deg,-2.1508378815*u.deg),length=40*u.arcsec,angle=25*u.deg, width=1)

    file_M0   = isotop['file'][:-5]+'_M0.fits'
    file_M1   = isotop['file'][:-5]+'_M1.fits'
    file_M2   = isotop['file'][:-5]+'_M2.fits'
    file_M8   = isotop['file'][:-5]+'_M8.fits'
    file_M9   = isotop['file'][:-5]+'_M9.fits'

    cubeMom_0 = imagecube(file_M0)
    cubeMom_1 = imagecube(file_M1)
    cubeMom_2 = imagecube(file_M2)
    cubeMom_8 = imagecube(file_M8)
    cubeMom_9= imagecube(file_M9)
    
    threshold = thresh * rms_mom0

    if cube == 0:
        cube = cubeMom_0
    elif cube == 1:
        cube = cubeMom_1
        cube.data =  np.where(cubeMom_0.data > threshold, cubeMom_1.data, np.nan)
    elif cube == 2:
        cube = cubeMom_2
        cube.data =  np.where(cubeMom_0.data > threshold, cubeMom_2.data, np.nan)
    elif cube == 8:
        cube = cubeMom_8
    elif cube == 9:
        cube= cubeMom_9
        cube.data =  np.where(cubeMom_0.data > threshold, cubeMom_9.data, np.nan)

    data, velax = bm.load_cube(isotop['file'])
    #data = cube.data
    velax=velax/1000
    #print("velocity axis", "velax")
    smoothed_data = bm.smooth_data(data=data, smooth=smooth, polyorder=polyorder)
    data = smoothed_data

    path = Path(path_coords,width=width)
    pv_diagram = extract_pv_slice(data, path,wcs=fit_wcs)
    

    crpix_hz = pv_diagram.header["CRPIX2"]
    crpix_deg = pv_diagram.header["CRPIX1"]
    restfreq_Hz=pv_diagram.header["RESTFRQ"]*u.Hz
    cdelt_Hz=pv_diagram.header["CDELT2"]* u.Hz
    cdelt_deg=pv_diagram.header["CDELT1"]*u.deg
    crval_Hz=pv_diagram.header["CRVAL2"]* u.Hz
    crval_deg=pv_diagram.header["CRVAL1"]*u.deg
    n_pix_hz = pv_diagram.header["NAXIS2"]
    n_pix_deg = pv_diagram.header["NAXIS1"]

    # --- Frequency axis in Hz ---
    pixels_freq = np.arange(n_pix_hz)
    freq_axis = crval_Hz + (pixels_freq - (crpix_hz - 1)) * cdelt_Hz

    # --- SpectralCoord: frequency → velocity (radio convention) ---
    spec_coord = SpectralCoord(freq_axis,
                            unit=u.Hz,
                            doppler_rest=restfreq_Hz,
                            doppler_convention='radio')
    vel_axis = spec_coord.to(u.km/u.s)  # Convert frequency to velocity (radio convention)

    #--spatial axis--

    pixels_deg = np.arange(n_pix_deg)
    deg_axis = crval_deg + (pixels_deg - (crpix_deg - 1)) * cdelt_deg

    # Plotting
    fig = plt.figure(figsize=(25, 15))  # Increase figure size if needed

    # Moment 0 map
    ax1 = fig.add_subplot(121, projection=fit_wcs, slices=('x', 'y', 0, 0))
    im0 = ax1.imshow(cube.data, origin='lower', cmap=cmap)
    path.show_on_axis(ax1, spacing=1, linestyle=':', edgecolor='r', linewidth=1.5)  # Slightly thicker line for path
    #mom0_cube.plot_beam(ax1, x0=0.1, y0=0.1, color='white')
    #mom1_cube.plot_beam(ax1, x0=70, y0=70, color='black')
    cbar = fig.colorbar(im0, ax=ax1, pad=0.01,shrink=shrink)  
    cbar.ax.tick_params(labelsize=20) 
    cbar.set_label(cbar_label, fontsize=22, labelpad=10) 


    # Add axis labels to the first plot
    ax1.set_xlabel('RA', fontsize=22)
    ax1.set_ylabel('Dec', fontsize=22)
    ax1.set_title(mom_title, fontsize=22)
    ax1.tick_params(axis='both', labelsize=18)  # Increase tick label size
    if moment_xlim is not None and moment_ylim is not None:
        ax1.set_xlim(moment_xlim[0],moment_xlim[1])
        ax1.set_ylim(moment_ylim[0],moment_ylim[1])
    else:
        pass

    # PV diagram with accurate velocity axis
    ax2 = fig.add_subplot(122)
    im = ax2.imshow(pv_diagram.data, origin='lower', cmap=cmocean.cm.ice)
    ax2.set_xlabel("Offset (degrees)", fontsize=22)  # x-axis in degrees
    ax2.set_ylabel("Velocity (km/s)", fontsize=22)  # y-axis for velocity
    ax2.set_yticks(np.linspace(0, n_pix_hz-1, 8))
    ax2.set_yticklabels([f"{v:.1f}" for v in np.linspace(vel_axis[0].value, vel_axis[-1].value, 8)], fontsize=18)
    ax2.tick_params(axis='x', labelsize=18)  # Increase x-tick label size

    # Update x-axis to show offsets in degrees (no unnecessary conversion)

    #offset_deg = (pixels - (crpix - 1)) * cdelt_deg  # Offset in degrees
    ax2.set_xticks(np.linspace(0, n_pix_deg-1, 10))
    ax2.set_xticklabels([f"{offset:.3f}" for offset in np.linspace(deg_axis[0].value, deg_axis[-1].value, 10)], fontsize=18)

    # Set y-limits (to zoom in the PV diagram as per your original code)


    # Create a properly sized colorbar

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.2)  # Wider colorbar with padding
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("", fontsize=22)  # Match other labels
    cb.ax.tick_params(labelsize=18)  # Increase colorbar tick label size
    cb.set_label("Jy/beam", fontsize=22) 
    if pv_xlim is not None:
        ax2.set_xlim(pv_xlim[0], pv_xlim[1])
    if pv_ylim is not None:
        ax2.set_ylim(pv_ylim[0],pv_ylim[1])
    else:
        pass
    # Adjust layout with more control
    plt.subplots_adjust(wspace=0.3)  # Add more space between subplots
    plt.tight_layout()

    # Show the figure
    plt.show()

    #return pv_diagram


    



'''mean intensity 1D profile'''
def mean_intensity_1D_profile(isotop, cube_FOV):
        file_M0   = isotop['file']+'_M0.fits'
        cubeMom_0 = imagecube(file_M0, FOV=cube_FOV)
        x, y, dy = cubeMom_0.radial_profile(inc=isotop['inc'], PA=isotop['PA'],x0=isotop['x0'], y0=isotop['y0'])
        fig, ax = plt.subplots(constrained_layout=True)
        ax.errorbar(x, y, dy)
        '''ax.plot range may change based on spec line'''
        ax.plot([0,12],[0,0],'--k')
        ax.set_yscale('linear')
        ax.set_xlabel('Radius (arcsec)')
        ax.set_ylabel('Integrated Intensity (mJy/beam m/s)')