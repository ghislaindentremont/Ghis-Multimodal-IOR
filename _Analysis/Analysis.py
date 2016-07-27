import os
import mne
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

os.chdir("/Volumes/Seagate Backup Plus Drive/Experiments/multimodal_ior/_Data/forMNE/BeforeSummer_ForAnalysis")


##############################################################################################
####                                Load Data                                             ####
##############################################################################################
# Participant
participant = "e12"

# for one participant: e12
raw = mne.io.read_raw_brainvision('multimodal_ior_%s.vhdr' % participant, preload = True)

# look at info
print(raw.info)

# look at channels in particular
print(raw.info['ch_names'])

# remove Aux
raw.drop_channels(['Aux1']) # AUX = microsensor 

# check that Aux was dropped
print(raw.info['ch_names'])

# look for bad channels 
raw.plot(n_channels = 64, scalings = dict(eeg = 100e-6), block = True)

# label bads for each P
if participant == "e12":
   raw.info['bads'] = []   

# channel names
ch_names = [
    # greens (1-32)
    'Fp1'
    , 'Fz'
    , 'F3'
    , 'F7'
    , 'FT9'
    , 'FC5'
    , 'FC1'
    , 'C3'
    , 'T7'
    , 'TP9'
    , 'CP5'
    , 'CP1'
    , 'Pz'
    , 'P3'
    , 'P7'
    , 'O1'
    , 'Oz'
    , 'O2'
    , 'P4'
    , 'P8'
    , 'TP10'
    , 'CP6'
    , 'CP2'
    , 'Cz'
    , 'C4'
    , 'T8'
    , 'FT10'
    , 'FC6'
    , 'FC2'
    , 'F4'
    , 'F8'
    , 'Fp2'
    # yellows (33-64)
    , 'AF7'
    , 'AF3'
    , 'AFz'
    , 'F1'
    , 'F5'
    , 'FT7'
    , 'FC3'
    , 'FCz'
    , 'C1'
    , 'C5'
    , 'TP7'
    , 'CP3'
    , 'P1'
    , 'P5'
    , 'PO7'
    , 'PO3'
    , 'POz'
    , 'PO4'
    , 'PO8'
    , 'P6'
    , 'P2'
    , 'CPz'
    , 'CP4'
    , 'TP8'
    , 'C6'
    , 'C2'
    , 'FC4'
    , 'FT8'
    , 'F6'
    , 'F2'
    , 'AF4'
    , 'AF8'
    ]

# get montage for interpolation and visualization
montage = mne.channels.read_montage(kind = "easycap-M1", ch_names = ch_names)
mne.viz.plot_montage(montage, show_names = True)

# rename channels in raw 
original_ch_names = raw.info['ch_names'][0:64]
ch_names_dict = dict(zip(original_ch_names, ch_names))
raw.rename_channels(mapping = ch_names_dict)

# add montage to raw instance 
raw.set_montage(montage)

# interpolate bads 
raw.interpolate_bads()

# see effect of interpolation 
raw.plot(n_channels = 64, scalings = dict(eeg = 100e-6), block = True)




##############################################################################################
####                                Re-Reference                                          ####
##############################################################################################
# events, LV/LV for re-reference visualization
events = mne.find_events(raw, stim_channel = 'STI 014', output = 'onset')

# visualize timeline of events 
mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp)

picks_all = mne.pick_types(raw.info, meg=False, eeg=True, eog=False, exclude='bads')
event_id, tmin, tmax = {'LV/LV': 20}, -0.2, 0.5
epochs_params_test = dict(picks = picks_all, events=events, event_id=event_id, tmin=tmin, tmax=tmax)


#------------------------- Show Effect of Reference on Evoked -------------------------------#
fig, ax = plt.subplots(2,1, sharex = True)

# no reference 
raw_no_ref, _ = mne.io.set_eeg_reference(raw, [])
evoked_no_ref = mne.Epochs(raw_no_ref, **epochs_params_test).average()
del raw_no_ref  # save memory

evoked_no_ref.plot(axes = ax[0], titles=dict(eeg='EEG Original reference'), show=False)

# Average reference 
raw_ref, _ = mne.io.set_eeg_reference(raw)
evoked_ref = mne.Epochs(raw_ref, **epochs_params_test).average()
del raw_ref  # save memory

evoked_ref.plot(axes = ax[1], titles=dict(eeg='EEG Average reference'))

# see addition of proj
print(raw.info['projs'])
#------------------------- Show Effect of Reference on Evoked -------------------------------#


#---------------------- Show Effect of Reference on Continuous ------------------------------#
raw.plot(n_channels = 64, scalings = dict(eeg = 100e-6), block = True)
#---------------------- Show Effect of Reference on Continuous ------------------------------#


#------------------------ Apply and Rid Reference Projector ---------------------------------#
raw.apply_proj()
raw.info['projs'] = []
#------------------------ Apply and Rid Reference Projector ---------------------------------#



##############################################################################################
####                                    Filters                                           ####
##############################################################################################

#---------------------------------- Butterworth Filter --------------------------------------#
# Pick a subset of channels 
picks = mne.pick_types(
    raw.info
    , meg = False
    , eeg = True
    , eog = False
    , selection = 
    ['Ch8', 'Ch25'     # somatosensory (left/right) 
    , 'Ch47', 'Ch51']  # occipital(left/right)
    )

# without filter
raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (1,0,0), show = False)

# apply filter 
iir_params = dict(order=2, ftype='butter')
raw.filter(
    .1
    , 50
    , picks=picks
    , method = 'iir'
    , iir_params = iir_params
    , n_jobs = 4
)
raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (0,1,0))
#---------------------------------- Butterworth Filter --------------------------------------#


#---------------------------------- Notch Filter --------------------------------------------#
# without filter
raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (1,0,0), show = False)

# apply notch
raw.notch_filter(
    np.arange(60,121,60)  # get 60 and 120 Hz 
    , picks=picks
    , n_jobs = 4
    )
raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (0,1,0))
#---------------------------------- Notch Filter --------------------------------------------#




##############################################################################################
####                                   Functions                                          ####
##############################################################################################

def get_evoked( raw, event_id, channel_name, tmin, tmax, reject_num, baseline = (-0.1,0), get_evoked = True ):
    reject = dict(eeg=reject_num)
    picks = mne.pick_types(
    raw.info
    , meg = False
    , eeg = True
    , eog = False
    , selection = channel_name
    )
    params = dict(
    picks =  picks  
    , events=events  # global variable
    , event_id=event_id
    , tmin=tmin
    , tmax=tmax
    )
    epochs= mne.Epochs(
    raw
    , **params
    , add_eeg_ref = False
    , baseline = baseline
    )
    # artifact rejection
    epochs.drop_bad(reject = reject)
    # percentage rejected by channel
    epochs.plot_drop_log()
    # visualize by channel, by epoch
    epochs.plot()
    
    if get_evoked:
        evoked = epochs.average()
        evoked.plot()
        # average channels 
        sums = np.zeros(len(evoked.data[0]))
        for i in range(0, len(evoked.data)):
            sums = evoked.data[i] + sums
        
        avg = np.array( [sums/len(evoked.data)] )
        info = mne.create_info(
            ch_names = [' '.join(channel_name)]
            , sfreq = raw.info['sfreq']
            , ch_types = 'eeg'
            )
        avg_wave = mne.EvokedArray(avg, info, tmin=tmin)
        avg_wave.plot()
        to_return = avg_wave
    else:
        to_return = epochs
    return to_return


def concatenate_epochs( epoch1, epoch2, ch_name, reversal = False ):
    epoch1_arr = epoch1.get_data()
    epoch2_arr = epoch2.get_data()
    comb_arr = np.concatenate( [epoch1_arr, epoch2_arr] )
    
    comb_arr_sum = np.zeros(len(comb_arr[0][0]))
    for idx in range(0, len(comb_arr)):
        temp = comb_arr[idx][0]
        comb_arr_sum = comb_arr_sum + temp
    
    comb_arr_avg = np.array( [comb_arr_sum/len(comb_arr)] ) # get type 2d array
    
    if reversal:
        comb_arr_avg = comb_arr_avg * (-1)
    
    info_comb = mne.create_info(
        ch_names = [ch_name]
        , sfreq = epoch1.info['sfreq']
        , ch_types = 'eeg'
        )
    
    evoked = mne.EvokedArray(comb_arr_avg, info_comb, tmin=tmin)
    evoked.plot()
    
    return evoked



# ##############################################################################################
# ####                          Visual Grand Average                                        ####
# ##############################################################################################

# vis_id = {
# 'LV/LV': 20
# , 'LV/RV': 22
# , 'LT/LV': 24
# , 'LT/RV': 26
# , 'RV/LV': 28
# , 'RV/RV': 30
# , 'RT/LV': 32
# , 'RT/RV': 34 
# }

# grand_avg_vis = get_evoked( raw, vis_id, ['Ch47','Ch51'], tmin, tmax, reject_num = 100e-6 )



# ##############################################################################################
# ####                          Tactile Grand Average                                       ####
# ##############################################################################################

# #---------------------------------- Left/Contra ---------------------------------------------#
# tact_LC_id = {
# 'LV/LT': 21
# , 'LT/LT': 25
# , 'RV/LT': 29
# , 'RT/LT': 33
# }
# epochs_LC_tact = get_evoked( raw, tact_LC_id, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
# #---------------------------------- Left/Contra ---------------------------------------------#


# #---------------------------------- Right/Contra --------------------------------------------#
# tact_RC_id = {
# 'LV/RT': 23
# , 'LT/RT': 27
# , 'RV/RT': 31
# , 'RT/RT': 35 
# }
# epochs_RC_tact = get_evoked( raw, tact_RC_id, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
# #---------------------------------- Right/Contra --------------------------------------------#


# #---------------------------------- Concatenate ---------------------------------------------#
# evoked_contra = concatenate_epochs(epochs_RC_tact, epochs_LC_tact, 'contra')
# #---------------------------------- Concatenate ---------------------------------------------#


# #---------------------------------- Left/Ipsi -----------------------------------------------#
# tact_LI_id = {
# 'LV/LT': 21
# , 'LT/LT': 25
# , 'RV/LT': 29
# , 'RT/LT': 33
# }
# epochs_LI_tact = get_evoked( raw, tact_LI_id, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
# #---------------------------------- Left/Ipsi -----------------------------------------------#


# #---------------------------------- Right/Ipsi ----------------------------------------------#
# tact_RI_id = {
# 'LV/RT': 23
# , 'LT/RT': 27
# , 'RV/RT': 31
# , 'RT/RT': 35 
# }
# epochs_RI_tact = get_evoked( raw, tact_RI_id, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
# #---------------------------------- Right/Ipsi ----------------------------------------------#


# #---------------------------------- Concatenate ---------------------------------------------#
# evoked_ipsi = concatenate_epochs(epochs_RI_tact, epochs_LI_tact, 'ipsi', reversal = True)
# #---------------------------------- Concatenate ---------------------------------------------#


# #-------------------------------- Join Ipsi/Contra ------------------------------------------#
# # average channels
# tact_avg = np.array( [(evoked_ipsi.data[0] + evoked_contra.data[0])/2] )  # get form (n_channel, n_times)
# info_tact = mne.create_info(
#     ch_names = ['Avg Ch25 & Ch8']
#     , sfreq = raw.info['sfreq']
#     , ch_types = 'eeg'
#     )
# grand_avg_tact = mne.EvokedArray(tact_avg, info_tact, tmin=tmin)
# grand_avg_tact.plot()
# #-------------------------------- Join Ipsi/Contra ------------------------------------------#




##############################################################################################
####                                    By Condition                                      ####
##############################################################################################


####################################### Tactile/Tactile ######################################

#---------------------------------- Cued/Contra/TT ------------------------------------------#
epochs_contra_LTLT = get_evoked( raw, {'LT/LT': 25}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RTRT = get_evoked( raw, {'RT/RT': 35}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_contra_TT = concatenate_epochs(epochs_contra_LTLT, epochs_contra_RTRT, 'cued contra TT')
cued_contra_TT = evoked_cued_contra_TT.data[0]
#---------------------------------- Cued/Contra/TT ------------------------------------------#


#---------------------------------- Uncued/Contra/TT ----------------------------------------#
epochs_contra_LTRT = get_evoked( raw, {'LT/RT': 27}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RTLT = get_evoked( raw, {'RT/LT': 33}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_contra_TT = concatenate_epochs(epochs_contra_LTRT, epochs_contra_RTLT, 'uncued contra TT')
uncued_contra_TT = evoked_uncued_contra_TT.data[0]
#---------------------------------- Uncued/Contra/TT ----------------------------------------#


#---------------------------------- Cued/Ipsi/TT ------------------------------------------#
epochs_ipsi_LTLT = get_evoked( raw, {'LT/LT': 25}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RTRT = get_evoked( raw, {'RT/RT': 35}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_ipsi_TT = concatenate_epochs(epochs_ipsi_LTLT, epochs_ipsi_RTRT, 'cued ipsi TT', reversal = True)
cued_ipsi_TT = evoked_cued_ipsi_TT.data[0]
#---------------------------------- Cued/Contra/TT ------------------------------------------#


#---------------------------------- Uncued/Ipsi/TT ----------------------------------------#
epochs_ipsi_LTRT = get_evoked( raw, {'LT/RT': 27}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RTLT = get_evoked( raw, {'RT/LT': 33}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_ipsi_TT = concatenate_epochs(epochs_ipsi_LTRT, epochs_ipsi_RTLT, 'uncued ipsi TT', reversal = True)
uncued_ipsi_TT = evoked_uncued_ipsi_TT.data[0]
#---------------------------------- Uncued/Contra/TT ----------------------------------------#



####################################### Visual/Tactile #######################################

#---------------------------------- Cued/Contra/VT ------------------------------------------#
epochs_contra_LVLT = get_evoked( raw, {'LV/LT': 21}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RVRT = get_evoked( raw, {'RV/RT': 31}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_contra_VT = concatenate_epochs(epochs_contra_LVLT, epochs_contra_RVRT, 'cued contra VT')
cued_contra_VT = evoked_cued_contra_VT.data[0]
#---------------------------------- Cued/Contra/VT ------------------------------------------#


#---------------------------------- Uncued/Contra/VT ----------------------------------------#
epochs_contra_LVRT = get_evoked( raw, {'LV/RT': 23}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RVLT = get_evoked( raw, {'RV/LT': 29}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_contra_VT = concatenate_epochs(epochs_contra_LVRT, epochs_contra_RVLT, 'uncued contra VT')
uncued_contra_VT = evoked_uncued_contra_VT.data[0]
#---------------------------------- Uncued/Contra/VT ----------------------------------------#


#---------------------------------- Cued/Ipsi/VT ------------------------------------------#
epochs_ipsi_LVLT = get_evoked( raw, {'LV/LT': 21}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RVRT = get_evoked( raw, {'RV/RT': 31}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_ipsi_VT = concatenate_epochs(epochs_ipsi_LVLT, epochs_ipsi_RVRT, 'cued ipsi VT', reversal = True)
cued_ipsi_VT = evoked_cued_ipsi_VT.data[0]
#---------------------------------- Cued/Contra/VT ------------------------------------------#


#---------------------------------- Uncued/Ipsi/VT ----------------------------------------#
epochs_ipsi_LVRT = get_evoked( raw, {'LT/RT': 23}, ['Ch25'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RVLT = get_evoked( raw, {'RT/LT': 29}, ['Ch8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_ipsi_VT = concatenate_epochs(epochs_ipsi_LVRT, epochs_ipsi_RVLT, 'uncued ipsi VT', reversal = True)
uncued_ipsi_VT = evoked_uncued_ipsi_VT.data[0]
#---------------------------------- Uncued/Contra/VT ----------------------------------------#



####################################### Tactile/Visual #######################################

#---------------------------------- Cued/Contra/TV ------------------------------------------#
epochs_contra_LTLV = get_evoked( raw, {'LT/LV': 24}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RTRV = get_evoked( raw, {'RT/RV': 34}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_contra_TV = concatenate_epochs(epochs_contra_LTLV, epochs_contra_RTRV, 'cued contra TV')
cued_contra_TV = evoked_cued_contra_TV.data[0]
#---------------------------------- Cued/Contra/TV ------------------------------------------#


#---------------------------------- Uncued/Contra/TV ----------------------------------------#
epochs_contra_LTRV = get_evoked( raw, {'LT/RV': 26}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RTLV = get_evoked( raw, {'RT/LV': 32}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_contra_TV = concatenate_epochs(epochs_contra_LTRV, epochs_contra_RTLV, 'uncued contra TV')
uncued_contra_TV = evoked_uncued_contra_TV.data[0]
#---------------------------------- Uncued/Contra/TV ----------------------------------------#


#---------------------------------- Cued/Ipsi/TV ------------------------------------------#
epochs_ipsi_LTLV = get_evoked( raw, {'LT/LV': 24}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RTRV = get_evoked( raw, {'RT/RV': 34}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_ipsi_TV = concatenate_epochs(epochs_ipsi_LTLV, epochs_ipsi_RTRV, 'cued ipsi TV')
cued_ipsi_TV = evoked_cued_ipsi_TV.data[0]
#---------------------------------- Cued/Contra/TV ------------------------------------------#


#---------------------------------- Uncued/Ipsi/TV ----------------------------------------#
epochs_ipsi_LTRV = get_evoked( raw, {'LT/RV': 26}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RTLV = get_evoked( raw, {'RT/LV': 32}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_ipsi_TV = concatenate_epochs(epochs_ipsi_LTRV, epochs_ipsi_RTLV, 'uncued ipsi TV')
uncued_ipsi_TV = evoked_uncued_ipsi_TV.data[0]
#---------------------------------- Uncued/Contra/TV ----------------------------------------#



####################################### Visual/Visual ########################################

#---------------------------------- Cued/Contra/VV ------------------------------------------#
epochs_contra_LVLV = get_evoked( raw, {'LV/LV': 20}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RVRV = get_evoked( raw, {'RV/RV': 30}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_contra_VV = concatenate_epochs(epochs_contra_LVLV, epochs_contra_RVRV, 'cued contra VV')
cued_contra_VV = evoked_cued_contra_VV.data[0]
#---------------------------------- Cued/Contra/VV ------------------------------------------#


#---------------------------------- Uncued/Contra/VV ----------------------------------------#
epochs_contra_LVRV = get_evoked( raw, {'LV/RV': 22}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_contra_RVLV = get_evoked( raw, {'RV/LV': 28}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_contra_VV = concatenate_epochs(epochs_contra_LVRV, epochs_contra_RVLV, 'uncued contra VV')
uncued_contra_VV = evoked_uncued_contra_VV.data[0]
#---------------------------------- Uncued/Contra/VV ----------------------------------------#


#---------------------------------- Cued/Ipsi/VV ------------------------------------------#
epochs_ipsi_LVLV = get_evoked( raw, {'LV/LV': 20}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RVRV = get_evoked( raw, {'RV/RV': 30}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_cued_ipsi_VV = concatenate_epochs(epochs_ipsi_LVLV, epochs_ipsi_RVRV, 'cued ipsi VV')
cued_ipsi_VV = evoked_cued_ipsi_VV.data[0]
#---------------------------------- Cued/Ipsi/VV ------------------------------------------#


#---------------------------------- Uncued/Ipsi/VV ----------------------------------------#
epochs_ipsi_LVRV = get_evoked( raw, {'LV/RV': 22}, ['Ch51'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
epochs_ipsi_RVLV = get_evoked( raw, {'RV/LV': 28}, ['Ch47'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
evoked_uncued_ipsi_VV = concatenate_epochs(epochs_ipsi_LVRV, epochs_ipsi_RVLV, 'uncued ipsi VV')
uncued_ipsi_VV = evoked_uncued_ipsi_VV.data[0]
#---------------------------------- Uncued/Ipsi/VV ----------------------------------------#




#------------------------------------ Plot Together -----------------------------------------#
X = np.linspace(-200, 500, 701)

f, ax = plt.subplots(2,4, sharex = True, sharey = True)

ax[0,0].plot(X, cued_contra_TT, label = 'cued')
ax[0,0].plot(X, uncued_contra_TT, label = 'uncued')
ax[0,0].axhline(y=0, color = 'black')
ax[0,0].axvline(x=0, linestyle='dashed', color = 'black')
ax[0,0].legend()

ax[0,1].plot(X, cued_ipsi_TT, label = 'cued')
ax[0,1].plot(X, uncued_ipsi_TT, label = 'uncued')
ax[0,1].axhline(y=0, color = 'black')
ax[0,1].axvline(x=0, linestyle='dashed', color = 'black')

ax[1,0].plot(X, cued_contra_VT, label = 'cued')
ax[1,0].plot(X, uncued_contra_VT, label = 'uncued')
ax[1,0].axhline(y=0, color = 'black')
ax[1,0].axvline(x=0, linestyle='dashed', color = 'black')

ax[1,1].plot(X, cued_ipsi_VT, label = 'cued')
ax[1,1].plot(X, uncued_ipsi_VT, label = 'uncued')
ax[1,1].axhline(y=0, color = 'black')
ax[1,1].axvline(x=0, linestyle='dashed', color = 'black')

ax[0,2].plot(X, cued_contra_TV, label = 'cued')
ax[0,2].plot(X, uncued_contra_TV, label = 'uncued')
ax[0,2].axhline(y=0, color = 'black')
ax[0,2].axvline(x=0, linestyle='dashed', color = 'black')

ax[0,3].plot(X, cued_ipsi_TV, label = 'cued')
ax[0,3].plot(X, uncued_ipsi_TV, label = 'uncued')
ax[0,3].axhline(y=0, color = 'black')
ax[0,3].axvline(x=0, linestyle='dashed', color = 'black')

ax[1,2].plot(X, cued_contra_VV, label = 'cued')
ax[1,2].plot(X, uncued_contra_VV, label = 'uncued')
ax[1,2].axhline(y=0, color = 'black')
ax[1,2].axvline(x=0, linestyle='dashed', color = 'black')

ax[1,3].plot(X, cued_ipsi_VV, label = 'cued')
ax[1,3].plot(X, uncued_ipsi_VV, label = 'uncued')
ax[1,3].axhline(y=0, color = 'black')
ax[1,3].axvline(x=0, linestyle='dashed', color = 'black')

plt.show()

# savefig( '/Users/ghislaindentremont/Documents/Multimodal_IOR/condition_averages/condition_averages_%s.png'%participant )
#------------------------------------ Plot Together -----------------------------------------#



##############################################################################################
####                                    Save Dataframe                                    ####
##############################################################################################

dataset = list(
    zip(
        X
        , cued_contra_TT.tolist()
        , uncued_contra_TT.tolist()
        , cued_ipsi_TT.tolist()
        , uncued_ipsi_TT.tolist()
        
        , cued_contra_VT.tolist()
        , uncued_contra_VT.tolist()
        , cued_ipsi_VT.tolist()
        , uncued_ipsi_VT.tolist()
        
        , cued_contra_TV.tolist()
        , uncued_contra_TV.tolist()
        , cued_ipsi_TV.tolist()
        , uncued_ipsi_TV.tolist()      

        , cued_contra_VV.tolist()
        , uncued_contra_VV.tolist()
        , cued_ipsi_VV.tolist()
        , uncued_ipsi_VV.tolist()            
        )
    )

df = pd.DataFrame(
    data = dataset
    , columns = 
        [
        'time (ms)'
        ,'CCTT'
        , 'UCTT'
        , 'CITT'
        , 'UITT'
        
        , 'CCVT'
        , 'UCVT'
        , 'CIVT'
        , 'UIVT'
        
        , 'CCTV'
        , 'UCTV'
        , 'CITV'
        , 'UITV'
        
        , 'CCVV'
        , 'UCVV'
        , 'CIVV'
        , 'UIVV'
        ]
 )

df.to_csv( '/Users/ghislaindentremont/Documents/Multimodal_IOR/condition_averages/condition_averages_%s.csv'%participant )
