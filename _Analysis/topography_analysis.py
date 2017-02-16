import os
import mne
from mne import io
import numpy as np
from matplotlib import pyplot as plt
import fnmatch



##############################################################################################
####                                   Functions                                          ####
##############################################################################################

def get_topo( raw, topo_ids, topo_times, average, tmin, tmax, reject_num, baseline=(-0.1, 0), AR=True ):
    topo_reject = dict(eeg=reject_num)
    topo_picks = mne.pick_types(
        raw.info
        , meg=False
        , eeg=True
        , eog=False
    )
    topo_params = dict(
        picks=topo_picks
        , events=events  # global variable
        , event_id=topo_ids
        , tmin=tmin
        , tmax=tmax
    )
    topo_epochs = mne.Epochs(
        raw
        , **topo_params
        , add_eeg_ref=False
        , baseline=baseline
    )
    # artifact rejection
    if AR:
        topo_epochs.drop_bad(reject=topo_reject)

    topo_evoked = topo_epochs.average()

    return (topo_evoked)



# Participant
participants = ['e02', 'e12', 'e16', 'e17', 'e20', 'e22', 'e27', 'e03', 'e04'
    , 'e05', 'e06', 'p06', 'e40', 'e41', 'e42', 'e44', 'e45', 'e46', 'e47', 'e48', 'e49']

os.chdir('/Volumes/LaCie/Experiments/MMIOR/Ghis/Topo')

dont_plot = True
first = True

for participant in participants:

    file = None  # to ensure that an error comes up later if a file cannot be found
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if fnmatch.fnmatch(name, '*%s_raw.fif'%participant):
                file = os.path.join(root, name)

    directory = '/Volumes/LaCie/Experiments/MMIOR/Ghis/P_analysis_topo/%s' % participant

    raw = io.read_raw_fif(file, add_eeg_ref = False, preload=True)


    # Visualize Events
    events = mne.find_events(raw, stim_channel = 'STI 014', output = 'onset')
    mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, show = False)
    plt.savefig('%s/P_topo/events_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()


    # define epoch window
    tmin, tmax = -0.2, 0.5



    ##############################################################################################
    ####                                    By Condition                                      ####
    ##############################################################################################


    ####################################### Tactile/Tactile ######################################

    evoked_left_cued_TT = get_topo( raw, {'LT/LT': 25}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    evoked_right_cued_TT = get_topo(raw, {'RT/RT': 35}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num=100e-6, AR=False)

    evoked_left_uncued_TT = get_topo( raw, {'RT/LT': 33}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    evoked_right_uncued_TT = get_topo( raw, {'LT/RT': 27}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )

    if first:
        sum_left_cued_TT = evoked_left_cued_TT
        sum_right_cued_TT = evoked_right_cued_TT

        sum_left_uncued_TT = evoked_left_uncued_TT
        sum_right_uncued_TT = evoked_right_uncued_TT

        first = False
    else:
        sum_left_cued_TT = mne.combine_evoked([sum_left_cued_TT, evoked_left_cued_TT], weights=[1, 1])  # sum
        sum_right_cued_TT = mne.combine_evoked([sum_right_cued_TT, evoked_right_cued_TT], weights=[1, 1])  # sum

        sum_left_uncued_TT =  mne.combine_evoked([sum_left_uncued_TT, evoked_left_uncued_TT], weights=[1, 1])  # sum
        sum_right_uncued_TT = mne.combine_evoked([sum_right_uncued_TT, evoked_right_uncued_TT], weights=[1, 1])  # sum




    # ####################################### Visual/Tactile #######################################
    #
    # #---------------------------------- Cued/VT ------------------------------------------#
    # evoked_cued_VT = get_topo( raw, {'LV/LT': 21, 'RV/RT': 31}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    #
    #
    #
    # #---------------------------------- Uncued/VT ----------------------------------------#
    # evoked_uncued_VT = get_topo( raw, {'LV/RT': 23, 'RV/LT': 29}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )





    # ####################################### Tactile/Visual #######################################
    #
    # #---------------------------------- Cued/Contra/TV ------------------------------------------#
    # epochs_contra_LTLV = get_topo( raw, {'LT/LV': 24}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_contra_RTRV = get_topo( raw, {'RT/RV': 34}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_cued_contra_TV = concatenate_epochs(epochs_contra_LTLV, epochs_contra_RTRV, 'cued contra TV')
    # cued_contra_TV = evoked_cued_contra_TV.data[0]
    # #---------------------------------- Cued/Contra/TV ------------------------------------------#
    #
    #
    # #---------------------------------- Uncued/Contra/TV ----------------------------------------#
    # epochs_contra_LTRV = get_topo( raw, {'LT/RV': 26}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_contra_RTLV = get_topo( raw, {'RT/LV': 32}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_uncued_contra_TV = concatenate_epochs(epochs_contra_LTRV, epochs_contra_RTLV, 'uncued contra TV')
    # uncued_contra_TV = evoked_uncued_contra_TV.data[0]
    # #---------------------------------- Uncued/Contra/TV ----------------------------------------#
    #
    #
    # #---------------------------------- Cued/Ipsi/TV ------------------------------------------#
    # epochs_ipsi_LTLV = get_topo( raw, {'LT/LV': 24}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_ipsi_RTRV = get_topo( raw, {'RT/RV': 34}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_cued_ipsi_TV = concatenate_epochs(epochs_ipsi_LTLV, epochs_ipsi_RTRV, 'cued ipsi TV')
    # cued_ipsi_TV = evoked_cued_ipsi_TV.data[0]
    # #---------------------------------- Cued/Contra/TV ------------------------------------------#
    #
    #
    # #---------------------------------- Uncued/Ipsi/TV ----------------------------------------#
    # epochs_ipsi_LTRV = get_topo( raw, {'LT/RV': 26}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_ipsi_RTLV = get_topo( raw, {'RT/LV': 32}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_uncued_ipsi_TV = concatenate_epochs(epochs_ipsi_LTRV, epochs_ipsi_RTLV, 'uncued ipsi TV')
    # uncued_ipsi_TV = evoked_uncued_ipsi_TV.data[0]
    # #---------------------------------- Uncued/Contra/TV ----------------------------------------#
    #
    #
    #
    # ####################################### Visual/Visual ########################################
    #
    # #---------------------------------- Cued/Contra/VV ------------------------------------------#
    # epochs_contra_LVLV = get_topo( raw, {'LV/LV': 20}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_contra_RVRV = get_topo( raw, {'RV/RV': 30}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_cued_contra_VV = concatenate_epochs(epochs_contra_LVLV, epochs_contra_RVRV, 'cued contra VV')
    # cued_contra_VV = evoked_cued_contra_VV.data[0]
    # #---------------------------------- Cued/Contra/VV ------------------------------------------#
    #
    #
    # #---------------------------------- Uncued/Contra/VV ----------------------------------------#
    # epochs_contra_LVRV = get_topo( raw, {'LV/RV': 22}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_contra_RVLV = get_topo( raw, {'RV/LV': 28}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_uncued_contra_VV = concatenate_epochs(epochs_contra_LVRV, epochs_contra_RVLV, 'uncued contra VV')
    # uncued_contra_VV = evoked_uncued_contra_VV.data[0]
    # #---------------------------------- Uncued/Contra/VV ----------------------------------------#
    #
    #
    # #---------------------------------- Cued/Ipsi/VV ------------------------------------------#
    # epochs_ipsi_LVLV = get_topo( raw, {'LV/LV': 20}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_ipsi_RVRV = get_topo( raw, {'RV/RV': 30}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_cued_ipsi_VV = concatenate_epochs(epochs_ipsi_LVLV, epochs_ipsi_RVRV, 'cued ipsi VV')
    # cued_ipsi_VV = evoked_cued_ipsi_VV.data[0]
    # #---------------------------------- Cued/Ipsi/VV ------------------------------------------#
    #
    #
    # #---------------------------------- Uncued/Ipsi/VV ----------------------------------------#
    # epochs_ipsi_LVRV = get_topo( raw, {'LV/RV': 22}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # epochs_ipsi_RVLV = get_topo( raw, {'RV/LV': 28}, np.arange(.070, .170, .010), 0.005, tmin, tmax, reject_num = 100e-6, AR = False )
    # evoked_uncued_ipsi_VV = concatenate_epochs(epochs_ipsi_LVRV, epochs_ipsi_RVLV, 'uncued ipsi VV')
    # uncued_ipsi_VV = evoked_uncued_ipsi_VV.data[0]
    # #---------------------------------- Uncued/Ipsi/VV ----------------------------------------#
    #
    #




##############################################################################################
####                                    Topographies                                      ####
##############################################################################################

directory = '/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Ghis/P_analysis_topo'


####################################### Tactile/Tactile ######################################

# get averages out of sums
avg_left_cued_TT = mne.combine_evoked([sum_left_cued_TT], weights=[1/len(participants)])
avg_right_cued_TT = mne.combine_evoked([sum_right_cued_TT], weights=[1/len(participants)])

avg_left_uncued_TT = mne.combine_evoked([sum_left_uncued_TT], weights=[1/len(participants)])
avg_right_uncued_TT = mne.combine_evoked([sum_right_uncued_TT], weights=[1/len(participants)])

# get ipsi and contra waveforms
avg_contra_cued_TT = (avg_left_cued_TT.data[24] + avg_right_cued_TT.data[7])/2  # average
avg_ipsi_cued_TT = -(avg_left_cued_TT.data[7] + avg_right_cued_TT.data[24])/2  # average and reverse for ipsi

avg_contra_uncued_TT = (avg_left_uncued_TT.data[24] + avg_right_uncued_TT.data[7])/2  # average
avg_ipsi_uncued_TT = -(avg_left_uncued_TT.data[7] + avg_right_uncued_TT.data[24])/2  # average and reverse for ipsi


#----------------------------------- Left Reduction ----------------------------------------#
avg_left_red_TT  = mne.combine_evoked([avg_left_cued_TT, avg_left_uncued_TT], weights=[1, -1])  # subtraction

avg_left_red_TT.plot_topomap([.090], ch_type='eeg', show_names=True, average= 0.030)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/left_red_TT_avg_topomap.png' % directory)
if dont_plot:
    plt.close()

avg_left_red_TT.plot_topomap(np.arange(.060, .130, .010), ch_type='eeg', show_names=True, average= 0.005)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/left_red_TT_timeseries_topomap.png' % directory)
if dont_plot:
    plt.close()

# INTERPOLATE
avg_left_red_TT.info['bads'] = ['T7']
avg_left_red_TT.interpolate_bads(reset_bads=True)

avg_left_red_TT.plot_topomap([.090], ch_type='eeg', show_names=True, average= 0.030)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/left_red_TT_interp_avg_topomap.png' % directory)
if dont_plot:
    plt.close()

avg_left_red_TT.plot_topomap(np.arange(.060, .130, .010), ch_type='eeg', show_names=True, average= 0.005)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/left_red_TT_interp_timeseries_topomap.png' % directory)
if dont_plot:
    plt.close()
#----------------------------------- Left Reduction ----------------------------------------#


#---------------------------------- Right Reduction ----------------------------------------#
avg_right_red_TT = mne.combine_evoked([avg_right_cued_TT, avg_right_uncued_TT], weights=[1, -1])  # subtraction

avg_right_red_TT.plot_topomap([.090], ch_type='eeg', show_names=True, average= 0.030)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/right_red_TT_avg_topomap.png' % directory)
if dont_plot:
    plt.close()

avg_right_red_TT.plot_topomap(np.arange(.060, .130, .010), ch_type='eeg', show_names=True, average= 0.005)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/right_red_TT_timeseries_topomap.png' % directory)
if dont_plot:
    plt.close()

# INTERPOLATE
avg_right_red_TT.info['bads'] = ['T7']
avg_right_red_TT.interpolate_bads(reset_bads=True)

avg_right_red_TT.plot_topomap([.090], ch_type='eeg', show_names=True, average= 0.030)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/right_red_TT_interp_avg_topomap.png' % directory)
if dont_plot:
    plt.close()

avg_right_red_TT.plot_topomap(np.arange(.060, .130, .010), ch_type='eeg', show_names=True, average= 0.005)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig('%s/right_red_TT_interp_timeseries_topomap.png' % directory)
if dont_plot:
    plt.close()
#---------------------------------- Right Reduction ----------------------------------------#




##############################################################################################
####                                               ERPs                                   ####
##############################################################################################

X = np.linspace(-200, 500, 701)

f, ax = plt.subplots(1, 2, sharex=True, sharey=True)

ax[0].plot(X, avg_contra_cued_TT, label='cued')
ax[0].plot(X, avg_contra_uncued_TT, label='uncued')
ax[0].axhline(y=0, color='black')
ax[0].axvline(x=0, linestyle='dashed', color='black')
ax[0].legend(prop={'size': 12})

ax[1].plot(X, avg_ipsi_cued_TT, label='cued')
ax[1].plot(X, avg_ipsi_uncued_TT, label='uncued')
ax[1].axhline(y=0, color='black')
ax[1].axvline(x=0, linestyle='dashed', color='black')

f.set_size_inches(18.5, 10.5)

plt.gca().invert_yaxis()

plt.savefig('%s/condition_averages.png' % directory)

if dont_plot:
    plt.close()



# evoked_left_cued_TT = get_topo(raw, {'LT/LT': 25}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num=100e-6,
#                                AR=False)
# evoked_right_cued_TT = get_topo(raw, {'RT/RT': 35}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num=100e-6,
#                                 AR=False)
#
# evoked_left_uncued_TT = get_topo(raw, {'RT/LT': 33}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num=100e-6,
#                                  AR=False)
# evoked_right_uncued_TT = get_topo(raw, {'LT/RT': 27}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num=100e-6,
#                                   AR=False)
#
# evoked_left_red_TT = mne.combine_evoked([evoked_left_cued_TT, evoked_left_uncued_TT], weights=[1, -1])  # subtraction
# evoked_right_red_TT = mne.combine_evoked([evoked_right_cued_TT, evoked_right_uncued_TT], weights=[1, -1])  # subtraction
#
# evoked_cued_TT = mne.combine_evoked([evoked_left_cued_TT, evoked_right_cued_TT], weights=[0.5, 0.5])  # average
# evoked_uncued_TT = mne.combine_evoked([evoked_left_uncued_TT, evoked_right_uncued_TT], weights=[0.5, 0.5])  # average
#
# evoked_red_TT = mne.combine_evoked([evoked_cued_TT, evoked_uncued_TT], weights=[1, -1])  # subtraction



# # left tactile
# evoked_cued_TT.info['ch_names'][7]
#
# # right tactile
# evoked_cued_TT.info['ch_names'][24]
#
# # left visual
# evoked_cued_TT.info['ch_names'][46]
#
# # right visual
# evoked_cued_TT.info['ch_names'][50]
#
# evoked_cued_m_uncued_TT.plot(np.array([7]))
#
# evoked_cued_m_uncued_TT.plot_topomap(np.arange(.060, .120, .010), ch_type='eeg', show_names=True, average= 0.005)
# plt.savefig('%s/P_topo/TT_reduction_topomap_%s.png' % (directory, participant))
# if dont_plot:
#     plt.close()
#
# evoked_cued_m_uncued_TT.plot_topo()
# plt.savefig('%s/P_topo/TT_reduction_topoplot_%s.png' % (directory, participant))
# if dont_plot:
#     plt.close()