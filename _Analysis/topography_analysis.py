import os
import mne
from mne import io
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import fnmatch



##############################################################################################
####                                   Functions                                          ####
##############################################################################################

def get_topo( raw, topo_ids, topo_times, average, tmin, tmax, reject_num, baseline=(-0.1, 0), AR=False ):
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

        f = open('%s/P_AR/%s_%s_AR_output_%s.txt' % (
        directory, list(topo_ids)[0][0:2], list(topo_ids)[0][3:5], participant), 'w')
        # get percentage of epochs dropped
        f.write("Total number of epochs: ")
        f.write(str(len(topo_epochs)))
        f.write("\nPercentage of epochs dropped: ")
        f.write(str(topo_epochs.drop_log_stats()))
        f.close()

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

    raw = io.read_raw_fif(file, add_eeg_ref = False, preload=True)

    directory = '/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Ghis/P_analysis_topo/%s' % participant
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists('%s/P_AR'%directory):
        os.makedirs('%s/P_AR'%directory)
    if not os.path.exists('%s/P_topo'%directory):
        os.makedirs('%s/P_topo'%directory)
    if not os.path.exists('%s/P_averages'%directory):
        os.makedirs('%s/P_averages'%directory)

    # Visualize Events
    events = mne.find_events(raw, stim_channel = 'STI 014', output = 'onset')
    mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, show = False)
    plt.savefig('%s/events_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()


    # define epoch window
    tmin, tmax = -0.2, 0.5



    ##############################################################################################
    ####                                    By Condition                                      ####
    ##############################################################################################


    ####################################### Tactile/Tactile ######################################

    evoked_left_cued_TT = get_topo( raw, {'LT/LT': 25}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6 )
    evoked_right_cued_TT = get_topo(raw, {'RT/RT': 35}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num=100e-6)

    evoked_left_uncued_TT = get_topo( raw, {'RT/LT': 33}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6 )
    evoked_right_uncued_TT = get_topo( raw, {'LT/RT': 27}, np.arange(.060, .120, .010), 0.005, tmin, tmax, reject_num = 100e-6 )

    # get individual effect plots
    ind_left_red_TT = mne.combine_evoked([evoked_left_cued_TT, evoked_left_uncued_TT], weights=[1, -1])  # subtraction
    ind_left_red_TT.plot_topomap([.090], ch_type='eeg', show_names=True, average=0.030)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_left_red_TT_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    # get individual effect plots
    ind_right_red_TT = mne.combine_evoked([evoked_right_cued_TT, evoked_right_uncued_TT],
                                         weights=[1, -1])  # subtraction
    ind_right_red_TT.plot_topomap([.090], ch_type='eeg', show_names=True, average=0.030)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_right_red_TT_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()


    if first:
        sum_left_cued_TT = evoked_left_cued_TT
        sum_right_cued_TT = evoked_right_cued_TT

        sum_left_uncued_TT = evoked_left_uncued_TT
        sum_right_uncued_TT = evoked_right_uncued_TT

    else:
        sum_left_cued_TT = mne.combine_evoked([sum_left_cued_TT, evoked_left_cued_TT], weights=[1, 1])  # sum
        sum_right_cued_TT = mne.combine_evoked([sum_right_cued_TT, evoked_right_cued_TT], weights=[1, 1])  # sum

        sum_left_uncued_TT =  mne.combine_evoked([sum_left_uncued_TT, evoked_left_uncued_TT], weights=[1, 1])  # sum
        sum_right_uncued_TT = mne.combine_evoked([sum_right_uncued_TT, evoked_right_uncued_TT], weights=[1, 1])  # sum


    # get ipsi and contra waveforms
    ERP_contra_cued_TT = (evoked_left_cued_TT.data[24] + evoked_right_cued_TT.data[7]) / 2  # average
    ERP_ipsi_cued_TT = -(evoked_left_cued_TT.data[7] + evoked_right_cued_TT.data[24]) / 2  # average and reverse for ipsi

    ERP_contra_uncued_TT = (evoked_left_uncued_TT.data[24] + evoked_right_uncued_TT.data[7]) / 2  # average
    ERP_ipsi_uncued_TT = -(evoked_left_uncued_TT.data[7] + evoked_right_uncued_TT.data[24]) / 2  # average and reverse for ipsi



    ######################################## Visual/Tactile #######################################

    evoked_left_cued_VT = get_topo(raw, {'LV/LT': 21}, np.arange(.060, .120, .010), 0.005, tmin, tmax,
                                   reject_num=100e-6)
    evoked_right_cued_VT = get_topo(raw, {'RV/RT': 31}, np.arange(.060, .120, .010), 0.005, tmin, tmax,
                                    reject_num=100e-6)

    evoked_left_uncued_VT = get_topo(raw, {'RV/LT': 29}, np.arange(.060, .120, .010), 0.005, tmin, tmax,
                                     reject_num=100e-6)
    evoked_right_uncued_VT = get_topo(raw, {'LV/RT': 23}, np.arange(.060, .120, .010), 0.005, tmin, tmax,
                                      reject_num=100e-6)

    # get individual effect plots
    ind_left_red_VT = mne.combine_evoked([evoked_left_cued_VT, evoked_left_uncued_VT],
                                         weights=[1, -1])  # subtraction
    ind_left_red_VT.plot_topomap([.090], ch_type='eeg', show_names=True, average=0.030)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_left_red_VT_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    # get individual effect plots
    ind_right_red_VT = mne.combine_evoked([evoked_right_cued_VT, evoked_right_uncued_VT],
                                          weights=[1, -1])  # subtraction
    ind_right_red_VT.plot_topomap([.090], ch_type='eeg', show_names=True, average=0.030)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_right_red_VT_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    if first:
        sum_left_cued_VT = evoked_left_cued_VT
        sum_right_cued_VT = evoked_right_cued_VT

        sum_left_uncued_VT = evoked_left_uncued_VT
        sum_right_uncued_VT = evoked_right_uncued_VT

    else:
        sum_left_cued_VT = mne.combine_evoked([sum_left_cued_VT, evoked_left_cued_VT], weights=[1, 1])  # sum
        sum_right_cued_VT = mne.combine_evoked([sum_right_cued_VT, evoked_right_cued_VT], weights=[1, 1])  # sum

        sum_left_uncued_VT = mne.combine_evoked([sum_left_uncued_VT, evoked_left_uncued_VT], weights=[1, 1])  # sum
        sum_right_uncued_VT = mne.combine_evoked([sum_right_uncued_VT, evoked_right_uncued_VT], weights=[1, 1])  # sum


    # get ipsi and contra waveforms
    ERP_contra_cued_VT = (evoked_left_cued_VT.data[24] + evoked_right_cued_VT.data[7]) / 2  # average
    ERP_ipsi_cued_VT = -(evoked_left_cued_VT.data[7] + evoked_right_cued_VT.data[24]) / 2  # average and reverse for ipsi

    ERP_contra_uncued_VT = (evoked_left_uncued_VT.data[24] + evoked_right_uncued_VT.data[7]) / 2  # average
    ERP_ipsi_uncued_VT = -(evoked_left_uncued_VT.data[7] + evoked_right_uncued_VT.data[24]) / 2  # average and reverse for ipsi



    # ####################################### Tactile/Visual #######################################

    evoked_left_cued_TV = get_topo(raw, {'LT/LV': 24}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                   reject_num=100e-6)
    evoked_right_cued_TV = get_topo(raw, {'RT/RV': 34}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                    reject_num=100e-6)

    evoked_left_uncued_TV = get_topo(raw, {'RT/LV': 32}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                     reject_num=100e-6)
    evoked_right_uncued_TV = get_topo(raw, {'LT/RV': 26}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                      reject_num=100e-6)

    # get individual effect plots
    ind_left_red_TV = mne.combine_evoked([evoked_left_cued_TV, evoked_left_uncued_TV],
                                         weights=[1, -1])  # subtraction
    ind_left_red_TV.plot_topomap([.120], ch_type='eeg', show_names=True, average=0.050)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_left_red_TV_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    # get individual effect plots
    ind_right_red_TV = mne.combine_evoked([evoked_right_cued_TV, evoked_right_uncued_TV],
                                          weights=[1, -1])  # subtraction
    ind_right_red_TV.plot_topomap([.120], ch_type='eeg', show_names=True, average=0.050)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_right_red_TV_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    if first:
        sum_left_cued_TV = evoked_left_cued_TV
        sum_right_cued_TV = evoked_right_cued_TV

        sum_left_uncued_TV = evoked_left_uncued_TV
        sum_right_uncued_TV = evoked_right_uncued_TV

    else:
        sum_left_cued_TV = mne.combine_evoked([sum_left_cued_TV, evoked_left_cued_TV], weights=[1, 1])  # sum
        sum_right_cued_TV = mne.combine_evoked([sum_right_cued_TV, evoked_right_cued_TV], weights=[1, 1])  # sum

        sum_left_uncued_TV = mne.combine_evoked([sum_left_uncued_TV, evoked_left_uncued_TV], weights=[1, 1])  # sum
        sum_right_uncued_TV = mne.combine_evoked([sum_right_uncued_TV, evoked_right_uncued_TV], weights=[1, 1])  # sum

    # get ipsi and contra waveforms
    ERP_contra_cued_TV = (evoked_left_cued_TV.data[50] + evoked_right_cued_TV.data[46]) / 2  # average
    ERP_ipsi_cued_TV = (evoked_left_cued_TV.data[46] + evoked_right_cued_TV.data[50]) / 2  # average

    ERP_contra_uncued_TV = (evoked_left_uncued_TV.data[50] + evoked_right_uncued_TV.data[46]) / 2  # average
    ERP_ipsi_uncued_TV = (evoked_left_uncued_TV.data[46] + evoked_right_uncued_TV.data[50]) / 2  # average



    # ####################################### Visual/Visual ########################################

    evoked_left_cued_VV = get_topo(raw, {'LV/LV': 20}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                   reject_num=100e-6)
    evoked_right_cued_VV = get_topo(raw, {'RV/RV': 30}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                    reject_num=100e-6)

    evoked_left_uncued_VV = get_topo(raw, {'RV/LV': 28}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                     reject_num=100e-6)
    evoked_right_uncued_VV = get_topo(raw, {'LV/RV': 22}, np.arange(.070, .180, .010), 0.005, tmin, tmax,
                                      reject_num=100e-6)

    # get individual effect plots
    ind_left_red_VV = mne.combine_evoked([evoked_left_cued_VV, evoked_left_uncued_VV],
                                         weights=[1, -1])  # subtraction
    ind_left_red_VV.plot_topomap([.120], ch_type='eeg', show_names=True, average=0.050)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_left_red_VV_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    # get individual effect plots
    ind_right_red_VV = mne.combine_evoked([evoked_right_cued_VV, evoked_right_uncued_VV],
                                          weights=[1, -1])  # subtraction
    ind_right_red_VV.plot_topomap([.120], ch_type='eeg', show_names=True, average=0.050)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig('%s/P_topo/%s_right_red_VV_topomap.png' % (directory, participant))
    if dont_plot:
        plt.close()

    if first:
        sum_left_cued_VV = evoked_left_cued_VV
        sum_right_cued_VV = evoked_right_cued_VV

        sum_left_uncued_VV = evoked_left_uncued_VV
        sum_right_uncued_VV = evoked_right_uncued_VV

        first = False
    else:
        sum_left_cued_VV = mne.combine_evoked([sum_left_cued_VV, evoked_left_cued_VV], weights=[1, 1])  # sum
        sum_right_cued_VV = mne.combine_evoked([sum_right_cued_VV, evoked_right_cued_VV], weights=[1, 1])  # sum

        sum_left_uncued_VV = mne.combine_evoked([sum_left_uncued_VV, evoked_left_uncued_VV], weights=[1, 1])  # sum
        sum_right_uncued_VV = mne.combine_evoked([sum_right_uncued_VV, evoked_right_uncued_VV], weights=[1, 1])  # sum

    # get ipsi and contra waveforms
    ERP_contra_cued_VV = (evoked_left_cued_VV.data[50] + evoked_right_cued_VV.data[46]) / 2  # average
    ERP_ipsi_cued_VV = (evoked_left_cued_VV.data[46] + evoked_right_cued_VV.data[50]) / 2  # average

    ERP_contra_uncued_VV = (evoked_left_uncued_VV.data[50] + evoked_right_uncued_VV.data[46]) / 2  # average
    ERP_ipsi_uncued_VV = (evoked_left_uncued_VV.data[46] + evoked_right_uncued_VV.data[50]) / 2  # average




    #------------------------------------ Plot Together -----------------------------------------#
    X = np.linspace(-200, 500, 701)

    f, ax = plt.subplots(2,4, sharex = True, sharey = True)

    ax111 = f.add_subplot(111, frameon = False)
    ax111.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax111.set_xlabel("time (ms)")
    ax111.set_ylabel("voltage (V)", labelpad = 15)

    # cue modality labels
    ax211 = f.add_subplot(211, frameon = False)
    ax211.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax211.set_ylabel("tactile cue", rotation = 270, labelpad = 20)
    ax211.yaxis.set_label_position("right")

    ax212 = f.add_subplot(212, frameon = False)
    ax212.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax212.set_ylabel("visual cue", rotation = 270, labelpad = 20)
    ax212.yaxis.set_label_position("right")

    # target modality labels
    ax121 = f.add_subplot(121, frameon = False)
    ax121.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax121.set_title("tactile target", y = 1.05)

    ax122 = f.add_subplot(122, frameon = False)
    ax122.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax122.set_title("visual target", y = 1.05)

    # contralaterality labels
    ax141 = f.add_subplot(141, frameon = False)
    ax141.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax141.set_title("contralateral")

    ax142 = f.add_subplot(142, frameon = False)
    ax142.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax142.set_title("ipsilateral")

    ax143 = f.add_subplot(143, frameon = False)
    ax143.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax143.set_title("contralateral")

    ax144 = f.add_subplot(144, frameon = False)
    ax144.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax144.set_title("ipsilateral")

    ax[0,0].plot(X, ERP_contra_cued_TT, label = 'cued')
    ax[0,0].plot(X, ERP_contra_uncued_TT, label = 'uncued')
    ax[0,0].axhline(y=0, color = 'black')
    ax[0,0].axvline(x=0, linestyle='dashed', color = 'black')
    ax[0,0].legend(prop = {'size':12})

    ax[0,1].plot(X, ERP_ipsi_cued_TT, label = 'cued')
    ax[0,1].plot(X, ERP_ipsi_uncued_TT, label = 'uncued')
    ax[0,1].axhline(y=0, color = 'black')
    ax[0,1].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,0].plot(X, ERP_contra_cued_VT, label = 'cued')
    ax[1,0].plot(X, ERP_contra_uncued_VT, label = 'uncued')
    ax[1,0].axhline(y=0, color = 'black')
    ax[1,0].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,1].plot(X, ERP_ipsi_cued_VT, label = 'cued')
    ax[1,1].plot(X, ERP_ipsi_uncued_VT, label = 'uncued')
    ax[1,1].axhline(y=0, color = 'black')
    ax[1,1].axvline(x=0, linestyle='dashed', color = 'black')

    ax[0,2].plot(X, ERP_contra_cued_TV, label = 'cued')
    ax[0,2].plot(X, ERP_contra_uncued_TV, label = 'uncued')
    ax[0,2].axhline(y=0, color = 'black')
    ax[0,2].axvline(x=0, linestyle='dashed', color = 'black')

    ax[0,3].plot(X, ERP_ipsi_cued_TV, label = 'cued')
    ax[0,3].plot(X, ERP_ipsi_uncued_TV, label = 'uncued')
    ax[0,3].axhline(y=0, color = 'black')
    ax[0,3].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,2].plot(X, ERP_contra_cued_VV, label = 'cued')
    ax[1,2].plot(X, ERP_contra_uncued_VV, label = 'uncued')
    ax[1,2].axhline(y=0, color = 'black')
    ax[1,2].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,3].plot(X, ERP_ipsi_cued_VV, label = 'cued')
    ax[1,3].plot(X, ERP_ipsi_uncued_VV, label = 'uncued')
    ax[1,3].axhline(y=0, color = 'black')
    ax[1,3].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,3].set_ylim(ax[1,3].get_ylim()[::-1])
    ax[1,3].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    f.set_size_inches(18.5, 10.5)

    plt.savefig( '%s/P_averages/condition_averages_%s.png'%(directory, participant) )

    if dont_plot:
        plt.close()
    #------------------------------------ Plot Together -----------------------------------------#




    ##############################################################################################
    ####                                    Save Dataframe                                    ####
    ##############################################################################################

    dataset = list(
        zip(
            X
            , ERP_contra_cued_TT.tolist()
            , ERP_contra_uncued_TT.tolist()
            , ERP_ipsi_cued_TT.tolist()
            , ERP_ipsi_uncued_TT.tolist()

            , ERP_contra_cued_VT.tolist()
            , ERP_contra_uncued_VT.tolist()
            , ERP_ipsi_cued_VT.tolist()
            , ERP_ipsi_uncued_VT.tolist()

            , ERP_contra_cued_TV.tolist()
            , ERP_contra_uncued_TV.tolist()
            , ERP_ipsi_cued_TV.tolist()
            , ERP_ipsi_uncued_TV.tolist()

            , ERP_contra_cued_VV.tolist()
            , ERP_contra_uncued_VV.tolist()
            , ERP_ipsi_cued_VV.tolist()
            , ERP_ipsi_uncued_VV.tolist()
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

    df.to_csv( '%s/P_averages/condition_averages_%s.csv'%(directory, participant) )





##############################################################################################
####                                       Save Files                                     ####
##############################################################################################

directory = '/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Ghis/P_analysis_topo'

if not os.path.exists('%s/grand_sums' % directory):
    os.makedirs('%s/grand_sums' % directory)

sum_left_cued_TT.save('%s/grand_sums/sum_left_cued_TT-ave.fif' % directory)
sum_right_cued_TT.save('%s/grand_sums/sum_right_cued_TT-ave.fif' % directory)
sum_left_uncued_TT.save('%s/grand_sums/sum_left_uncued_TT-ave.fif' % directory)
sum_right_uncued_TT.save('%s/grand_sums/sum_right_uncued_TT-ave.fif'% directory)

sum_left_cued_VT.save('%s/grand_sums/sum_left_cued_VT-ave.fif' % directory)
sum_right_cued_VT.save('%s/grand_sums/sum_right_cued_VT-ave.fif' % directory)
sum_left_uncued_VT.save('%s/grand_sums/sum_left_uncued_VT-ave.fif' % directory)
sum_right_uncued_VT.save('%s/grand_sums/sum_right_uncued_VT-ave.fif' % directory)

sum_left_cued_TV.save('%s/grand_sums/sum_left_cued_TV-ave.fif' % directory)
sum_right_cued_TV.save('%s/grand_sums/sum_right_cued_TV-ave.fif' % directory)
sum_left_uncued_TV.save('%s/grand_sums/sum_left_uncued_TV-ave.fif' % directory)
sum_right_uncued_TV.save('%s/grand_sums/sum_right_uncued_TV-ave.fif' % directory)

sum_left_cued_VV.save('%s/grand_sums/sum_left_cued_VV-ave.fif' % directory)
sum_right_cued_VV.save('%s/grand_sums/sum_right_cued_VV-ave.fif' % directory)
sum_left_uncued_VV.save('%s/grand_sums/sum_left_uncued_VV-ave.fif' % directory)
sum_right_uncued_VV.save('%s/grand_sums/sum_right_uncued_VV-ave.fif' % directory)





