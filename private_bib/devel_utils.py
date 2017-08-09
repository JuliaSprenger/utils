import os
import shutil
import subprocess
import quantities as pq
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from private_bib import developmentio as DIO
from private_bib import utils as ut
from private_bib import spikesort

import elephant


import neo


################################################

def block_to_movie(block,moviedir,movie_name,fps=10,deceleration_factor=4,framewidth=500*pq.ms,dpi=100,format='mp4',ffmpeg_options=None,remove_frames=None,subplot_labels=None):
    """
    SPECIALIZED VERSION OF NEO_VIS.PY FOR PLOTTING DEVEL DATA
    converting analogsignalarrays of first segment of a neo block into a movie using ffmpeg. Different
    analogsignals are centered around their 'electrode_id' annotation value. EpochArrays are visualized by
    shading and EventArrays by vertical lines. Zscoring and subsequent amplitude decrease by a factor of 3 is
    advised for good visibility of multiple analogsignal(array)s.

    Requires: ffmpeg

    Arguments:
        block:      neo block object containing AnalogSignalArrays
        moviedir:   string, directory in which to save the movie
        movie_name: string, filename of the movie
        fps:        integer, number of frames per second to be used in the movie.
                    Default: 10
        deceleration_factor: integer, deceleration factor between realtime and
                    recording time. Default: 4
        framewidth: Quantity, width of the recording time presented per frame.
                    Default: 500*pq.s
        dpi:        integer, resolution of the pngs used for movie creation.
                    Default: 100
        format:     string, output format of the movie. For a list of possible
                    formats check your local ffmpeg installation (ffmpeg -formats).
                    Default: 'mp4'
        ffmpeg_options: string, which is directly passed to ffmpeg.
                    If format='mp4' and ffmpeg=None additional options are used
                    for better quality during movie generation. To avoid this use
                    ffmpeg_options=''.
                    Default None
        remove_frames: string, can be None, 'all', 'png' or 'pdf'. Removes frames
                    of given category after movie generation (for memory reasons)
                    Default: None

        Note:   *This code only supports up to 10^5 frames per movie.
                *If you want to include your movie in a power point presentation
                 it is best to generate an 'avi' movie with
                 '-crf 20 -qscale 10  -vcodec wmv2' option. To adapt file size and
                 quality to your needs change the qscale parameter in the range
                 from 0 to 100 (high quality to low quality).
    """

    # def get(st):
    #     # if 'sp_win_extract' in st.annotations['extraction_params']:
    #     #     t_before, t_after = st.annotations['extraction_params']['sp_win_extract']
    #     # elif 'extract_s_before' in st.annotations['extraction_params']:
    #     #     before = st.annotations['extraction_params']['extract_s_before']
    #     #     after = st.annotations['extraction_params']['extract_s_after']
    #     #     t_before = (-1) * before / st.sampling_rate
    #     #     t_after = after / st.sampling_rate
    #     # else:
    #     #     raise ValueError('Could not extract waveform times from spike train with annotations %s'%st.annotations)
    #     #
    #     # return t_before, t_after
    #
    #     return st.annotations['left_sweep'], st.annotations['right_sweep']

    for s,segment in enumerate(block.segments):
        movie_name_s = movie_name + '_%i'%s
        # Setting up folder structure and paths
        png_folder = moviedir + '/%s'%movie_name_s
        if not os.path.exists(png_folder):
            os.makedirs(png_folder)

        # waveforms = False
        # if len(segment.spiketrains)>0 and sum([len(s.waveforms) for s in segment.spiketrains])>0:
        #     waveforms = True

        # Extracting Analogsignalarrays and basic parameters
        anasigs = segment.analogsignalarrays
        spiketrains = segment.spiketrains
        # tmin = segment.t_start
        # tmax = segment.t_stop
        tmin = min([a.t_start for a in anasigs])
        tmax = max([a.t_stop for a in anasigs])

        # Basic frame parameters
        frames =  int((fps*pq.Hz * (tmax-tmin-framewidth) * deceleration_factor).rescale('dimensionless'))
        overlapping_frames = int((framewidth * frames / (tmax-tmin-framewidth)).rescale('dimensionless'))

        if frames < 0:
            warnings.warn('Framewidth is longer than provided time segment. Only generating single frame.')
            frames = 1

        # Basic settings
        frame_format = 'png'
        num_colors = 17
        #Anasig electrode mapping to visualization
        yorder = range(15,-1,-1)  + range(32,16,-1)
        # el_id  -> baseline location in plot, ls_id
        el_mapping = dict(zip(range(len(yorder)),yorder))
        el_mapping_reverse = {v:k for k,v in el_mapping.iteritems()}

        # Figure generation and basic layout settings
        # fig, axarr = plt.subplots(2, facecolor='white',figsize=(16,9),dpi=dpi,sharex=True)
        fig = plt.figure(facecolor='white',figsize=(16,9),dpi=dpi)

        stim_ax = plt.subplot2grid((10,10), (0,0), colspan=9)
        lfp_ax = plt.subplot2grid((10,10), (1,0), colspan=9, rowspan=9)
        txt_ax = plt.subplot2grid((10,10), (0,9))
        wf_ax = plt.subplot2grid((10,10), (1,9), rowspan=9)

        plt.hold(True)

        stim_ax.set_ylim(-1,2)
        stim_ax.set_ylabel('')
        stim_ax.set_xticks([0,1])
        stim_ax.set_xticklabels(['On/Off','Power'],rotation=90)
        stim_ax.set_yticklabels([])
        stim_ax.grid(which='major', axis='y',zorder=-1)
        if subplot_labels!=None and 'stim_title' in subplot_labels:
            stim_ax.set_title(subplot_labels['stim_title'],**{'fontname':'Arial', 'size':'10'})
        lfp_ax.set_ylim(-2,max([a for a in el_mapping.itervalues()])+2)
        lfp_ax.set_ylabel('HP          Electrode ID          PFC')
        lfp_ax.set_xlabel('Time [ms]')
        lfp_ax.set_yticks([a for a in el_mapping.itervalues()])
        lfp_ax.set_yticklabels(['%i'%a for a in el_mapping.iterkeys()])
        lfp_ax.axhline(16,color='k')
        lfp_ax.grid(which='major', axis='both',zorder=-1)
        if subplot_labels!=None and 'lfp_title' in subplot_labels:
            lfp_ax.set_title(subplot_labels['lfp_title'],**{'fontname':'Arial', 'size':'10'})
        # txt_ax.spines['right'].set_visible(False)
        # txt_ax.spines['top'].set_visible(False)
        txt_ax.set_xticks([])
        txt_ax.set_yticks([])
        wf_ax.set_ylim(-2,max([a for a in el_mapping.itervalues()])+2)
        wf_ax.set_xticks(range(0,31,1))
        wf_ax.set_yticks([a for a in el_mapping.itervalues()])
        wf_ax.set_yticklabels([])
        wf_ax.grid(which='major', axis='both',zorder=-1)
        if subplot_labels!=None and 'waveform_title' in subplot_labels:
            wf_ax.set_title(subplot_labels['waveform_title'],**{'fontname':'Arial', 'size':'10'})
        plt.subplots_adjust(hspace = 0.5)

        #get waveform length
        if len(segment.spiketrains)>0 and sum([len(train.waveforms) for train in segment.spiketrains if train.waveforms!=None])>0:
            wmin = min([train.left_sweep.rescale('ms').magnitude for train in segment.spiketrains])
            wmax = max([train.right_sweep.rescale('ms').magnitude for train in segment.spiketrains])
        else:
            wmin, wmax = -1, 1
        wf_ax.set_xlim(wmin,wmax)
        wf_ax.set_xlabel('Time [ms]')

        # Generating colormap
        colormap = plt.cm.gist_ncar
        lfp_ax.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_colors)])
        unit_colors = [colormap(i) for i in np.linspace(0, 0.9, 8)]

        # Pre-generation of plots
        ls = []
        ws = []
        ss = []
        for i in range(max([u for u in el_mapping.itervalues()])+1):
            l, = lfp_ax.plot([], [], linewidth=0.5)
            ls.append(l)
            w, = wf_ax.plot([],[], linewidth=0.5,zorder=5)
            ws.append(w)
        for i in range(2):
            s, = stim_ax.plot([],[], linewidth=0.5,zorder=4-i)
            ss.append(s)

        # plt.subplots_adjust(left=0.04, bottom=0.06, right=0.97, top=0.97, wspace=0, hspace=0)

        # Global wf scaling factor
        waveforms = False
        if len(segment.spiketrains)>0 and any([s.waveforms.shape[0]>0 for s in segment.spiketrains]):
            waveforms=True
            frac = 0.5
            wf_yscale = (np.std(np.vstack([st.waveforms.rescale('mV').magnitude for st in segment.spiketrains if st.waveforms!=None])))**(-1)*frac
            wf_ax.vlines(-0.3,ymin=-1.9,ymax=-0.9,linewidth=4,color='k')
            wf_ax.text(-0.2, -1.9+0.5*wf_yscale**(-1), '%.3f$mV$'%wf_yscale**(-1), fontdict=None)

        # Adding Spikes
        for st in segment.spiketrains:
            if len(st.times)==0:
                continue
            lfp_ax.plot(st.rescale('ms'),np.ones(len(st))*el_mapping[st.annotations['electrode_id']],'ro',zorder=-1)#,marker=st.annotations['unit_id'])

            # Adding Waveforms
            if st.waveforms!=None and st.waveforms.shape[0]>0:
                wf_times = np.linspace(st.left_sweep.rescale('ms').magnitude,
                                       st.right_sweep.magnitude,
                                       st.waveforms.shape[1])

                for wf_id in range(st.waveforms.shape[0]):
                    wf_ax.plot(wf_times,st.waveforms[wf_id,:].rescale('mV').magnitude*wf_yscale + el_mapping[st.annotations['electrode_id']],color=lighter(unit_colors[st.annotations['unit_id']],0.5),alpha=0.1,zorder=-1,linewidth=0.5)

        # Adding epochs
        for epocharray in segment.epocharrays:
            if 'Stimulation' in epocharray.name:
                for i in range(len(epocharray.times)):
                    lfp_ax.axvspan(epocharray.times[i].rescale('ms'),
                                      (epocharray.times[i]+epocharray.durations[i]).rescale('ms'),
                                      alpha=0.3, color='k')
                    stim_ax.axvspan(epocharray.times[i].rescale('ms'),
                                      (epocharray.times[i]+epocharray.durations[i]).rescale('ms'),
                                      alpha=0.3, color='k')
            elif 'Spindle' in epocharray.name:
                spindle_electrode = el_mapping[epocharray.annotations['channel_id']-1] #pythonic index shift...
                for i in range(len(epocharray.times)):
                    # Spindle box style
                    p_fancy = FancyBboxPatch(((epocharray.times[i]).rescale('ms').magnitude,
                                             spindle_electrode-0.5),
                                             epocharray.durations[i].rescale('ms').magnitude, 1,
                                             boxstyle="round,pad=0.1",
                                             fc=(1., .8, 1.),
                                             ec=(1., 0.5, 1.),
                                             alpha=0.3)
                    lfp_ax.add_patch(p_fancy)

        # Adding Events
        for events in segment.eventarrays:
            for i in range(len(events.times)):
                lfp_ax.axvline(x=events.times[i].rescale('ms'), linewidth=3, color='b')



        # Generating individual frames
        print 'Starting to save individual frames'
        for frame in range(frames):
            # Checking if this frame already exists or needs to be generated
            if os.path.exists(png_folder + "/frame%05i.%s"%(frame,frame_format)):
                print 'Frame %i already exists. Skipping this frame.'%frame
                continue

            #Generating new frame
            print 'Plotting frame %i/%i'%(frame+1,frames)
            plot_present=False

            # Calculating border times of current frame and adjusting frame
            tleft = (float(frame)/frames * (tmax-tmin-framewidth) + tmin).rescale('ms')
            tright = (tleft + framewidth).rescale('ms')
            lfp_ax.set_xlim([tleft,tright])
            lfp_ax.set_xbound([tleft,tright])
            stim_ax.set_xlim([tleft,tright])
            stim_ax.set_xbound([tleft,tright])

            # Plotting all existing analogsignals
            for i,anasig in enumerate(anasigs):
                if tleft > anasig.t_start or tright < anasig.t_stop:
                    sig = anasig.time_slice(tleft,tright)
                else:
                    sig = anasig
                if isinstance(sig,neo.core.AnalogSignalArray):
                    # recording signals
                    if anasig.annotations['electrode_id'] < len(el_mapping):
                        plot_present = True
                        ls[el_mapping[sig.annotations['electrode_id']]].set_data(sig.times.rescale('ms'),
                                                                                 np.asarray(sig) + el_mapping[sig.annotations['electrode_id']])
                    #stimulation signals
                    else:
                        if anasig.annotations['electrode_id']==32:
                            id = 0
                        elif anasig.annotations['electrode_id']==35:
                            id = 1
                        else:
                            raise ValueError('Unknown electrode_id %i'%anasig.annotations['electrode_id'])
                        ss[id].set_data(sig.times.rescale('ms'), np.asarray(sig) + id)

            for i,st in enumerate(spiketrains):
                if st.waveforms!=None and st.waveforms.shape[0]>0:
                    current_ids = np.where(np.logical_and(st>tleft,st<tright))[0]

                    wf_times = np.linspace(st.left_sweep.rescale('ms').magnitude,
                                           st.right_sweep.rescale('ms').magnitude,
                                           st.waveforms.shape[1])
                    el_id = el_mapping[st.annotations['electrode_id']]
                    if len(current_ids)>0:
                        ws[el_id].set_data(wf_times,st.waveforms[current_ids[-1],:].rescale('mV').magnitude * wf_yscale  + el_id)
                        ws[el_id].set_color(unit_colors[st.annotations['unit_id']])
                    # else:
                    #     ws[el_id].set_data([],[])

            # Cancel plotting if no analogsignal is present in current time frame
            if plot_present==False:
                break

            # Saving each frame as png
            fig.savefig(png_folder + "/frame%05i.%s"%(frame,frame_format), facecolor = fig.get_facecolor(),
                        transparent=True, dpi=dpi,  format=frame_format)
            #saving non-overlapping data also as pdf
            if frame%overlapping_frames==0 or frame==frames-1:
                fig.savefig(png_folder + "/frame%05i.%s"%(frame,'pdf'), facecolor = fig.get_facecolor(),
                            bbox_inches='tight', transparent=True,  format='pdf')

        # Using ffmpeg to convert pngs to mp4 or other format movie
        print 'Converting frames to video...'

        # Using advanced ffmpeg options for better mp4 quality
        if format == 'mp4' and ffmpeg_options==None:
                ffmpeg_options = '-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p'
        ffmpeg_options = '' if ffmpeg_options==None else ffmpeg_options

        command =\
        """
        cd %s
        ffmpeg -framerate %i -i frame%%05d.%s %s %s.%s
        """%(moviedir + '/' + movie_name_s + '/',fps, frame_format,ffmpeg_options, '../' + movie_name_s,format)
        print command
        subprocess.call(command,shell=True)

        # Removing frames to free memory
        if remove_frames != None:
            print 'Removing %s source frames.'%remove_frames
            if remove_frames == 'all':
                # deleting frame source folder
                shutil.rmtree(os.path.join(moviedir,movie_name_s))
            else:
                # selecting files for deletion
                files = [os.path.join(moviedir,movie_name_s,f) for f in os.listdir(os.path.join(moviedir,movie_name_s))
                          if os.path.isfile(os.path.join(moviedir,movie_name_s,f)) and f.endswith('.' + remove_frames)]
                # deleting files
                for f in files:
                    os.remove(f)


def make_movie(sessiondir,savefigdir,mdatadir=None,t_starts=None,t_stops=None,**kwargs):

    session = sessiondir.split('/')[-1]

    IO = DIO.DevelopmentIO(sessiondir=sessiondir,mdatadir=mdatadir,cachedir=ut.get_cache_location(session),use_cache='datesize',print_diagnostic=True)

    if 'stim_ids' in kwargs and t_starts==None and t_stops==None:
        print 'Loading times based on stimulation epochs requested'
        pre,post = kwargs['ep_shift'] if 'ep_shift' in kwargs else [0*pq.s,0*pq.s]
        stim_ids = kwargs['stim_ids']
        t1s, t2s = [],[]
        if IO.odML_doc == None:
            stim_sections = None
            raise ValueError('This session has no associated odml. No stimulation times available.')
        else:
            stim_sections = [s for s in IO.odML_doc.find_related('Stimulations',type='recordings').sections if s.name.startswith('Stimulation_')]

        for stim_sec in stim_sections:
            if stim_ids == [] or int(stim_sec.name.replace('Stimulation_','')) in stim_ids:
                t1_val = stim_sec.properties['StimulationPeriodStart'].value
                t1s.append(pq.Quantity(t1_val.data,t1_val.unit) - IO.parameters_global['t_start'] + pre)
                t2_val = stim_sec.properties['StimulationPeriodEnd'].value
                t2s.append(pq.Quantity(t2_val.data,t1_val.unit) - IO.parameters_global['t_start'] + post)

        t_starts, t_stops = t1s, t2s
        kwargs.pop('stim_ids')
        if 'ep_shift'in kwargs:
            kwargs.pop('ep_shift')
        movie_prefix = session + ('_stims%s'%(stim_ids)).replace(' ','')

    elif 'spindle_ids' in kwargs and t_starts==None and t_stops==None:
        print 'Loading times based on spindle epochs requested'
        pre,post = kwargs['ep_shift'] if 'ep_shift' in kwargs else [0*pq.s,0*pq.s]
        channel_ids, stim_ids = kwargs['spindle_ids']
        t1s, t2s = [],[]
        stim_sections = [s for s in IO.odML_doc.itersections() if s.name.startswith('Spindle_')]
        for stim_sec in stim_sections:
            if (stim_ids == [] or int(stim_sec.name.replace('Spindle_','')) in stim_ids) and \
                    (channel_ids == [] or int(stim_sec.parent.name.replace('Channel_','')) in channel_ids):
                t1_val = stim_sec.properties['SpindleStart'].value
                t1s.append(pq.Quantity(t1_val.data,t1_val.unit) - IO.parameters_global['t_start'] + pre)
                t2_val = stim_sec.properties['SpindleEnd'].value
                t2s.append(pq.Quantity(t2_val.data,t1_val.unit) - IO.parameters_global['t_start'] + post)

        t_starts, t_stops = t1s, t2s
        kwargs.pop('spindle_ids')
        if 'ep_shift'in kwargs:
            kwargs.pop('ep_shift')
        movie_prefix = session + ('_spindle%s'%(stim_ids)).replace(' ','')
    else:
        movie_prefix = session + ('_%s-%s'%(t_starts,t_stops)).replace(' ','')


    print 'Loading data...'
    block = IO.read_block(t_starts=t_starts,t_stops=t_stops,electrode_list=[],
                          unit_list=None)
    print 'Completed data loading.'

    if 'spikesort' in kwargs:
        if kwargs['spikesort']:
            # loading spike sorting data
            try:
                print 'Trying to load sorted spikes.'
                print('WARNING, not using saved data. Correct this here!')
                raise ValueError()
                block = spikesort.load_spikesort(block,sessiondir=sessiondir,
                                                 sortdir= ut.get_sort_location(session),
                                                 mdatadir= ut.get_metadata_location(session),
                                                 cachedir= ut.get_cache_location(session))
            except ValueError:
                print 'Generating new spike extraction file...'
                spikesort.spikesorting(sessiondir=sessiondir)
                block = spikesort.load_spikesort(block,sessiondir=sessiondir,
                                                 sortdir= ut.get_sort_location(session),
                                                 mdatadir= ut.get_metadata_location(session),
                                                 cachedir= ut.get_cache_location(session))

            # print 'Spike sorting..'
            # jelephant.analysis.spike_sorting.generate_spiketrains(block,'k_means_plus',waveforms=True,sort=True,extraction_dict={},sorting_dict={})
            # # jelephant.analysis.spike_sorting.spikesorting2hdf5(savefigdir,block,'k_means_plus',waveforms=False,sort=True,extraction_dict={},sorting_dict={})
        kwargs.pop('spikesort')

    post_label = '_'

    subplot_labels ={'lfp_title': 'LFP','stim_title':'Stimulation / Laser Control Signals','waveform_title':'Waveforms'}
    if 'highpass' in kwargs and 'lowpass' in kwargs and (kwargs['highpass']!=None or kwargs['lowpass']!=None):
        subplot_labels['lfp_title'] +=  ' %s - %s'%(kwargs['highpass'],kwargs['lowpass'])

    for s in range(len(block.segments)):
        anasigs = block.segments[s].analogsignalarrays
        rcgs = [a.recordingchannelgroup for a in anasigs]
        if 'highpass' in kwargs and 'lowpass' in kwargs and (kwargs['highpass']!=None or kwargs['lowpass']!=None):
            print 'Filtering signals'
            post_label += 'filter%s-%s'%(kwargs['highpass'],kwargs['lowpass'])
            for i in range(len(anasigs)):
                # skipping analog signals from laser
                if anasigs[i].annotations['electrode_id']>=32:
                    continue
                anasigs[i] = elephant.signal_processing.butter(anasigs[i], highpass_freq=kwargs['highpass'], lowpass_freq=kwargs['lowpass'], order=3, filter_function='filtfilt', fs=1.0, axis=-1)

        print 'Zscoring data...'
        for i in range(len(anasigs)):
            if anasigs[i].annotations['electrode_id']<32:
                elephant.signal_processing.zscore(anasigs[i])
                anasigs[i] = anasigs[i]/3.
            else:
                rescaleda = anasigs[i]-np.mean(anasigs[i])
                anasigs[i] = (rescaleda)/np.max(np.abs(rescaleda))

            anasigs[i].segment = block.segments[s]
            anasigs[i].recordingchannelgroup = rcgs[i]

    for kw in ['highpass','lowpass']:
        if kw in kwargs:
            kwargs.pop(kw)


    movie_name = ("%s_fps%i_decfactor%i_dpi%i_frame%s%s"%(movie_prefix,kwargs['fps'],kwargs['deceleration_factor'],kwargs['dpi'],kwargs['framewidth'],post_label)).replace(' ','')
    print 'Generating movie...'
    block_to_movie(block,savefigdir,movie_name, subplot_labels = subplot_labels,**kwargs)

######### supplementory functions ################################
def lighter(color, percent):
    '''assumes color is rgb between (0, 0, 0) and (255, 255, 255)'''
    color = np.array(color)
    white = np.array([1, 1, 1, 1])
    vector = white-color
    vector[-1] = 0
    return color + vector * percent



