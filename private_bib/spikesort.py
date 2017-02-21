import os
import quantities as pq
import numpy as np
import neo
import elephant
import private_bib.developmentio as DIO


set = {'sessiondir':'/home/julia/data/SPP/data/2015-02-03_14-22-50',
       # 'sessiondir': '/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/data/2015-02-20_10-19-03',
       'sortdir': '/home/julia/projects/SPP1665/analysis/spikesortings',
       'mdatadir': '/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/metadata',
       'cachedir': '/home/julia/cache'}

software = 'phy'
method = ''
waveforms = True
sort = True
parameter_dict = {'keep_temp_files':True,'num_cpus': 2,'threshold_strong_std_factor':6,'threshold_weak_std_factor':3}#{'threshold': 5.0}
all_channels_at_once = True
starts = 0*pq.s
stops = 1*pq.min



def spikesorting(**kwargs):

    set.update(kwargs)

    session = set['sessiondir'].split('/')[-1]

    IO = DIO.DevelopmentIO(sessiondir=set['sessiondir'], mdatadir=set['mdatadir'], cachedir=set['cachedir'],
                           use_cache='datesize', print_diagnostic=True)

    # Name generation for spike sorting file
    sorting_prefix = '%s_%s' % (session, software)
    sorting_prefix += str(parameter_dict)
    if waveforms:
        sorting_prefix += '_wf'

    block = None
    channel_index = None

    all_at_once = all_channels_at_once
    t_start = starts
    t_stop = stops

    if 'all_at_once' in kwargs:
        all_at_once = kwargs['all_at_once']
        kwargs.pop('all_at_once')
    if 't_start' in kwargs:
        t_start = kwargs['t_start']
        kwargs.pop('t_start')
    if 't_stop' in kwargs:
        t_stop = kwargs['t_stop']
        kwargs.pop('t_start')


    if all_at_once:
        iterations = 1
        electrode_list = range(32)
    else:
        iterations = 32
        electrode_list = None
    # Loading only one anasig at a time to save memory
    for electrode in range(iterations):
        print 'Loading data of electrode %i' % electrode
        t_starts = t_start
        t_stops = t_stop
        electrode_list = electrode if electrode_list==None else electrode_list
        anasig_block = IO.read_block(t_starts=t_starts, t_stops=t_stops,
                                     electrode_list=[electrode], units=None)
        print 'Completed data loading.'

        if block is None:
            block = anasig_block
            if len(block.segments[0].analogsignalarrays) == 0:
                continue
            channel_index = anasig_block.segments[0].analogsignalarrays[0].channel_index
        else:
            for s in range(len(block.segments)):
                block.segments[s].analogsignalarrays = anasig_block.segments[s].analogsignalarrays
                for a in range(len(block.segments[s].analogsignalarrays)):
                    block.segments[s].analogsignalarrays[a].segment = block.segments[s]
                    block.segments[s].analogsignalarrays[a].channel_index = channel_index

        elephant.spike_sorting.generate_spiketrains(block,'phy',
                                                              waveforms=waveforms,
                                                              sort=sort,
                                                              parameter_dict=parameter_dict)
    save_spikesorting(os.path.join(set['sortdir'], session), block)





def save_spikesorting(savedir,block):
    """
    :param savedir:
    :param block:
    :return:
    """

    filename = savedir + '_spikesorting.hdf5'
    print filename
    nix_file = neo.nixio.NixIO(filename,'rw')

    chidx = neo.ChannelIndex(0, name='spike sorting')
    chidxannotations = {}
    # Extracting only spikes with fitting spike extraction and sorting parameters as the first unit and attach them
    # channelindex
    units = block.list_units
    for unit in units:
        if unit.annotations['parameters'] == units[0].annotations['parameters'] \
            and unit.annotations['sorted'] == units[0].annotations['sorted']:

            for key in ['parameters','sorted']:
                if key not in chidxannotations:
                    chidxannotations[key] = unit.annotations[key]
                elif chidxannotations[key] != unit.annotations[key]:
                    raise ValueError('Block contains units of different sortings. Can not be saved into one recording channel group')

            chidx.units.append(unit)

    chidx.annotations.update(chidxannotations)

    #check if there is already a block with the same annotation and units
    nix_block = nix_file.read_block()

    if nix_block is None:
        nix_block = neo.Block(name='spike sorting block')
    for chidx in nix_block.channel_indexes:
        if chidx.annotations == chidxannotations:
            # check if units are overlapping
            overlapping_units = np.intersect1d(chidx.channel_indexes,
                                               [u.spiketrains[0].annotations[
                                                    'channel_index'] for u in
                                                chidx.units])
            if overlapping_units:
                raise ValueError('There is already a block with these '
                                 'annotations and units %s saved in the '
                                 'hdf5. Not saving new spikesorting '
                                 '(Annotations: %s)'%(overlapping_units,
                                                      chidxannotations))
            else:
                chidx.units.extend(chidx.units)
                rcg = chidx

    nix_file.file.close()


def load_spikesort(block,session,sorting_dir,software,parameter_dict,
                   sort=True, electrode_list='all'):
    """
    :param load_from:
    :param block:
    :return:
    """

    corr_annotations = {'sorted':sort,'parameters':elephant.spike_sorting.get_updated_parameters(software, parameter_dict)}
    # Warning this needs to be reset to line above!!!!!!!!!!!!!!!!!!!!!!!!!!
    # corr_annotations ={'sorted': False, 'parameters': {'filter_low': None,
    #                                                    'n_post': 10,
    #                                                    'filter_high':np.array(400.0) * pq.Hz, 'threshold': -3.5, 'n_pre': -10, 'alignment': 'min'}}
        # {'extraction_params': {'edge': 'falling',
        #       'filter': [np.array(500.0) * pq.Hz, None],
        #       'filter_order': 4,
        #       'sp_win_align': [np.array(-1.0) * pq.ms, np.array(1.0) * pq.ms],
        #       'sp_win_extract': [np.array(-0.5) * pq.ms, np.array(1.5) * pq.ms],
        #       'threshold': 5.0},
        #      'sorted': True,
        #      'sorting_params': {'method': 'k_means_plus', 'ncomps': 2, 'num_units': 3}}

    filename = os.path.join(sorting_dir,session + '_spikesorting.hdf5')
    if not os.path.exists(filename):
        raise IOError('File does not exist!')
    print 'Loading extracted spikes from %s'%filename
    nix_file = neo.io.NixIO(filename,'ro')

    chidx = None
    nix_block = nix_file.read_block()
    nix_file.file.close()
    channel_indexes = []
    if nix_block is not None:
        channel_indexes = nix_block.channel_indexes
    for curr_chidx in channel_indexes:
        # print rcg_x.annotations
        if curr_chidx.annotations == corr_annotations:
            chidx = curr_chidx
            break

    if chidx == None:
        raise ValueError('No spike sorting with parameters %s found.'%corr_annotations)

    block.channel_indexes.append(chidx)
    for seg in block.segments:
        tstart, tstop = seg.t_start, seg.t_stop
        selected_sts = []
        for unit in chidx.units:
            if (electrode_list == 'all' or unit.spiketrains[0].annotations[
                'channel_index'] in electrode_list):
                for spiketrain in unit.spiketrains:
                    # loading only spikes in current segment
                    st = spiketrain.time_slice(tstart,tstop)
                    # if len(st) > 0:
                    selected_sts.append(st)
                        #seg.spiketrains.append(st)

        if electrode_list != 'all':
            missing_els = np.setdiff1d(electrode_list,
                                       [st.annotations['channel_index'] for
                                        st in selected_sts])
            if len(missing_els):
                raise ValueError('Can not load spiketrain for electrode ids '
                                 '%s'%missing_els)

        seg.spiketrains.extend(selected_sts)

    return block


# deprecated
# def spikesorting2hdf5(savedir,block):
#     """
#     :param savedir:
#     :param block:
#     :return:
#     """
#
#     filename = savedir + '_spikesorting.hdf5'
#     print filename
#     Nhdf5 = neo.io.NeoHdf5IO(filename)
#
#     rcg = neo.core.RecordingChannelGroup(name='spike sorting')
#     rcgannotations = {}
#     # Extracting only spikes with fitting spike extraction and sorting parameters as the first unit and attach them
#     # recordingchannelgroup
#     units = block.list_units
#     for unit in units:
#         if unit.annotations['parameters'] == units[0].annotations['parameters'] \
#             and unit.annotations['sorted'] == units[0].annotations['sorted']:
#
#             for key in ['parameters','sorted']:
#                 if key not in rcgannotations:
#                     rcgannotations[key] = unit.annotations[key]
#                 elif rcgannotations[key] != unit.annotations[key]:
#                     raise ValueError('Block contains units of different sortings. Can not be saved into one recording channel group')
#
#             rcg.units.append(unit)
#
#     rcg.channel_indexes = np.unique(np.array([st.annotations['channel_index']
#                                               for u in rcg.units for st in
#                                               u.spiketrains]))
#     rcg.annotations.update(rcgannotations)
#
#     #check if there is already a block with the same annotation and units
#     for rcg_id in range(Nhdf5.get_info()['RecordingChannelGroup']):
#         rcg_x = Nhdf5.get('/RecordingChannelGroup_%i'%rcg_id)
#         if rcg_x.annotations == rcgannotations:
#             # check if units are overlapping
#             overlapping_units = np.intersect1d(rcg_x.channel_indexes,
#                                                [u.spiketrains[0].annotations['channel_index'] for u in rcg.units])
#             if overlapping_units:
#                 raise ValueError('There is already a block with these '
#                                  'annotations and units %s saved in the '
#                                  'hdf5. Not saving new spikesorting '
#                                  '(Annotations: %s)'%(overlapping_units,
#                                                       rcgannotations))
#             else:
#                 rcg_x.units.extend(rcg.units)
#                 rcg = rcg_x
#
#
#
#
#
#     # # checking if rcg with this sorting annotations already exists
#     # for rcg_id in range(Nhdf5.get_info()['RecordingChannelGroup']):
#     #     rcg_old = Nhdf5.get('/RecordingChannelGroup_%i'%rcg_id)
#     #     if rcg_old.annotations == rcgannotations:
#     #         warnings.warn('Not saving sorting for session %s. '
#     #                       'RCG with same annotation already exists!'%(savedir.split('/')[-1]))
#     #         Nhdf5.close()
#     #         return
#
#     Nhdf5.save(rcg)
#     Nhdf5.close()




# def load_spikesort(block,session,sorting_dir,software,parameter_dict,
#                    sort=True, electrode_list='all'):
#     """
#     :param load_from:
#     :param block:
#     :return:
#     """
#
#     corr_annotations = {'sorted':sort,'parameters':elephant.spike_sorting.get_updated_parameters(software, parameter_dict)}
#     # Warning this needs to be reset to line above!!!!!!!!!!!!!!!!!!!!!!!!!!
#     # corr_annotations ={'sorted': False, 'parameters': {'filter_low': None,
#     #                                                    'n_post': 10,
#     #                                                    'filter_high':np.array(400.0) * pq.Hz, 'threshold': -3.5, 'n_pre': -10, 'alignment': 'min'}}
#         # {'extraction_params': {'edge': 'falling',
#         #       'filter': [np.array(500.0) * pq.Hz, None],
#         #       'filter_order': 4,
#         #       'sp_win_align': [np.array(-1.0) * pq.ms, np.array(1.0) * pq.ms],
#         #       'sp_win_extract': [np.array(-0.5) * pq.ms, np.array(1.5) * pq.ms],
#         #       'threshold': 5.0},
#         #      'sorted': True,
#         #      'sorting_params': {'method': 'k_means_plus', 'ncomps': 2, 'num_units': 3}}
#
#     filename = os.path.join(sorting_dir,session + '_spikesorting.hdf5')
#     if not os.path.exists(filename):
#         raise IOError('File does not exist!')
#     print 'Loading extracted spikes from %s'%filename
#     Nhdf5 = neo.io.NeoHdf5IO(filename)
#
#     rcg = None
#     for rcg_id in range(Nhdf5.get_info()['RecordingChannelGroup']):
#         rcg_x = Nhdf5.get('/RecordingChannelGroup_%i'%rcg_id)
#         # print rcg_x.annotations
#         if rcg_x.annotations == corr_annotations:
#             rcg = rcg_x
#             break
#
#     if rcg == None:
#         raise ValueError('No spike sorting with parameters %s found.'%corr_annotations)
#
#     block.recordingchannelgroups.append(rcg)
#     for seg in block.segments:
#         tstart,tstop = seg.t_start, seg.t_stop
#         selected_sts = []
#         for unit in rcg.units:
#             if (electrode_list == 'all' or unit.spiketrains[0].annotations[
#                 'channel_index'] in electrode_list):
#                 for spiketrain in unit.spiketrains:
#                     # loading only spikes in current segment
#                     st = spiketrain.time_slice(tstart,tstop)
#                     # if len(st) > 0:
#                     selected_sts.append(st)
#                         #seg.spiketrains.append(st)
#
#         if electrode_list != 'all':
#             missing_els = np.setdiff1d(electrode_list,
#                                        [st.annotations['channel_index'] for
#                                         st in selected_sts])
#             if len(missing_els):
#                 raise ValueError('Can not load spiketrain for electrode ids '
#                                  '%s'%missing_els)
#
#         seg.spiketrains.extend(selected_sts)
#
#
#
#
#
#     Nhdf5.close()
#     return block



if __name__ == '__main__':

    spikesorting()

    pass