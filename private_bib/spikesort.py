import os
import sys
import copy
import warnings
import quantities as pq
import numpy as np
import neo
import elephant.spike_sorting
from private_bib import neo_utils
from private_bib.developmentio import DevelopmentIO


# settings = {'sessiondir':'/home/julia/data/SPP/data/2015-02-03_14-22-50',
#        # 'sessiondir': '/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/data/2015-02-20_10-19-03',
#        'sortdir': '/home/julia/projects/SPP1665/analysis/spikesortings',
#        'mdatadir': '/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/metadata',
#        'cachedir': '/home/julia/cache'}
#
# software = 'phy'
# method = ''
# waveforms = True
# sort = True
# # parameter_dict = {'keep_temp_files':True , 'num_cpus': 2,
# #                   'threshold_strong_std_factor':6,
# #                   'threshold_weak_std_factor':3}#{'threshold': 5.0}
# all_channels_at_once = True
# starts = 0*pq.s
# stops = 1*pq.min



def spikesort_session(session=None):
    from private_bib.spikesort_configuration import (cache_directory,
                                                     extraction_parameters,
                                                     sorting_directory)
    if session is None:
        try:
            from private_bib.spikesort_configuration import session
        except:
            from private_bib.spikesort_configuration import sessions as session

    # listify session
    if not isinstance(session, list):
        session = [session]
    sessions = session


    for session in sessions:
        print(session)
        session_name = os.path.splitext(os.path.basename(session))[0]
        print('Sorting session {}'.format(session_name))
        IO = DevelopmentIO(session,
            cachedir=cache_directory,
            use_cache='datesize',
            print_diagnostic=False)

        chids = range(1,33) #TODO: This should not be hardcoded here! Better
        # use available channel ids provided by IO
        for chid in chids:
            print('Ch{}'.format(chid))
            block = IO.read_block(t_starts=None, t_stops=None,
                                  analogsignals=True,
                                  electrode_list=[chid], unit_list=None)

            # loading spikesorting
            sorter = elephant.spike_sorting.SpikeExtractor(**extraction_parameters)
            sorter.sort_block(block)
            save_spikesorting(os.path.join(sorting_directory, session_name),
                              block, sorting_hash=sorter.sorting_hash)


    # settings .update(kwargs)
    #
    # session = settings ['sessiondir'].split('/')[-1]
    #
    # IO = DIO.DevelopmentIO(sessiondir=settings['sessiondir'],
    #                        mdatadir=settings['mdatadir'],
    #                        cachedir=settings['cachedir'],
    #                        use_cache='datesize',
    #                        print_diagnostic=True)
    #
    # # Name generation for spike sorting file
    # sorting_prefix = '%s_%s' % (session, software)
    # sorting_prefix += str(parameter_dict)
    # if waveforms:
    #     sorting_prefix += '_wf'
    #
    # block = None
    # channel_index = None

        # electrode_list = electrode if electrode_list==None else electrode_list
        # anasig_block = IO.read_block(t_starts=t_starts, t_stops=t_stops,
        #                              electrode_list=[electrode], units=None)
        # print('Completed data loading.')
        #
        # if block is None:
        #     block = anasig_block
        #     if len(block.segments[0].analogsignalarrays) == 0:
        #         continue
        #     channel_index = anasig_block.segments[0].analogsignalarrays[0].channel_index
        # else:
        #     for s in range(len(block.segments)):
        #         block.segments[s].analogsignalarrays = anasig_block.segments[s].analogsignalarrays
        #         for a in range(len(block.segments[s].analogsignalarrays)):
        #             block.segments[s].analogsignalarrays[a].segment = block.segments[s]
        #             block.segments[s].analogsignalarrays[a].channel_index = channel_index

    #     extractor = elephant.spike_sorting.SpikeExtractor(parameter_dict)
    #     extractor.sort_block(block)
    #     # elephant.spike_sorting.generate_spiketrains(block,'phy',
    #     #                                           waveforms=waveforms,
    #     #                                           sort=sort,
    #     #                                           parameter_dict=parameter_dict)
    # save_spikesorting(os.path.join(settings ['sortdir'], session), block,
    #                   sorting_hash=extractor.sorting_hash)


def save_spikesorting(sorting_file, block, sorting_hash=None,
                      parameter_dict=None):
    """
    :param savedir:
    :param block:
    :return:
    """

    if parameter_dict is None and sorting_hash is None:
        raise ValueError('Please provide either a parameter dictionary or a '
                         'sorting hash to specify the sorting you want to '
                         'load.')
    elif parameter_dict is not None:
        sorting_hash = elephant.spike_sorting.SpikeSorter.get_sorting_hash(parameter_dict)

    filename = sorting_file + '_spikesorting.hdf5'
    print('Saving sorted spikes (sorting hash {}) at {}'
          ''.format(sorting_hash, filename))

    sorting_chidx = [i for i in block.channel_indexes if
                     ('sorting_hash' in i.annotations
                      and i.annotations['sorting_hash']==sorting_hash)]

    if len(sorting_chidx) == 0:
        warnings.warn('No channel_index selected for saving spikesorting with '
                      'hash "{}". Not saving anything.'.format(sorting_hash))
        return

    elif len(sorting_chidx) > 1:
        warnings.warn('Multiple channels with sorting hash "{}" exist. Not '
                      'saving anything.'.format(sorting_hash))
        return
    sorting_chidx = sorting_chidx[0]

    with neo.nixio.NixIO(filename,'rw') as nix_file:
        # check if there is already a block with the same channel_index annotation
        nix_blocks = nix_file.read_all_blocks()
        new_block = False
        if len(nix_blocks) is 0:
            new_block = True
            nix_blocks = [neo.Block(name='spike sorting block')]
            nix_blocks[0].segments.append(neo.Segment(name='spike sorting '
                                                          'segment'))

        assert len(nix_blocks[0].segments) == 1
        nix_block = nix_blocks[0]
        seg = nix_blocks[0].segments[0]

        # check if channel_index for this sorting already exists
        chidx = None
        for chidx_i in nix_block.channel_indexes:
            if chidx_i.annotations['sorting_hash'] == sorting_hash:
                chidx = chidx_i
                break
        if chidx is None:
            chidx = neo.ChannelIndex([-1], name='spike sorting',
                                     sorting_hash=sorting_hash)
            chidx.block = nix_block
            nix_block.channel_indexes.append(chidx)

        duplicated_units = [copy.deepcopy(u) for u in sorting_chidx.units]
        for uid, unit in enumerate(duplicated_units):
            chidx.units.append(unit)
            unit.channel_index = chidx
            # spiketrain is automatically duplicated when rescaled
            original_sts = sorting_chidx.units[uid].spiketrains
            unit.spiketrains = []
            for o_st in original_sts:
                # TODO: Remove this rescaling quickfix once nixpy is fixed
                # but duplicate spiketrain instead
                duplicated_st = o_st.rescale('s')
                unit.spiketrains.append(duplicated_st)

            # duplicating and decoupling of spiketrains from original neo structure
            for st in unit.spiketrains:
                ######################
                # # TODO: Remove this quickfix once nixpy is fixed
                # # Quickfix
                # st_new = st.rescale('s')
                # unit.spiketrains.remove(st)
                # unit.spiketrains.append(st_new)
                # st_new.segment = seg

                # TODO: Remove this quickfix once neo nixio is fixed
                if st.left_sweep:
                    st.left_sweep = [st.left_sweep.rescale('s').magnitude]*pq.s
                st.sampling_rate = st.sampling_rate.rescale('Hz')
                if ('invalid_waveforms' in st.annotations
                    and st.annotations['invalid_waveforms'] == []):
                    st.annotations['invalid_waveforms'] = 0
                ########################

                st.unit = unit
                st.segment = seg

            seg.spiketrains.extend(unit.spiketrains)

        nix_block.create_relationship()
        nix_file.write_block(nix_block)


def load_spikesorting(block, sorting_file, parameter_dict=None,
                      sorting_hash=None, electrode_list='all'):
    """
    :param load_from:
    :param block:
    :return:
    """

    if parameter_dict is None and sorting_hash is None:
        raise ValueError('Please provide either a parameter dictionary or a '
                         'sorting hash to specify the sorting you want to '
                         'load.')
    elif parameter_dict is not None:
        sorting_hash = elephant.spike_sorting.SpikeSorter.get_sorting_hash(parameter_dict)

    filename = sorting_file +'_spikesorting.hdf5'
    if not os.path.exists(filename):
        raise IOError('File does not exist!')

    print('Attempting to load extracted spikes from {} for sorting {}...'
          ''.format(filename, sorting_hash))
    with neo.io.NixIO(filename, 'ro') as nix_file:
        nix_block = nix_file.read_block()
    channel_indexes = []
    if nix_block is not None:
        channel_indexes = [c for c in nix_block.channel_indexes
                           if ('sorting_hash' in c.annotations
                               and c.annotations['sorting_hash'] == sorting_hash)]

    if len(channel_indexes) != 1:
        raise ValueError('No spike sorting found with sorting hash {} ('
                         'parameters {}).'.format(sorting_hash, parameter_dict))
    nix_sorted_chidx = channel_indexes[0]

    sorted_chidx = [c for c in block.channel_indexes
                    if ('sorting_hash' in c.annotations
                        and c.annotations['sorting_hash'] == sorting_hash)]

    if len(sorted_chidx)>1:
        raise ValueError('Multiple channel indexes for same sorting present')
    elif len(sorted_chidx)==1:
        sorted_chidx = sorted_chidx[0]
    else:
        sorted_chidx = copy.deepcopy(nix_sorted_chidx)
        sorted_chidx.block = block
        block.channel_indexes.append(sorted_chidx)

    for seg in block.segments:
        tstart, tstop = seg.t_start, seg.t_stop
        duplicated_units = []
        for unit in nix_sorted_chidx.units:
            # skipping non-selected units (and spiketrains)
            if not (electrode_list == 'all' or unit.annotations['channel_id'] in electrode_list):
                continue

            duplicated_unit = copy.deepcopy(unit)
            duplicated_unit.channel_index = sorted_chidx
            sorted_chidx.units.append(duplicated_unit)
            duplicated_units.append(duplicated_unit)

            for spiketrain in unit.spiketrains:
                duplicated_st = spiketrain.time_slice(tstart, tstop)
                duplicated_st.unit = duplicated_unit
                duplicated_unit.spiketrains.append(duplicated_st)
                seg.spiketrains.append(duplicated_st)
                duplicated_st.segment = seg

                ############ Fix for issue
                # https://github.com/NeuralEnsemble/python-neo/issues/373
                if duplicated_st.waveforms.shape == (0,):
                    duplicated_st.waveforms = \
                        duplicated_st.waveforms.reshape((0,0,0))
                ################


        if electrode_list != 'all':
            missing_els = np.setdiff1d(electrode_list,
                                       [u.annotations['channel_id'] for
                                        u in duplicated_units])
            if len(missing_els):
                raise ValueError('Can not load spiketrain for electrode ids '
                                 '%s'%missing_els)
    print('Loading successful')
    return block


def get_sorting(IO, block, sorter, sorting_dir, ellist):
    # try to load spiketrains for each electrode individually
    filename = os.path.basename(block.file_origin)
    for elid in ellist:
        try:
            load_spikesorting(block=block,
                              session=filename,
                              sorting_dir=sorting_dir,
                              parameter_dict=sorter.parameter_dict,
                              electrode_list=[elid])
        except (ValueError, IOError) as e:
            messages = ['File does not exist!',
                        'No spike sorting found',
                        'Can not load spiketrain for electrode ids']

            if any([msg in e.message for msg in messages]):
                generate_new_sorting(IO, block, elid, sorter,
                                     sorting_dir, filename)
            else:
                raise e


def generate_new_sorting(IO, block, electrode_id, sorter, sorting_dir,
                        filename):
    print('Generating new spikes for electrode {}'.format(electrode_id))
    print('Reading full data block...')
    # global block
    # loading electrode-wise for memory reasons
    t_starts, t_stops = [], []
    for seg_id, seg in enumerate(block.segments):
        if seg.t_start == seg.t_stop:
            print('WARNING: segment {} ({}/{}) has same start and stop time {' \
                  '}'.format(seg, seg_id, len(block.segments)-1, seg.t_start))
            continue
        t_starts.append(seg.t_start)
        t_stops.append(seg.t_stop)
    print(t_starts)
    print(t_stops)
    anasig_block = IO.read_block(t_starts=t_starts, t_stops=t_stops,
                                 analogsignals=True,
                                 electrode_list=[electrode_id],
                                 unit_list=None)

    neo_utils.check_neo_compliant(anasig_block)

    print('Finished reading analogsignal block')
    sorter.sort_block(anasig_block)

    # neo_utils.check_neo_compliant(anasig_block)

    save_spike_dir = os.path.join(sorting_dir, filename)
    sorting_hash = sorter.sorting_hash
    save_spikesorting(save_spike_dir, anasig_block, sorting_hash)
    neo_utils.check_neo_compliant(anasig_block)

    # loading generated spikes
    load_spikesorting(block=block, session=filename,
                      sorting_dir=sorting_dir,
                      parameter_dict=sorter.parameter_dict,
                      electrode_list=[electrode_id])
    neo_utils.check_neo_compliant(block)


if __name__ == '__main__':
    print('Main!')
    session = None
    if len(sys.argv)>1:
        session = sys.argv[1]

    print(sys.path)

    spikesort_session(session)
