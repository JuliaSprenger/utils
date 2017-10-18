import unittest
import os.path
import numpy as np
import quantities as pq
from neo import Block, SpikeTrain, ChannelIndex, Unit, Segment
import elephant.spike_sorting
from private_bib.spikesort import load_spikesorting, save_spikesorting

def default_data(block=None, n_chidx=1, n_units=1):

    # generate new block if none provided, otherwise attach to provided block
    if block is None:
        block = Block()
    if not block.segments:
        block.segments = [Segment()]

    for id in range(n_chidx):
        sorting_hash = elephant.spike_sorting.SpikeSorter.get_sorting_hash({
            'channel_index':id,
            'random annotation':np.random.randint(0, 10**10)})
        chidx = ChannelIndex([], sorting_hash=sorting_hash)
        chidx.block = block
        block.channel_indexes.append(chidx)

    for chid, chidx in enumerate(block.channel_indexes):
        for id in range(n_units):
            unit = Unit(unit_id=id, channel_id=chid)
            chidx.units.append(unit)
            unit.channel_index = chidx

            for st_id in range(id):
                st = SpikeTrain(np.random.uniform(0, st_id, 1) * pq.s,
                                t_start=0 * pq.s,
                                t_stop=st_id * pq.s, spiketrain_id=st_id)
                unit.spiketrains.append(st)
                block.segments[0].spiketrains.append(st)
                st.unit = unit
                st.segment = block.segments[0]

    block.create_relationship()
    return block

class SpikeSaveLoadTestCase(unittest.TestCase):

    def setUp(self):
        self.block = default_data(n_chidx=1,
                                  n_units=1)

        sorting_file = 'testdata'
        if os.path.exists(sorting_file + '_spikesorting.hdf5'):
            os.remove(sorting_file + '_spikesorting.hdf5')

        self.sorting_hash = self.block.channel_indexes[0].annotations[
            'sorting_hash']

        save_spikesorting(sorting_file, self.block,
                          sorting_hash=self.sorting_hash)

        self.new_block = Block(type='loaded block')
        load_spikesorting(self.new_block, sorting_file='testdata',
                          sorting_hash=self.sorting_hash)

        self.object_classes = ['ChannelIndex', 'Unit', 'SpikeTrain', 'Segment',
                               'AnalogSignal']

    def test_data_exist(self):
        for obj_class in self.object_classes:
            old_objs = self.block.list_children_by_class(obj_class)
            new_objs = self.new_block.list_children_by_class(obj_class)
            self.assertEqual(len(old_objs), len(new_objs))

    def test_for_annotations(self):
        for obj_class in self.object_classes:
            old_objs = self.block.list_children_by_class(obj_class)
            new_objs = self.new_block.list_children_by_class(obj_class)
            for id in range(len(old_objs)):
                d2 = old_objs[id].annotations
                d1 = new_objs[id].annotations
                self.assertTrue(set(d2.items()).issubset(set(d1.items())))

    def test_for_data(self):
        # This test will fail until neuralensemble issue #410 is solved
        for obj_class in self.object_classes:
            old_objs = self.block.list_children_by_class(obj_class)
            new_objs = self.new_block.list_children_by_class(obj_class)
            for id in range(len(old_objs)):
                if hasattr(old_objs[id],'index'):
                    np.testing.assert_array_equal(old_objs[id].index,
                                                  new_objs[id].index)
                if hasattr(old_objs[id],'times'):
                    np.testing.assert_array_equal(old_objs[id].times,
                                                  new_objs[id].times)



class SpikeSaveLoadComplexTestCase(unittest.TestCase):

    def setUp(self):
        self.block = default_data(n_chidx=3,
                                  n_units=4)

        sorting_file = 'testdata'
        if os.path.exists(sorting_file + '_spikesorting.hdf5'):
            os.remove(sorting_file + '_spikesorting.hdf5')

        self.sorting_hash = self.block.channel_indexes[2].annotations[
            'sorting_hash']

        save_spikesorting(sorting_file, self.block,
                          sorting_hash=self.sorting_hash)

        self.new_block = Block(type='loaded block')
        load_spikesorting(self.new_block, sorting_file='testdata',
                          sorting_hash=self.sorting_hash)

        self.object_classes = ['ChannelIndex', 'Unit', 'SpikeTrain', 'Segment',
                               'AnalogSignal']

    def test_data_exist(self):
        for obj_class in self.object_classes:
            old_objs = self.block.channel_indexes[2].list_children_by_class(
                obj_class)
            new_objs = self.new_block.channel_indexes[0].list_children_by_class(
                obj_class)
            self.assertEqual(len(old_objs), len(new_objs), obj_class)

    def test_for_data(self):
        # This test will fail until neuralensemble issue #410 is solved
        for obj_class in self.object_classes:
            old_objs = self.block.list_children_by_class(obj_class)
            new_objs = self.new_block.list_children_by_class(obj_class)
            for id in range(len(old_objs)):
                if hasattr(old_objs[id], 'index'):
                    np.testing.assert_array_equal(old_objs[id].index,
                                                  new_objs[id].index)
                if hasattr(old_objs[id], 'times'):
                    np.testing.assert_array_equal(old_objs[id].times,
                                                  new_objs[id].times)


    #
    # def test_traincount(self):
    #     self.default_sorting()
    #     self.assertEqual(len(self.block.segments[0].spiketrains),
    #                      sum([a.shape[-1] for a in
    #                           self.block.segments[0].analogsignals]))
    #
    # def test_spiketrain_relations(self):
    #     self.default_sorting()
    #     self.assertEqual(len(self.block.channel_indexes),1)
    #     sts = self.block.segments[0].spiketrains
    #     chidx = self.block.channel_indexes[0]
    #     anasigs = self.block.segments[0].analogsignals
    #     self.assertEqual(len(chidx.units),
    #                      len(self.block.segments[0].analogsignals))
    #     for unit_idx in range(len(chidx.units)):
    #         self.assertEqual(len(chidx.units[unit_idx].spiketrains),
    #                          anasigs[unit_idx].shape[-1])
    #     for st in sts:
    #         self.assertTrue(st.unit is not None)