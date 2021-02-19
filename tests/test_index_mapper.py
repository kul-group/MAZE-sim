from unittest import TestCase
from maze.index_mapper import IndexMapper

class TestIndexMapper(TestCase):
    def test_init(self):
        atom_indices = [i for i in range(10, 20)]
        index_mapper = IndexMapper(atom_indices)
        with self.subTest(msg='test main index keys'):
            self.assertCountEqual(index_mapper.main_index.keys(), range(0, len(atom_indices)))

        with self.subTest(msg='test main index values'):
            for value, parent_index in zip(index_mapper.main_index.values(), range(10, 20)):
                tmp_dict = {'parent': parent_index}
                self.assertDictEqual(tmp_dict, value)

    def test_get_id(self):
        atom_indices = [i for i in range(0, 10)]
        index_mapper = IndexMapper(atom_indices)
        with self.subTest(msg='test that id is a string'):
            self.assertIsInstance(index_mapper.get_id(), str)

        with self.subTest(msg='ensure call 1 < call 2'):
            id1 = index_mapper.get_id()
            id2 = index_mapper.get_id()
            self.assertLess(int(id1),int(id2))

        with self.subTest(msg='ensure two seperate index mappers produce different indices'):
            index_mapper_2 = IndexMapper(atom_indices)
            id3 = index_mapper.get_id()
            id4 = index_mapper_2.get_id()
            self.assertLess(int(id3), int(id4))


    def test_get_unique_name(self):
        atom_indices = [i for i in range(0, 10)]
        index_mapper = IndexMapper(atom_indices)
        name1 = index_mapper.get_unique_name('zeotype')
        name2 = index_mapper.get_unique_name('zeotype')
        name3 = index_mapper.get_unique_name('fish')
        with self.subTest(msg='testing that names are different'):
            self.assertNotEqual(name1, name2)
        with self.subTest(msg='testing that names contain input string'):
            self.assertIn('zeotype', name1)
            self.assertIn('fish', name3)


    def test_get_reverse_main_index(self):
        index_mapper = IndexMapper([i for i in range(10, 20)])
        reverse_main_index = index_mapper.get_reverse_main_index('parent')
        with self.subTest(msg='testing keys'):
            self.assertCountEqual(reverse_main_index.keys(), range(10,20))
        with self.subTest(msg='testing values'):
            self.assertCountEqual(reverse_main_index.values(), range(0, 10))

    def test_get_name1_to_name2_map(self):
        index_mapper = IndexMapper([i for i in range(10, 20)])
        old_to_new_map = dict(zip(range(10, 20), range(20, 30)))
        index_mapper.register('parent', 'fish', old_to_new_map)
        n1_to_n2_map = index_mapper.get_name1_to_name2_map('parent', 'fish')
        self.assertCountEqual(n1_to_n2_map.keys(), range(10, 20))
        self.assertCountEqual(n1_to_n2_map.values(), range(20, 30))


    def test_register_with_main(self):
        parent_indices = [i for i in range(10, 20)]
        fish_indices = [i for i in range(20, 30)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(range(0,len(parent_indices)), fish_indices))
        index_mapper.register_with_main('fish', old_to_new_map)
        with self.subTest(msg='test main index keys'):
            self.assertCountEqual(range(0, len(parent_indices)), index_mapper.main_index.keys())

        with self.subTest(msg='test main index values'):
            for value, parent_index, fish_index in zip(index_mapper.main_index.values(), parent_indices, fish_indices):
                tmp_dict = {'parent': parent_index, 'fish': fish_index}
                self.assertDictEqual(tmp_dict, value)


    def test_register(self):
        parent_indices = [i for i in range(10, 20)]
        fish_indices = [i for i in range(20, 30)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(parent_indices, fish_indices))
        index_mapper.register('parent', 'fish', old_to_new_map)
        with self.subTest(msg='test main index keys'):
            self.assertCountEqual(range(0, len(parent_indices)), index_mapper.main_index.keys())

        with self.subTest(msg='test main index values'):
            for value, parent_index, fish_index in zip(index_mapper.main_index.values(), parent_indices, fish_indices):
                tmp_dict = {'parent': parent_index, 'fish': fish_index}
                self.assertDictEqual(tmp_dict, value)

    def test_add_name(self):
        parent_indices = [i for i in range(0, 10)]
        fish_indices = [i for i in range(30, 20, -1)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(parent_indices, fish_indices))
        index_mapper.add_name('fish', 'parent', old_to_new_map)
        with self.subTest(msg='test main index keys'):
            self.assertCountEqual(range(0, len(parent_indices)), index_mapper.main_index.keys())

        with self.subTest(msg='test main index values'):
            for value, parent_index, fish_index in zip(index_mapper.main_index.values(), parent_indices, fish_indices):
                tmp_dict = {'parent': parent_index, 'fish': fish_index}
                self.assertDictEqual(tmp_dict, value)

    def test__make_none_dict(self):
        parent_indices = [i for i in range(10, 20)]
        fish_indices = [i for i in range(20, 30)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(parent_indices, fish_indices))
        index_mapper.register('parent', 'fish', old_to_new_map)
        none_dict = index_mapper._make_none_dict()
        none_dict_desired = {'parent': None, 'fish': None}
        self.assertDictEqual(none_dict, none_dict_desired)

    def test_add_atoms(self):
        parent_indices = [i for i in range(10, 20)]
        fish_indices = [i for i in range(20, 30)]
        index_mapper = IndexMapper(parent_indices)
        index_mapper.add_atoms('fish', fish_indices)
        with self.subTest(msg='test main index keys'):
            self.assertCountEqual(range(0, len(parent_indices) + len(fish_indices)), index_mapper.main_index.keys())

        parent_iter = iter(parent_indices)
        fish_iter = iter(fish_indices)
        with self.subTest(msg='test main index values'):
            for value in index_mapper.main_index.values():
                try:
                    parent_index = next(parent_iter)
                    fish_index = None
                except StopIteration:  # called too many times but okay because this is a test
                    parent_index = None
                    fish_index = next(fish_iter)

                tmp_dict = {'parent': parent_index, 'fish': fish_index}
                self.assertDictEqual(tmp_dict, value)

    def test_delete_atoms(self):
        parent_indices = [i for i in range(0, 10)]
        fish_indices = [i for i in range(30, 20, -1)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(parent_indices, fish_indices))
        index_mapper.add_name('fish', 'parent', old_to_new_map)
        indices_to_delete = [27, 28, 29]
        index_mapper.delete_atoms('fish', indices_to_delete)
        with self.subTest(msg='test main index keys'):
            self.assertCountEqual(range(0, len(parent_indices)), index_mapper.main_index.keys())

        with self.subTest(msg='test main index values'):
            for value, parent_index, fish_index in zip(index_mapper.main_index.values(), parent_indices, fish_indices):
                tmp_dict = {'parent': parent_index, 'fish': fish_index}
                if fish_index in indices_to_delete:
                    tmp_dict['fish'] = None
                self.assertDictEqual(tmp_dict, value)

    def test_get_index(self):
        parent_indices = [i for i in range(0, 10)]
        fish_indices = [i for i in range(30, 20, -1)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(parent_indices, fish_indices))
        index_mapper.add_name('fish', 'parent', old_to_new_map)
        with self.subTest(msg='test mapping between indices'):
            self.assertEqual(30, index_mapper.get_index('parent', 'fish', 0))

    def test_delete_name(self):
        parent_indices = [i for i in range(0, 10)]
        fish_indices = [i for i in range(30, 20, -1)]
        index_mapper = IndexMapper(parent_indices)
        old_to_new_map = dict(zip(parent_indices, fish_indices))
        index_mapper.add_name('fish', 'parent', old_to_new_map)
        self.assertIn('fish', index_mapper.names)
        index_mapper.delete_name('fish')
        self.assertNotIn('fish', index_mapper.names)
        for values in index_mapper.main_index.values():
            self.assertNotIn('fish', values.keys())

