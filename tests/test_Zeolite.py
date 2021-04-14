from unittest import TestCase
from maze.zeolite import Zeolite, PerfectZeolite
from maze.index_mapper import IndexMapper
from ase import Atoms
import numpy as np
from ase.visualize import view
from maze.adsorbate import Adsorbate
import copy


class TestImperfectZeotype(TestCase):
    def test_build_cap_atoms(self):
        with self.subTest(msg="testing simple case"):
            cap_atoms_dict = {'H': [np.array([0, 0, 0])],
                              'O': [np.array([1, 2, 3])],
                              'N': [np.array([1, 3, 5])]}
            cap_atoms = Zeolite.build_cap_atoms(cap_atoms_dict)
            for atom in cap_atoms:
                self.assertIn(atom.symbol, cap_atoms_dict)
                self.assertCountEqual(cap_atoms_dict[atom.symbol][0], atom.position)
        with self.subTest(msg="testing multi atoms of same type"):
            cap_atoms_dict = {'H': [np.array([0, 0, 0]), np.array([1, 1, 1])],
                              'O': [np.array([3, 3, 3]), np.array([3, 4, 5])]}
            cap_atoms = Zeolite.build_cap_atoms(cap_atoms_dict)
            # this complicated for loop checks that the position is in the
            for atom in cap_atoms:
                array_in = False
                for pos in cap_atoms_dict[atom.symbol]:
                    all_equal = True
                    for dict_el, atom_el in zip(pos, atom.position):
                        if dict_el != atom_el:
                            all_equal = False
                            break
                    if all_equal:
                        array_in = True
                        break
                self.assertTrue(array_in)


    def test_cap_atom(self):
        with self.subTest(msg="test hydrogen capping"):
            iz = Zeolite(Atoms('OSiOSi', positions=[[0, 0, 0], [0, 0, -1], [0, 0, 1], [1,1,1]]))
            iz = iz.delete_atoms(2)  # delete Si
            iz = iz.cap_atoms()
            x = 0
            #self.assertEqual(len(iz), 3)
            for atom in iz:
                if atom.symbol == "H":
                    pass
                    #self.assertTrue(np.all(atom.position == np.array([0, 0, 1])))
                if atom.symbol == "O":
                    pass
                    #self.assertTrue(np.all(atom.position == np.array([0, 0, -1])))

        with self.subTest('oxygen and hydrogen capping'):
            cha = Zeolite.make('CHA')
            cha = cha.delete_atoms([i for i in range(0, 10)])
            for atom in cha:
                if atom.symbol == 'O':
                    atom.symbol = 'Po'

            cha = cha.cap_atoms()
            view(cha)

    def test_get_cluster(self):
        iz = Zeolite.make('BEA')
        cluster, od = iz.get_cluster(174)

        with self.subTest(msg='testing length is valid'):
            self.assertEqual(len(cluster) + len(od), len(iz))

        with self.subTest(msg='assert ztype correct'):
            self.assertEqual(iz.ztype, 'Zeolite')
            self.assertEqual(cluster.ztype, 'Cluster')
            self.assertEqual(od.ztype, 'Open Defect')

        with self.subTest(msg='assert name is correct'):
            self.assertIn('Cluster', cluster.name)
            self.assertIn('Open Defect', od.name)

        with self.subTest(msg='test names registered with index mapper'):
            # check names are in
            self.assertIn(cluster.name, iz.index_mapper.names)
            self.assertIn(od.name, iz.index_mapper.names)
            self.assertIn(iz.name, iz.index_mapper.names)

        with self.subTest(msg='test indices are in index mapper'):
            self.assertTrue(iz.index_mapper.get_reverse_main_index(iz.name))
            self.assertTrue(cluster.index_mapper.get_reverse_main_index(iz.name))
            self.assertTrue(od.index_mapper.get_reverse_main_index(iz.name))


    def test_integrate(self):
        iz = Zeolite.make('BEA')
        cluster, od = iz.get_cluster(174)
        with self.subTest(msg='integrate other zeolite'):
            new_iz = cluster.integrate(od)
            self.assertEqual(len(new_iz), len(iz))
            # this should be fixed eventually
                #new_iz_index = iz.index_mapper.get_index(iz.name, new_iz.name, atom.index)
                #print(new_iz_index)
                #print(atom.symbol, new_iz[new_iz_index].symbol)
                #print()
                #self.assertEqual(atom.symbol, new_iz[new_iz_index].symbol)
                #self.assertEqual(atom.position, new_iz[new_iz_index].position)


    def test_integrate_adsorbate(self):
        iz = Zeolite.make('BEA')
        n3 = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
        n3_ads = Adsorbate(n3, name='n3')
        iz2, n3_ads2 = iz.integrate_adsorbate(n3_ads)
        with self.subTest(msg='test index map integration'):
            self.assertIn(n3_ads2.name, iz2.index_mapper.names)
        with self.subTest(msg='test lengths make sense'):
            self.assertEqual(len(iz2), len(iz) + len(n3_ads))
        with self.subTest(msg='Test that adsorbate name is in additions adsorbate list'):
            self.assertIn(n3_ads2.name, list(iz2.additions['adsorbate']))

    def test_remove_adsorbate(self):
        iz = Zeolite.make('BEA')
        n3 = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
        n3_ads = Adsorbate(n3, name='n3')
        iz2, n3_ads2 = iz.integrate_adsorbate(n3_ads)
        iz3 = iz2.remove_adsorbate(n3_ads2)
        with self.subTest(msg='test that adsorbate is still in index map'):
            self.assertIn(n3_ads2.name, iz3.index_mapper.names)
        with self.subTest(msg='test length back to original length'):
            self.assertEqual(len(iz3), len(iz))
        with self.subTest(msg='Test that adsorbate name is not in additions adsorbate list'):
            self.assertNotIn(n3_ads2.name, list(iz3.additions['adsorbate']))

    def test_change_atom_properties(self):
        self.fail()
    def test_build_h_atoms_cap_dict(self):
        self.fail()

    def test_build_o_atoms_cap_dict(self):
       pass


    def test_build_all_atoms_cap_dict(self):
        self.fail()

    def test_get_h_pos(self):
        self.fail()

    def test_find_missing_atom(self):
        self.fail()

    def test_delete_atoms(self):
        iz = Zeolite.make('BEA')
        sites_to_delete = iz.site_to_atom_indices['T1']
        iz2 = iz.delete_atoms(sites_to_delete)
        with self.subTest(msg='Test that the lengths are correct'):
            self.assertEqual(len(iz2) + len(sites_to_delete), len(iz))
        with self.subTest(msg='Test that the T sites are removed'):
            self.assertFalse(iz2.site_to_atom_indices['T1'])
        with self.subTest(msg='test that all sites in iz2 are mappable to iz1'):
            for atom in iz2:
                iz_index = iz2.index_mapper.get_index(iz2.name, iz.name, atom.index)
                self.assertIsNotNone(iz_index)
                for p1, p2 in zip(atom.position, iz[iz_index].position):
                    self.assertEqual(p1, p2)
                self.assertEqual(atom.symbol, iz[iz_index].symbol)
                self.assertEqual(atom.tag, iz[iz_index].tag)
        with self.subTest(msg='test that T sites map to None'):
            for T1_site in iz.site_to_atom_indices['T1']:
                self.assertIsNone(iz.index_mapper.get_index(iz.name, iz2.name, T1_site))

    def test_set_attrs_source(self):
        self.fail()

    def test__change_atoms(self):
        self.fail()

    def test_add_atoms(self):
        iz = Zeolite.make('BEA')
        iz_name_index = int( str(iz.name.split('_')[-1]))
        n3 = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
        addition_type = 'fish'
        addition_name = addition_type + '_' + str(iz_name_index + 1)
        iz2 = iz.add_atoms(n3, addition_type)
        with self.subTest(msg='test that added atom is in the index map'):
            self.assertIn(addition_name, iz2.index_mapper.names) # this should be the second addition
        with self.subTest(msg='test that added atoms is in additions'):
            self.assertIn(addition_name, iz2.additions[addition_type])
        with self.subTest(msg='test that the length is correct'):
            self.assertEqual(len(iz2), len(iz) + len(n3))
        with self.subTest(msg='test mapping between n3 and iz'):
            for atom in n3:
                iz_index = iz2.index_mapper.get_index(addition_name, iz2.name, atom.index)
                self.assertIsNotNone(iz_index)
                for p1, p2 in zip(atom.position, iz2[iz_index].position):
                    self.assertEqual(p1, p2)
                self.assertEqual(atom.symbol, iz2[iz_index].symbol)
                self.assertEqual(atom.tag, iz2[iz_index].tag)

        with self.subTest(msg='test that all of the indices of the main z are there'):
            for key, value in iz2.index_mapper.main_index.items():
                self.assertIsNotNone(value[iz2.name])



    def test_needs_cap(self):
        self.fail()

    def test_get_oxygen_cap_pos(self):
        cha = Zeolite.make('CHA')
        cha = cha.delete_atoms([i for i in range(0, 10)])
        atoms_to_cap = []
        for atom in cha:
            if atom.symbol == 'O':
                if cha.needs_cap(atom.index, 2):
                    atoms_to_cap.append(atom.index)
            if atom.symbol == 'Si':
                if cha.needs_cap(atom.index, 4):
                    atoms_to_cap.append(atom.index)

        for index in atoms_to_cap:
            oxygen_cap_pos = cha.get_oxygen_cap_pos(index)

    def test_get_hydrogen_cap_pos_simple(self):
        self.fail()

    def test_retag_self(self):
        with self.subTest(msg="retag doesn't throw errrors"):
            cha = Zeolite.make('CHA')
            cha = cha.retag_self()
        with self.subTest(msg='test main index matches tags'):
            reverse_index_map = cha.index_mapper.get_reverse_main_index(cha.name)
            for atom in cha:
                self.assertEqual(atom.tag, reverse_index_map[atom.index])

    def test_build_additions_map(self):
        z = Zeolite.make('CHA')
        water = Adsorbate(Atoms('H2O', positions=[[-1,-1,0], [0,0,0], [1,1,0]]))
        z, ads = z.integrate_adsorbate(water)
        z, ads = z.integrate_adsorbate(water)
        z, ads = z.integrate_adsorbate(water)

        print(z.build_additions_map())

    def test_register_with_parent(self):
        with self.subTest(msg='test simple case, no additions'):
            parent = PerfectZeolite.make('BEA')
            child = Zeolite.make('BEA')
            child = child.retag_self()
            child.register_with_parent(parent)
            self.assertEqual(parent.index_mapper, child.index_mapper)
            self.assertEqual(parent, child.parent_zeotype)
            self.assertIn(child.name, child.index_mapper.names)

        with self.subTest(msg='test additions'):
            parent = PerfectZeolite.make('BEA')
            additions_dict = copy.deepcopy(parent.additions)
            child = Zeolite.make('BEA')
            child = child.add_atoms(Atoms("H2", positions=[[0,0,0], [1,1,1]]), 'H2').delete_atoms([1, 2, 3, 4])
            ads_map = child.build_additions_map()
            child.additions = additions_dict
            child = child.retag_self()
            child.register_with_parent(parent, ads_map)
            self.assertEqual(parent.index_mapper, child.index_mapper)
            self.assertEqual(parent, child.parent_zeotype)
            self.assertIn(child.name, child.index_mapper.names)
            name_list = []
            for key in ads_map:
                name_list.extend(list(ads_map[key]))
            for name in name_list:
                self.assertIn(name, child.index_mapper.names)

    def test_get_self_to_main_index_map(self):
        # TODO: This function isn't needed and can be removed
        import pandas as pd
        z = Zeolite.make('BEA')
        z = z.add_atoms(Atoms("H2",  positions=[[0,0,0], [1,1,1]]), 'H2').delete_atoms([1, 2, 3, 4]).cap_atoms()
        for atom in z:
            if atom.index == 194:
                atom.symbol = 'He'
        view(z)
        print(z.additions['o_caps'])
        #z = z.add_atoms(Atoms("H2", positions=[[0,0,0], [1,1,1]]), 'H2') #.cap_atoms()
        df = pd.DataFrame(z.index_mapper.main_index)
        print(df.T.to_string())
        #print(z.get_self_to_main_index_map())

    def test__check_unique_positions(self):
        overlapping = Atoms('H2')
        non_overlapping = Atoms('H2', positions=[[0, 0, 0], [1, 1, 1]])
        with self.subTest('check if non_overlapping pass'):
            try:
                Zeolite._check_unique_positions(non_overlapping.get_positions())
                self.assertTrue(True)
            except ValueError:
                self.fail('error thrown in wrong case')
        with self.subTest('check if overlapping fail'):
            try:
                Zeolite._check_unique_positions(overlapping.get_positions())
                self.fail('error not thrown')
            except ValueError:
                self.assertTrue(True)

    def test_remove_caps(self):
        bea = Zeolite.make('BEA')
        cluster, od = bea.get_cluster(154)
        with self.subTest('testing removing caps without arguments'):
            capped_cluster = cluster.cap_atoms()
            uncapped_cluster = capped_cluster.remove_caps()
            for cap_type in uncapped_cluster.additions.keys():
                if 'cap' not in cap_type:
                    continue
                self.assertFalse(uncapped_cluster.additions[cap_type])  # assert empty
            self.assertEqual(len(uncapped_cluster), len(cluster), msg='num atoms in uncapped = num atoms in original')
            for atom in uncapped_cluster:
                self.assertFalse(atom.symbol=='H', msg='no hydrogen in uncapped cluster')  # check for no hydrogens

        with self.subTest('testing oxygen and hydrogen caps'):
            cluster_wo_o = cluster.delete_atoms([1])
            capped_cluster = cluster_wo_o.cap_atoms()
            uncapped_cluster = capped_cluster.remove_caps()
            for cap_type in uncapped_cluster.additions.keys():
                if 'cap' not in cap_type:
                    continue
                self.assertFalse(uncapped_cluster.additions[cap_type])  # assert empty
            self.assertEqual(len(uncapped_cluster), len(cluster_wo_o), 'num atoms in uncapped = num atoms in original')
            for atom in uncapped_cluster:
                self.assertFalse(atom.symbol=='H', msg='no hydrogen in uncapped cluster')  # check for no hydrogens

        with self.subTest('remove only H caps'):
            cluster_wo_o = cluster.delete_atoms([1])
            capped_cluster = cluster_wo_o.cap_atoms()
            uncapped_cluster = capped_cluster.remove_caps(cap_type='h_cap')
            for cap_type in uncapped_cluster.additions.keys():
                if 'h_cap' not in cap_type:
                    continue
                self.assertFalse(uncapped_cluster.additions[cap_type])  # assert empty
            self.assertEqual(len(uncapped_cluster), len(cluster),
                             msg='num atoms in uncapped = num atoms in original')
            for atom in uncapped_cluster:
                self.assertFalse(atom.symbol == 'H',
                                 msg='no hydrogen in uncapped cluster')  # check for no hydrogens

        view(uncapped_cluster)