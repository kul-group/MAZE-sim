import json
from copy import deepcopy
from pathlib import Path
from typing import Dict

import ase
from ase.db import connect

from maze.zeotypes import Zeotype, Cluster, OpenDefect, ImperfectZeotype


class ZeotypeDatabase:
    zeotype_dict = {'Zeotype': Zeotype, 'Cluster': Cluster, 'OpenDefect': OpenDefect,
                    'ImperfectZeotype': ImperfectZeotype}

    def __init__(self, ase_db_name: str, datajson_name: str):
        """
        Initializes a ZeotypeDatabase object. If ase_db_name and datajson_name exist, they are loaded. If they do not
        exist they are created.
        :param ase_db_name: filepath (including name) of the ase database (must end in db)
        :type ase_db_name: str
        :param datajson_name: name of the additional json database
        :type datajson_name: str
        """

        self.parent_zeotype_dict: Dict[str, int] = {}
        self.data_dict = {}
        self.ase_db = None
        self.json_filepath = datajson_name if 'json' in datajson_name else datajson_name + '.json'
        self.ase_db = connect(ase_db_name)
        json_path = Path(datajson_name)
        if json_path.is_file():
            with open(self.json_filepath, 'r') as f:
                self.data_dict = json.load(f)
                for key in self.data_dict.keys():
                    if self.data_dict[key]['name'] == 'parent':
                        self.parent_zeotype_dict[self.data_dict[key]['GUID']] = int(key)

        else:
            self.data_dict = {}

    def add_data(self, key: int, new_data: Dict) -> None:
        """
        Add data to the json database
        :param key: db id of Zeotype object data belongs to
        :type key: int
        :param new_data: dict of new data to add to database
        :type new_data: Dict
        :return: None
        :rtype: None
        """

        self.data_dict[key] = new_data
        self.write_data()

    def _get(self, key: int, default=None) -> Dict:
        """
        get an atoms object from the ase database
        :param key: id of atoms object in database
        :type key: int
        :param default: default value
        :type default: Optional[Union[str, int]]
        :return: Data dict
        :rtype: Dict
        """

        if default is None:
            row = self.ase_db.get(key)  # ase_db.get cannot except default=None as arg
        else:
            row = self.ase_db.get(key, default=default)

        atoms = row.toatoms()
        if key is default:
            return default
        data = deepcopy(self.data_dict[key])
        data['atoms'] = atoms
        return data

    def _write(self, zeotype: Zeotype, data: Dict) -> int:
        """
        Write a Zeotype to the ase db and jsondata dict
        :param zeotype: Zeotype object to write to the database
        :type zeotype: Zeotype
        :param data: data corresponding to Zeotype
        :type data: Dict
        :return: id of added Zeotype in both dbs
        :rtype: int
        """

        key = self.ase_db.write(zeotype)
        self.add_data(key, data)
        return key

    def write_data(self) -> None:
        """
        Write data dict to the json database
        :return: None
        :rtype: None
        """

        with open(self.json_filepath, 'w') as f:
            json.dump(self.data_dict, f, indent=4)


    def make_parent_dict(self, zeotype: Zeotype) -> Dict:
        """
        Make a dictionary for the parent zeotype
        :param zeotype: Parent zeotype to make dictionary from
        :type zeotype: Zeotype
        :return: Dict of zeotype's data
        :rtype: Dict
        """

        data = {'name': zeotype.name,
                'GUID': zeotype.unique_id,
                'type': 'Zeotype',
                'parent_index': 0,  # hardcode to 0 when self
                'additions': dict(zeotype.additions),
                'atom_indices_to_site': zeotype.atom_indices_to_site,
                'site_to_atom_indices': zeotype.site_to_atom_indices,
                'main_index': zeotype.index_mapper.main_index,
                'id': zeotype.index_mapper.id,
                'names': zeotype.index_mapper.names}

        return data

    def write_parent(self, zeotype: Zeotype) -> int:
        """
        Write parent zeotype to both databases
        :param zeotype: parent zeotype to write to database
        :type zeotype: Zeotype
        :return: id of added zeotype
        :rtype: int
        """

        data = self.make_parent_dict(zeotype)
        db_index = self._write(zeotype, data)
        self.parent_zeotype_dict[zeotype.unique_id] = db_index
        return db_index

    def update_parent(self, zeotype: Zeotype, key: int) -> None:
        """
        Update the parent zeotype
        :param zeotype: zeotype to add to the
        :type zeotype:
        :param key:
        :type key:
        :return:
        :rtype:
        """
        data = self.make_parent_dict(zeotype)
        self.add_data(key, data)

    def write(self, zeotype: Zeotype):
        if zeotype.parent_zeotype is zeotype:  # writing parent zeotype
            if zeotype.unique_id in self.parent_zeotype_dict.keys():
                db_index = self.parent_zeotype_dict[zeotype.unique_id]
                self.update_parent(zeotype, db_index)
            else:
                db_index = self.write_parent(zeotype)
                self.parent_zeotype_dict[zeotype.unique_id] = db_index
                return db_index

        else:
            if zeotype.parent_zeotype.unique_id in self.parent_zeotype_dict.keys():
                db_index = self.parent_zeotype_dict[zeotype.parent_zeotype.unique_id]
                self.update_parent(zeotype.parent_zeotype, db_index)
            else:
                db_index = self.write_parent(zeotype.parent_zeotype)
                self.parent_zeotype_dict[zeotype.parent_zeotype.unique_id] = db_index

            if type(zeotype).__name__ not in self.zeotype_dict.keys():
                raise ValueError(f'Type {type(zeotype).__name__} not in self.zeotype_dict')
            data = {'name': zeotype.name,
                    'type': type(zeotype).__name__,
                    'parent_index': db_index,
                    'additions': zeotype.additions}

            db_index = self._write(zeotype, data)
            return db_index

    def build_parent_zeotype(self, data):
        assert data['type'] == 'Zeotype', 'only use build_parent_zeotype for building Zeotypes'
        zeotype = self.zeotype_dict[data['type']](data['atoms'])
        zeotype.additions = data['additions']
        zeotype.index_mapper.main_index = data['main_index']
        zeotype.index_mapper.names = data['names']
        zeotype.index_mapper.id = data['id']
        zeotype.atom_indices_to_site = data['atom_indices_to_site']
        zeotype.site_to_atom_indices = data['site_to_atom_indices']

        return zeotype

    def get(self, key, default=None):
        data = self._get(key, default=None)
        if data is None:
            return default
        p_index = data['parent_index']
        if p_index == 0:
            return self.build_parent_zeotype(data)
        else:
            parent_zeotype = self.build_parent_zeotype(self._get(p_index))
            zeotype = self.zeotype_dict[data['type']](data['atoms'])
            zeotype.parent_zeotype = parent_zeotype
            zeotype.index_mapper = parent_zeotype.index_mapper
            zeotype.additions = data['additions']
            zeotype.name = data['name']
            return zeotype

    @classmethod
    def connect(cls, db_name: str) -> "ZeotypeDatabase":
        path = Path(db_name)
        if path.suffix != '.db':
            raise ValueError('db_name must end in .db')
        ase_db_path = db_name
        jsondata_path = path.with_suffix('json')
        return ZeotypeDatabase(ase_db_path, jsondata_path)


if __name__ == "__main__":
    import pandas as pd
    db = ZeotypeDatabase('test.db', 'test.json')
    bea = ImperfectZeotype.make('BEA')
    cluster = bea.get_cluster(cluster_indices=[0,2,13])
    index = db.write(bea)
    print(index)
    bea_2 = db.get(index)
    print(bea_2.parent_zeotype)
    print(pd.DataFrame(bea_2.index_mapper.main_index).T.head())
    print(bea_2)
    print(db.parent_zeotype_dict)