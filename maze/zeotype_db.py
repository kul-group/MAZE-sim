from ase.db import connect
from typing import Union, Dict
from ase.db.core import Database
from maze.zeotypes import Zeotype, Cluster, OpenDefect, ImperfectZeotype
import json
from copy import deepcopy
import ase
from pathlib import Path

class ZeotypeDatabase:
    zeotype_dict = {'Zeotype': Zeotype, 'Cluster': Cluster, 'OpenDefect': OpenDefect,
                    'ImperfectZeotype': ImperfectZeotype}

    def __init__(self, ase_db_name: str, datajson_name: str):
        self.parent_zeotype_dict: Dict[str, int] = {}
        self.data_dict = {}
        self.ase_db = None
        self.json_filepath = datajson_name if 'json' in datajson_name else datajson_name + '.json'
        self.ase_db = connect(ase_db_name)
        json_path = Path(datajson_name)
        if json_path.is_file():
            with open(self.json_filepath, 'r') as f:
                self.data_dict = json.load(f)
        else:
            self.data_dict = {}

    def add_data(self, key, new_data):
        self.data_dict[key] = new_data
        self.write_data()

    def _get(self, key, default=None):
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

    def _write(self, zeotype: Zeotype, data: Dict):
        key = self.ase_db.write(zeotype)
        self.add_data(key, data)
        return key

    def write_data(self):
        with open(self.json_filepath, 'w') as f:
            json.dump(self.data_dict, f, indent=4)


    def make_parent_dict(self, zeotype: Zeotype) -> Dict:
        data = {'name': zeotype.name,
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
        data = self.make_parent_dict(zeotype)
        db_index = self._write(zeotype, data)
        self.parent_zeotype_dict[zeotype.unique_id] = db_index
        return db_index

    def update_parent(self, zeotype: Zeotype, key: int) -> None:
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