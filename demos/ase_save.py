import copy
import glob
from ase.io import read, write
from maze.perfect_zeotype import Zeotype
import os
from pathlib import Path
import json
from source.index_mapper import IndexMapper
from typing import List
import shutil


def make_folder_to_zeo(folder_path):
    shutil.make_archive(folder_path, 'zip', folder_path)
    shutil.move(folder_path + '.zip', folder_path + '.zeo')


def delete_folder(folder_path):
    shutil.rmtree(folder_path)


def unpack_zeo_file(filename) -> str:
    file_dir = Path(filename).parents[0]
    file_stem = Path(filename).stem
    output_path = os.path.join(file_dir, file_stem)
    shutil.unpack_archive(filename, output_path, 'zip')
    return output_path


def save_zeotypes(folder_path: str, zeotype_list: List[Zeotype], ase_ext: str = '.traj'):
    assert '.' not in folder_path, 'do not add file extension when saving zeolites'
    my_path = Path(folder_path)
    my_path.mkdir(parents=True, exist_ok=True)
    name_list = []
    for z in zeotype_list:
        name_list.append(z.name)
        for names in z.additions.values():
            name_list.extend(names)

    name_list = [z.name for z in zeotype_list]
    assert "parent" in name_list, 'parent must be in zeotype list'
    new_index_mapper = copy.deepcopy(zeotype_list[0].index_mapper)

    for z in zeotype_list:
        assert z.index_mapper is zeotype_list[0].index_mapper, 'index mappers must match'

    # save index mapper
    index_mapper_path = os.path.join(my_path, 'index_mapper.json')
    for name in new_index_mapper.names:
        if name not in name_list:
            new_index_mapper.delete_name(name)
    save_index_mapper(index_mapper_path, new_index_mapper)

    for z in zeotype_list:
        zeotype_folder = os.path.join(my_path, z.name)
        Path(zeotype_folder).mkdir(parents=True, exist_ok=True)
        dict_json = {'name': z.name,
                     'type': str(type(z)),
                     'additions': z.additions}
        if z.name == 'parent':
            additional_params = {'site_to_atom_indices': z._site_to_atom_indices,
                                 'atom_indices_to_site': z._atom_indices_to_site}
            dict_json.update(additional_params)

        # TODO: impliment save index mapper
        # save things
        binary_path = os.path.join(zeotype_folder, z.name + ase_ext)
        dict_path = os.path.join(zeotype_folder, z.name + '.json')
        write(binary_path, z)
        with open(dict_path, 'w') as f:
            json.dump(dict_json, f, indent=4, ensure_ascii=True)

        make_folder_to_zeo(folder_path)
        delete_folder(folder_path)


def save_index_mapper(filepath, index_mapper: IndexMapper) -> None:
    json_dict = {'main_index': index_mapper.main_index,
                 'id': index_mapper.id,
                 'names': index_mapper.names}
    with open(filepath, 'w') as f:
        json.dump(json_dict, f, indent=4, ensure_ascii=True)


def load_index_mapper(filepath) -> IndexMapper:
    with open(filepath, 'r') as f:
        index_mapper_json = json.load(f)
    new_im = IndexMapper([0, 1, 2])
    new_im.main_index = index_mapper_json['main_index']
    if new_im.id < index_mapper_json['id']:
        new_im.id = index_mapper_json['id'] + 1
    new_im.names = index_mapper_json['names']
    return new_im


def read_zeotypes(file_path, str_ext: str = '.traj'):
    if '.' not in file_path:
        file_path = file_path + '.zeo'
    folder_path = unpack_zeo_file(file_path)
    zeotype_dict = {}
    folder_list = glob.glob(os.path.join(folder_path, '*/'))
    for folder in folder_list:
        name = Path(folder).stem
        json_path = os.path.join(folder, name + '.json')
        binary_path = os.path.join(folder, name + str_ext)
        with open(json_path, 'r') as f:
            my_zeotype = Zeotype(read(binary_path))
            attr_dict = json.load(f)
            my_zeotype.name = name
            my_zeotype.additions = attr_dict['additions']

            if name == 'parent':
                my_zeotype._atom_indices_to_site = attr_dict['atom_indices_to_site']
                my_zeotype._site_to_atom_indices = attr_dict['site_to_atom_indices']

            zeotype_dict[name] = my_zeotype

    for name, z in zeotype_dict.items():
        z.parent_zeotype = zeotype_dict['parent']
        z.index_mapper = zeotype_dict['parent'].index_mapper

    delete_folder(folder_path)
    return zeotype_dict


if __name__ == "__main__":
    data_path = '/Users/dda/Code/MAZE-sim/data/'
    cif_path = os.path.join(data_path, 'BEA.cif')
    output_path = os.path.join(data_path, 'my_first_zeotype')
    my_zeotype = Zeotype.build_from_cif_with_labels(cif_path)
    cluster, od = my_zeotype.get_cluster(1,10,10)
    # data_dir = os.path.join(Path(os.getcwd()).parent, 'data')
    # output_traj = os.path.join(data_dir, 'test.traj')
    # cif_dir = os.path.join(data_dir, 'BEA.cif')
    # zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    save_zeotypes(output_path, [my_zeotype, cluster, od])
    new_dict = read_zeotypes(output_path + '.zeo')
    print(new_dict)
    print(new_dict['parent'])
    print(new_dict['parent'].index_mapper.main_index)
    # Trajectory(output_traj,'w', zeolite)
