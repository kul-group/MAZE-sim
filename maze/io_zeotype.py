import copy
import glob
from ase.io import read, write
from maze.zeotypes import Zeotype
import os
from pathlib import Path
import json
from maze.index_mapper import IndexMapper
from typing import List, Iterable, Dict
import shutil


def make_folder_to_zeo(folder_path: str) -> None:
    """
    Creates a zip archive with a .zeo extension
    :param folder_path: path of folder to be zipped
    :type folder_path: str
    :return: None
    :rtype: None
    """

    shutil.make_archive(folder_path, 'zip', folder_path)
    shutil.move(folder_path + '.zip', folder_path + '.zeo')


def delete_folder(folder_path: str) -> None:
    """
    Deletes a folder
    :param folder_path: folder path to delete
    :type folder_path: str
    :return: None
    :rtype: None
    """

    shutil.rmtree(folder_path)


def unpack_zeo_file(filename) -> str:
    """
    unpacks a zeo file
    :param filename: name of zeo file with .zeo extension
    :type filename: str
    :return: path to unpacked file
    :rtype: str
    """

    file_dir = Path(filename).parents[0]
    file_stem = Path(filename).stem
    output_path = os.path.join(file_dir, file_stem)
    shutil.unpack_archive(filename, output_path, 'zip')
    return output_path


def save_zeotypes(folder_path: str, zeotype_list: Iterable[Zeotype], ase_ext: str = '.traj') -> None:
    """
    This saves a list of Zeotypes to a .zeo file
    :param folder_path: path and name of the zeo file without the zeo extension
    :type folder_path: str
    :param zeotype_list: an iterable collection of Zeotype objects (or subclasses)
    :type zeotype_list: Iterable[Zeotype]
    :param ase_ext: extension to save the Zeotype files (default .traj)
    :type ase_ext: str
    :return: None
    :rtype: None
    """

    assert '.' not in folder_path, 'do not add file extension when saving zeolites'
    my_path = Path(folder_path)
    my_path.mkdir(parents=True, exist_ok=True)
    name_list = []
    for z in zeotype_list:
        name_list.append(z.name)
        for names in z.additions.values():
            name_list.extend(names)

    print(name_list)
    assert "parent" in name_list, 'parent must be in zeotype list'
    new_index_mapper = copy.deepcopy(zeotype_list[0].index_mapper)

    for z in zeotype_list:
        assert z.index_mapper is zeotype_list[0].index_mapper, 'index mappers must match'

    # save index mapper
    index_mapper_path = os.path.join(my_path, 'index_mapper.json')
    for name in zeotype_list[0].index_mapper.names:
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

        binary_path = os.path.join(zeotype_folder, z.name + ase_ext)
        dict_path = os.path.join(zeotype_folder, z.name + '.json')
        write(binary_path, z)
        with open(dict_path, 'w') as f:
            json.dump(dict_json, f, indent=4, ensure_ascii=True)

    make_folder_to_zeo(folder_path)
    delete_folder(folder_path)


def save_index_mapper(filepath, index_mapper: IndexMapper) -> None:
    """
    Saves an index mapper to json file format
    :param filepath: json filepath
    :type filepath: str
    :param index_mapper: IndexMapper object to save
    :type index_mapper: IndexMapper
    :return: None
    :rtype: None
    """

    json_dict = {'main_index': index_mapper.main_index,
                 'id': index_mapper.id,
                 'names': index_mapper.names}
    with open(filepath, 'w') as f:
        json.dump(json_dict, f, indent=4, ensure_ascii=True)


def load_index_mapper(filepath) -> IndexMapper:
    """
    Loads an IndexMapper object from a json file
    :param filepath: filepath to load IndexMapper from
    :type filepath: str
    :return: loaded IndexMapper
    :rtype: IndexMapper
    """

    with open(filepath, 'r') as f:
        index_mapper_json = json.load(f)
    new_main_dict = {}
    for key, name_dict in index_mapper_json['main_index'].items():
        new_name_dict = {}
        for inner_key, value in name_dict.items():
            if value is None:
                new_name_dict[inner_key] = None
            else:
                new_name_dict[inner_key] = int(value)
        new_main_dict[int(key)] = new_name_dict

    new_im = IndexMapper([0, 1, 2])
    new_im.main_index = new_main_dict
    if new_im.id < index_mapper_json['id']:
        new_im.id = index_mapper_json['id'] + 1
    new_im.names = index_mapper_json['names']
    return new_im


def read_zeotypes(file_path: str, str_ext: str = '.traj') -> Dict[str, Zeotype]:
    """
    Read the zeotypes from a .zeo file
    :param file_path: path to .zeo file with or without .zeo extension
    :type file_path: str
    :param str_ext: type of files in .zeo zip file
    :type str_ext: str
    :return: Dictionary of Zeotypes loaded from the files
    :rtype: Dict[str, Zeotype]
    """

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
    index_mapper = load_index_mapper(os.path.join(folder_path, 'index_mapper.json'))
    for name, z in zeotype_dict.items():
        z.parent_zeotype = zeotype_dict['parent'].index_mapper
        z.index_mapper = index_mapper

    delete_folder(folder_path)
    return zeotype_dict


