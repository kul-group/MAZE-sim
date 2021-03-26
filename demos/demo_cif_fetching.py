from maze import zeolite as zeo
import copy


if __name__ == '__main__':
    cha = zeo.Zeolite.make('CHA')
    cha4 = zeo.Zeolite(cha)
    print(cha.index_mapper.main_index)
    cha2 = copy.copy(cha)
    print(cha2.index_mapper.main_index)
    print(cha.index_mapper.main_index)
    print(cha.index_mapper==cha2.index_mapper)
