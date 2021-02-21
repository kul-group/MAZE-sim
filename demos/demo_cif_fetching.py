from maze import zeotypes as zeo
import copy


if __name__ == '__main__':
    cha = zeo.ImperfectZeotype.make('CHA')
    cha4 = zeo.ImperfectZeotype(cha)
    print(cha.index_mapper.main_index)
    cha2 = copy.copy(cha)
    print(cha2.index_mapper.main_index)
    print(cha.index_mapper.main_index)
    print(cha.index_mapper==cha2.index_mapper)
