from maze import download_cif, build_zeolite_from_code


def main():
    my_dir = "/Users/dda/Code/zeotype/data"
    download_cif("OFF", data_dir=my_dir)
    download_cif("GOO", data_dir=my_dir)
    zeolite = build_zeolite_from_code("BEA", data_dir=my_dir)
    print(zeolite, zeolite.get_tags(), sep='\n')


if __name__ == '__main__':
    main()
