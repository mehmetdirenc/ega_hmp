


def compare_kirsten(our_csv, db_csv, to_write):
    db_dict = {}
    with open(db_csv, 'r') as db:
        for line in db:
            split = line.split(';')
            id = split[0]
            if id == "id":
                continue
            new_id = id[id.index('_'):][1:]
            db_dict[new_id] = line
    x = 0
    with open(our_csv, 'r') as ours:
        with open(to_write, 'w') as w:
            for line in ours:
                if line == ';;;;;;;;;;;;\n' or line == ';;;;;;;;;;;;':
                    continue
                to_search = line.split(';')[2] + "_" + line.split(';')[4]
                if to_search not in db_dict.keys():
                    x += 1
                    w.write(line)
                    # print(line)

    print(x)

if __name__ == '__main__':
    our_csv = "/home/direnc/Downloads/Kopie von dual_LM_final.csv"
    db_csv = "/home/direnc/Downloads/Zuber_dual_mouse_library.csv"
    to_write = "/home/direnc/Downloads/kristen_not_found.csv"
    compare_kirsten(our_csv, db_csv, to_write)