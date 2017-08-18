import numpy as np
import sys
import itertools

# star, all edges point inside to center

def main():
    filename = "C-elegans-frontal.txt"
    n = 131

    size = int(sys.argv[1])

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    v_from = {}
    for line in content[4:]:
        pair = line.split()
        if int(pair[1]) not in v_from:
            v_from[int(pair[1])]=[]
        v_from[int(pair[1])].append(int(pair[0]))


    l = []
    for item in v_from:
        # print(item, v_from[item])
        els = [list(x) for x in itertools.combinations(v_from[item], size)]
        # print els

        for it in els:
            mat = np.zeros(n)
            mat[item] = 1
            for i in range(len(it)):
                mat[it[i]] = 1
            l.append(mat)


    # print(len(l))
    for row in l:
        for ele in row:
            print ele,
        print ""
    # l list of arrays, each array is length n


if __name__ == '__main__':
    main()
