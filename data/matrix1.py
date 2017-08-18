import numpy as np
import sys
import itertools

# star, all edges point outside from center

def main():
    filename = "C-elegans-frontal.txt"
    n = 131

    size = int(sys.argv[1])

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    v_to = {}
    for line in content[4:]:
        pair = line.split()
        if int(pair[0]) not in v_to:
            v_to[int(pair[0])]=[]
        v_to[int(pair[0])].append(int(pair[1]))


    l = []
    for item in v_to:
        print(item, v_to[item])
        els = [list(x) for x in itertools.combinations(v_to[item], size)]
        # print els

        for it in els:
            mat = np.zeros(n)
            mat[item] = 1
            for i in range(len(it)):
                mat[it[i]] = 1
            l.append(mat)

    #print(len(l))
    #print(l[:6])
    for row in l:
        for ele in row:
            print ele,
        print ""
    # l list of arrays, each array is length n


if __name__ == '__main__':
    main()
