import numpy as np
import sys
import itertools

# motif M6

def main():
    filename = "Florida-bay.txt"
    n = 128

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    v_to = {}
    for line in content[5:]:
        pair = line.split()
        if int(pair[0]) not in v_to:
            v_to[int(pair[0])]=[]
        v_to[int(pair[0])].append(int(pair[1]))


    l = []
    for item in v_to:
        # print(item, v_to[item])
        els = [list(x) for x in itertools.combinations(v_to[item], 2)]
        for it in els:
            if it[1] in v_to and it[1] in v_to:
                if it[0] in v_to[it[1]] and it[1] in v_to[it[0]]:
                    mat = np.zeros(n)
                    mat[item] = 1
                    mat[it[0]] = 1
                    mat[it[1]] = 1
                    l.append(mat)

    # print len(l)
    # print l[:10]
    for row in l:
        for ele in row:
            print ele,
        print ""

if __name__ == '__main__':
    main()
