import numpy as np
import sys
import itertools
import csv

# star, all edges point inside to center

def main():
    filename = "toy1.txt"
    # the value of n is provided by the filename.txt
    n = 16
    # the value of num_comment is the line of comment at the beginning of the filename
    num_comment = 4


    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]

    v_from = {}

    for line in content[4:]:
        pair = line.split()
        if int(pair[1]) not in v_from:
            v_from[int(pair[1])]=[]
        v_from[int(pair[1])].append(int(pair[0]))


    v_to = {}

    for line in content[4:]:
        pair = line.split()
        if int(pair[0]) not in v_to:
            v_to[int(pair[0])]=[]
        v_to[int(pair[0])].append(int(pair[1]))


    l = []

    del_v_to = []
    for size in list(reversed(list(range(3,21)))):
    #remove list that is already used once
        if(len(del_v_to) != 0):
            #print(del_v_to)
            for item in del_v_to:
                #print(item)
                #print(v_to)
                del v_to[item]
                del_v_to.remove(item)

        for item in v_to:

            els = [list(x) for x in itertools.combinations(v_to[item], size)]
            #double counting for stars
            if(len(els) != 0): del_v_to.append(item)
            
                #print('els1')
            for it in els:
                mat = np.zeros(n)
                mat[item] = 1
                for i in range(len(it)):
                    mat[it[i]] = 1

                l.append(mat)


    print(l)
    with open('toy1_output_s.csv','w') as f:
        writer = csv.writer(f)
        writer.writerows(l)


if __name__=='__main__':
    main()
