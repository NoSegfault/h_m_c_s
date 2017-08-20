import numpy as np
import sys
import itertools
import csv

# motif M6

def main():
    filename = "Florida-bay.txt"
    n = 128

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    v_to = {}
    for i in range(n):
        v_to[i] = []
    for line in content[5:]:
        pair = line.split()
        v_to[int(pair[0])].append(int(pair[1]))

    l = []
    mark = []
    for item in v_to:
        #print(item, v_to[item])
        els = [list(x) for x in itertools.combinations(v_to[item], 2)]
        for it in els:
            #print item,it[0],it[1]
            if it[0] in v_to[it[1]] and it[1] in v_to[it[0]] and item not in v_to[it[0]] and item not in v_to[it[1]]:
                mat = np.zeros(n)
                mat[item] = 1
                mat[it[0]] = 1
                mat[it[1]] = 1
                l.append(mat)

                sub = [item, it[0], it[1]]
                mark.append(sub)

    #print len(l)
    #print mark

    set_list = []
    for row in mark:
        find = 0
        for s in set_list:
            for it in row:
                if it in s:
                    s.update(row)
                    find = 1
                    break
            if find == 1:
                break
        if find == 0:
            #need a new set
            #print("need a new set")
            #print(row)
            d = set()
            d.update(row)
            set_list.append(d)

    #print('set_list before merge is: ')
    #print(set_list)

    for i in range(0, len(set_list)):
        for j in range(i+1, len(set_list)):
            if len(set_list[i].intersection(set_list[j])) != 0 :
                set_list[i] =  set_list[i].union(set_list[j])
                set_list[j] = set()

    set_list.sort(key=len, reverse=True)
    #print set_list
    new_n = len(set_list[0])
    print('size of lcc is {}'.format(new_n))
    lcc = []
    lcc.append(list(set_list[0]))
    #print lcc
    for i in range(len(l)):
        if mark[i][0] in set_list[0]:
            #print mark[i]
            row = [l[i][j] for j in set_list[0]]
            lcc.append(row)
    #print lcc
    print('number of hyperedges is {}'.format(len(lcc)-1))


    with open('Florida-m6.csv','w') as f:
        writer = csv.writer(f)
        writer.writerows(lcc);
        f.close()



if __name__ == '__main__':
    main()
