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
    d_c_l = {}
    fst = True
    #print('v_to')
    #print(v_to)


    for size in list(reversed(list(range(3,21)))):
        d_c_l[size] = []
        for item in v_to:
            #print('size is {}'.format(size))
            #print(l)
            #print(d_c_l)
            print()
            if(fst):

                    els = [list(x) for x in itertools.combinations(v_to[item], size)]
                    #print('els1')
                    for it in els:
                        mat = np.zeros(n)
                        mat[item] = 1
                        for i in range(len(it)):
                            mat[it[i]] = 1

                        l.append(mat)
                        #d_c_l[size] = list(map(set,list(mat)))
                        d_c_l[size].append( (np.where(mat != 0)[0]).tolist() )
 
            else:
                #filter key value that are less or equal to the current size
                s_l_1 = [x for x in d_c_l]
                s_l_2 = [x > size for x in s_l_1]
                s_l_1 = np.array(s_l_1)
                s_l_2 = np.array(s_l_2)
                s_l = (s_l_1 * s_l_2).tolist()
                s_l = list(filter(lambda x: x != 0, s_l))

                print('works here ?')    

                els = [list(x) for x in itertools.combinations(v_to[item], size)]
                #filter els
                #print('els2')
                #print(els)
                els_new = []
                for it1 in els:
                    #print('how about here?')
                    flag = 0
                    for it2 in s_l:
                        '''
                                if(set(it1) in d_c_l[it2]):
                                flag = 1
                        '''
                        for i in range(len(d_c_l[it2])):
                            print(it1)
                            print(d_c_l)
                            print(d_c_l[it2][i])
                            if(set(it1).issubset(d_c_l[it2][i])):
                                flag = 1
    
                    if(not(flag)):
                        els_new.append(it1)
    
                
                            
                    for it in els_new:
                        mat = np.zeros(n)
                        mat[item] = 1
                        for i in range(len(it)):
                            mat[it[i]] = 1
    
                        l.append(mat)
                        d_c_l[size].append( (np.where(mat != 0)[0]).tolist() )
                    
            
            fst = 0

    print(l)
    with open('toy1_output.csv','w') as f:
        writer = csv.writer(f)
        writer.writerows(l)


if __name__=='__main__':
    main()
