import numpy as np
import sys

# star, all edges point inside to center

def main():
    filename = "Florida-bay.txt"
    global n
    n = 128

    size = int(sys.argv[1])

    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    v_from = {}
    for line in content[5:]:
        pair = line.split()
        if int(pair[1]) not in v_from:
            v_from[int(pair[1])]=[]
        v_from[int(pair[1])].append(int(pair[0]))



    global l
    l = []
    for item in v_from:
        print(item, v_from[item])
        if size == 3:
            three(item, v_from[item])
        elif size == 4:
            four(item, v_from[item])
        elif size == 5:
            five(item, v_from[item])
        elif size == 6:
            six(item, v_from[item])
        elif size == 7:
            seven(item, v_from[item])
        elif size == 8:
            eight(item, v_to[item])
        elif size == 9:
            nine(item, v_to[item])
        elif size == 10:
            ten(item, v_to[item])
        elif size == 11:
            eleven(item, v_to[item])
        elif size == 12:
            twelve(item, v_to[item])
        elif size == 13:
            thirteen(item, v_to[item])
        elif size == 14:
            fourteen(item, v_to[item])
        elif size == 15:
            fifteen(item, v_to[item])

    print(len(l))
    #for row in l:
        #for ele in row:
            #print ele,
        #print ""
    # l list of arrays, each array is length n

def three(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                mat = np.zeros(n)
                mat[item] = 1
                mat[item_list[i]] = 1
                mat[item_list[j]] = 1
                mat[item_list[k]] = 1
                l.append(mat)

def four(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                for ii in range(k+1, len(item_list)):
                    mat = np.zeros(n)
                    mat[item] = 1
                    mat[item_list[i]] = 1
                    mat[item_list[j]] = 1
                    mat[item_list[k]] = 1
                    mat[item_list[ii]] = 1
                    l.append(mat)

def five(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                for ii in range(k+1, len(item_list)):
                    for jj in range(ii+1, len(item_list)):
                        mat = np.zeros(n)
                        mat[item] = 1
                        mat[item_list[i]] = 1
                        mat[item_list[j]] = 1
                        mat[item_list[k]] = 1
                        mat[item_list[ii]] = 1
                        mat[item_list[jj]] = 1
                        l.append(mat)

def six(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                for ii in range(k+1, len(item_list)):
                    for jj in range(ii+1, len(item_list)):
                        for kk in range(jj+1, len(item_list)):
                            mat = np.zeros(n)
                            mat[item] = 1
                            mat[item_list[i]] = 1
                            mat[item_list[j]] = 1
                            mat[item_list[k]] = 1
                            mat[item_list[ii]] = 1
                            mat[item_list[jj]] = 1
                            mat[item_list[kk]] = 1
                            l.append(mat)

def seven(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                for ii in range(k+1, len(item_list)):
                    for jj in range(ii+1, len(item_list)):
                        for kk in range(jj+1, len(item_list)):
                            for i3 in range(kk+1, len(item_list)):
                                mat = np.zeros(n)
                                mat[item] = 1
                                mat[item_list[i]] = 1
                                mat[item_list[j]] = 1
                                mat[item_list[k]] = 1
                                mat[item_list[ii]] = 1
                                mat[item_list[jj]] = 1
                                mat[item_list[kk]] = 1
                                mat[item_list[i3]] = 1
                                l.append(mat)

def eight(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                for ii in range(k+1, len(item_list)):
                    for jj in range(ii+1, len(item_list)):
                        for kk in range(jj+1, len(item_list)):
                            for i3 in range(kk+1, len(item_list)):
                                for j3 in range(i3+1, len(item_list)):
                                    mat = np.zeros(n)
                                    mat[item] = 1
                                    mat[item_list[i]] = 1
                                    mat[item_list[j]] = 1
                                    mat[item_list[k]] = 1
                                    mat[item_list[ii]] = 1
                                    mat[item_list[jj]] = 1
                                    mat[item_list[kk]] = 1
                                    mat[item_list[i3]] = 1
                                    mat[item_list[j3]] = 1
                                    l.append(mat)

def nine(item, item_list):
    for i in range(0, len(item_list)):
        for j in range(i+1, len(item_list)):
            for k in range(j+1, len(item_list)):
                for ii in range(k+1, len(item_list)):
                    for jj in range(ii+1, len(item_list)):
                        for kk in range(jj+1, len(item_list)):
                            for i3 in range(kk+1, len(item_list)):
                                for j3 in range(i3+1, len(item_list)):
                                    for k3 in range(j3+1, len(item_list)):
                                        mat = np.zeros(n)
                                        mat[item] = 1
                                        mat[item_list[i]] = 1
                                        mat[item_list[j]] = 1
                                        mat[item_list[k]] = 1
                                        mat[item_list[ii]] = 1
                                        mat[item_list[jj]] = 1
                                        mat[item_list[kk]] = 1
                                        mat[item_list[i3]] = 1
                                        mat[item_list[j3]] = 1
                                        mat[item_list[k3]] = 1
                                        l.append(mat)

def ten(item, item_list):
    end = len(item_list)
    for i in range(0, end):
        for j in range(i+1, end):
            for k in range(j+1, end):
                for ii in range(k+1, end):
                    for jj in range(ii+1, end):
                        for kk in range(jj+1, end):
                            for i3 in range(kk+1, end):
                                for j3 in range(i3+1, end):
                                    for k3 in range(j3+1, end):
                                        for i4 in range(k3+1, end):
                                            mat = np.zeros(n)
                                            mat[item] = 1
                                            mat[item_list[i]] = 1
                                            mat[item_list[j]] = 1
                                            mat[item_list[k]] = 1
                                            mat[item_list[ii]] = 1
                                            mat[item_list[jj]] = 1
                                            mat[item_list[kk]] = 1
                                            mat[item_list[i3]] = 1
                                            mat[item_list[j3]] = 1
                                            mat[item_list[k3]] = 1
                                            mat[item_list[i4]] = 1
                                            l.append(mat)

def eleven(item, item_list):
    end = len(item_list)
    for i in range(0, end):
        for j in range(i+1, end):
            for k in range(j+1, end):
                for ii in range(k+1, end):
                    for jj in range(ii+1, end):
                        for kk in range(jj+1, end):
                            for i3 in range(kk+1, end):
                                for j3 in range(i3+1, end):
                                    for k3 in range(j3+1, end):
                                        for i4 in range(k3+1, end):
                                            for j4 in range(i4+1, end):
                                                mat = np.zeros(n)
                                                mat[item] = 1
                                                mat[item_list[i]] = 1
                                                mat[item_list[j]] = 1
                                                mat[item_list[k]] = 1
                                                mat[item_list[ii]] = 1
                                                mat[item_list[jj]] = 1
                                                mat[item_list[kk]] = 1
                                                mat[item_list[i3]] = 1
                                                mat[item_list[j3]] = 1
                                                mat[item_list[k3]] = 1
                                                mat[item_list[i4]] = 1
                                                mat[item_list[j4]] = 1
                                                l.append(mat)

def twelve(item, item_list):
    end = len(item_list)
    for i in range(0, end):
        for j in range(i+1, end):
            for k in range(j+1, end):
                for ii in range(k+1, end):
                    for jj in range(ii+1, end):
                        for kk in range(jj+1, end):
                            for i3 in range(kk+1, end):
                                for j3 in range(i3+1, end):
                                    for k3 in range(j3+1, end):
                                        for i4 in range(k3+1, end):
                                            for j4 in range(i4+1, end):
                                                for k4 in range(j4+1, end):
                                                    mat = np.zeros(n)
                                                    mat[item] = 1
                                                    mat[item_list[i]] = 1
                                                    mat[item_list[j]] = 1
                                                    mat[item_list[k]] = 1
                                                    mat[item_list[ii]] = 1
                                                    mat[item_list[jj]] = 1
                                                    mat[item_list[kk]] = 1
                                                    mat[item_list[i3]] = 1
                                                    mat[item_list[j3]] = 1
                                                    mat[item_list[k3]] = 1
                                                    mat[item_list[i4]] = 1
                                                    mat[item_list[j4]] = 1
                                                    mat[item_list[k4]] = 1
                                                    l.append(mat)

def thirteen(item, item_list):
    end = len(item_list)
    for i in range(0, end):
        for j in range(i+1, end):
            for k in range(j+1, end):
                for ii in range(k+1, end):
                    for jj in range(ii+1, end):
                        for kk in range(jj+1, end):
                            for i3 in range(kk+1, end):
                                for j3 in range(i3+1, end):
                                    for k3 in range(j3+1, end):
                                        for i4 in range(k3+1, end):
                                            for j4 in range(i4+1, end):
                                                for k4 in range(j4+1, end):
                                                    for i5 in range(k4+1, end):
                                                        mat = np.zeros(n)
                                                        mat[item] = 1
                                                        mat[item_list[i]] = 1
                                                        mat[item_list[j]] = 1
                                                        mat[item_list[k]] = 1
                                                        mat[item_list[ii]] = 1
                                                        mat[item_list[jj]] = 1
                                                        mat[item_list[kk]] = 1
                                                        mat[item_list[i3]] = 1
                                                        mat[item_list[j3]] = 1
                                                        mat[item_list[k3]] = 1
                                                        mat[item_list[i4]] = 1
                                                        mat[item_list[j4]] = 1
                                                        mat[item_list[k4]] = 1
                                                        mat[item_list[i5]] = 1
                                                        l.append(mat)



def fourteen(item, item_list):
    end = len(item_list)
    for i in range(0, end):
        for j in range(i+1, end):
            for k in range(j+1, end):
                for ii in range(k+1, end):
                    for jj in range(ii+1, end):
                        for kk in range(jj+1, end):
                            for i3 in range(kk+1, end):
                                for j3 in range(i3+1, end):
                                    for k3 in range(j3+1, end):
                                        for i4 in range(k3+1, end):
                                            for j4 in range(i4+1, end):
                                                for k4 in range(j4+1, end):
                                                    for i5 in range(k4+1, end):
                                                        for j5 in range(i5+1, end):
                                                            mat = np.zeros(n)
                                                            mat[item] = 1
                                                            mat[item_list[i]] = 1
                                                            mat[item_list[j]] = 1
                                                            mat[item_list[k]] = 1
                                                            mat[item_list[ii]] = 1
                                                            mat[item_list[jj]] = 1
                                                            mat[item_list[kk]] = 1
                                                            mat[item_list[i3]] = 1
                                                            mat[item_list[j3]] = 1
                                                            mat[item_list[k3]] = 1
                                                            mat[item_list[i4]] = 1
                                                            mat[item_list[j4]] = 1
                                                            mat[item_list[k4]] = 1
                                                            mat[item_list[i5]] = 1
                                                            mat[item_list[j5]] = 1
                                                            l.append(mat)

def fifteen(item, item_list):
    end = len(item_list)
    for i in range(0, end):
        for j in range(i+1, end):
            for k in range(j+1, end):
                for ii in range(k+1, end):
                    for jj in range(ii+1, end):
                        for kk in range(jj+1, end):
                            for i3 in range(kk+1, end):
                                for j3 in range(i3+1, end):
                                    for k3 in range(j3+1, end):
                                        for i4 in range(k3+1, end):
                                            for j4 in range(i4+1, end):
                                                for k4 in range(j4+1, end):
                                                    for i5 in range(k4+1, end):
                                                        for j5 in range(i5+1, end):
                                                            for k5 in range(j5+1, end):
                                                                mat = np.zeros(n)
                                                                mat[item] = 1
                                                                mat[item_list[i]] = 1
                                                                mat[item_list[j]] = 1
                                                                mat[item_list[k]] = 1
                                                                mat[item_list[ii]] = 1
                                                                mat[item_list[jj]] = 1
                                                                mat[item_list[kk]] = 1
                                                                mat[item_list[i3]] = 1
                                                                mat[item_list[j3]] = 1
                                                                mat[item_list[k3]] = 1
                                                                mat[item_list[i4]] = 1
                                                                mat[item_list[j4]] = 1
                                                                mat[item_list[k4]] = 1
                                                                mat[item_list[i5]] = 1
                                                                mat[item_list[j5]] = 1
                                                                mat[item_list[k5]] = 1
                                                                l.append(mat)



if __name__ == '__main__':
    main()
