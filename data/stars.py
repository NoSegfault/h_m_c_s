import numpy as np
import sys
import itertools

# star, all edges point inside to center

def double_counting_list()
def main():
    filename = "C-elegans-frontal.txt"
    # the value of n is provided by the filename.txt
    n = 131
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
    bool fst = 1
    for item in v_from:
	for size in list(reversed(list(range(3,8)))):
		'''
        	els = [list(x) for x in itertools.combinations(v_from[item], size)]

        	for it in els:
            		mat = np.zeros(n)
            		mat[item] = 1
            		for i in range(len(it)):
                		mat[it[i]] = 1

		l.append(mat)
		'''
		if(fst):
        		els = [list(x) for x in itertools.combinations(v_from[item], size)]

        		for it in els:
            			mat = np.zeros(n)
            			mat[item] = 1
            			for i in range(len(it)):
               		 		mat[it[i]] = 1

			l.append(mat)
			d_c_l[size] = list(map(set,list(mat)))

		else:
			#filter key value that are less or equal to the current size
			s_l_1 = [x for x in d_c_l]
			s_l_2 = [x > size for x in s_l_1]
			s_l_1 = np.array(s_l_1)
			s_l_2 = np.array(s_l_2)
			s_l = list(s_l_1 * s_l_2)
			s_l = list(filter(lambda x: x != 0, s_l))

			els = [list(x) for x in itertools.combinations(v_from[item], size)]
			#filter els
			els_new = []
			for it1 in els:
				bool flag = 0
				for it2 in s_l:
					if(set(it1) in d_c_l[it2]):
						flag = 1

				if(!flag):
					els_new.append(it1)

			
						
        		for it in els_new:
            			mat = np.zeros(n)
            			mat[item] = 1
            			for i in range(len(it)):
               		 		mat[it[i]] = 1

			l.append(mat)
			d_c_l[size] = list(map(set,list(mat)))
				
		
		fst = 0


    # print(len(l))
    for row in l:
        for ele in row:
            print ele,
        print ""

    # l list of arrays, each array is length n


if __name__ == '__main__':
    main()
