import numpy as np
from IDPP.mgraph import MGraph

def get_duplicates(items):
    ''' 
    returns the indices of all duplicates in a list
    @param items a 1D list
    '''
    locs = []

    for i in range(len(items)):
        loc_tmp = [i]
        item = items[i]
        start_at = i+1
        while True:
            try:
                loc = items.index(item,start_at)
            except ValueError:
                break
            else:
                loc_tmp.append(loc)
                start_at = loc + 1
        
        if len(loc_tmp)>1:
            locs.append(loc_tmp)
            
    #print("locs: ", locs)
    return(locs)

def remove_duplicates(m_prime, c_prime):
    '''
    remove any duplicate columns and merge associated features
    @param m_prime M' matrix
    @param c_prime features list corresponding to columns of M'
    '''
    #print(np.rot90(m_prime)[::-1])
    # add list()
    m_prime_tmp = list(map(lambda x: '.'.join(map(str,x)),np.rot90(m_prime)[::-1]))
    #print(m_prime_tmp)
    dups = get_duplicates(m_prime_tmp)

    n_dup = len(dups)
    to_del = []
    if n_dup > 0:
        for idx,dup in enumerate(dups):
            to_del.extend(dup[1:])
            c_prime[dup[0]] = b'_'.join(c_prime[dup])

    m_prime = np.delete(m_prime,to_del,axis=1)
    c_prime = np.delete(c_prime,to_del)

    return(m_prime, c_prime)

def remove_S_semi_universal(m3,c3):
    '''
    remove any columns that have no 0 entries
    @param m3 the M matrix
    @param c3 column features corresponding to the M matrix
    '''
    mi = np.empty(0,dtype='int')
    ncol = len(m3[0])
    idxs_to_delete = [i for i in range(ncol) if not np.any(m3[:,i]==0)]
    m3 = np.delete(m3,idxs_to_delete,axis=1)
    c3 = np.delete(c3,idxs_to_delete)
    return(m3,c3)

def m_column_sort(m_prime, c_prime):
    m_prime = np.rot90(m_prime)
    
    # count binary score of columns
    binary_strings = []
    for col in m_prime:
        # change all -1 to 0
        col = np.array([ci if ci>0 else 0  for ci in col])
        col_string = '0b'+''.join(map(str,col))
        binary_strings.append(int(col_string,2))
        
    # sort by binary score
    order = np.argsort(binary_strings)[::-1]
    m_prime = m_prime[order] 
    m_prime = np.rot90(m_prime)[::-1] #rotate again
    c_order = (len(c_prime) - 1) - order #translate order of rotated matrix to order of columns
    c_prime = c_prime[c_order]
    return m_prime,c_prime

def m_initial(m,c):
    '''
    Construct M' matrix, remove S_semi_universal c and duplicates, sort the column
    '''
    m_prime, c_prime = remove_S_semi_universal(m,c)
    m_prime, c_prime = remove_duplicates(m_prime, c_prime)
    m_prime, c_prime = m_column_sort(m_prime, c_prime)
    
    return m_prime, c_prime

def get_k1_matrix(m_prime,features):
    '''
    Generate k1 from m' matrix
    Allows for checking of perfect phylogeny
    @param mp the m prime matrix
    @param nodes corresponding to the matrix
    '''
    ncol = len(m_prime[0])
    k = np.empty( [0,ncol], dtype='|S15' )

    for m in m_prime:
        row_feats = features[m!=0] #features in the row
        mrow = np.zeros(ncol,dtype='|S15')
        mrow.fill('0')

        for idx,feature in enumerate(row_feats):
            mrow[idx] = feature

        n_feat = len(row_feats)    
        if n_feat < ncol: 
            mrow[n_feat]='#'

        k = np.append(k,[mrow],axis=0)

    return(k)

def perfect_phylogeny_exists(k1,features):
    '''
    Determine whether perfect phylogeny exists from a k1 matrix    
    @param k1 the k1 matrix (output from get_k1_matrix)
    @param features the column features
    '''
    locations = []
    for feature in features:
        present_at = set([])
        for k_i in k1:
            [ present_at.add(loc_list) for loc_list in list(np.where(k_i==feature)[0]) ]
        #print("present at",present_at)
        locations.append(present_at)
    
    print("!!!!!!!locations:",locations)
    loc_test = np.array([len(loc_list)>1 for loc_list in locations])
    if np.any(loc_test):
        print('No phylogeny found!')
        return(False)
    else:    
        print('Success! Found phylogeny!\nK1 matrix:')
        print(k1)
        return(True)

def perfect_phylogeny_exists_new(k1,features):
    '''
    Determine whether perfect phylogeny exists from a k1 matrix    
    @param k1 the k1 matrix (output from get_k1_matrix)
    @param features the column features
    '''
    #locations = []
    for feature in features:
        present_at = []
        row_at = []
        for row, k_i in enumerate(k1):
            if len(list(np.where(k_i==feature)[0])) != 0:
                present_at.append(np.where(k_i==feature)[0][0])
                row_at.append(row)
            elif len(list(np.where(k_i==feature)[0])) > 1:
                print('No phylogeny found!')
                return False

        if len(set(present_at)) > 1:
            print('No phylogeny found!')
            return False
        elif len(present_at) != 0:
            col = present_at[0]
            if col == 0:
                pass
            str_set = set([])
            for r in row_at:
                str_r = ""
                for c in range(col):
                    str_r+=str(k1[r][c].decode("utf-8"))
                str_set.add(str_r)
            print("str_set: ", str_set)
            if len(str_set) > 1:
                print('No phylogeny found!')
                return False
     
    print('Success! Found phylogeny!\nK1 matrix:')
    print(k1)
    return(True)

def get_connected_componet(k,q,m_graph):
    '''
    return the first K vector with E[K] >= 1
    @param k holds the connection vector, initialise with set() 
    @param q chain of connected vertices
    @param m_graph graph object containing matrix connections
    '''
    # find connected components
    if not q:
        return(k)
    else:
        q1 = q.pop()
        k.add(q1) 
        pairs = m_graph.get_pairs_containing(q1)
        pairs = set([node for pair in pairs for node in pair]) #flatten set of elements
        for node in pairs:
            if node not in k:
                q.add(node)
        return(get_connected_componet(k,q,m_graph))

def get_S_prime_semi_universal(m3,s,c3,k):
    s_prime   = set(s).intersection(k)    # s in connected component K
    s_indexes = np.array([np.where(s_i==s)[0][0] for s_i in s_prime])
    m_tmp     = m3.copy()[s_indexes]

    print('k:')
    print(k)
    print("S':")
    print(s_prime)
    
    u = []
    c_in_k  = set(c3).intersection(k)   # c in connected component K
    for c_i in c_in_k:
        c_index = np.where(c_i==c3)[0][0]
        m_col = m_tmp[:,c_index:c_index+1]
        if np.all(m_col!=0):      # c column not all zero ==> s prime semi-universal
            cm = np.where(c_i==c3)[0]
            u.append(c_i)
    
    return s_prime, u
