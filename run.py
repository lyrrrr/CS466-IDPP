#!/usr/bin/env python
import argparse
import numpy as np
from IDPP.utils import *
from IDPP import draw_tree

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument('input_file',help='Matrix file')
parser.add_argument('outname',nargs='?',default='tree',help='Output name')
parser.add_argument('--plot',dest='plot',action='store_true',help='Plot output (if valid phylogeny). Requires Graphviz.')
args = parser.parse_args()

def get_value_from_file(in_file):
    """
    get name of species, charaters and matrix values
    """
    m_txt = np.genfromtxt(in_file,delimiter='\t',dtype=None,invalid_raise=False)
    print(m_txt)
    s = np.array(m_txt[1:,0],dtype='S100') 
    c = np.array(m_txt[0,1:],dtype='S100') 
    m = np.array(m_txt[1:,1:],dtype=int)   

    return m, s, c

def IDPP_solution(m,s,c):
    m3, c3 = m_initial(m,c)
    m_graph = MGraph()
    m_graph.build_graph(m3,s,c3)
    m_pairs = m_graph.get_edge_pairs()
    tree_set = [set(s)] #initialise a tree
    u_set = []

    while len(m_pairs) > 1:
        # find connected component k
        k = get_connected_componet(set(),set(m_pairs[0]),m_graph)  
        if len(k) < 3:
            for n in k:
                m_graph.delNode(n)
            m_pairs = m_graph.get_edge_pairs()
        else:
            s_prime, u = get_S_prime_semi_universal(m3,s,c3,k)
            
            if not u: # empty set.. no phylogeny tree
                break
            else:
                u_set.append(u)
                print('u:')
                print(u)
                tree_set.append(s_prime)
                for n in u:
                    m_graph.delNode(n)
                m_pairs = m_graph.get_edge_pairs()
                print('Tree:')
                print(tree_set)

    m_new = m3.copy()
    
    for s_set_i, s_set in enumerate(tree_set[1:]):
        #print("s_set_i: ",s_set_i)
        #print("###### s_set ", s_set)
        s_prime = set(s).intersection(s_set)
        #print("###### s_prime ",s_prime)
        s_indexes = np.array([np.where(s_i==s)[0][0] for s_i in s_prime])
        m_tmp = m3.copy()[s_indexes]

        current_u = u_set[s_set_i]
        #print("current_u: ", current_u)
        for c_i in current_u:
            c_idx = np.where(c_i==c3)[0][0]
            m_col = m_tmp[:,c_idx:c_idx+1]
            if np.all(m_col!=0) and np.any(m_col==-1):
                cm = np.where(c_i==c3)[0]
                m_col[m_col==-1] = 1     # set -1 as 1   in the tree and s_prime semi-supervised
                m_new[s_indexes,c_idx:c_idx+1] = m_col

    for i in range(len(c3)):
        tcol = m_new[:,i:i+1]
        tcol[tcol==-1] = 0     # set rest -1 as 0, other set 0
        m_new[:,i:i+1] = tcol
    
    print("New ori M':")
    print(m_new)
    print(c3)
    m_new, c_prime = m_initial(m_new,c3)
    k1 = get_k1_matrix(m_new,c_prime)
    print("New M':")
    print(m_new)
    
    return(m_new, c_prime, k1, tree_set)

if __name__ == '__main__':
    m,s,c = get_value_from_file(args.input_file)
    print(m,s,c)
    # s = np.array(['s1', 's2', 's3', 's4', 's5']
    #              ,dtype='S100') 
    # c = np.array(['c1', 'c2', 'c3', 'c4', 'c5'],dtype='S100') 
    # m = np.array([[0, -1, -1, -1, 0],
    #               [0, 0, -1, 1, 1],
    #               [1, 1, -1, 0, -1],
    #               [-1, 0, 0, -1, -1],
    #               [-1, -1, -1, 1, 0]]
    #              ,dtype='int')   
    
    m_new, c_prime, k1, tree = IDPP_solution(m,s,c)

    if perfect_phylogeny_exists_new(k1,c_prime) and args.plot:
        print('Plotting tree...')
        draw_tree.dot_doc(k1,s,args.outname)

