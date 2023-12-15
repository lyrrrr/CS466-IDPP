import numpy as np

def dot_doc(k1,s,outname):
    '''    
    takes a k1 matrix and writes a dot source file of the node 
    connections and edges, then converts the dot to a postscript file. 
    @param k1 matrix output from perfect phylogeny checking
    @param s samples (rows)
    '''
    dot_lines = []

    for r_idx,k_row in enumerate(k1):
        for c_idx,k_i in enumerate(k_row):
            if k_i=='#': break
            if c_idx == len(k_row) - 1: break

            k_next = k_row[c_idx+1]
            k_i = k_i.decode("utf-8")
            k_next = k_next.decode("utf-8")
            if c_idx==0:
                dot_lines.append('\troot [label=""];\n')
                dot_lines.append('\troot -> node_%s [label="%s"];\n' % (k_i,k_i))
            if k_next=='#':
                dot_lines.append('\tnode_%s [label=""];\n' % k_i)
                dot_lines.append('\tnode_%s -> %s;\n' % (k_i,s[r_idx].decode("utf-8")))
                break                
            dot_lines.append('\tnode_%s [label=""];\n\tnode_%s [label=""];\n' % (k_i,k_next))
            dot_lines.append('\tnode_%s -> node_%s [label="%s"];\n' % (k_i,k_next,k_next))
   
    dot_lines = np.unique(np.array(dot_lines))
    with open('%s.dot'%outname,'w') as fout:
        fout.write('digraph {\n')
        fout.write('\tgraph[size="7.75,10.25"]\n')
        for line in dot_lines:
            fout.write(line)
        fout.write('}\n')