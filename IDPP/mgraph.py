import numpy as np

class Node:
    def __init__(self,key):
        self.id = key
        self.connectedTo = {}

    def addNeighbor(self,nbr):
        self.connectedTo[nbr] = 0

    def __str__(self):
        return str(self.id) + ' connectedTo: ' + str([x.id for x in self.connectedTo])

    def getConnections(self):
        return self.connectedTo.keys()

    def delConnection(self,key):
        to_del = [n for n in self.connectedTo if n.getId()==key]
        if len(to_del)>0:            
            del self.connectedTo[to_del[0]]
            
    def getId(self):
        return self.id

class Graph:
    def __init__(self):
        self.nodeList = {}
        self.numNodes = 0

    def addNode(self,key):
        self.numNodes = self.numNodes + 1
        newNode = Node(key)
        self.nodeList[key] = newNode
        return newNode

    def getNode(self,n):
        if n in self.nodeList:
            return self.nodeList[n]
        else:
            return None

    def delNode(self,key):
        try:
            node_to_del = self.getNode(key)
            for n in self.nodeList:
                node = self.getNode(n)
                node.delConnection(key)
            del self.nodeList[key]
            self.numNodes = self.numNodes - 1
            return None
        except KeyError:
            raise Exception('Node %s does not exist' % key)        
        
    def __contains__(self,n):
        return n in self.nodeList

    def addEdge(self,f,t,cost=0):
        if f not in self.nodeList:
            nv = self.addNode(f)
        if t not in self.nodeList:
            nv = self.addNode(t)
        self.nodeList[f].addNeighbor(self.nodeList[t])

    def getNodes(self):
        return self.nodeList.keys()

    def __iter__(self):
        return iter(self.nodeList.values()) 

'''
Author: Marek Cmero
Graph implementation for building incomplete phylogeny _M_ Graphs
'''
class MGraph(Graph):    
    def build_graph(self,m,s,c):
        '''
        take the _M_ matrix with its correspoinding samples (s) and
        feature columns (c), returning the corresponding connection graph    
        '''
        for si in s:
            self.addNode(si)
        for ci in c:
            self.addNode(ci)

        for i in range(len(m)):
            nc = c[np.where(m[i]==1)]
            for ci in nc:
                self.addEdge(s[i],ci)
        return None

    def get_edge_pairs(self):
        '''
        get all the pairwise connections of the graph
        '''
        pairs = []
        for n in self:
            for w in n.getConnections():
                pairs.append([n.getId(),w.getId()])
        return pairs

    def get_pairs_containing(self,x):
        '''
        get all pairs of conections that contain the element ID x
        '''
        pairs = []
        for pair in self.get_edge_pairs():
            if pair[0]==x or pair[1]==x:
                pairs.append(pair)
        return pairs

