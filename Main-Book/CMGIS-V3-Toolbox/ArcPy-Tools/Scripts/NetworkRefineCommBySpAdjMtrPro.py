# **************************************************************************************************************
# Date: December 8, 2020
# Author: Changzhen Wang
# Description: This script implements the ScLouvain and ScLeiden methods
# **************************************************************************************************************


import networkx as nx
import copy
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

state1 = -1 # no flows to adjacent community


# Compute the modularity of a graph partition
def modularity(partition, graph):
    if type(graph) != nx.Graph:
        raise TypeError("Bad graph type, use only non directed graph")

    inc = dict([]) # sum of weights of all edges within each community c
    deg = dict([]) # sum of weights of all edges connected with a node in each community c
    links = graph.size(weight='weight') # return total of all edge weights.
    if links == 0:
        raise ValueError("A graph without link has an undefined modularity")

    for node in graph:
        if bool(partition.get(node)):
            com = partition[node]
            deg[com] = deg.get(com, 0.) + graph.degree(node, weight='weight')
            for neighbor, datas in graph[node].items():
                weight = datas.get("weight", 1)
                if bool(partition.get(neighbor)):
                    if partition[neighbor] == com:
                        if neighbor == node:  
                            inc[com] = inc.get(com, 0.) + float(weight) #inc[com] is the sum of weights of all edges within the community com.
                        else:
                            inc[com] = inc.get(com, 0.) + float(weight) / 2.

    res = 0.
    for com in set(partition.values()):
        res += (inc.get(com, 0.) / links) - ((deg.get(com, 0.) / (2. * links)) ** 2)  # Q=∑c[Σin2m−(Σtot2m)2]; annotation should be ∑(in/links-r*(ds/2*links)**2) by Wang 12-15-2018
    return res


# split each community by spatial adjacency matrix and update the status
def __refine_graph(current_graph, com2node_dict, status, adjmx, graph):
    com2node_dict_temp = copy.deepcopy(com2node_dict)
    
    for node in current_graph.nodes():
        #com_node = status.node2com[node]
        node_original_list = com2node_dict_temp[node]
        sub_coms = dict()
        if len(node_original_list) > 1:
            sub_coms = split_com_by_adj(node_original_list, adjmx) 
        if len(sub_coms) > 1: # exist multiple communities, then remove and extend the current com2node_dict_temp
            for i in range(1,len(sub_coms)):
                temp_pop = 0
                for sub_com_node in sub_coms[i]:
                    com2node_dict_temp[node].remove(sub_com_node)  
                    temp_pop +=  graph.nodes[sub_com_node]['pop']                               
                ext_com = {max(com2node_dict_temp.keys()) + 1 : sub_coms[i]}
                com2node_dict_temp.update(ext_com)
                status.com2node.update({max(status.com2node.keys()) + 1 : len(status.com2node)})
                status.node2com.update({max(status.node2com.keys()) + 1 : node})
                status.com_attr1.update({max(status.com_attr1.keys()) + 1: temp_pop})
                status.com_attr1[node] -= temp_pop
    
    partition_init = __com2node_partition(com2node_dict_temp)
    current_graph = induced_graph(partition_init, graph, None, status)
    return current_graph, com2node_dict_temp


# split a list of nodes by spatial adjacent matrix
def split_com_by_adj(nodeslist, adjmx):
    sort_nodeslist = sorted(nodeslist)
    n = len(sort_nodeslist)
    nodes_dict = {k : sort_nodeslist[k] for k in range(0, len(sort_nodeslist))}
    adj_new = np.array([[0]*n for x in range(n)])  
    rows_bl = np.isin(adjmx[0], sort_nodeslist)
    col_bl = np.isin(adjmx[1], sort_nodeslist)

    rows_index = np.nonzero(rows_bl)
    cols_index = np.nonzero(col_bl)
    comm_index = np.intersect1d(rows_index, cols_index)

    dd1 = np.take(adjmx[0], comm_index)
    dd2 = np.take(adjmx[1], comm_index)

    for lvalue, rvalue in zip(dd1, dd2):
        lindex = sort_nodeslist.index(lvalue)
        rindex = sort_nodeslist.index(rvalue)
        adj_new[lindex][rindex] = adj_new[rindex][lindex] = 1
    mlabel = adj2cluster(adj_new)
    sub_partition = { k: [] for k in mlabel.keys()}
    for key, value in mlabel.items():
        for lower_value in value:
            sub_partition[key].append(nodes_dict[lower_value])
    return sub_partition


# It is a very important and efficient function to construct cluster/island
def adj2cluster(A):
    # symmetrize adjacency matrix
    S = A + A.transpose()
    #print(S)
    # transform dense matrix into sparse matrix
    G_sparse = csr_matrix(S)
    #print(G_sparse)
    # Reverse Cuthill-McKee ordering
    r = np.flip(reverse_cuthill_mckee(G_sparse))
    #print(r)
    # Get the clusters
    mlabel = dict()
    clusterNum = 0
    mlabel[clusterNum] = [r[0]]
    for i in range(1, len(r)):
        if S[mlabel[clusterNum], r[i]].any():
            c = np.append(mlabel[clusterNum],r[i])
            mlabel[clusterNum] = c
        else:
            clusterNum = clusterNum + 1
            mlabel[clusterNum] = [r[i]]
    return mlabel


# calculate LI, average travel time and the number of destination nodes
# partition is the final partition
# part_init is the initial partition
def calculateindices(digraph, partition, part_init):
    internalvols = dict([])
    outgoingvols = dict([])
    sum_weight = dict()
    sum_weighttime = dict()
    dest_nodes = dict()

    links = digraph.size(weight = "weight")
    if links == 0:
        raise ValueError("A directed graph without link is invalid!")
    
    alledgeattr = list(list(digraph.edges(data=True))[0][-1].keys())
    for node in digraph.nodes():        
        com = partition[node]
        com_init = part_init[node]
        if com_init != state1:
            if bool(dest_nodes.get(com)):
                dest_nodes[com] += 1
            else:
                dest_nodes[com] = 1
        for neighbor, datas in digraph[node].items():        
            odweight = 0
            odesttime = 0
            if datas is not None:
                odweight = datas["weight"]                  
                if "estTime" in alledgeattr:
                    odesttime = datas["estTime"]
                sum_weight[com] = sum_weight.get(com, 0) + odweight
                sum_weighttime[com] = sum_weighttime.get(com, 0) + odweight * odesttime
            if partition[neighbor] == com:  # node and its neighbor are in the sam community                
                internalvols[com] = internalvols.get(com, 0) + float(odweight)
                outgoingvols[com] = outgoingvols.get(com, 0)
            else:
                outgoingvols[com] = outgoingvols.get(com, 0) + float(odweight)
                internalvols[com] = internalvols.get(com, 0)
 
    values = sorted(set(partition.values()))
    res = {k: [] for k in list(values)}

    for com in values:
        if not bool(internalvols.get(com)):
            internalvols[com] = 0
        if not bool(outgoingvols.get(com)):
            outgoingvols[com] = 0
        if not bool(sum_weight.get(com)):
            sum_weight[com] = 0
        internalvol = internalvols[com]
        LI = -1
        if internalvol != links:
            if internalvol + outgoingvols[com] > 0:
                LI = internalvol/(internalvol + outgoingvols[com])
        res[com].append(LI)
        weightedtime = 0        
        if sum_weight[com] != 0:
            weightedtime = sum_weighttime[com] / sum_weight[com]
        res[com].append(weightedtime)
        if bool(dest_nodes.get(com)):
            res[com].append(dest_nodes[com])
        else:
            res[com].append(0)
    return res    

# make a new graph where discovered communities become its nodes,
# links between communities become its edges (weights are summed).
def induced_graph(partition, graph, com_old_vs_com_new, status):
    # com_old_vs_com_new is a dict that stores the relationship between old community id vs. new id after renumbering
    # status is the current status
    ret = nx.Graph()
    ret.add_nodes_from(partition.values(), pop =0) # partition.values() has the community number of each node after renumbering
    #ret.add_nodes_from(partition.values())
    # add and update the attribute1 (e.g., pop) of the new communities after renumbered
    if status is not None:
        if com_old_vs_com_new is not None:
            for com_old in com_old_vs_com_new.keys():
                com_new = com_old_vs_com_new[com_old]
                com_attr1 = status.com_attr1[com_old]
                ret.nodes[com_new]['pop'] = com_attr1                
        else:
            for com, popvalue in status.com_attr1.items():
                if ret.has_node(com):
                    ret.nodes[com]['pop'] = popvalue
    else:
        sum_pop1 = 0
        for node in graph.nodes():
            com = partition[node]
            if ret.has_node(com):
                ret.nodes[com]['pop'] += graph.nodes[node]['pop'] 
                sum_pop1 += graph.nodes[node]['pop']     
               
    for node1, node2, datas in graph.edges(data=True): #  If data=True, then return edge attribute dict in 3-tuple (u,v,ddict)
        weight = datas.get("weight", 1) # an edge weight should have a column name 'weight'!!!!! if 'weight' column is not found, then return 1 for this variable       
        com1 = partition[node1] # community of node 1 of an edge
        com2 = partition[node2] # community of node 2 of an edge  
        w_prec = ret.get_edge_data(com1, com2, {"weight": 0}).get("weight", 1) # Return the attribute dictionary associated with edge (u,v); {"weight": 0} is the value to return if the edge (com1,com2) is not found in the network ret.     
        ret.add_edge(com1, com2, weight = w_prec + weight)

    return ret


# Renumber the values of the dictionary from 1 to n
# It is primarily used to renumber the community id in the partition result list
# Also, it updates the community id in the global dict com2node_globaldict after its renumbering
def __renumber(partition_dict, com2node_dict):
    count = 0 # new community id after renumber starts from 0
    re_part = partition_dict.copy()
    new_values = dict([]) # a dict that stores the relationship between old and new community id, see below
    com2node_tempdict = dict([]) # a temporal dict of com2node

    for key in partition_dict.keys():
        value = partition_dict[key]
        new_value = new_values.get(value, -1)
        if new_value == -1:
            new_values[value] = count # new_values is a dict that stores {526: 1, 486: 2...} where 526 is the
            # original community name (i.e., value) and 1 is the community name after renumbering (i.e., count)
            com2node_tempdict[count] = com2node_dict[value] # get ids that a community contains
            new_value = count   # pass value to the outside of the loop
            count += 1
        re_part[key] = new_value

    return re_part, com2node_tempdict, new_values


# has to change as the gain (even negative) is useless, merge the population to the o
def __second_level_pop(graph, status, attr1_limit, resolution, com2node_dict, isloop, adjmx):
    """
        merge all communities (i.e., new nodes in the current graph) whose attribute
        (e.g., total pop) are less than the predefined threshold to their neighbors
    """
    com2node_dict_temp = copy.deepcopy(com2node_dict)
    for node in graph.nodes():
        com_node = status.node2com[node]
        com_node_attr1 = status.com_attr1[com_node]
        if com_node_attr1 < attr1_limit:
            degc_totw = status.gdegrees.get(node, 0.) / (status.total_weight * 2.) 
            edgeneigh_communities, toponeigh_communities = __toponeighcom(node, graph, status, com2node_dict, adjmx)            
            __remove(node, com_node, edgeneigh_communities.get(com_node, 0.), status) # remove the current node from its original community          
            best_com = com_node            
            best_increase = 0
            if isloop is False:
                if len(toponeigh_communities) == 1:
                    best_com = list(toponeigh_communities.keys())[0]
                else:
                    if all(value == 0 for value in toponeigh_communities.values()):
                        com_pop_dict = {key: status.com_attr1[key] for key in toponeigh_communities.keys()}
                        com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1]))
                        neigh_coms = {key : value for (key, value) in com_pop_dict_sort.items() if value >= attr1_limit}
                        if len(neigh_coms) == 0:
                            best_com = list(com_pop_dict_sort.keys())[-1]   
                        else:
                            best_com = [key for (key, value) in neigh_coms.items() if value == min(neigh_coms.values())][0]                                                             
                    else:
                        incr = 0                       
                        for neigh_com, dnc in toponeigh_communities.items():
                            incr = dnc - resolution * status.degrees.get(neigh_com, 0.) * degc_totw
                            if incr > best_increase:
                                best_increase = incr
                                best_com = neigh_com
            else:
                com_pop_dict = {key: status.com_attr1[key] for key in toponeigh_communities.keys()}
                com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1]))
                neigh_coms = {key : value for (key, value) in com_pop_dict_sort.items() if value >= attr1_limit}
                if len(neigh_coms) == 0:
                    best_com = list(com_pop_dict_sort.keys())[-1]   
                else:
                    best_com = [key for (key, value) in neigh_coms.items() if value == min(neigh_coms.values())][0]  
            __insert(node, best_com, edgeneigh_communities.get(best_com, 0.), status) # assign this node to the community with the max deltaQ
            if best_com != com_node:
                node_origin_list = com2node_dict[node][:]  # only pass value not reference (another way is to use deepcopy)
                for nd in node_origin_list:
                    com2node_dict_temp[com_node].remove(nd)  # remove
                    com2node_dict_temp[best_com].append(nd)  # update the com2node_dict               
    com2node_dict = copy.deepcopy(com2node_dict_temp)
    return com2node_dict


def __ensure_one_more_destnodes(graph, status, resolution, com2node_dict, com_no_dnodes, isloop, adjmx):
    com2node_dict_temp = copy.deepcopy(com2node_dict)

    for node in graph.nodes():
        com_node = status.node2com[node]
        if com_node in com_no_dnodes:
            degc_totw = status.gdegrees.get(node, 0.) / (status.total_weight * 2.)
            edgeneigh_communities, toponeigh_communities = __toponeighcom(node, graph, status, com2node_dict, adjmx)
            __remove(node, com_node, edgeneigh_communities.get(com_node, 0.), status) # remove the current node from its original community
            best_com = com_node
            best_increase = 0
            if isloop is False:
                if len(toponeigh_communities) == 1:
                    best_com = list(toponeigh_communities.keys())[0]
                else:
                    if all(value == 0 for value in toponeigh_communities.values()): # no flows to topocommunity
                        com_has_dnodes = [key for key in toponeigh_communities.keys() if key not in com_no_dnodes] # find topocommunity with dnodes
                        if len(com_has_dnodes) > 0 :
                            com_pop_dict = {key: status.com_attr1[key] for key in com_has_dnodes} # find the population
                            com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1])) # ascending order
                            best_com = list(com_pop_dict_sort.keys())[0] # get com with the smallest population size                        
                    else:
                        incr = 0
                        for neigh_com, dnc in toponeigh_communities.items():
                            if neigh_com not in com_no_dnodes:
                                incr = dnc - resolution * status.degrees.get(neigh_com, 0.) * degc_totw
                                if incr > best_increase: 
                                    best_increase = incr
                                    best_com = neigh_com                                                          
            else:
                # just in case the dead loop and no optimial topocommunity is found from above scripts
                com_has_dnodes = [key for key in toponeigh_communities.keys() if key not in com_no_dnodes] # find topocommunity with dnodes
                if len(com_has_dnodes) > 0 :
                    com_pop_dict = {key: status.com_attr1[key] for key in com_has_dnodes} # find the population
                    com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1])) # ascending order
                    best_com = list(com_pop_dict_sort.keys())[0] # get com with the smallest population size
            __insert(node, best_com, edgeneigh_communities.get(best_com, 0.), status) # assign this node to the community with the max deltaQ
            if best_com != com_node:
                node_original_list = com2node_dict[node][:]  # only pass value not reference (another way is to use deepcopy)
                for nd in node_original_list:
                    com2node_dict_temp[com_node].remove(nd)  # remove
                    com2node_dict_temp[best_com].append(nd)  # update the com2node_dict
    
    com2node_dict = copy.deepcopy(com2node_dict_temp)
    return com2node_dict

        

# ensure all community has at least one destination node
def __ensure_nonzero_LI(graph, status, resolution, com2node_dict, isloop, adjmx):
    com2node_dict_temp = copy.deepcopy(com2node_dict)
    for node in graph.nodes():
        com_node = status.node2com[node]
        incflows = status.internals[com_node] # get internal flows
        if incflows == 0:
            degc_totw = status.gdegrees.get(node, 0.) / (status.total_weight * 2.)
            edgeneigh_communities, toponeigh_communities = __toponeighcom(node, graph, status, com2node_dict, adjmx)
            __remove(node, com_node, edgeneigh_communities.get(com_node, 0.), status) # remove the current node from its original community
            best_com = com_node
            best_increase = 0                    
            if isloop is False:
                if len(toponeigh_communities) == 1:
                    best_com = list(toponeigh_communities.keys())[0]
                else:
                    if all(value == 0 for value in toponeigh_communities.values()):
                        com_pop_dict = {key: status.com_attr1[key] for key in toponeigh_communities.keys() if status.internals[key] != 0 }
                        com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1]))
                        best_com = list(com_pop_dict_sort.keys())[0]
                    else:
                        incr = 0
                        for neigh_com, dnc in toponeigh_communities.items():
                            incr = dnc - resolution * status.degrees.get(neigh_com, 0.) * degc_totw
                            if incr > best_increase:
                                best_increase = incr
                                best_com = neigh_com
            else:
                com_pop_dict = {key: status.com_attr1[key] for key in toponeigh_communities.keys() if status.internals[key] != 0}
                com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1]))
                if len(com_pop_dict_sort) > 0:
                    best_com = list(com_pop_dict_sort.keys())[0]

            __insert(node, best_com, edgeneigh_communities.get(best_com, 0.), status) # assign this node to the community with the max deltaQ
            if best_com != com_node:
                node_original_list = com2node_dict[node][:]  # only pass value not reference (another way is to use deepcopy)
                for nd in node_original_list:
                    com2node_dict_temp[com_node].remove(nd)  # remove
                    com2node_dict_temp[best_com].append(nd)  # update the com2node_dict

    com2node_dict = copy.deepcopy(com2node_dict_temp)    
    return com2node_dict


# convert com2node to partition
def __com2node_partition(com2node):
    partition = {}
    for key, values in com2node.items():
        if len(values) > 0:
            for value in values:
                partition[value] = key
    return partition

# convert partition to com2node
def partition_com2node(partition):
    com2node_dict = {k: [] for k in list(set(partition.values()))}
    for node, com in partition.items():
        if bool(com2node_dict.get(com)):
            com2node_dict[com].append(node)
        else:
            com2node_dict[com] = [node]
    return com2node_dict

# part_init is the partition from network optimization method without any spatial constraints
# part_init2 is the partition based on the destination node
def refined_partition_network(graph, adjmx, attr1_limit, resolution, part_init, part_init2):
    if type(graph) != nx.Graph:
        raise TypeError("Bad graph type, use only non directed graph")

    status = Status()
    status.init(graph, part_init)
    com2node_globaldict = copy.deepcopy(status.com2node)
    current_graph = induced_graph(part_init, graph, None, status)
    # will decompose community by spatial adjacency matrix first
    current_graph, com2node_globaldict = __refine_graph(current_graph, com2node_globaldict, status, adjmx, graph) 
    status.init(current_graph)

    # the following codes examine the attribute (e.g., total population) in each community, and merge the
    # community that is below the predefined threshold to its neighbor having the most modularity gain.
    minpop = min(status.com_attr1.values()) # the minimum pop of all communities of current status
    if attr1_limit is not None and attr1_limit > 0:
        if minpop < attr1_limit:
            isloop = False
            while True:
                num_coms = len(com2node_globaldict)
                com2node_globaldict = __second_level_pop(current_graph, status, attr1_limit, resolution, com2node_globaldict, isloop, adjmx)                  
                partition, com2node_globaldict, com_old_vs_com_new = __renumber(status.node2com, com2node_globaldict)
                current_graph = induced_graph(partition, current_graph, com_old_vs_com_new, status)
                status.init(current_graph)
                minpop = min(status.com_attr1.values())
                if minpop >= attr1_limit:
                    break
                else:
                    if num_coms == len(com2node_globaldict):
                        isloop = True # means infinite loop will begin, so merge strategy should be changed


    # find all destination nodes
    allvalues = [value for value in part_init2.values() if value != state1]
    # find com without destination nodes
    com_no_destnodes = [com for com, nodes in com2node_globaldict.items() if len(set(nodes).intersection(set(allvalues))) == 0]
    if len(com_no_destnodes) > 0:
        isloop = False
        while True:
            num_coms = len(com2node_globaldict)
            com2node_globaldict =__ensure_one_more_destnodes(current_graph, status, resolution, com2node_globaldict, com_no_destnodes, isloop, adjmx)
            partition, com2node_globaldict, com_old_vs_com_new = __renumber(status.node2com, com2node_globaldict)
            current_graph = induced_graph(partition, current_graph, com_old_vs_com_new, status)
            status.init(current_graph)
            com_no_destnodes = [com for com, nodes in com2node_globaldict.items() if len(set(nodes).intersection(set(allvalues))) == 0]
            if len(com_no_destnodes) == 0:
                break
            else:
                if num_coms == len(com2node_globaldict):
                    isloop = True

    # processing orphan nodes (a) no flows to all, (b) no internal flows
    miniinnflows = min(status.internals.values()) # ensure all LI >0
    if miniinnflows == 0:
        isloop = False
        while True:
            num_coms = len(com2node_globaldict)               
            com2node_globaldict = __ensure_nonzero_LI(current_graph, status, resolution, com2node_globaldict, isloop, adjmx)
            partition, com2node_globaldict, com_old_vs_com_new = __renumber(status.node2com, com2node_globaldict)
            current_graph = induced_graph(partition, current_graph, com_old_vs_com_new, status)
            status.init(current_graph)
            miniinnflows = min(status.internals.values())
            if miniinnflows > 0:
                break
            else:
                if num_coms == len(com2node_globaldict):
                    isloop = True 

    partition = __com2node_partition(com2node_globaldict)    
    return partition


class Status(object):
    """
    To handle several data in one struct.

    Could be replaced by named tuple, but don't want to depend on python 2.6
    """
    node2com = {} # the community of each node
    com2node = {} # ids of nodes in each community, added by hu
    com_attr1 = {} # one attribute of each community, e.g., total population, added by hu
    node_attr1 = {} # one attribute of each node, e.g., total population, added by hu
    total_weight = 0 # sum of weights of all edges
    internals = {} # sum of weights of all edges within a community
    degrees = {} # degree list for all communities
    gdegrees = {} # degree list for all nodes in graph

    def __init__(self):
        self.node2com = dict([])
        self.com2node = dict([])
        self.com_attr1 = dict() # its value of each key is a number, not a list

        self.node_attr1 = dict()
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([]) # edge weights list for all communities
        self.loops = dict([]) # edge weights list for all nodes in graph

    def copy(self):
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = copy.deepcopy(self.node2com)
        new_status.com2node = copy.deepcopy(self.com2node)
        new_status.com_attr1 = copy.deepcopy(self.com_attr1)
        new_status.node_attr1 = copy.deepcopy(self.node_attr1)
        new_status.internals = copy.deepcopy(self.internals)
        new_status.degrees = copy.deepcopy(self.degrees)
        new_status.gdegrees = copy.deepcopy(self.gdegrees)
        new_status.loops = copy.deepcopy(self.loops)
        new_status.total_weight = self.total_weight

        return new_status

    def init(self, graph, part=None):
        """Initialize the status of a graph with every node in one community"""
        count = 0 # new community id starts from 0
        self.node2com = dict([])
        self.com2node = dict([])
        self.com_attr1 = dict()
        self.node_attr1 = dict()
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.total_weight = graph.size(weight='weight') # .size function returns the number of edges or sum of edge weights in the graph

        if part is None:
            self.com2node = {k: [] for k in range(0, graph.number_of_nodes())} # initialize a dictionary of empty list
            for node in graph.nodes(): # node is the node name
                self.node2com[node] = count # if no partition provided, each node belongs to one community             
                self.node_attr1[node] = graph.nodes[node]['pop']
                if bool(self.com2node.get(count)): # if the same community exists, append nodes and plus attributes, changed by Wang 12-15-2018
                    self.com2node[count].append(node) # node is the node name. Other tips: g.node.keys()[i] return
                    # the ith key in the dict
                    self.com_attr1[count] += graph.nodes[node]['pop'] # may need to change the field name for other attributes
                else: # if no same community exists, append nodes and set attributes, changed by Wang 12-15-2018 
                    self.com2node[count].append(node)
                    self.com_attr1[count] = graph.nodes[node]['pop']
                deg = float(graph.degree(node, weight='weight')) # calculate the degree of each node
                if deg < 0:
                    error = "Bad graph type ({})".format(type(graph))
                    raise ValueError(error)
                self.degrees[count] = deg
                self.gdegrees[node] = deg
                edge_data = graph.get_edge_data(node, node, {"weight": 0})
                self.loops[node] = float(edge_data.get("weight", 1))
                self.internals[count] = self.loops[node]                
                count += 1
        else: # when a partition is provided as input
            #self.com2node = {k: [] for k in range(min(list(set(part.values()))), max(list(set(part.values()))) + 1 )}
            self.com2node = {k: [] for k in list(set(part.values()))}
            for node in graph.nodes():
                com = part[node]
                self.node2com[node] = com
                self.node_attr1[node] = graph.nodes[node]['pop']
                if bool(self.com2node.get(com)):
                    self.com2node[com].append(node)
                    self.com_attr1[com] += graph.nodes[node]['pop']
                else:
                    self.com2node[com].append(node)
                    self.com_attr1[com] = graph.nodes[node]['pop']
                deg = float(graph.degree(node, weight='weight'))
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.gdegrees[node] = deg

                inc = 0.
                for neighbor, datas in graph[node].items():
                    edge_weight = datas.get("weight", 1)              
                    if edge_weight <= 0:
                        raise ValueError("Bad graph type, use positive weights")
                    if part[neighbor] == com: # if one of this node's neighbors is also in the same community
                        if neighbor == node: # if this neighbor is itself (i.e., a loop)
                            inc += float(edge_weight) # inc indicates the sum of edge weights inside a community
                        else:
                            inc += float(edge_weight) / 2

                self.internals[com] = self.internals.get(com, 0) + inc


def __isadj(nodelist, neighborlist, adjmx):
    rows_bl = np.isin(adjmx[0], nodelist)
    col_bl = np.isin(adjmx[1], neighborlist)

    rows_index = np.nonzero(rows_bl)
    cols_index = np.nonzero(col_bl)
    comm_index = np.intersect1d(rows_index, cols_index)
    if(len(comm_index) > 0):
        return True
    else:
        return False

'''
    Compute the communities of the neighborood of the node in the graph given
    with the decomposition node2com; note that the neighborhood of the node is defined
    more accurately by considering the spatial adjacency matrix 'adjmx' calculated
    from the polygon features. Note that the adjmx starts index from 0 rather than 1.
    The input com2node_dict (i.e., com2node_globaldict) is used to
    find the original node names. com_node is the community id of this node.

    In other words, a neighbor of one node i in the graph is a node j that has an
    edge ij between them as well as polygon i and polygon j are spatially touched.
'''
def __toponeighcom(node, graph, status, com2node_dict, adjmx):
    weights_EDGEneighbor = {} # it stores results for neighbors defined based on edges only
    weights_TOPOneighbor = {} # it stores results for neighbors defined based on topology only
    # sometimes, a node is only connected with nontopological neighbors
    nodelist = com2node_dict[node]
    neighcom_list = list()
    for neighbor, datas in graph[node].items(): # graph[node].items() returns the node name
        # and edge weight of each edge that has one end being this node
        if neighbor != node: # if one of this node's neighbors is not itself
            weight = datas.get("weight", 1)
            neighborcom = status.node2com[neighbor] # the community of its neighbor            
            neighborlist = com2node_dict[neighborcom]
            weights_EDGEneighbor[neighborcom] = weights_EDGEneighbor.get(neighborcom, 0) + weight                               
            isadj = __isadj(nodelist, neighborlist, adjmx) 
            if isadj:
                #isConnected_toponeigh = True                    
                weights_TOPOneighbor[neighborcom] = weights_TOPOneighbor.get(neighborcom, 0) + weight
                neighcom_list.append(neighborcom)
    # get spatial adjacent neighbors 
    sortnodelist = sorted(nodelist)
    rows_bl = np.isin(adjmx[0], sortnodelist)
    rows_index = np.nonzero(rows_bl)
    node_toponbr_list = np.take(adjmx[1], rows_index)[0]
    diff_list = set(node_toponbr_list) - set(sortnodelist)
    toponbr_lists = set(diff_list)

    '''                
    for node_i in nodelist:
        node_toponbr_list = np.nonzero(adjmx[node_i-1])[0] 
        node_toponbr_list = [value + 1 for value in node_toponbr_list]
        diff_list = set(node_toponbr_list) - set(nodelist)
        toponbr_lists = list(set(toponbr_lists) | set(diff_list))
    '''
    toponeighcom_list = set([key for key, value in com2node_dict.items() if len(list(set(value).intersection(set(toponbr_lists))))]) - set(neighcom_list) - set([node])
    for neighborcom in toponeighcom_list:
        edge_data = graph.get_edge_data(node, neighborcom)
        weight = 0
        if edge_data is not None:
            weight = edge_data.get("weight")
        weights_TOPOneighbor[neighborcom] = weights_TOPOneighbor.get(neighborcom, 0) + weight 
        #weights_EDGEneighbor[neighborcom] = 0             
    return weights_EDGEneighbor, weights_TOPOneighbor # it is a dictionary which stores the community of its neighbor and the weight of edge between them


def __remove(node, com, weight, status):
    """ Remove node from community com and modify status"""
    status.degrees[com] = (status.degrees.get(com, 0.)
                           - status.gdegrees.get(node, 0.))
    status.internals[com] = float(status.internals.get(com, 0.) -
                                  weight - status.loops.get(node, 0.))
    status.node2com[node] = -1 # then this node has no community assigned
    status.com2node[com].remove(node) # remove this node ('node' is the node name here not the index)
    # from the original community
    status.com_attr1[com] -= status.node_attr1[node] # update the attr1 of this community


def __insert(node, com, weight, status):
    """ Insert node into community and modify status"""
    status.node2com[node] = com
    if bool(status.com2node.get(com)): # if there is at least a value of key=com in the dictionary
        status.com2node[com].append(node)
        status.com_attr1[com] += status.node_attr1[node]
    else:
        # status.com2node[com] = [node] # make it a list
        status.com2node[com].append(node)
        status.com_attr1[com] = status.node_attr1[node]
    status.degrees[com] = (status.degrees.get(com, 0.) +
                           status.gdegrees.get(node, 0.))
    status.internals[com] = float(status.internals.get(com, 0.) +
                                  weight + status.loops.get(node, 0.))

# Fast compute the modularity of the partition of the graph using status precomputed
def __modularity(status, resolution):
    links = float(status.total_weight)
    result = 0.
    for community in set(status.node2com.values()):
        in_degree = status.internals.get(community, 0.) # Σin, i.e., sum of weights of edges within a community
        degree = status.degrees.get(community, 0.) # Σtot, i.e., sum of weights of edges that have one end node inside this community
        if links > 0:
            result += in_degree / links - resolution * ((degree / (2. * links)) ** 2)
    return result
