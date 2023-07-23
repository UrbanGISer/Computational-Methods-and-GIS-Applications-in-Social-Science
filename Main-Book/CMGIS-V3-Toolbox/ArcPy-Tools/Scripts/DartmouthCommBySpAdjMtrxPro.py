# **************************************************************************************************************
# Date: December 13, 2020
# Author: Changzhen Wang
# Description: This script implements the plurality rule and provides functions for the refined Dartmouth method
# **************************************************************************************************************

import networkx as nx
import copy
import numpy as np


# digraph is the original directed graph
# partition is the final partition
# part_init is the initial partition
def calculateindices(digraph, partition, part_init):
    internalvols = dict()
    outgoingvols = dict()
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
        if com_init != -node:
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
            if partition[neighbor] == com:  # node and its neighbor are in the same community                
                internalvols[com] = internalvols.get(com, 0) + float(odweight)
                outgoingvols[com] = outgoingvols.get(com, 0)
            else:
                outgoingvols[com] = outgoingvols.get(com, 0) + float(odweight)
                internalvols[com] = internalvols.get(com, 0)            
    
    values = sorted(set(partition.values()))
    res = {k: [] for k in list(values)}
    
    # accommodate indices for orphan nodes
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
def induced_graph(partition, digraph, com_old_vs_com_new, status):
    # com_old_vs_com_new is a dict that stores the relationship between old community id vs. new id after renumbering
    # status is the current status
    diret = nx.DiGraph()
    diret.add_nodes_from(partition.values(), pop = 0)
    # add and update the attribute1 (e.g., pop) of the new communities after renumbered
    if status is not None:
        if com_old_vs_com_new is not None:
            for com_old in com_old_vs_com_new.keys():
                com_new = com_old_vs_com_new[com_old]
                com_attr1 = status.com_attr1[com_old]
                diret.nodes[com_new]['pop'] = com_attr1
        else:
            for com, popvalue in status.com_attr1.items():
                if diret.has_node(com):
                    diret.nodes[com]['pop'] = popvalue
    else:        
        sum_pop2 = 0
        for node in digraph.nodes():
            com = partition[node]
            if diret.has_node(com):
                diret.nodes[com]['pop'] += digraph.nodes[node]['pop']
                sum_pop2 += digraph.nodes[node]['pop']

    for node1, node2, datas in digraph.edges(data = True):
        weight = datas.get("weight", 1)
        com1 = partition[node1] # community of node 1 of an edge
        com2 = partition[node2] # community of node 2 of an edge
        w_prec = diret.get_edge_data(com1, com2, {"weight": 0}).get("weight", 1)
        diret.add_edge(com1, com2, weight = w_prec + weight)
            
    return diret


# Renumber the values of the dictionary from 1 to n
# It is primarily used to renumber the community id in the partition result list
# Also, it updates the community id in the global dict com2node_globaldict after its renumbering
def __renumber(partition_dict, com2node_dict):
    count = 1 # new community id after renumber starts from 1
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

# convert com2node to partition
def __com2node_partition(com2node):
    partition = {}
    for key, values in com2node.items():
        if len(values) > 0:
            for value in values:
                partition[value] = key
    return partition

# four options:
# (1), no constraints
# (2), only population
# (3), only LI
# (4), population and LI combined
def find_partition_Dartmouth(digraph, adjmx, part_init, attr1_limit, localizationindex):
    if type(digraph) != nx.DiGraph:
        raise TypeError("Bad graph type, use one directed graph")
    dicurrent_graph = digraph.copy()
    status = Status()
    status.init(dicurrent_graph, part_init)   
    com2node_globaldict = __one_level_Dartmouth(dicurrent_graph, status, adjmx)
    # it is required to renumber the HSA for defining HRRs (build spatial adjacency matrix)
    partition = __com2node_partition(com2node_globaldict)    
    partition, com2node_globaldict, com_old_vs_com_new = __renumber(partition, com2node_globaldict)
    dicurrent_graph = induced_graph(partition, dicurrent_graph, com_old_vs_com_new, None)
    status.init(dicurrent_graph)

    # when define population threshold
    if attr1_limit is not None and attr1_limit > 0:      
        minpop = min(status.com_attr1.values())
        if minpop < attr1_limit:
            isloop = False
            while True:
                num_coms = len(com2node_globaldict)
                upadjmx = __update_adjmx(adjmx, com2node_globaldict)
                com2node_globaldict = __second_level_pop_Dartmouth(dicurrent_graph, status, com2node_globaldict, upadjmx, attr1_limit, isloop)
                partition, com2node_globaldict, com_old_vs_com_new = __renumber(status.node2com, com2node_globaldict)
                dicurrent_graph = induced_graph(partition, dicurrent_graph, com_old_vs_com_new, status)
                status.init(dicurrent_graph)
                minpop = min(status.com_attr1.values())
                if minpop >= attr1_limit:
                    break
                else:
                    if num_coms == len(com2node_globaldict):
                        isloop = True        

    # when define Li threshold               
    if localizationindex is not None and localizationindex > 0:               
        minli = min(status.li.values())
        if minli < localizationindex:            
            isloop = False
            while True:
                num_coms = len(com2node_globaldict)
                upadjmx = __update_adjmx(adjmx, com2node_globaldict)
                com2node_globaldict = __second_level_li_Dartmouth(dicurrent_graph, status, com2node_globaldict, upadjmx, localizationindex, isloop)
                partition, com2node_globaldict, com_old_vs_com_new = __renumber(status.node2com, com2node_globaldict)
                dicurrent_graph = induced_graph(partition, dicurrent_graph, com_old_vs_com_new, status)
                status.init(dicurrent_graph)
                minli = min(status.li.values())
                if minli >= localizationindex:
                    break
                else:
                    if num_coms == len(com2node_globaldict):
                        isloop = True

    partition = __com2node_partition(com2node_globaldict) # a required step to get original node and community  
    print(min(partition.keys()), max(partition.keys()), min(partition.values()), max(partition.values()))
    return partition


# update spatial adjacency matrix to current level
def __update_adjmx(adjmx, com2node_dict):
    n = len(com2node_dict.keys())
    adj = np.zeros((n, n), dtype=int) # define an zero matrix
    for i in range(n):
        origin_node_list = com2node_dict[i + 1]
        origin_node_list = [(key-1) for key in origin_node_list]
        rows_bl = np.isin(adjmx[0], origin_node_list)
        rows_index = np.nonzero(rows_bl)
        node_toponbr_list = np.take(adjmx[1], rows_index)[0]
        toponbr_lists = list(set(node_toponbr_list) - set(origin_node_list))
        toponbr_lists = [(key + 1) for key in toponbr_lists]
        toponbr_com_lists = [com for (com, nodes) in com2node_dict.items() if bool(set(nodes) & set(toponbr_lists)) == True]
        for toponbr_com in toponbr_com_lists:            
            if i != toponbr_com:
                adj[i][toponbr_com -1] = adj[toponbr_com -1][i] = 1
    return adj

# calculate the three-step refined Dartmouth method
def __one_level_Dartmouth(digraph, status, adjmx):    
    nodes_move0 = copy.deepcopy(list(status.node2com.keys()))
    nodes_move = [nd for nd in nodes_move0 if digraph.degree(nd) > 0]
    new_unassignedNodes = len(nodes_move)
    unassignedNodes = new_unassignedNodes
    while True:        
        for node in digraph.nodes():
            if node in nodes_move:
                com_node = status.node2com[node]
                out_edgeneigh_communities, in_edgeneigh_communities, out_toponeigh_communities = __toponeighcom1(node, digraph, status, adjmx)
                __remove(node, com_node, out_edgeneigh_communities.get(com_node, 0.), in_edgeneigh_communities.get(com_node, 0.), status)
                best_com = com_node
                out_neighborcom_mxwght = [key for (key, value) in out_edgeneigh_communities.items() if value == max(out_edgeneigh_communities.values()) and value != 0]
                out_toponeighborcom_mxwght = [key for (key, value) in out_toponeigh_communities.items() if value == max(out_toponeigh_communities.values()) and value !=0]
                comm_neighborcom_mxwght = list(set(out_neighborcom_mxwght) & set(out_toponeighborcom_mxwght))
                if len(comm_neighborcom_mxwght) == 1:
                    best_com = comm_neighborcom_mxwght[0]
                elif len(comm_neighborcom_mxwght) > 1:
                    com_pop_dict = {key: status.com_attr1[key] for key in comm_neighborcom_mxwght}
                    com_pop_list_min = [key for (key, value) in com_pop_dict.items() if value == min(com_pop_dict.values())]
                    best_com = com_pop_list_min[0]
                else:
                    pass
                __insert(node, best_com, out_edgeneigh_communities.get(best_com, 0.), in_edgeneigh_communities.get(best_com, 0.), status)
                if best_com != com_node:                       
                    if node in nodes_move:
                        nodes_move.remove(node)
                else:
                    if len(comm_neighborcom_mxwght) > 0:
                        nodes_move.remove(node)        
        new_unassignedNodes = len(nodes_move)
        if new_unassignedNodes == unassignedNodes:
            break
        unassignedNodes = new_unassignedNodes

    # step 3
    while True:
        for node in digraph.nodes():
            if node in nodes_move:
                com_node = status.node2com[node]
                out_edgeneigh_communities, in_edgeneigh_communities, out_toponeigh_communities = __toponeighcom1(node, digraph, status, adjmx)
                __remove(node, com_node, out_edgeneigh_communities.get(com_node, 0.), in_edgeneigh_communities.get(com_node, 0.), status)
                best_com = com_node                
                if all(value == 0 for value in out_toponeigh_communities.values()):
                    if com_node != node:
                        com_pop_dict = {key: status.com_attr1[key] for key in out_toponeigh_communities.keys() if key > 0}
                        com_pop_list_min = [key for (key, value) in com_pop_dict.items() if value == min(com_pop_dict.values())]
                        if len(com_pop_list_min) > 0:
                            best_com = com_pop_list_min[0] 
                else:
                    out_toponeighborcom_mxwght = [key for (key, value) in out_toponeigh_communities.items() if value == max(out_toponeigh_communities.values())]
                    if len(out_toponeighborcom_mxwght) == 1:
                        best_com = out_toponeighborcom_mxwght[0]
                    else:
                        com_pop_dict = {key: status.com_attr1[key] for key in out_toponeighborcom_mxwght}
                        com_pop_list_min = [key for (key, value) in com_pop_dict.items() if value == min(com_pop_dict.values())]
                        best_com = com_pop_list_min[0]
                __insert(node, best_com, out_edgeneigh_communities.get(best_com, 0.), in_edgeneigh_communities.get(best_com, 0.), status)
                if best_com != com_node:
                    if node in nodes_move:
                        nodes_move.remove(node)
                else:
                    if len(comm_neighborcom_mxwght) > 0:
                        nodes_move.remove(node)
        new_unassignedNodes = len(nodes_move)
        if new_unassignedNodes == unassignedNodes or new_unassignedNodes == 0:            
            break
        unassignedNodes = new_unassignedNodes

    # to ensure the spatial adjacency and at least one hospital
    com2node_dict_temp = copy.deepcopy(status.com2node)
    com2node_dict_temp = {key: value for (key, value) in status.com2node.items() if len(value) > 0}
    com2node_dict = copy.deepcopy(com2node_dict_temp)
    for com in com2node_dict_temp.keys():
        nodes = com2node_dict_temp[com][:]
        if com not in nodes:
            com_lists = [key for key, value in com2node_dict_temp.items() if com in value]
            if len(com_lists) > 0:
                best_com = com_lists[0]
                for node in nodes:
                    com2node_dict[com].remove(node)
                    com2node_dict[best_com].append(node)
    com2node_dict = {key: value for (key, value) in com2node_dict.items() if len(value) > 0}
    return com2node_dict

# refined Dartmouth method to account for population constraint
def __second_level_pop_Dartmouth(digraph, status, com2node_dict, adjmx, attr1_limit, isloop):
    com2node_dict_temp = copy.deepcopy(com2node_dict)
    for node in digraph.nodes():
        com_node = status.node2com[node]
        com_node_attr1 = status.com_attr1[com_node]
        if com_node_attr1 < attr1_limit:
            out_edgeneigh_communities, in_edgeneigh_communities, out_toponeigh_communities = __toponeighcom2(node, digraph, status, adjmx)
            __remove(node, com_node, out_edgeneigh_communities.get(com_node, 0.), in_edgeneigh_communities.get(com_node, 0.), status)
            best_com = com_node
            out_neighborcom_mxwght = [key for (key, value) in out_edgeneigh_communities.items() if value == max(out_edgeneigh_communities.values())]
            out_toponeighborcom_mxwght = [key for (key, value) in out_toponeigh_communities.items() if value == max(out_toponeigh_communities.values())]
            comm_neighborcom_mxwght = list(set(out_neighborcom_mxwght) & set(out_toponeighborcom_mxwght))
            if isloop is False:
                if len(comm_neighborcom_mxwght) == 1:
                    best_com = comm_neighborcom_mxwght[0]
                elif len(comm_neighborcom_mxwght) > 1:
                    com_pop_dict = {key: status.com_attr1[key] for key in comm_neighborcom_mxwght}
                    com_pop_list_min = [key for (key, value) in com_pop_dict.items() if value == min(com_pop_dict.values())]
                    best_com = com_pop_list_min[0]
                else:
                    pass
            else:                
                com_pop_dict = {key: status.com_attr1[key] for key in out_toponeighborcom_mxwght}
                com_pop_dict_sort = dict(sorted(com_pop_dict.items(), key=lambda value: value[1]))
                neigh_coms = {key : value for (key, value) in com_pop_dict_sort.items() if value >= attr1_limit}
                if len(neigh_coms) == 0:
                    best_com = list(com_pop_dict_sort.keys())[-1]   
                else:
                    best_com = [key for (key, value) in neigh_coms.items() if value == min(neigh_coms.values())][0]
            __insert(node, best_com, out_edgeneigh_communities.get(best_com, 0.), in_edgeneigh_communities.get(best_com, 0.), status)
            if best_com != com_node:
                node_origin_list = com2node_dict[com_node][:]
                for nd in node_origin_list:
                    com2node_dict_temp[com_node].remove(nd)  # remove
                    com2node_dict_temp[best_com].append(nd)  # update the com2node_dict

    com2node_dict = copy.deepcopy(com2node_dict_temp)            
    return com2node_dict


# refined Dartmouth method to account for LI constraint
def __second_level_li_Dartmouth(digraph, status, com2node_dict, adjmx, li_limit, isloop):
    com2node_dict_temp = copy.deepcopy(com2node_dict)
    for node in digraph.nodes():
        com_node = status.node2com[node]
        com_node_li = status.li[com_node]
        if com_node_li < li_limit:
            out_edgeneigh_communities, in_edgeneigh_communities, out_toponeigh_communities = __toponeighcom2(node, digraph, status, adjmx)
            __remove(node, com_node, out_edgeneigh_communities.get(com_node, 0.), in_edgeneigh_communities.get(com_node, 0.), status)
            best_com = com_node            
            # find the one neighbor community with the maximal weight
            out_neighborcom_mxwght = [key for (key, value) in out_edgeneigh_communities.items() if value == max(out_edgeneigh_communities.values())]
            out_toponeighborcom_mxwght = [key for (key, value) in out_toponeigh_communities.items() if value == max(out_toponeigh_communities.values())]
            comm_neighborcom_mxwght = list(set(out_neighborcom_mxwght) & set(out_toponeighborcom_mxwght))
            if isloop is False:
                if len(comm_neighborcom_mxwght) == 1:
                    best_com = comm_neighborcom_mxwght[0]
                elif len(comm_neighborcom_mxwght) > 1: # select the one with minimal LI
                    com_li_dict = {key: status.li[key] for key in comm_neighborcom_mxwght}
                    com_li_list_min = [key for (key, value) in com_li_dict.items() if value == min(com_li_dict.values())]
                    if len(com_li_list_min) == 1:
                        best_com = com_li_list_min[0]
                    elif len(com_li_list_min) > 1: # if multiple communities have the same li, select the one with minimal population
                        com_pop_dict = {key: status.com_attr1[key] for key in com_li_list_min}
                        com_pop_list_min = [key for (key, value) in com_pop_dict.items() if value == min(com_pop_dict.values())]
                        best_com = com_pop_list_min[0] # if multiple mini population, select the first one
                    else:
                        print("bugs in the code, please check!")
                else:
                    pass
                    #print("the neighbor and toponeighbor to the community with maximal weight is different")
            else:
                if len(out_toponeighborcom_mxwght) == 1:
                    best_com = out_toponeighborcom_mxwght[0]
                elif len(out_toponeighborcom_mxwght) > 1:
                    com_li_dict = {key: status.li[key] for key in out_toponeighborcom_mxwght}
                    com_li_list_min = [key for (key, value) in com_li_dict.items() if value == min(com_li_dict.values())]
                    if len(com_li_list_min) == 1:
                        best_com = com_li_list_min[0]
                    elif len(com_li_list_min) > 1: # if multiple communities have the same li, select the one with minimal population
                        com_pop_dict = {key: status.com_attr1[key] for key in com_li_list_min}
                        com_pop_list_min = [key for (key, value) in com_pop_dict.items() if value == min(com_pop_dict.values())]
                        best_com = com_pop_list_min[0] # if multiple mini population, select the first one
                    else:
                        print("bugs in the code, please check!")
            __insert(node, best_com, out_edgeneigh_communities.get(best_com, 0.), in_edgeneigh_communities.get(best_com, 0.), status)
            if best_com != com_node:
                node_origin_list = com2node_dict[com_node][:] # keep consistent with status
                for nd in node_origin_list:
                    com2node_dict_temp[com_node].remove(nd)  # remove, it is okay to remove here
                    com2node_dict_temp[best_com].append(nd)  # update the com2node_dict
    com2node_dict = copy.deepcopy(com2node_dict_temp)            
    return com2node_dict
        

class Status(object):
    """
    To handle several data in one struct.
    """
    node2com = {} # the community of each node
    com2node = {} # ids of nodes in each community
    com_attr1 = {} # one attribute of each community, e.g., total population
    node_attr1 = {} # one attribute of each node, e.g., total population
    total_weight = 0 # sum of weights of all edges
    internals = {} # sum of weights of all edges within a community
    degrees = {} # degree list for all communities    
    out_degrees = {} # outdegree list for all communities
    in_degrees = {} # indegree list for all communities
    gdegrees = {} # degree list for all nodes in graph    
    out_gdegrees = {} # outdegree list for all nodes in graph
    in_gdegrees = {} # indegree list for all nodes in graph
    li = {}

    def __init__(self):
        self.node2com = dict([])
        self.com2node = dict([])
        self.com_attr1 = dict() # its value of each key is a number, not a list

        self.node_attr1 = dict()
        self.total_weight = 0
        self.degrees = dict([])
        self.out_degrees = dict([])  
        self.in_degrees = dict([])              
        self.gdegrees = dict([])        
        self.out_gdegrees = dict([])
        self.in_gdegrees = dict([])
        self.internals = dict([]) # edge weights list for all communities
        self.loops = dict([]) # edge weights list for all nodes in graph
        self.li = dict() # localilzation index all communities

    def copy(self):
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = copy.deepcopy(self.node2com)
        new_status.com2node = copy.deepcopy(self.com2node)
        new_status.com_attr1 = copy.deepcopy(self.com_attr1)
        new_status.node_attr1 = copy.deepcopy(self.node_attr1)
        new_status.internals = copy.deepcopy(self.internals)
        new_status.degrees = copy.deepcopy(self.degrees)
        new_status.out_degrees = copy.deepcopy(self.out_degrees)
        new_status.in_degrees = copy.deepcopy(self.in_degrees)        
        new_status.gdegrees = copy.deepcopy(self.gdegrees)
        new_status.out_gdegrees = copy.deepcopy(self.out_gdegrees)
        new_status.in_gdegrees = copy.deepcopy(self.in_gdegrees)        
        new_status.loops = copy.deepcopy(self.loops)
        new_status.total_weight = self.total_weight
        new_status.li = copy.deepcopy(self.li)

        return new_status

    def init(self, digraph, part=None):
        """Initialize the status of a graph with every node in one community"""
        count = 1 # new community id starts from 0
        self.node2com = dict([])
        self.com2node = dict([])
        self.com_attr1 = dict()
        self.node_attr1 = dict()
        self.total_weight = 0
        self.degrees = dict([])
        self.out_degrees = dict([]) 
        self.in_degrees = dict([])               
        self.gdegrees = dict([])
        self.out_gdegrees = dict([])
        self.in_gdegrees = dict([])        
        self.internals = dict([])
        self.total_weight = digraph.size(weight='weight') # .size function returns the number of edges or sum of edge weights in the graph
        self.li = dict()

        if part is None:
            self.com2node = {k: [] for k in range(1, digraph.number_of_nodes() + 1)} # initialize a dictionary of empty list
            for node in digraph.nodes(): # node is the node name
                self.node2com[node] = count # if no partition provided, each node belongs to one community             
                self.node_attr1[node] = digraph.nodes[node]['pop']                
                if bool(self.com2node.get(count)):
                    self.com2node[count].append(node)
                    self.com_attr1[count] += digraph.nodes[node]['pop']
                else:
                    self.com2node[count]= [node]
                    self.com_attr1[count] = digraph.nodes[node]['pop']
                deg = float(digraph.degree(node, weight='weight')) # calculate the degree of each node
                out_deg = float(digraph.out_degree(node, weight = 'weight'))
                in_deg = float(digraph.in_degree(node, weight = 'weight'))
                if deg < 0 or in_deg < 0 or out_deg < 0:
                    error = "Bad directed graph type ({})".format(type(digraph))
                    raise ValueError(error)
                self.degrees[count] = deg
                self.out_degrees[count] = out_deg
                self.in_degrees[count] = in_deg                
                self.gdegrees[node] = deg
                self.out_gdegrees[node] = out_deg
                self.in_gdegrees[node] = in_deg
                edge_data = digraph.get_edge_data(node, node, {"weight": 0})
                self.loops[node] = float(edge_data.get("weight", 1))
                self.internals[count] = self.loops[node]
                self.li[count] = -1
                if self.out_degrees[count] > 0:
                    self.li[count] = self.internals[count]/self.out_degrees[count] 
                count += 1
        else: # when a partition is provided as input
            self.com2node = {k: [] for k in list(set(part.values()))}
            for node in digraph.nodes():
                com = part[node]
                self.node2com[node] = com
                self.node_attr1[node] = digraph.nodes[node]['pop']                
                if bool(self.com2node.get(com)):
                    self.com2node[com].append(node)                    
                    self.com_attr1[com] += digraph.nodes[node]['pop']
                else:                    
                    self.com2node[com]= [node]
                    self.com_attr1[com] = digraph.nodes[node]['pop']
                deg = float(digraph.degree(node, weight='weight'))
                out_deg = float(digraph.out_degree(node, weight = 'weight'))
                in_deg = float(digraph.in_degree(node, weight = 'weight'))
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.out_degrees[com] = self.out_degrees.get(com, 0) + out_deg
                self.in_degrees[com] = self.in_degrees.get(com, 0) + in_deg
                self.gdegrees[node] = deg
                self.out_gdegrees[node] = out_deg
                self.in_gdegrees[node] = in_deg              

                inc = 0.
                for neighbor, datas in digraph[node].items():                    
                    edge_weight = datas.get("weight", 1)
                    if edge_weight <= 0:
                        raise ValueError("Bad directed graph type, use positive weights")
                    if part[neighbor] == com: # if one of this node's neighbors is also in the same community
                        inc += float(edge_weight) # inc indicates the sum of edge weights inside a community
                self.internals[com] = self.internals.get(com, 0) + inc 
            
            for com in self.com2node.keys():
                self.li[com] = -1
                if self.out_degrees[com] > 0:
                    self.li[com] = self.internals[com]/self.out_degrees[com]

# check if two lists of polygons are adjacent
def __isadj(nodelist, neighborlist, adjmx):
    nodelist = [(node-1) for node in nodelist]
    neighborlist = [(neighbor -1) for neighbor in neighborlist]
    rows_bl = np.isin(adjmx[0], nodelist)
    col_bl = np.isin(adjmx[1], neighborlist)

    rows_index = np.nonzero(rows_bl)
    cols_index = np.nonzero(col_bl)
    comm_index = np.intersect1d(rows_index, cols_index)
    if(len(comm_index) > 0):
        return True
    else:        
        return False

# find the topology and flow adjacent communities for directed network
def __toponeighcom1(node, digraph, status, adjmx):
    out_weights_EDGEneighbor = {}
    in_weights_EDGEneighbor = {}
    out_weights_TOPOneighbor = {}
    out_neighcom_list = list()
    out_neighbors = list(digraph.successors(node))
    in_neighbors = list(digraph.predecessors(node))
    if node in in_neighbors:
        in_neighbors.remove(node)

    for out_neighbor in out_neighbors:
        out_weight = digraph.get_edge_data(node, out_neighbor, {"weight": 0}).get("weight", 0)
        out_neighborcom = status.node2com[out_neighbor]                
        out_weights_EDGEneighbor[out_neighborcom] = out_weights_EDGEneighbor.get(out_neighborcom, 0) + out_weight

        out_neighbors2 = status.com2node[out_neighborcom]
        isadj = __isadj([node], out_neighbors2, adjmx) # the node and the community (all nodes)
        if node in out_neighbors2:
            isadj = True        
        if isadj:
            out_weights_TOPOneighbor[out_neighborcom] = out_weights_TOPOneighbor.get(out_neighborcom, 0) + out_weight
            out_neighcom_list.append(out_neighborcom)
    
    for in_neighbor in in_neighbors:
        in_weight = digraph.get_edge_data(in_neighbor, node, {"weight": 0}).get("weight", 0)
        in_neighborcom = status.node2com[in_neighbor]        
        in_weights_EDGEneighbor[in_neighborcom] = in_weights_EDGEneighbor.get(in_neighborcom, 0) + in_weight
    
    rows_bl = np.isin(adjmx[0], [node-1])
    rows_index = np.nonzero(rows_bl)
    node_toponbr = np.take(adjmx[1], rows_index)[0]
    node_toponbr = [(value + 1) for value in node_toponbr]

    node_toponbr_coms = [value for (key, value) in status.node2com.items() if key in node_toponbr]
    out_toponbr_coms = set(node_toponbr_coms) - set(out_neighcom_list)

    if len(out_toponbr_coms) > 0:
        for out_toponbr_com in out_toponbr_coms:
            out_node_toponbrs = status.com2node[out_toponbr_com]
            for out_node_toponbr in out_node_toponbrs:
                if out_node_toponbr != node:
                    out_edge_data = digraph.get_edge_data(node, out_node_toponbr)
                    out_weight = 0
                    if out_edge_data is not None:
                        out_weight = out_edge_data.get("weight")
                    out_weights_TOPOneighbor[out_toponbr_com] = out_weights_TOPOneighbor.get(out_toponbr_com, 0) + out_weight
    return out_weights_EDGEneighbor, in_weights_EDGEneighbor, out_weights_TOPOneighbor

# find the topology and flow adjacent communities for directed network, used in population and LI constraints
def __toponeighcom2(node, graph, status, adjmx):
    out_weights_EDGEneighbor = {}
    in_weights_EDGEneighbor = {}
    out_weights_TOPOneighbor = {}
    com_node = status.node2com[node]
    out_neighcom_list = list()

    out_neighbors = list(graph.successors(node))
    if len(out_neighbors) > 0:
        if node in out_neighbors:
            out_neighbors.remove(node)

    in_neighbors = list(graph.predecessors(node))
    if len(in_neighbors) > 0:
        if node in in_neighbors:
            in_neighbors.remove(node)

    node_toponbr_com = np.where(adjmx[node -1] == 1)[0] # get the toponeigbor communities
    node_toponbr_com = [(key +1) for key in node_toponbr_com]
    
    for out_neighbor in out_neighbors:
        out_weight = graph.get_edge_data(node, out_neighbor, {"weight": 0}).get("weight", 0)
        out_neighborcom = status.node2com[out_neighbor]
        out_weights_EDGEneighbor[out_neighborcom] = out_weights_EDGEneighbor.get(out_neighborcom, 0) + out_weight

        out_neighbors2 = status.com2node[out_neighborcom]
        comm_out_neighbors2 = list(set(node_toponbr_com) & set(out_neighbors2))
        if len(comm_out_neighbors2) > 0:
            out_weights_TOPOneighbor[out_neighborcom] = out_weights_TOPOneighbor.get(out_neighborcom, 0) + out_weight
            out_neighcom_list.append(out_neighborcom)

    for in_neighbor in in_neighbors:
        in_weight = graph.get_edge_data(in_neighbor, node, {"weight": 0}).get("weight", 1)
        in_neighborcom = status.node2com[in_neighbor]
        in_weights_EDGEneighbor[in_neighborcom] = in_weights_EDGEneighbor.get(in_neighborcom, 0) + in_weight

    out_node_toponbr_com = set(node_toponbr_com) - set(out_neighcom_list)

    if len(out_node_toponbr_com) > 0:
        for out_neighborcom in out_node_toponbr_com:
            out_edge_data = graph.get_edge_data(com_node, out_neighborcom)
            out_weight = 0
            if out_edge_data is not None:
                out_weight = out_edge_data.get("weight")
            out_weights_TOPOneighbor[out_neighborcom] = out_weights_TOPOneighbor.get(out_neighborcom, 0) + out_weight
    return out_weights_EDGEneighbor, in_weights_EDGEneighbor, out_weights_TOPOneighbor

# remove a node from a community
def __remove(node, com, out_weight, in_weight, status):
    """ Remove node from community com and modify status"""
    status.degrees[com] = (status.degrees.get(com, 0.) - status.gdegrees.get(node, 0.))
    status.out_degrees[com] = (status.out_degrees.get(com, 0.) - status.out_gdegrees.get(node, 0.))
    status.in_degrees[com] = (status.in_degrees.get(com, 0.) - status.in_gdegrees.get(node, 0.))
    status.internals[com] = float(status.internals.get(com, 0.) - out_weight - in_weight - status.loops.get(node, 0.))  
    status.node2com[node] = 0 # then this node has no community assigned
    status.com2node[com].remove(node) # remove this node ('node' is the node name here not the index)
    status.li[com] = 0
    if status.out_degrees[com] > 0:
        status.li[com] = status.internals[com]/status.out_degrees[com]
    # from the original community
    status.com_attr1[com] -= status.node_attr1[node] # update the attr1 of this community

# insert a node into a community
def __insert(node, com, out_weight, in_weight, status):
    """ Insert node into community and modify status"""
    status.node2com[node] = com
    if bool(status.com2node.get(com)): # if there is at least a value of key=com in the dictionary
        status.com2node[com].append(node)
        status.com_attr1[com] += status.node_attr1[node]
    else:
        # status.com2node[com] = [node] # make it a list
        status.com2node[com].append(node)
        status.com_attr1[com] = status.node_attr1[node]
    #print(status.gdegrees.get(node, 0.)-status.out_gdegrees.get(node, 0.) - status.in_gdegrees.get(node, 0.))
    status.degrees[com] = (status.degrees.get(com, 0.) + status.gdegrees.get(node, 0.))
    status.out_degrees[com] = (status.out_degrees.get(com, 0.) + status.out_gdegrees.get(node, 0.))
    status.in_degrees[com] = (status.in_degrees.get(com, 0.) + status.in_gdegrees.get(node, 0.))
    status.internals[com] = float(status.internals.get(com, 0.) + out_weight + in_weight + status.loops.get(node, 0.))
    status.li[com] = 0
    if status.out_degrees[com] > 0:
        status.li[com] = status.internals[com]/status.out_degrees[com]


