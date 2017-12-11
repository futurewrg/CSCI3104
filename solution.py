#!/usr/bin/env python2
#-*- coding: UTF-8 -*-
#Data
import random
import math
import networkx as nx
import matplotlib.pyplot as plt
import collections

#problem2_c
def BFS(graph, s, t, parent):
    '''Returns true if there is a path from source 's' to sink 't' in
    residual graph. Also fills parent[] to store the path '''

    # Mark all the vertices as not visited
    nodes = graph.nodes()
    visited = {node: False for node in nodes}

    # Create a queue for BFS
    queue = collections.deque()

    # Mark the source node as visited and enqueue it
    queue.append(s)
    visited[s] = True

    # Standard BFS Loop
    while queue:
        u = queue.popleft()

        # Get all adjacent vertices's of the dequeued vertex u
        # If a adjacent has not been visited, then mark it
        # visited and enqueue it
        for ind in graph.adj[u]:
            val = graph[u][ind]['capacity']
            if visited[ind] == False and val > 0 :
                queue.append(ind)
                visited[ind] = True
                parent[ind] = u
 
        # If we reached sink in BFS starting from source, then return
        # true, else false
    return visited[t]

# Returns the maximum flow from s to t in the given graph
def FordFulkerson(graph, source, sink):

    # This array is filled by BFS and to store path
    parent = dict()

    max_flow = 0 # There is no flow initially

    # Augment the flow while there is path from source to sink
    while BFS(graph, source, sink, parent) :

        # Find minimum residual capacity of the edges along the
        # path filled by BFS. Or we can say find the maximum flow
        # through the path found.
        path_flow = float("Inf")
        s = sink
        while s !=  source:
            path_flow = min (path_flow, graph[parent[s]][s]['capacity'])
            s = parent[s]

        # Add path flow to overall flow
        max_flow +=  path_flow

        # update residual capacities of the edges and reverse edges
        # along the path
        v_node = sink
        while v_node !=  source:
            u_node = parent[v_node]
            graph[u_node][v_node]['capacity'] -= path_flow
            graph[v_node][u_node]['capacity'] += path_flow
            v_node = parent[v_node]

    return max_flow

#problem2_c
#input a graph with 
def problem2_c(graph, issue_param, person_param):
    """
    convert the graph G{P ∪ I , E} to G^ = { V^, E^}
    
    """
    new_graph = nx.DiGraph()

    #--------------Part1-----------------------------
    #convert the graph G{P ∪ I , E} to G' = { V', E'}
    # V' = {budget, survay}∪ P ∪ I
    # E' = {(budget, p)| p ∈ P} ∪ {(i, survay) | i ∈ I} ∪ E}

    #Add edges E
    new_graph.add_edges_from([(edge[0], edge[1], {'capacity':1}) for edge in graph.edges])

    # Add edges  {(budget, p)| p ∈ P} , 
    new_graph.add_edges_from([('budget', person, {'capacity': bp_tp[1] - bp_tp[0]}) for person, bp_tp in person_param.items()])
    # Add edges {(i, survay) | i ∈ I} ∪ E} 
    new_graph.add_edges_from([(issue, 'survay', {'capacity': li_ui[1] - li_ui[0]}) for issue, li_ui in issue_param.items()])

    
    #---------------part2----------------------------
    #convert the graph G' = {V' , E'} to G^ = { V^, E^}
    # V^ = {s,t}∪ V'
    # E^ = {(s, p), ('budget',t)| ('budget', p) ∈ E' }
    #      ∪{(s, 'survay') , (i, t) | (i, survay) ∈ E' }
    #      ∪E'

    # add edges (s,p) and ('budget',t)
    # capacity (s,p) and capacity ('budget', t) both are b_p
    new_graph.add_edges_from([('s', person, {'capacity': bp_tp[0]}) for person, bp_tp in person_param.items()])
    new_graph.add_edge('survay', 'budget', capacity = float('Inf'))
    cap = 0
    for person, bp_tp in person_param.items():
        cap += bp_tp[0]
    new_graph.add_edge('budget', 't', capacity = cap)

    # add edges (i,t) and (s, 'survay')
    # capacity (s, 'survay') and capacity(i,t) both are l_i
    new_graph.add_edges_from([(issue, 't', {'capacity': li_ui[0]}) for issue, li_ui in issue_param.items()])
    
    cap = 0
    for issue, li_ui in issue_param.items():
        cap += li_ui[0]
    new_graph.add_edge('s', 'survay', capacity = cap)
    return new_graph

def problem2_d(graph):
    residual_graph = nx.DiGraph()

    residual_graph.add_edges_from([(edge[1], edge[0], {'capacity':0}) for edge in graph.edges()])
    residual_graph.add_edges_from([(edge[0], edge[1], {'capacity': graph[edge[0]][edge[1]]['capacity']}) for edge in graph.edges()])

    max_flow = FordFulkerson(residual_graph, 's', 't')
    full_flow = graph['s']['survay']['capacity']
    if max_flow == full_flow:
        print "Success, the parameters are feasible."
    else:
        print "Failure, the parameters are unfeasible."
    
def problem2_e():
    """Generate test case for problem2_c
    Return:
    G --A graph G = {P ∪ I , E}
    person_param --A dict for {issue_id : (l_i, u_i)}
    issure_param --A dict for {person_id : (b_p, t_p) }
    """
    G = nx.DiGraph()
    issue_size = 10  #n = 10
    issue_id = ['i'+str(i) for i in range(0,issue_size)]
    issue_param = dict() # A dict for {issure_id:(l_i,u_i)}
    
    person_size = 1000  #m = 1000
    person_id = ['p'+str(i) for i in range(0,person_size)]
    person_param = dict() # A dict for {person_id: (b_p, t_p)}
    
    for person in person_id:
        t_p = 0
        for issue in issue_id:
            # a person has a probability of 50% to have a opition aboult the issue
           if random.randint(0, 1) == 0:
               G.add_edge(person, issue)
               t_p += 1
        b_p =  t_p / 2
        person_param [person] = (b_p, t_p)

    for issue in issue_id:
        l_i = random.randint(300,400)   # l_i is drawn unifromly from [300,400]
        u_i = random.randint(500,700)   # u_i is drawn uniformly from [500,700]
        issue_param[issue] = (l_i, u_i)
    return G, issue_param, person_param

    
graph, issue_param, person_param = problem2_e()
graph = problem2_c(graph, issue_param, person_param)

""" if you want to test problem2, uncomment next line."""
#problem2_d(graph)




    
#problem4
#input: a graph
#return: the number of all of the simple paths from Spider to Fly
def problem4(graph):
    edges = graph.edges()
    v_reach = set()            #store the reached node
    v = set()                  #store the unreached node
    prev_nodes = dict()        #store the previous node set of every node 
    simple_paths = dict()      #store the number of simple paths from Spider to node
    for edge in edges:
        #init the prev_nodes
        if edge[1] in prev_nodes:
            prev_nodes[edge[1]].add(edge[0])
        else:
            prev_nodes[edge[1]] = set()
            prev_nodes[edge[1]].add(edge[0])

        #init the node set V
        v.add(edge[0])
        v.add(edge[1])

    #init the reached node set V_reach
    v_reach.add('Spider')
    v.remove('Spider')
    simple_paths['Spider'] = 1
    
    while 'Fly' not in v_reach and len(v) > 0:
        for ver in v:
            #if ver is not reached and all of the previous node of ver is reached
            #put the ver into reached set V_reach
            if ver not in v_reach and prev_nodes[ver].issubset(v_reach):
                simple_paths[ver] = 0
                v_reach.add(ver)
                for u in prev_nodes[ver]:#caculate the number of simple paths from 'Fly' to ver
                    simple_paths[ver] += simple_paths[u]
        v = v - v_reach      # remove the reached ver from v;
    return simple_paths['Fly']
                
        
#probelm4_test
#  generater the graph for probelm4
def problem4_test():
    #create a directed graph
    G = nx.DiGraph()

    #adding an edge also adds the node
    G.add_edge('Spider', 'A', weight=1.0)
    G.add_edge('Spider', 'H', weight=1.0)
    G.add_edge('Spider', 'J', weight=1.0)
    
    G.add_edge('H', 'G', weight=1.0)
    G.add_edge('H', 'K', weight=1.0)
    
    G.add_edge('G', 'L', weight=1.0)
    G.add_edge('G', 'F', weight=1.0)

    G.add_edge('F', 'E', weight=1.0)

    G.add_edge('E', 'Fly', weight=1.0)

    G.add_edge('J', 'S', weight=1.0)
    G.add_edge('J', 'K', weight=1.0)

    G.add_edge('K', 'L', weight=1.0)
    G.add_edge('L', 'M', weight=1.0)
    G.add_edge('M', 'N', weight=1.0)
    G.add_edge('M', 'F', weight=1.0)

    G.add_edge('N', 'O', weight=1.0)
    G.add_edge('N', 'E', weight=1.0)
    
    G.add_edge('O', 'Fly', weight=1.0)

    G.add_edge('A', 'S', weight=1.0)
    G.add_edge('A', 'B', weight=1.0)

    G.add_edge('B', 'R', weight=1.0)
    G.add_edge('B', 'C', weight=1.0)

    G.add_edge('S', 'R', weight=1.0)
    G.add_edge('R', 'Q', weight=1.0)
    
    G.add_edge('Q', 'C', weight=1.0)
    G.add_edge('Q', 'P', weight=1.0)

    G.add_edge('C', 'D', weight=1.0)
    G.add_edge('D', 'Fly', weight=1.0)
    G.add_edge('P', 'D', weight=1.0)
    G.add_edge('P', 'O', weight=1.0)
    G.add_edge('O', 'Fly', weight=1.0)

    G.add_edge('T', 'Q', weight=1.0)
    G.add_edge('T', 'P', weight=1.0)
    G.add_edge('T', 'O', weight=1.0)
    G.add_edge('T', 'N', weight=1.0)
    G.add_edge('T', 'M', weight=1.0)

    G.add_edge('R', 'T', weight=1.0)
    G.add_edge('S', 'T', weight=1.0)
    G.add_edge('J', 'T', weight=1.0)
    G.add_edge('K', 'T', weight=1.0)
    G.add_edge('L', 'T', weight=1.0)
    return G        
    #end probelm4_test

def problem5(dis_list, num_person):
    """problem 5

    Input :

    dis_list --the distance tuple list
        we use a distance tuple (personA, personB, dis) to discrible the distance between personA and personB 

    num_person --the total number of people

    Return:
    True  --at least one person not hit
    False --everyone is hit
    """    

    neighbor_hit = set();

    person_hurl = set();
    
    dis_list = sorted(dis_list, cmp = lambda x,y: cmp(x[2],y[2])) 
    
    for dis in dis_list:
        personA = dis[0]
        personB = dis[1]

        #personA have the nearest neighbor personB, so hurl balloon on personB
        if personA not in person_hurl:
            neighbor_hit.add(personB)
            person_hurl.add(personA)

        #personB have the nearest neighbor personA, so hurl balloon on personA
        if personB not in person_hurl:
            neighbor_hit.add(personA)
            person_hurl.add(personB)

    #Exist people not hit by balloon, return true
    if len(neighbor_hit) != num_person:
        print "True : there always remains at least one person not hit by a balloon"
    else:
        print "False : everyone hit by a balloon"
    return neighbor_hit

def problem5_test():
    """
    generater the test data for problem 5.
    ensure that everyone have a unique nearest neighbor.
    """
    num_person = 2 * int( random.uniform(1, 20)) + 1
    dis_list = []
    
    min_dis = []
    pos_list = []
    for personA in range(0, num_person):
        unique = False

        # randomly generater the position of personA
        #if not satisfy everyone have unique nearest neighbor, randomly generater again.
        while (not unique):
            personA_x = random.uniform(1,1000)
            personA_y = random.uniform(1,1000)
            personA_set = dict()

            minA = [float('Inf'), 1]

            #check whether everyone have unique nearest neighbor after add the personA. 
            for pos in pos_list:
                x = pos[0] - personA_x
                y = pos[1] - personA_y
                personB = pos[2]
                dis = x * x + y * y

                #personB has two nearest neighbors, regenerator the position of personA
                if(min_dis[personB] == pos):
                    minA[1] = 2
                    break
                
                if dis == minA[0]:
                    minA[1] += 1
                elif dis < minA[0]:
                    minA= [dis,1]
            #personA has two nearest neighbor, regenerator the position of personA
            if(minA[1] == 1):
                break
            #end for
        #end while

        #print personA, personA_x, personA_y
        pos_list.append([personA_x, personA_y, personA])
        min_dis.append(minA[0])
        for idB in range(0, personA):
            personB = pos_list[idB]
            x = personA_x - personB[0]
            y = personA_y - personB[1]
            dis = x*x + y * y
            dis_list.append([personA, idB, dis])
            if dis < min_dis[idB]:
                min_dis[idB] = dis;
    #end for
    return dis_list, num_person


#problem5_test()
#result = problem5(dis_list, num_person)    
#print problem4(problem4_test())

    
    
