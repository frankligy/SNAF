#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle
import plotly.graph_objects as go
import networkx as nx

def disentangle_uid(uid,add_stable_id=False,add_symbol=False):
    fg,bg = uid.split('|')
    gene_symbol = fg.split(':')[0]
    stable_id = fg.split(':')[1]
    exons = ':'.join(fg.split(':')[2:])
    output = exons
    if add_stable_id:
        output = stable_id + ':' + output
    if add_symbol:
        output = gene_symbol + ':' + output
    return output


network = pd.read_csv('../altanalyze_output/leucegene_inference/bbsr_weight1/network_aug.txt',sep='\t',index_col=0)
all_sfs = network.columns.tolist()
network = network.stack().reset_index()
network.columns = ['event','sf','weight']
network.sort_values(by='weight',ascending=False,inplace=True)
network.set_index(keys=np.arange(network.shape[0]),inplace=True)
network['percentile'] = 1- network.index.values / network.shape[0]
# no query for sf, but for events
psi = pd.read_csv('/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.U2AF1-Q157_variants_vs_Others.txt',sep='\t',index_col=0)
events = [disentangle_uid(item,add_symbol=True) for item in psi.index]
network = network.loc[network['event'].isin(events),:]
# queries = ['U2AF1','U2AF2','SF3B1','HNRNPK','FUS']
# if queries != 'all':
#     formula_string = ''
#     for q in queries:
#         if q not in all_sfs:
#             sys.exit('{} not in all sfs'.format(q))
#         else:
#             formula_string += '(network[\'sf\']==\'{}\')|'.format(q)
#     formula_string = formula_string.rstrip('|')
#     network = network.loc[eval(formula_string),:]
network = network.loc[network['weight']>0.5,:]
# network = network.sample(frac=0.1)
print(network)


# original
G = nx.from_pandas_edgelist(network,source='event',target='sf',edge_attr='weight')
nx.set_node_attributes(G,nx.spring_layout(G,seed=42),'pos')
edge_x = []
edge_y = []
middle_node_x = []
middle_node_y = []
middle_node_text = []
for edge in G.edges(data=True):
    x0,y0 = G.nodes[edge[0]]['pos']
    x1,y1 = G.nodes[edge[1]]['pos']
    x_middle = (x0 + x1)/2
    y_middle = (y0 + y1)/2
    edge_x.append(x0)
    edge_y.append(y0)
    edge_x.append(x1)
    edge_y.append(y1)
    edge_x.append(None)
    edge_y.append(None)
    middle_node_x.append(x_middle)
    middle_node_y.append(y_middle)
    middle_node_text.append('->'.join(['{}'.format(edge[0]),'{}'.format(edge[1]), str(round(edge[2]['weight'],4))]))

# # complicated
# G = nx.from_pandas_edgelist(network,source='event',target='sf',edge_attr='weight')
# nx.set_node_attributes(G,nx.spring_layout(G,seed=42),'pos')
# edge_x = []
# edge_y = []
# both_middle_node_x = []
# both_middle_node_y = []
# both_middle_node_text = []
# func_middle_node_x = []
# func_middle_node_y = []
# func_middle_node_text = []
# bind_middle_node_x = []
# bind_middle_node_y = []
# bind_middle_node_text = []
# for edge in G.edges(data=True):
#     x0,y0 = G.nodes[edge[0]]['pos']
#     x1,y1 = G.nodes[edge[1]]['pos']
#     x_middle = (x0 + x1)/2
#     y_middle = (y0 + y1)/2
#     edge_x.append(x0)
#     edge_y.append(y0)
#     edge_x.append(x1)
#     edge_y.append(y1)
#     edge_x.append(None)
#     edge_y.append(None)
#     if edge[2]['weight'] == 1:
#         both_middle_node_x.append(x_middle)
#         both_middle_node_y.append(y_middle)
#         both_middle_node_text.append('->'.join(['{}'.format(edge[0]),'{}'.format(edge[1]), str(round(edge[2]['weight'],4))]))
#     elif edge[2]['weight'] == 0.51:
#         bind_middle_node_x.append(x_middle)
#         bind_middle_node_y.append(y_middle)
#         bind_middle_node_text.append('->'.join(['{}'.format(edge[0]),'{}'.format(edge[1]), str(round(edge[2]['weight'],4))])) 
#     elif edge[2]['weight'] == 0.49:
#         func_middle_node_x.append(x_middle)
#         func_middle_node_y.append(y_middle)
#         func_middle_node_text.append('->'.join(['{}'.format(edge[0]),'{}'.format(edge[1]), str(round(edge[2]['weight'],4))]))        

node_x = []
node_y = []
node_name = []
for node in G.nodes(data=True):
    x,y = node[1]['pos']
    node_x.append(x)
    node_y.append(y)
    node_name.append('{} -- Degree {}'.format(node[0],G.degree(node[0])))
node_color = ['red' if ':' in item else 'blue' for item in node_name]
node_size = [4 if ':' in item else 15 for item in node_name]

# original
edge_trace = go.Scatter(x=edge_x,y=edge_y,mode='lines',line={'width':0.1,'color':'black'})
middle_trace = go.Scatter(x=middle_node_x,y=middle_node_y,mode='markers',hoverinfo='text',text=middle_node_text,marker={'color':'green','size':2})
node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',hoverinfo='text',text=node_name,marker={'color':node_color,'size':node_size})
fig = go.Figure(data=[edge_trace,node_trace,middle_trace],layout=go.Layout(showlegend=False))
fig.write_html('../altanalyze_output/leucegene_inference/bbsr_weight1/network.html',include_plotlyjs='cdn')

# # more complicated 
# edge_trace = go.Scatter(x=edge_x,y=edge_y,mode='lines',line={'width':0.1,'color':'black'})
# both_middle_trace = go.Scatter(x=both_middle_node_x,y=both_middle_node_y,mode='markers',hoverinfo='text',text=both_middle_node_text,marker={'color':'green','size':2})
# func_middle_trace = go.Scatter(x=func_middle_node_x,y=func_middle_node_y,mode='markers',hoverinfo='text',text=func_middle_node_text,marker={'color':'orange','size':2})
# bind_middle_trace = go.Scatter(x=bind_middle_node_x,y=bind_middle_node_y,mode='markers',hoverinfo='text',text=bind_middle_node_text,marker={'color':'magenta','size':2})
# node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',hoverinfo='text',text=node_name,marker={'color':node_color,'size':node_size})
# fig = go.Figure(data=[edge_trace,node_trace,both_middle_trace,func_middle_trace,bind_middle_trace],layout=go.Layout(showlegend=False,title='_'.join(queries)))
# fig.write_html('../altanalyze_output/prior/whole_prior/shRNA_K562_rule.html',include_plotlyjs='cdn')





