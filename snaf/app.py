import pandas as pd
from sklearn.decomposition import PCA
from umap import UMAP
import numpy as np
import os,sys
from ast import literal_eval
import dash
from dash import dcc,html
import plotly.graph_objects as go
import plotly.express as px
from dash.dependencies import Input,Output,State



def run_dash_app():
    # global preprocessing

    ### df
    df = pd.read_csv('/Users/ligk2e/Desktop/learn_dash/frequency_stage3.txt',sep='\t',index_col=0)
    data = []
    for item in df.index:
        data.append(item.split(','))
    n_total = len(list(zip(*data)))
    if n_total == 5:
        df.index = list(zip(*data))[0]
        df['HLA'] = list(zip(*data))[1]
        df['extra'] = list(zip(*data))[2]
        df['n_first'] = list(zip(*data))[3]
        df['uid'] = list(zip(*data))[4]
    df = df.loc[np.logical_not(df.index.duplicated()),:]
    df['length'] = [len(item) for item in df.index]
    df['samples'] = [literal_eval(item) for item in df['samples']]  # from str to list


    ### umap
    after_pca = np.loadtxt('/Users/ligk2e/Desktop/learn_dash/after_pca.txt')
    def aaindex(peptide,after_pca):
        amino = 'ARNDCQEGHILKMFPSTWYV-'
        matrix = np.transpose(after_pca)   # [12,21]
        encoded = np.empty([len(peptide), 12])  # (seq_len,12)
        for i in range(len(peptide)):
            query = peptide[i]
            if query == 'X': query = '-'
            query = query.upper()
            encoded[i, :] = matrix[:, amino.index(query)]
        return encoded.flatten()
    df9 = df.loc[df['length']==9,:]
    df10 = df.loc[df['length']==10,:]
    embed = []
    for df_s,l_s in zip([df9,df10],[9,10]):
        mer_encoded = np.empty([df_s.shape[0],l_s*12])
        for i,pep in enumerate(df_s.index):
            mer_encoded[i,:] = aaindex(pep,after_pca)
        model = PCA()
        model.fit(mer_encoded)
        pca_scoring = model.transform(mer_encoded)[:,:55]
        reducer = UMAP(random_state=42,min_dist=0.5)
        embedding = reducer.fit_transform(pca_scoring)
        embed.append(embedding)

    # building app

    app = dash.Dash(__name__)
    dropdown_cols = list(df.columns)
    dropdown_cols.remove('samples')
    dropdown_cols.remove('length')
    dropdown_cols.remove('uid')
    app.layout = html.Div([
        # heading
        html.Div(html.H1('SNAF Neoantigen Viewer')),
        # selection
        html.Div([
            # dropdown label
            html.Label('Metadata to display'),
            dcc.Dropdown(
                id = 'metadata_dropdown',
                options = [{'label':col,'value':col} for col in dropdown_cols],
                value = 'HLA'),
            # length radioitem label
            html.Label('Length'),
            dcc.RadioItems(id='length_radioitem',options=[{'label':9,'value':9},{'label':10,'value':10}],value=9),
            # button to submit
            html.Button(id='submit_button',n_clicks=0,children='Submit')

        ]),
        # graph
        html.Div([
            html.H2('Scatter Plot'),
            dcc.Graph(id='scatter_figure')]),
        # table for display    
        html.Div([
            html.H2('Selected Neoantigen'),
            dcc.Graph(id='display_table')]) 

    ])

    # app callback
    @app.callback(
        Output('scatter_figure','figure'),
        State('metadata_dropdown','value'),
        State('length_radioitem','value'),
        Input('submit_button','n_clicks'))
    def scatter_figure(dropdown_value,length_value,n_clicks):
        nonlocal df
        # filter length
        plot_df = df.copy()
        plot_df = plot_df.loc[plot_df['length']==length_value]  
        embedding = embed[0] if length_value==9 else embed[1]
        # plot
        fig = px.scatter(x=embedding[:,0],y=embedding[:,1],
                        color=plot_df[dropdown_value],text=plot_df.index.values)
        fig.update_traces(mode='markers',customdata=plot_df['uid'].values,hovertemplate='%{text}<br>%{customdata}')
        fig.update_layout(title='test',margin=dict(l=300,r=300),plot_bgcolor='rgba(0,0,0,0)',hovermode='closest')
        fig.update_xaxes(title='umap_x',type='linear')
        fig.update_yaxes(title='umap_y',type='linear')
        return fig

    @app.callback(
        Output('display_table','figure'),
        Input('scatter_figure','selectedData'))
    def display_df(selectedData):
        nonlocal df
        selected_index = []
        for p in selectedData['points']:
            selected_index.append(p['text'])
        df = df.loc[selected_index,:]
        fig = go.Figure(data=[go.Table(header=dict(values=['neoantigen','uid'],fill_color='paleturquoise',align='left'),
                                    cells=dict(values=[df.index,df.uid],fill_color='lavender',align='left'))])
        #fig.update_layout(width=600)
        return fig

    # run app
    app.run_server()







