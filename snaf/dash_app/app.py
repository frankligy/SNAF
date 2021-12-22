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
from .pweblogo import run_pweblogo
from datetime import date,datetime
import atexit


@atexit.register
def clear_assets():
    imgs = os.listdir('assets')
    for img in imgs:
        os.remove(os.path.join('assets',img))


def run_dash_app(intpath,remove_cols=None,host='127.0.0.1',port='8050'):
    # since this module heavily relies on relative path
    os.chdir(os.path.dirname(__file__))
    print('changed working directory to {}'.format(os.getcwd()))
    ### df
    df = pd.read_csv(intpath,sep='\t',index_col=0)
    df['length'] = [len(item) for item in df.index]
    ### umap
    after_pca = np.loadtxt(os.path.join('..','deepimmuno','data','after_pca.txt'))
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
    print('computing embedding, it may take a while')
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
    print('embedding ready, app will be built soon')

    # building app

    app = dash.Dash(__name__)
    dropdown_cols = list(df.columns)
    for col in remove_cols:
        dropdown_cols.remove(col)
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
            dcc.Graph(id='display_table')]),
        # weblogo
        html.Div([html.H2('Selected Weblogo'),
                  html.Img(alt='weblogo',id='display_weblogo')])
        

    ])

    # app callback
    @app.callback(
        Output('scatter_figure','figure'),
        State('metadata_dropdown','value'),
        State('length_radioitem','value'),
        Input('submit_button','n_clicks'))
    def scatter_figure(dropdown_value,length_value,n_clicks):
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
        display_df = df.copy()
        selected_index = []
        for p in selectedData['points']:
            selected_index.append(p['text'])
        display_df = display_df.loc[selected_index,:]
        fig = go.Figure(data=[go.Table(header=dict(values=['neoantigen','uid'],fill_color='paleturquoise',align='left'),
                                    cells=dict(values=[display_df.index,display_df.uid],fill_color='lavender',align='left'))])
        #fig.update_layout(width=600)
        return fig

    @app.callback(
        Output('display_weblogo','src'),
        Input('scatter_figure','selectedData'))
    def display_weblogo(selectedData):
        selected_index = []
        for p in selectedData['points']:
            selected_index.append(p['text'])
        suffix = str(date.today()) + datetime.now().strftime('%H-%M-%S-%f')
        run_pweblogo(selected_index,'./assets/pweblogo_{}.png'.format(suffix))
        return app.get_asset_url('pweblogo_{}.png'.format(suffix))

    # run app
    app.run_server(host=host,port=port)







