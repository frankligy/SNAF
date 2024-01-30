import pandas as pd
from sklearn.decomposition import PCA
from umap import UMAP
import numpy as np
import os,sys
from ast import literal_eval
import dash
from dash import dcc,html,dash_table
import plotly.graph_objects as go
import plotly.express as px
from dash.dependencies import Input,Output,State
from .pweblogo import run_pweblogo
from datetime import date,datetime
import subprocess
import atexit


@atexit.register
def clear_assets():
    imgs = os.listdir('assets')
    for img in imgs:
        os.remove(os.path.join('assets',img))


def run_dash_T_antigen(input_abs_path,remove_cols=['uid'],host=None,port='8050',output_abs_path=None):
    '''
    run the dash T antigen viewer

    :param input_abs_path: string, the absolute path to the input df
    :param remove_cols: list, the column name to remove from the input df
    :param host: string or None, if None, program will run hostname to automatically detect
    :param port: string, default is 8050
    :param output_abs_path: string or None, if you want to have the umap embedding, specify the output path and name

    Example::

        snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/shared_vs_unique_neoantigen_all.txt')
    '''
    # since this module heavily relies on relative path
    os.chdir(os.path.dirname(__file__))
    print('changed working directory to {}'.format(os.getcwd()))
    ### df
    df = pd.read_csv(input_abs_path,sep='\t',index_col=0)
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
    if output_abs_path is not None:
        # write the embedding for other usages
        df9['umap_x'] = embed[0][:,0]
        df9['umap_y'] = embed[0][:,1]
        df10['umap_x'] = embed[1][:,0]
        df10['umap_y'] = embed[1][:,1]
        df9.to_csv(os.path.join(output_abs_path,'mer9_umap_embed_df.txt'),sep='\t')
        df10.to_csv(os.path.join(output_abs_path,'mer10_umap_embed_df.txt'),sep='\t')

    # building app
    app = dash.Dash(__name__)
    dropdown_cols = list(df.columns)
    for col in remove_cols:
        dropdown_cols.remove(col)
    app.layout = html.Div([
        html.Div(html.H1('SNAF T-antigen Viewer'),style={'text-align':'center'}),
        html.Div([html.Label('Metadata to display',style={'font-weight':'bold'}),dcc.Dropdown(id = 'metadata_dropdown',options = [{'label':col,'value':col} for col in dropdown_cols],value = 'identity'),
                  html.Br(),
                  html.Label('Length',style={'font-weight':'bold'}),dcc.RadioItems(id='length_radioitem',options=[{'label':9,'value':9},{'label':10,'value':10}],value=9),
                  html.Br(),
                  html.Button(id='submit_button',n_clicks=0,children='Submit')],style={'width':'30%','float':'left'}),
        html.Div([dcc.Graph(id='scatter_figure')],style={'width':'60%','float':'right','margin-top':'100px'}),
        html.Div([html.Br(),html.H2('Selected Weblogo'),html.Img(alt='weblogo',id='display_weblogo',width='95%',height='80%',style={'border-style':'dashed'})],style={'clear':'left','width':'30%'}), 
        html.Div([html.Br(),html.H2('Selected Neoantigen'),dash_table.DataTable(id='display_table',columns=[{'name':column,'id':column} for column in ['neoantigen','uid','mean_percent_samples_junction_present','actual_percent_samples_neoantigen_present','identity','length']],page_size=10)],style={'width':'100%','clear':'both'})     
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
        fig = px.scatter(x=embedding[:,0],y=embedding[:,1],color=plot_df[dropdown_value],text=plot_df.index.values)
        fig.update_traces(mode='markers',customdata=plot_df['uid'].values,hovertemplate='%{text}<br>%{customdata}')
        fig.update_layout(title='Embedded based on physiochemical properties',margin=dict(t=30,l=0,r=0),plot_bgcolor='rgba(0,0,0,0)',hovermode='closest')
        fig.update_xaxes(title='umap_x',type='linear')
        fig.update_yaxes(title='umap_y',type='linear')
        return fig

    @app.callback(
        Output('display_table','data'),
        Input('scatter_figure','selectedData'))
    def display_df(selectedData):
        display_df = df.copy()
        selected_index = []
        for p in selectedData['points']:
            selected_index.append(p['text'])
        display_df = display_df.loc[selected_index,:]
        data_table = []
        for row in display_df.itertuples():
            data_table.append({'neoantigen':row.Index,'uid':row.uid,'mean_percent_samples_junction_present':row.mean_percent_samples_junction_present,'actual_percent_samples_neoantigen_present':row.actual_percent_samples_neoantigen_present,'identity':row.identity,'length':row.length})
        return data_table

    @app.callback(
        Output('display_weblogo','src'),
        Input('scatter_figure','selectedData'))
    def display_weblogo(selectedData):
        selected_index = []
        for p in selectedData['points']:
            selected_index.append(p['text'])
        suffix = str(date.today()) + datetime.now().strftime('%H-%M-%S-%f')
        if not os.path.exists('./assets'):
            os.mkdir('./assets')
        run_pweblogo(selected_index,'./assets/pweblogo_{}.png'.format(suffix))
        print(app.get_asset_url('pweblogo_{}.png'.format(suffix)))
        return app.get_asset_url('pweblogo_{}.png'.format(suffix))

    # run app
    if host is None:
        host = subprocess.run(['hostname'],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0]
    app.run_server(host=host,port=port)







