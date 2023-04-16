#!/data/salomonis-archive/FASTQs/NCI-R01/encode_shrna_kd/K562/inferelator/inferelator_env/bin/python3.6

import pandas as pd
import numpy as np
from inferelator import inferelator_workflow
import os,sys
from inferelator import MPControl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
from sklearn.metrics import precision_recall_curve,roc_curve,auc,confusion_matrix

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


def draw_PR(y_true,y_pred,outdir,name):
    try:
        total = len(y_true)
        n_true = np.count_nonzero(y_true.values) / len(y_true)
        precision,recall,_ = precision_recall_curve(y_true,y_pred,pos_label=1)
    except:
        print('{} either all pos or neg in gs'.format(name))
        return None,None,total
    else:
        area_PR = auc(recall,precision)
        baseline = np.sum(np.array(y_true) == 1) / len(y_true)

        plt.figure()
        lw = 2
        plt.plot(recall,precision, color='darkorange',
                lw=lw, label='PR (area = %0.2f)' % area_PR)
        plt.plot([0, 1], [baseline, baseline], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 0.2])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('PR curve {}'.format(name))
        plt.legend(loc="lower right")
        plt.savefig(os.path.join(outdir,'{}_aupr.pdf'.format(name)),bbox_inches='tight')
        plt.close()
    return area_PR,n_true,total



if __name__ == '__main__':
    # set up multiprocessing
    MPControl.set_multiprocess_engine("multiprocessing")
    MPControl.client.processes = 30
    MPControl.connect()   
    # run bbsr
    outdir = '../altanalyze_output/HepG2_inference/bbsr_weight1.5'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    worker = inferelator_workflow(regression="bbsr", workflow="tfa")
    worker.set_file_paths(input_dir='../altanalyze_output/HepG2_inference',
                        output_dir=outdir,
                        expression_matrix_file="expr.tsv",   # sample * events
                        tf_names_file="sf_name.tsv",
                        priors_file="prior.tsv",)
#                        gold_standard_file="gs.tsv")
    worker.set_run_parameters(num_bootstraps=50, random_seed=42)
    worker.set_regression_parameters(prior_weight=1.5)
    worker.set_network_data_flags(use_no_gold_standard=True)
    network_result = worker.run()

    # # augment
    # expr = pd.read_csv('../altanalyze_output/leucegene_inference/expr.tsv',sep='\t',index_col=0)
    # network = pd.read_csv('../altanalyze_output/leucegene_inference/bbsr_weight1/network.tsv',sep='\t')
    # network = network.groupby(by=['target','regulator'])['combined_confidences'].mean().unstack(fill_value=0)
    # missing_events = list(set(expr.columns).difference(set(network.index)))
    # supp_net = pd.DataFrame(data=np.full((len(missing_events),network.shape[1]),0),index=missing_events,columns=network.columns)
    # network = pd.concat([network,supp_net],axis=0)
    # network.to_csv('../altanalyze_output/leucegene_inference/bbsr_weight1/network_aug.txt',sep='\t')

    # # TFA
    # ## 1. derive sfa
    # network = pd.read_csv('../altanalyze_output/leucegene_inference/prior.tsv',sep='\t',index_col=0)
    # expr = pd.read_csv('../altanalyze_output/leucegene_inference/expr.tsv',sep='\t',index_col=0).T
    # missing_events = list(set(expr.index).difference(set(network.index)))
    # common_events = list(set(expr.index).intersection(set(network.index)))
    # supp_net = pd.DataFrame(data=np.full((len(missing_events),network.shape[1]),0),index=missing_events,columns=network.columns)
    # network = pd.concat([network.loc[network.index.isin(common_events),:],supp_net],axis=0).loc[expr.index,:]
    # print(network)
    # print(expr)
    # sfa = np.linalg.pinv(network.values)@expr.values
    # sfa = pd.DataFrame(data=sfa,columns=expr.columns,index=network.columns)
    # ## 2. attach meta data
    # # meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
    # # meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
    # # meta = pd.concat([meta_k,meta_c])
    # # meta.columns = [item.replace(' ','_') for item in meta.columns]
    # # meta['Experiment_target'].fillna(value='control',inplace=True)
    # # meta['Experiment_target'] = [item.split('-')[0] if item != 'control' else item for item in meta['Experiment_target']]
    # # pair1_to_target = meta['Experiment_target'].to_dict()
    # # mi_array = [tuple([pair1_to_target[item.split('_')[0]] for item in sfa.columns]),tuple(sfa.columns.tolist())]
    # # mi = pd.MultiIndex.from_arrays(arrays=mi_array,names=('target','name'))
    # # sfa.columns = mi
    # meta = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/AML/Leucegene/groups.expanded-U2AF1-CV-all-patients.txt',sep='\t',header=None,index_col=0)
    # other = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/AML/Leucegene/other.csv',index_col=0)
    # srs_to_meta = meta[1].to_dict()
    # srs_to_other = [other[item].to_dict() for item in other.columns]
    # srs_to_all = [srs_to_meta] + srs_to_other
    # mi_array = []
    # for dic in srs_to_all:
    #     mi_array.append(tuple([dic[item.split('_')[0]] for item in sfa.columns]))
    # mi_array.append(tuple(sfa.columns.tolist()))
    # mi = pd.MultiIndex.from_arrays(arrays=mi_array,names= ['group'] + other.columns.tolist() + ['name'])
    # sfa.columns = mi
    # sfa.to_csv('../altanalyze_output/leucegene_inference/bbsr_weight1/sfa_prior.txt',sep='\t')




    # aupr
    # network = pd.read_csv('../altanalyze_output/inference/bbsr_weight1/network_aug.txt',sep='\t',index_col=0)
    # prior = pd.read_csv('../altanalyze_output/prior/whole_prior/shRNA_K562_rule.txt',sep='\t',index_col=0)
    # gs = pd.read_csv('../altanalyze_output/prior/whole_prior/shRNA_HepG2_rule.txt',sep='\t',index_col=0)
    # all_sfs = gs.columns
    # missing_events = list(set(network.index).difference(set(gs.index)))
    # supp = pd.DataFrame(data=np.full((len(missing_events),gs.shape[1]),0),index=missing_events,columns=gs.columns)
    # gs = pd.concat([gs,supp],axis=0)

    # network = network.stack().to_frame(name='weight')
    # network.index = [event + ',' + sf for event,sf in network.index.tolist()]
    # prior = prior.stack().to_frame(name='K562_prior')
    # prior.index = [event + ',' + sf for event,sf in prior.index.tolist()]
    # gs = gs.stack().to_frame(name='HepG2_prior')
    # gs.index = [event + ',' + sf for event,sf in gs.index.tolist()]
    # merged = network.join(gs,how='inner').join(prior,how='left')
    # merged['gs'] = [1 if item > 0 else 0 for item in merged['HepG2_prior']]
    # merged.sort_values(by='weight',ascending=False,inplace=True)
    # merged['sf'] = [item.split(',')[1] for item in merged.index]
    # merged.to_csv('../altanalyze_output/inference/bbsr_weight1/aupr_eval/correspond.txt',sep='\t')

    # # using network
    # draw_PR(merged['gs'],merged['weight'],outdir='../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_network',name='all')
    # aupr_dic = {}
    # for sf in tqdm(all_sfs):
    #     merged_sub = merged.loc[merged['sf']==sf,:]
    #     aupr,n_true,total = draw_PR(merged_sub['gs'],merged_sub['weight'],outdir='../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_network',name=sf)
    #     aupr_dic[sf] = [aupr,n_true,total]
    # pd.DataFrame(data=aupr_dic,index=['aupr','n_true','total']).T.to_csv('../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_network/aupr_per_sf.txt',sep='\t')

    # just using K562 prior
    # merged['K562_prior'].fillna(value=0,inplace=True)
    # draw_PR(merged['gs'],merged['K562_prior'],outdir='../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_prior',name='all')
    # aupr_dic = {}
    # for sf in tqdm(all_sfs):
    #     merged_sub = merged.loc[merged['sf']==sf,:]
    #     aupr,n_true,total = draw_PR(merged_sub['gs'],merged_sub['K562_prior'],outdir='../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_prior',name=sf)
    #     aupr_dic[sf] = [aupr,n_true,total]
    # pd.DataFrame(data=aupr_dic,index=['aupr','n_true','total']).T.to_csv('../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_prior/aupr_per_sf.txt',sep='\t')

    # comparison between two modes
    # b = pd.read_csv('../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_network/aupr_per_sf.txt',sep='\t',index_col=0)
    # p = pd.read_csv('../altanalyze_output/inference/bbsr_weight1/aupr_eval/aupr_prior/aupr_per_sf.txt',sep='\t',index_col=0)
    # bp = b.join(p,lsuffix='_b',rsuffix='_p',how='inner')
    # bp = bp.loc[bp['aupr_b'].notna(),:].loc[bp['aupr_p'].notna(),:].sort_values(by='aupr_b')
    # print(bp)
    # import plotly.graph_objects as go
    # node_b_x = []
    # node_b_y = []
    # node_p_x = []
    # node_p_y = []
    # node_text = []
    # for i,row in enumerate(bp.itertuples()):
    #     node_b_x.append(i)
    #     node_b_y.append(row.aupr_b)
    #     node_p_x.append(i)
    #     node_p_y.append(row.aupr_p)
    #     node_text.append(row.Index)
    # node_b_trace = go.Scatter(x=node_b_x,y=node_b_y,mode='markers',hoverinfo='text',text=node_text,marker={'color':'green','size':4})
    # node_p_trace = go.Scatter(x=node_p_x,y=node_p_y,mode='markers',hoverinfo='text',text=node_text,marker={'color':'red','size':4})
    # fig_layout = go.Layout(showlegend=False,title='comparison')
    # fig = go.Figure(data=[node_b_trace,node_p_trace],layout=fig_layout)
    # fig.write_html('../altanalyze_output/inference/bbsr_weight1/aupr_eval/comparison_stack.html',include_plotlyjs='cdn')

    # node_x = []
    # node_y = []
    # node_text = []
    # for i,row in enumerate(bp.itertuples()):
    #     node_x.append(row.aupr_p)
    #     node_y.append(row.aupr_b)
    #     node_text.append(row.Index)
    # node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',hoverinfo='text',text=node_text,marker={'color':'black','size':4})
    # edge_x = [0,1]
    # edge_y = [0,1]
    # edge_trace = go.Scatter(x=edge_x,y=edge_y,mode='lines',line={'width':1,'color':'black','dash':'dash'})
    # fig_layout = go.Layout(showlegend=False,title='comparison',xaxis=dict(title_text='only prior (aupr)',range=[0,1]),yaxis=dict(title_text='using network (aupr)',range=[0,1]),margin=dict(l=300,r=300))
    # fig = go.Figure(data=[node_trace,edge_trace],layout=fig_layout)
    # fig.write_html('../altanalyze_output/inference/bbsr_weight1/aupr_eval/comparison_diagonal.html',include_plotlyjs='cdn')




















