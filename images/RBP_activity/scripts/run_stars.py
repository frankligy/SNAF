#!/data/miraldiLab/team/Frank/SRN/stars_alasso/glmnet_python_env/bin/python3.6

import os,sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from scipy.special import comb
import pingouin as pg
import math
import glmnet_python
from glmnet import glmnet
import multiprocessing as mp
import pickle
import copy
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from sklearn.metrics import precision_recall_curve,roc_curve,auc,confusion_matrix

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


def plot_total_instability(instability_for_plot,lambdau,events_to_draw=None,save=True,outdir='./',name='total_instability_plot'):
    # instabiligy_for_plot is a dict {event:[0.0,0.03,...0.0], the value is a list, of length=len(lambdau), order is descending order of lambda}
    if events_to_draw is not None:
        fig,ax = plt.subplots() 
        ax.set_xscale('log')
        for event,total_instability_list in instability_for_plot.items():
            if event in events_to_draw:
                ax.plot(lambdau,total_instability_list,label=event)
    else:
        ave_instability = np.array(list(instability_for_plot.values())).mean(axis=0)
        fig,ax = plt.subplots()
        ax.set_xscale('log')
        ax.plot(lambdau,ave_instability)
    if save:
        plt.savefig(os.path.join(outdir,'{}.pdf'.format(name)),bbox_inches='tight')

def get_sfa_and_full_prior(expr,prior):
    n_event_expr, n_sample = expr.shape
    n_event_prior, n_sf = prior.shape
    common_events = list(set(expr.index).intersection(set(prior.index)))
    missing_events = list(set(expr.index).difference(set(prior.index)))  
    enlarge_prior = pd.DataFrame(data=np.zeros((len(missing_events),n_sf)),columns=prior.columns,index=missing_events) 
    full_prior = pd.concat([prior.loc[common_events,:],enlarge_prior],axis=0).loc[expr.index,:]  # n_event_expr (or n_event), n_sf
    sfa_data = np.linalg.pinv(full_prior.values)@expr.values
    sfa = pd.DataFrame(data=sfa_data,columns=expr.columns,index=full_prior.columns) # n_sf * n_sample
    return sfa, full_prior


def subsampling(sfa,expr,b_subsample=None,n_subsample=None):
    n_event, n_sample = expr.shape
    n_sf, n_sample = sfa.shape
    combine = pd.concat([expr,sfa],axis=0)  # n_event+n_sf, n_sample
    if b_subsample is None:
        b_subsample = math.floor(0.63 * n_sample)
    if n_subsample is None:
        n_subsample = 50
    subsamples = []
    for i in range(n_subsample):
        subsample_indices = np.random.choice(np.arange(combine.shape[1]),size=b_subsample,replace=False)
        subsamples.append(combine.iloc[:,subsample_indices])
    stack = np.empty([n_subsample,n_event+n_sf,b_subsample])
    for i,subsample in enumerate(subsamples):
        stack[i,:,:] = subsamples[i].values
    return stack, (n_event,n_sf)

def monotonize(total_instability_list):
    # D bar is the supremum of D hat, see either StARs paper or bStARs paper
    # also, since lambdau is descending order, the total_instability_list is acending order
    max_value = max(total_instability_list)
    max_value_index = total_instability_list.index(max_value)
    new_list = copy.deepcopy(total_instability_list)
    new_list[max_value_index:] = (len(total_instability_list) - max_value_index) * [max_value]  
    return new_list,max_value



def single_task_regression(event_slice,sfs_slice,event_prior,bias=1.0,base=1,lambdau=None,step=10,cutoff=0.05,max_iter=20,tol=1e-4,sign=False):
    # event_slice: n_subsample * 1 * b_subsample
    # sfs_slice: n_subsample * n_sf * b_subsample
    # event_prior: n_sf (1d)
    n_subsample, _ , b_subsample = event_slice.shape
    n_subsample, n_sf, b_subsample = sfs_slice.shape
    if lambdau is None:
        lambdau = np.logspace(-3,1,200)
        lambdau = np.flip(np.sort(lambdau))  # make sure it is in descending order
    else:
        lambdau = np.flip(np.sort(lambdau))  # make sure it is in descending order
    event_prior_weighted = np.where(event_prior==0,base,bias*base)
    store = dict()  # storing beta values for each lambda (key is the index of lambda in lambdau), and for each lambda, will be a stacked array n_subsample * n_sf
    for i,item in enumerate(lambdau):
        store[i] = np.empty(shape=(n_subsample,n_sf))
    for i in range(n_subsample):
        Y = event_slice[i,:,:].T  # b_subsample * 1
        X = sfs_slice[i,:,:].T   # b_subsample * n_sf
        fit = glmnet(x=X,y=Y,alpha=1,lambdau=lambdau,penalty_factor=event_prior_weighted)
        path = fit['lambdau']  # large to small, 1d array
        beta = fit['beta']   # n_sf * n_lambdau
        for j,item in enumerate(path):
            store[j][i,:] = beta[:,j]
    total_instability_list = []  # length of n_lambdau
    v_list = []  # length of n_lambdau, each item means a n_subsample * n_sf array
    for k,v in store.items():
        v_ = np.where(v==0,0,1)
        p = v_.mean(axis=0)
        variance = p*(1-p)
        instability = 2*variance
        total_instability = instability.mean()
        total_instability_list.append(total_instability)
        v_list.append(v)
    # monotonize the total_instability_list
    total_instability_list,max_instability = monotonize(total_instability_list)
    # begin to search for supremum {total_instability(lambdau)<0.05}, start from the largest lambda, then descrease, should be parabolic curve
    if max_instability < cutoff:  # the peak of parabolic curve is still below cutoff, supremum is just the x-coord under the peak
        current_index = total_instability_list.index(max_instability)
    else: # need to use bisection method to find the supremum
        current_index = 0
        start_total_instability = total_instability_list[current_index]
        if start_total_instability >= cutoff:
            raise Exception('considering increasing the lambdau upper limit')
        current_total_instability = start_total_instability
        i = 0
        while i < max_iter and cutoff - current_total_instability > tol:
            next_index = current_index+step
            if next_index < len(total_instability_list):
                next_total_instability = total_instability_list[next_index]
                if next_total_instability >= cutoff:
                    if step > 1:
                        step = math.ceil(step/2)
                    else:
                        break
                else:
                    current_index = next_index  
                    current_total_instability = next_total_instability  
                    i += 1
            else:
                if step > 1:
                    step = math.ceil(step/2)
                else:
                    break
    final_lambda = lambdau[current_index]
    final_total_instability = total_instability_list[current_index]
    final_v = v_list[current_index]  # n_subsample * n_sf
    # compute confidence for the edges (event->sfs)
    if not sign:
        nonzero_time = np.count_nonzero(final_v,axis=0)
    else:
        positive_time = np.sum(final_v>0,axis=0)
        negative_time = np.sum(final_v<0,axis=0)
        zero_time = np.sum(final_v==0,axis=0)
        nonzero_time = positive_time + negative_time
    return final_lambda,final_total_instability,nonzero_time,total_instability_list


def wrapper_single_task_regression(event_chunk,sfs_slice,full_prior_chunk):
    print('{} events are running on {} for regression'.format(len(event_chunk),os.getpid()))
    network = np.empty([len(event_chunk),sfs_slice.shape[1]])
    event_total_instability_dict = {}
    for i,event in tqdm(enumerate(event_chunk),total=len(event_chunk)):
        event_slice = stack[:,i,:][:,np.newaxis,:]   # n_subsample * 1 * b_subsample
        event_prior = full_prior_chunk.loc[event,:].values # n_sf
        final_lambda, final_total_instability, nonzero_time, total_instability_list = single_task_regression(event_slice,sfs_slice,event_prior)
        network[i,:] = nonzero_time
        event_total_instability_dict[event] = total_instability_list
    return network, event_total_instability_dict

def wrapper_partial_correlation(event_chunk,sfs,df_chunk):
    print('{} events are running on {} for partial correlation'.format(len(event_chunk),os.getpid()),datetime.now())
    pc_mat = np.empty(shape=(len(event_chunk),len(sfs)))
    for i,event in enumerate(event_chunk):
        print('starting event{} on {}'.format(i,os.getpid()),datetime.now())
        for j,sf in enumerate(sfs):
            other_sfs = [sfs[k] for k in range(n_sf) if k != j]
            pc = pg.partial_corr(data=df_chunk,x=sf,y=event,y_covar=other_sfs).loc['pearson','r']
            pc_mat[i,j] = pc
    return pc_mat


if __name__ == '__main__':
    # ##### running the stars-lasso algorithm ######
    # expr = pd.read_csv('../altanalyze_output/HepG2_inference/expr.tsv',sep='\t',index_col=0)   # sample * event
    # prior = pd.read_csv('../altanalyze_output/HepG2_inference/prior.tsv',sep='\t',index_col=0)
    # outdir = '../altanalyze_output/HepG2_inference/stars_b1.0'

    # if not os.path.exists(outdir):
    #     os.mkdir(outdir)

    # # other parameters
    # cores = 50
    
    # # step0: get the dimension and some preprocessing
    # expr = expr.T
    # n_event, n_sample = expr.shape
    # _, n_sf = prior.shape

    # # step1: get sfa and full_prior
    # sfa, full_prior = get_sfa_and_full_prior(expr,prior)

    # # step2: subsamping
    # stack, event_sf_ratio = subsampling(sfa,expr)

    # # step3: regression
    # events = expr.index.values
    # sfs = sfa.index.values
    # sfs_slice = stack[:,event_sf_ratio[0]:,:]  # n_subsample * n_sf * b_subsample
    # print('{} cores deployable for regression'.format(cores))
    # event_chunks = np.array_split(events,cores)
    # full_prior_chunks = [full_prior.loc[event_chunk,:] for event_chunk in event_chunks]
    # pool = mp.Pool(processes=cores)
    # r = [pool.apply_async(func=wrapper_single_task_regression,args=(event_chunk,sfs_slice,full_prior_chunk,)) for event_chunk,full_prior_chunk in zip(event_chunks,full_prior_chunks)]
    # pool.close()
    # pool.join()
    # networks = []
    # instability_for_plot = {}
    # for collect in r:
    #     result = collect.get()
    #     networks.append(result[0])  # each result contains two, first is (n_event_chunks,n_sf), and the order is kept, second is event_total_instability_dict
    #     instability_for_plot.update(result[1])
    # network = np.concatenate(networks,axis=0)  # n_events, n_sf


    # # step4: print out and plot
    # final_network = network 
    # final_df = pd.DataFrame(data=final_network,columns=sfs,index=events)
    # final_df.to_csv(os.path.join(outdir,'./stars_alasso_network.txt'),sep='\t')
    # with open(os.path.join(outdir,'./instability_for_plot.p'),'wb') as f:
    #     pickle.dump(instability_for_plot,f)
    # plot_total_instability(instability_for_plot,lambdau=np.logspace(1,-3,200),events_to_draw=None,outdir=outdir,name='instability_all')
    # plot_total_instability(instability_for_plot,lambdau=np.logspace(1,-3,200),events_to_draw=np.random.choice(events,30,replace=False),outdir=outdir,name='instability_random')

    # ##### aupr #########
    # network = pd.read_csv('../altanalyze_output/HepG2_inference/stars_b1.0/stars_alasso_network.txt',sep='\t',index_col=0)
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
    # merged.to_csv('../altanalyze_output/HepG2_inference/stars_b1.0/aupr_eval/correspond.txt',sep='\t')

    # # using network
    # draw_PR(merged['gs'],merged['weight'],outdir='../altanalyze_output/HepG2_inference/stars_b1.0/aupr_eval/aupr_network',name='all')
    # aupr_dic = {}
    # for sf in tqdm(all_sfs):
    #     merged_sub = merged.loc[merged['sf']==sf,:]
    #     aupr,n_true,total = draw_PR(merged_sub['gs'],merged_sub['weight'],outdir='../altanalyze_output/HepG2_inference/stars_b1.0/aupr_eval/aupr_network',name=sf)
    #     aupr_dic[sf] = [aupr,n_true,total]
    # pd.DataFrame(data=aupr_dic,index=['aupr','n_true','total']).T.to_csv('../altanalyze_output/HepG2_inference/stars_b1.0/aupr_eval/aupr_network/aupr_per_sf.txt',sep='\t')

    # # just using K562 prior
    # merged['K562_prior'].fillna(value=0,inplace=True)
    # draw_PR(merged['gs'],merged['K562_prior'],outdir='../altanalyze_output/HepG2_inference/stars_b0.5/aupr_eval/aupr_prior',name='all')
    # aupr_dic = {}
    # for sf in tqdm(all_sfs):
    #     merged_sub = merged.loc[merged['sf']==sf,:]
    #     aupr,n_true,total = draw_PR(merged_sub['gs'],merged_sub['K562_prior'],outdir='../altanalyze_output/HepG2_inference/stars_b0.5/aupr_eval/aupr_prior',name=sf)
    #     aupr_dic[sf] = [aupr,n_true,total]
    # pd.DataFrame(data=aupr_dic,index=['aupr','n_true','total']).T.to_csv('../altanalyze_output/HepG2_inference/stars_b0.5/aupr_eval/aupr_prior/aupr_per_sf.txt',sep='\t')

    # comparison between two modes
    b = pd.read_csv('../altanalyze_output/HepG2_inference/stars_b1.0/aupr_eval/aupr_network/aupr_per_sf.txt',sep='\t',index_col=0)
    p = pd.read_csv('../altanalyze_output/HepG2_inference/stars_b0.5/aupr_eval/aupr_network/aupr_per_sf.txt',sep='\t',index_col=0)
    bp = b.join(p,lsuffix='_b',rsuffix='_p',how='inner')
    bp = bp.loc[bp['aupr_b'].notna(),:].loc[bp['aupr_p'].notna(),:].sort_values(by='aupr_b')
    print(bp)
    import plotly.graph_objects as go
    node_x = []
    node_y = []
    node_text = []
    for i,row in enumerate(bp.itertuples()):
        node_x.append(row.aupr_p)
        node_y.append(row.aupr_b)
        node_text.append(row.Index)
    node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',hoverinfo='text',text=node_text,marker={'color':'black','size':4})
    edge_x = [0,1]
    edge_y = [0,1]
    edge_trace = go.Scatter(x=edge_x,y=edge_y,mode='lines',line={'width':1,'color':'black','dash':'dash'})
    fig_layout = go.Layout(showlegend=False,title='comparison',xaxis=dict(title_text='with prior stars (aupr)',range=[0,1]),yaxis=dict(title_text='no prior stars (aupr)',range=[0,1]),margin=dict(l=300,r=300))
    fig = go.Figure(data=[node_trace,edge_trace],layout=fig_layout)
    fig.write_html('../altanalyze_output/HepG2_inference/stars_b1.0/aupr_eval/comparison_diagonal_cross_with0.5_without1.0_prior.html',include_plotlyjs='cdn')


    # #### TFA #######
    # ## 1. derive sfa
    # network = pd.read_csv('../altanalyze_output/leucegene_inference/stars_b0.5/stars_alasso_network.txt',sep='\t',index_col=0)
    # expr = pd.read_csv('../altanalyze_output/leucegene_inference/expr.tsv',sep='\t',index_col=0).T
    # missing_events = list(set(expr.index).difference(set(network.index)))
    # common_events = list(set(expr.index).intersection(set(network.index)))
    # supp_net = pd.DataFrame(data=np.full((len(missing_events),network.shape[1]),0),index=missing_events,columns=network.columns)
    # network = pd.concat([network.loc[network.index.isin(common_events),:],supp_net],axis=0).loc[expr.index,:]
    # print(network)
    # print(expr)
    # sfa = np.linalg.pinv(network.values)@expr.values
    # sfa = pd.DataFrame(data=sfa,columns=expr.columns,index=network.columns)
    # # 2. attach meta data
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
    # sfa.to_csv('../altanalyze_output/leucegene_inference/stars_b0.5/sfa_network.txt',sep='\t')


















