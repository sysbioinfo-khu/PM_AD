'''
This script processes spatial transcriptomics data using the cell2location model. It includes functions to prepare visium and reference data, train a regression model, and classify cell types in spatial transcriptomics data.
Functions:
    prepare_visum_data(filename):
        Prepares visium data by reading the input file, copying raw counts, identifying and removing mitochondrial genes, and returning the processed AnnData object.
    prepare_reference_data(filename):
        Prepares reference data by reading the input file, copying raw counts, filtering genes based on specified cutoffs, and returning the processed AnnData object.
Variables:
    visium_mode (bool): Flag to indicate whether to load a pre-trained model or train a new one.
    results_folder (str): Path to the folder where results will be saved.
    ref_run_name (str): Path and name for the reference regression model results.
    run_name (str): Path and name for the cell2location model results.
Main Workflow:
    - If visium_mode is True, loads a pre-trained regression model.
    - If visium_mode is False, prepares reference data, trains a new regression model, and saves the model and results.
    - Processes a list of healthy, lesion, and non-lesion files using the cell2location model.
    - For each file, prepares visium data, finds shared genes, subsets data, trains the cell2location model, and saves the results.
'''

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tqdm as notebook_tqdm
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

# helper function to prepare visum data

def prepare_visum_data(filename):
    
    adata_vis = sc.read_h5ad(filename)
    adata_vis.layers["imputed"] = adata_vis.X.copy()
    adata_vis.X = adata_vis.layers["raw_count"].copy()
    adata_vis.var[adata_vis.var['highly_variable']]

    # Plotting
    #sc.pl.spatial(adata_vis, color='FCER1G', gene_symbols='gene_ids')
    
    # are there mitochondrial genes in the data?
    mito_genes = adata_vis.var_names.str.startswith('MT-')                           
    mito_genes.sum()
    # find mitochondria-encoded (MT) genes
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['gene_ids']]

    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
    
    return adata_vis


def prepare_reference_data(filename):
    # Load reference data
    adata_ref = sc.read_h5ad(filename)
    
    # Get raw counts
    adata_ref.layers['imputed'] = adata_ref.X.copy()
    adata_ref.X = adata_ref.layers['counts']

    # Filter genes 
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata_ref = adata_ref[:, selected].copy()
    return adata_ref

visium_mode = False

results_folder = '/home/yang/Spacial_Transcriptomics/Cell2Location_analysis'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures_AD_merged_healthy'           # adjust if wanted
run_name = f'{results_folder}/cell2location_map'                                    # same here

if visium_mode: # was not used in this project, because of errors
    print('Visium mode: loading ref model')
    # Loading trained and saved model
    adata_file = f"{ref_run_name}/sc.h5ad"
    print(f'\nReference File  {adata_file} loading\n')
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
    print(f'{adata_file} loaded')
    
else:
    print('scRNA-seq mode: training ref model')
    ref_file = '/home/yang/Spacial_Transcriptomics/AD_scRNA_seq_merged.h5ad'
    adata_ref = prepare_reference_data(ref_file)
        
    # get excema + lesion data
    lesion_ref = adata_ref[adata_ref.obs['Status'] == 'Eczema']
    lesion_ref = lesion_ref[lesion_ref.obs['Site'] == 'lesion']

    # get healthy skin data
    healthy_ref = adata_ref[adata_ref.obs['Status'] == 'Healthy']

    # get AD non-lesion data
    nonlesion_ref = adata_ref[adata_ref.obs['Status'] == 'Eczema']
    nonlesion_ref = nonlesion_ref[nonlesion_ref.obs['Site'] == 'non_lesion']
        
    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=healthy_ref,       # TODO: select reference data
                            # 10X reaction / sample / batch
                            batch_key='donor_id',
                            # cell type, covariate used for constructing signatures
                            labels_key='AD_clustering',
                            # multiplicative technical effects (platform, 3' vs 5', donor effect)
                            categorical_covariate_keys=None
                        )

    # create the regression model
    mod = RegressionModel(healthy_ref)                                          # TODO: select reference data

    # view anndata_setup as a sanity check
    mod.view_anndata_setup()

    #  Trining sc model
    mod.train(max_epochs=250)#, use_gpu=True)

    mod.plot_history(20)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}#, 'use_gpu': True}
    )

    # Save model
    mod.save(f"{ref_run_name}", overwrite=True)

    # Save anndata object with results
    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref.write(adata_file)
    adata_file

    adata_ref = mod.export_posterior(
        adata_ref, use_quantiles=True,
        # choose quantiles
        add_to_varm=["q05","q50", "q95", "q0001"],
        sample_kwargs={'batch_size': 2500}#, 'use_gpu': True}
    )

""" # Loading trained and saved model

adata_file = f"{ref_run_name}/lesion_sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)  """

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']




# ------------------------------------- Cell2location model classification loop -------------------------------------
lesion_files = [#'/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907082_AD_6_LS_matrix_processed.h5ad',
                '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907077_AD_1_LS_matrix_processed.h5ad',
                '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907078_AD_2_LS_matrix_processed.h5ad',
                #'/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907079_AD_3_LS_matrix_processed.h5ad',
                '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907080_AD_4_LS_matrix_processed.h5ad',
                '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907081_AD_5_LS_matrix_processed.h5ad',
                '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907083_AD_7_LS_matrix_processed.h5ad'
                ]

nonlesion_files = [#'/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907089_AD_6_NL_matrix_processed.h5ad',
                   '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907084_AD_1_NL_matrix_processed.h5ad',
                   '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907085_AD_2_NL_matrix_processed.h5ad',
                   '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907086_AD_3_NL_matrix_processed.h5ad',
                   '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907088_AD_5_NL_matrix_processed.h5ad',
                   '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907090_AD_7_NL_matrix_processed.h5ad'
                   ]

healthy_files = [#'/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907095_HE_5_matrix_processed.h5ad',
                 '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907091_HE_1_matrix_processed.h5ad',
                 '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907092_HE_2_matrix_processed.h5ad',
                 '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907093_HE_3_matrix_processed.h5ad',
                 '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907094_HE_4_matrix_processed.h5ad',
                 '/home/yang/Spacial_Transcriptomics/downloaded_files/GSM5907096_HE_6_matrix_processed.h5ad'
                ]
for file in healthy_files:  # TODO: change to lesion_files or nonlesion_files
        
    # Getting the visum data
    vis_file = file
    sample = vis_file.split('/')[-1].split('_matrix')[0]
    print(f'Processing {sample}')
    adata_vis = prepare_visum_data(vis_file)

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    
    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=30,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )
    mod.view_anndata_setup()

    mod.train(max_epochs=30000,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1#,
            #use_gpu=True,
            )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}#, 'use_gpu': True}
    )

    # Save model
    mod.save(f"{run_name}_{sample}", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # Save anndata object with results
    adata_file = f"{run_name}_{sample}/_sp.h5ad"
    adata_vis.write(adata_file)
    adata_file
    
"""
    # saving model
    adata_file = f"{run_name}/{sample}_sp.h5ad"
    adata_vis = sc.read_h5ad(adata_file)
    mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

    # select one slide
 
    slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

    # plot in spatial coordinates
    with mpl.rc_context({'axes.facecolor':  'black',
                        'figure.figsize': [4.5, 5]}):

        sc.pl.spatial(slide, cmap='magma',
                    # show first 8 cell types
                    color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC',
                            'B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],
                    ncols=4, size=1.3,
                    img_key='hires',
                    # limit color scale at 99.2% quantile of cell abundance
                    vmin=0, vmax='p99.2'
                    ) """
                    