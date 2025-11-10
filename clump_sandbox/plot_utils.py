import sys
sys.dont_write_bytecode = True
import pandas as pd
import numpy as np

from misc.misc_utils import ExtractParams
#====================================================================
def PlotSimul(ax, 
              file_path,
              band_type,
              replacement_val,
              x_name='GFP genomes (scaled)',
              y_name='GFP IU', 
              color='green',
              marker='^',
              s=20,
              scatter=True,
              return_gen_well=False):
    
    df_simul = pd.read_csv(file_path)
    PARAM_DICT = ExtractParams(df_simul)

    num_simulations = 0
    for key in df_simul.keys():
        if ("IU run" in key):
            num_simulations += 1

    cell_count     = PARAM_DICT['cell_count']

    GEN_WELL_SIMUL = df_simul.loc[:, x_name].values
    GEN_CELL_SIMUL = GEN_WELL_SIMUL / cell_count
    INF_WELL_SIMUL = np.empty((num_simulations, df_simul.shape[0]), float)

    for i in range(num_simulations):
        INF_WELL_SIMUL[i, :] = df_simul.loc[:, y_name+' run='+str(i)]

    INF_CELL_SIMUL_MIN  = np.min(INF_WELL_SIMUL, axis=0) / cell_count
    INF_CELL_SIMUL_MAX  = np.max(INF_WELL_SIMUL, axis=0) / cell_count
    INF_CELL_SIMUL_MEAN = np.mean(INF_WELL_SIMUL, axis=0) / cell_count

    print(" >>", np.shape(INF_WELL_SIMUL), np.shape(INF_CELL_SIMUL_MEAN))

    #--------------------------------------------------------------------
    for i in range(np.shape(INF_CELL_SIMUL_MIN)[0]):
        if (INF_CELL_SIMUL_MIN[i] < replacement_val):
            INF_CELL_SIMUL_MIN[i] = replacement_val
        
        if (INF_CELL_SIMUL_MEAN[i] < replacement_val):
            INF_CELL_SIMUL_MEAN[i] = replacement_val
    #--------------------------------------------------------------------
    if (num_simulations == 1):
        ax.scatter(GEN_CELL_SIMUL, INF_WELL_SIMUL / cell_count, s=s, facecolors='none', edgecolors=color, marker=marker)

    elif (num_simulations > 1):
        ax.fill_between(GEN_CELL_SIMUL, INF_CELL_SIMUL_MIN, INF_CELL_SIMUL_MAX, color='black', alpha=.3)
        
        ax.plot(GEN_CELL_SIMUL, INF_CELL_SIMUL_MEAN, color='black', linestyle='-', linewidth=.5)
    
    if (scatter == True):
        ax.scatter(GEN_CELL_SIMUL, INF_CELL_SIMUL_MEAN, s=s, color=color, marker=marker)
    #--------------------------------------------------------------------
    if not return_gen_well:
        return GEN_CELL_SIMUL, INF_CELL_SIMUL_MEAN
    else:
        return GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL_MEAN

#====================================================================
def PlotText(ax, PARAM_DICT, xMin, xMax, yMin, yMax):
    remove_str = ''
    if (PARAM_DICT['remove'] == 1 or PARAM_DICT['remove'] == True):
        remove_str = 'Successful virions\nremoved'
    else:
        remove_str = 'Successful virions\nleft in'

    params_str = "Scale: 1/" + str(PARAM_DICT['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT['gamma'])+'\nvMax = '+str(PARAM_DICT['vMax'])+', cells = '+str(PARAM_DICT['cell_count'])
    
    if (PARAM_DICT['fixed_mean'] == True or PARAM_DICT['fixed_mean'] == 'True' or PARAM_DICT['fixed_mean'] == 1 or PARAM_DICT['fixed_mean'] == '1'):
        params_str += '\nFixed distribution='+str(PARAM_DICT['fixed_mean'])
        
        if (PARAM_DICT['distribution'] in ['geometric', '1inf-geo']):
            params_str += '\n' + PARAM_DICT['distribution'] + ', p='+str(PARAM_DICT['mean'])
            if (PARAM_DICT['distribution'] == '1inf-geo'):
                params_str += ', f_1='+str(PARAM_DICT['f_1'])

        elif (PARAM_DICT['distribution'] == 'poisson'):
            params_str += '\n' + PARAM_DICT['distribution'] + ', '+r'$\lambda=$'+str(PARAM_DICT['mean'])
    else:
        params_str += '\n' + PARAM_DICT['distribution'] + ', vMax_c='+str(PARAM_DICT['vMax_c'])
        if (PARAM_DICT['distribution'] == '1inf-geo'):
                params_str += ', f_1='+str(PARAM_DICT['f_1'])
    
    if ('clump' in PARAM_DICT['simul_name']):
        if (PARAM_DICT['simul_name'] == 'clump'):
            ax.text(1.1 * xMin, .01 * yMax, params_str + "\n"+remove_str)
        
        elif (PARAM_DICT['simul_name'] == 'clump_comp'):
            ax.text(1.1 * xMin, .01 * yMax, params_str + '\n' + r'$\kappa=$'+str(PARAM_DICT['kappa'])+"\n"+remove_str)
        
        elif (PARAM_DICT['simul_name'] == 'clump_acc_dam'):
            ax.text(1.1 * xMin, .01 * yMax, params_str + '\n' + r'$\beta=$'+str(PARAM_DICT['beta'])+"\n"+remove_str)