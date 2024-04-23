import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit
import numpy.random
from scipy.stats import norm, lognorm, skewnorm
from scipy.optimize import fsolve
import pandas as pd
#====================================================================
def Trapezoid(x_data, y_data):
    # Check that the lengths of x_data and y_data are equal
    assert len(x_data) == len(y_data), "x_data and y_data must have the same length\nGot len(x)="+str(len(x_data))+', len(y)='+str(len(y_data))

    integral = 0
    # Loop over the data points
    for i in range(len(x_data) - 1):
        # Calculate the width of the subinterval
        h = x_data[i + 1] - x_data[i]
        # Calculate the area of the trapezoid
        area = h * (y_data[i] + y_data[i + 1]) / 2
        # Add the area to the integral value
        integral += area
    # Return the integral value
    return integral
#====================================================================
def CreateMuSigAmp(num_gauss, model_type,
                   min_mu=100, max_mu=5000, init_mu=None,
                   min_std=0, max_std=200, init_std=None, # 150 good
                   min_amp=0.01, max_amp=2, init_amp=None,
                   min_skew=-100, max_skew=100, init_skew=None):
    assert (type(num_gauss) == int) and (num_gauss >= 1)
    assert model_type in ['normal', 'lognormal', 'skewnormal'], "'model_type' must be 'norma', 'lognormal', or 'skewnormal'"
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (init_mu == None):
        init_mu = .5 * (max_mu - min_mu)
    if (init_std == None):
        init_std = .5 * (max_std - min_std)
    if (init_amp == None):
        init_amp = .5 * (max_amp - min_amp)
    if (init_skew == None):
        init_skew = .5 * (max_amp - min_amp)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    params = Parameters()

    params.add('mu_1', value=init_mu, min=min_mu, max=max_mu)

    expr = 'mu_1'
    for i in range(1, num_gauss+1):
        params.add('std_'+str(i), value=init_std, min=min_std, max=max_std)
        params.add('amp_'+str(i), value=init_amp, min=min_amp, max=max_amp)

        if i > 1:
            params.add('offset_'+str(i-1), value=max_mu/2, min=0, max=max_mu)
            expr += '+offset_'+str(i-1)
            params.add('mu_'+str(i), expr=expr)

    if (model_type=='skewnormal'):
        for i in range(1, num_gauss+1):
            params.add('skew_'+str(i), value=init_skew, min=min_skew, max=max_skew)

    return params
#====================================================================
def NormalPDF(x, mu, sigma, amp):
    return (amp / (sigma * np.sqrt(2*np.pi))) * (np.exp(-.5 * np.power((x - mu) / sigma, 2)))
#====================================================================
def LogNormalPDF(x, mu, sigma, amp):
    return (amp / (x * sigma * np.sqrt(2*np.pi))) * (np.exp(-(1/np.sqrt(2)) * np.power((np.log(x) - mu) / sigma, 2)))
#====================================================================
def SkewNormalPDF(x, mu, sigma, amp, skew):
    return amp * skewnorm.pdf(x, skew, mu, sigma)
#====================================================================
def Model(x, params, model_type='normal'):
    assert model_type in ['normal', 'lognormal', 'skewnormal'], "Valid model types: normal, lognormal"
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    y = np.zeros_like(x)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    mus  = []
    stds = []
    amps = []
    skews = []

    for key in params.keys():
        if ('mu' in key):
            mus.append(key)
        elif ('std' in key):
            stds.append(key)
        elif ('amp' in key):
            amps.append(key)
        elif ('skew' in key):
            skews.append(key)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (model_type == 'normal'):
        for i in range(len(mus)):
            # y += gaussian(x, amplitude=params[amps[i]].value,
            #                  center=params[mus[i]].value,
            #                  sigma=params[stds[i]].value)
            y += NormalPDF(x, params[mus[i]].value, params[stds[i]].value, params[amps[i]].value)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elif (model_type == 'lognormal'):
        for i in range(len(mus)):
            # y += lognormal(x, amplitude=params[amps[i]].value,
            #                   center=params[mus[i]].value,
            #                   sigma=params[stds[i]].value)
            y += LogNormalPDF(x, params[mus[i]].value, params[stds[i]].value, params[amps[i]].value)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elif (model_type == 'skewnormal'):
        for i in range(len(mus)):
            # y += lognormal(x, amplitude=params[amps[i]].value,
            #                   center=params[mus[i]].value,
            #                   sigma=params[stds[i]].value)
            y += SkewNormalPDF(x, params[mus[i]].value, params[stds[i]].value, params[amps[i]].value, params[skews[i]].value)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return y
#====================================================================
def Residuals(params, x, y_data, model_type):
    y_model = Model(x, params, model_type)

    return abs(y_data - y_model)
#====================================================================
def NegLogLike(params, x, y_data, model_type):
    new_x_data = []
    new_y_data = []

    for i in range(len(y_data)):    
        if ((y_data[i] != 0) and (x[i] != 0)):
            new_x_data.append(x[i])
            new_y_data.append(y_data[i])

    new_x_data = np.asarray(new_x_data)
    new_y_data = np.asarray(new_y_data)
    #print("New x:", new_x_data, ",new y:", new_y_data)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    y_model = Model(new_x_data, params, model_type)

    nll = 0
    epsilon = 10**-10
    for i in range(len(new_y_data)):

        if (y_model[i] < epsilon):
            y_model[i] = epsilon

        if (y_model[i] > 1 - epsilon):
            y_model[i] = 1 - epsilon

        nll += new_y_data[i] * np.log(y_model[i]) + (1 - new_y_data[i]) * np.log(1 - y_model[i])

    return -nll
#====================================================================
def FlattenMeans(x_data, y_data, cutoff, replacement_value=0):
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)

    assert cutoff in x_data, str(cutoff) +" must be in x_data."

    idx = np.where(x_data == cutoff)[0][0]
    for i in range(idx + 1):
        y_data[i] = replacement_value

    return y_data 
#====================================================================
def GenerateClumpDiameter(num_virions, mean, lb, ub, scheme='linear', dist='uniform', target_x=None, target_prob=None):
    assert (scheme=='linear' or scheme=='regular_polygon'), "scheme must be 'linear' or 'regular_polygon'."
    assert (dist=='uniform' or dist=='normal'), "dist must be 'uniform' or 'normal'."
    #----------------------------------------------------------------
    ''' Calculate standard dev. necesary to put lb or ub at target_prob. '''
    if (dist == 'normal'):
        std_lb = FindStdev(mean, lb, target_prob)
        std_ub = FindStdev(mean, ub, target_prob)
        #std2 = FindStdev2(diameter, mean_of_means, mean, lb, ub) # Also gives 70
    #----------------------------------------------------------------
    if (scheme == 'linear'):
        if (dist == 'uniform'):
            diameter = sum(np.random.uniform(low=lb, high=ub, size=num_virions))
        elif (dist == 'normal'):
            # Best values I've gotten: mean=230, sigma=70
            diameter = sum(np.random.normal(loc=mean, scale=std_ub, size=num_virions))
    #----------------------------------------------------------------
    elif (scheme == 'regular_polygon'):
        if (dist == 'uniform'):
            virion_diam = np.random.uniform(low=lb, high=ub, size=1)

        elif (dist == 'normal'):
            virion_diam = np.random.normal(loc=mean, scale=std_ub, size=1)

        diameter = virion_diam / np.sin(np.pi / num_virions)
    #----------------------------------------------------------------
    return diameter
#====================================================================
def FindStdev(mean, target_x, target_prob):

    f = lambda sigma: (1 / (np.sqrt(2*np.pi)*sigma)) * np.exp(-.5 * ((target_x - mean) / sigma)**2) - target_prob
    return fsolve(f, [10])
#====================================================================
def FindStdev2(x, y, mean, lb, ub):
    new_x = []
    new_y = []
    for i in range(len(x)):
        if (x[i] >= lb and x[i] <= ub):
            new_x.append(x[i])
            new_y.append(y[i])

    params = Parameters()
    params.add('sigma', value=70, min=50, max=200)
    params.add('amp', value=1, min=.5, max=1.5)
    
    def Res(params, new_x, new_y):
        y = params['amp'].value * norm.pdf(new_x, mean, params['sigma'].value)
        return abs(y - new_y)

    result = minimize(Res, params, method='least_squares', args=(new_x, new_y),)
    print("++++++++++++++++++++++++++")
    report_fit(result)
    print("++++++++++++++++++++++++++")
    return result.params['sigma'].value, result.params['amp'].value
#====================================================================
def ExtractParams(DataFrame):
    COL = DataFrame.loc[:, 'Parameters']
    ParamDict = {}
    
    index = 0

    for string in COL:
        if not pd.isna(string):
            for i in range(len(string)):   
                if ((string[i] == '=')):
                    index = i
                    break
            try:
                ParamDict[string[:index]] = float(string[(index+1):])
                #print("Is float:", string[:index], string[(index+1):])
            except ValueError:
                ParamDict[string[:index]] = string[(index+1):]
                #print("Is string:", string[:index], string[(index+1):])
                #print("No value")
    
    return ParamDict
#====================================================================
def Replace(input_string, target, replacement):
    if (type(target) == list):
        for i in range(len(target)):
            input_string.replace(target[i], replacement)
    elif (type(target)==str):
        input_string.replace(target, replacement)

    return input_string