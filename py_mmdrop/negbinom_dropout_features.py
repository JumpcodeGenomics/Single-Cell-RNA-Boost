from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import numpy as np
import statsmodels.api as sm
import pandas as pd

def negbinom_dropout_features(fit, ntop=None, method="fdr_bh", qval_thresh=0.01, suppress_plot=True):

    def negbinom_dispvsexp(fit, suppress_plot=True):
        vals = fit['vals']
        size_g = np.array(fit['sizes'])

        forfit = (fit['sizes'] < max(size_g)) & (vals['tjs'] > 0) & (size_g > 0)
        higher = np.log(vals['tjs'] / vals['nc'])/np.log(2) > 4

        if np.sum(higher) > 2000:
            forfit = forfit & higher
            
        X = np.log((vals['tjs'] / vals['nc'])[forfit])
        X = sm.add_constant(X)
            
        y = np.log(size_g[forfit])
            
        model = sm.OLS(y, X).fit()

        if not suppress_plot:
            plt.scatter(X[:, 1], y)
            plt.xlabel("Log Mean Expression")
            plt.ylabel("Log Size")
            plt.plot(X[:, 1], model.predict(X), color='red')
            plt.show()
            
        return model.params


    vals = fit['vals']
    coeffs = negbinom_dispvsexp(fit, suppress_plot=True)
    exp_size = np.exp(coeffs[0] + coeffs[1] * np.log(vals['tjs'] / vals['nc']))

    droprate_exp = np.zeros(vals['ng'].astype(int))
    droprate_exp_err = np.zeros(vals['ng'].astype(int))

    for i in range(vals['ng'].astype(int)):
        mu_is = vals['tjs'][i] * vals['tis'] / vals['total']
        p_is = (1 + mu_is / exp_size[i]) ** (-exp_size[i])
        p_var_is = p_is * (1 - p_is)
        droprate_exp[i] = np.sum(p_is) / vals['nc']
        droprate_exp_err[i] = np.sqrt(np.sum(p_var_is) / (vals['nc'] ** 2))

    droprate_exp[droprate_exp < 1 / vals['nc']] = 1 / vals['nc']

    droprate_obs = vals['djs'] / vals['nc']
    droprate_obs_err = np.sqrt(droprate_obs * (1 - droprate_obs) / vals['nc'])

    diff = droprate_obs - droprate_exp
    combined_err = np.sqrt(droprate_exp_err**2 + droprate_obs_err**2)

    Zed = diff / combined_err
    p_value = (1 - norm.cdf(Zed))
    p_value = pd.DataFrame({'p_val':p_value, "inv_diff" : droprate_obs - droprate_exp}, index=pd.array(fit['gene_names']))
    p_value = p_value.sort_values(by=['p_val', 'inv_diff'], ascending=[True,False])

    # Adjust p-values for multiple testing
    reject, q_val, _, _ = multipletests(np.array(p_value['p_val']), method=method)
    p_value['q_val'] = q_val

    if ntop is not None:
        top_genes = np.argsort(list(p_value_dict.values()))[:ntop]
    else:
        top_genes = [gene for i, gene in enumerate(p_value.index) if p_value.loc[gene]['q_val'] < qval_thresh]

    p_value['highly_variable'] =  p_value.index.isin(top_genes)

    # You'll need to implement NBumiFitDispVsMean function or use an equivalent.
    # Make sure to import necessary libraries for that function.

    return p_value