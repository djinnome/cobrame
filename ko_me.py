

import argparse
import pickle
import cobrame
from cobrame.io.json import load_reduced_json_me_model, load_json_me_model
import os

# In[8]:


with open('/home/meuser/me_models/iJL1678b.pickle', 'rb') as f:
    me = pickle.load(f)
    





def solve_me_model(me, max_mu, precision=1e-6, min_mu=0, using_soplex=True,
                  compiled_expressions=None):
    if using_soplex:
        from cobrame.solve.algorithms import binary_search
        binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True, mu_accuracy=precision,
                      compiled_expressions=compiled_expressions)
    else:
        from qminospy.me1 import ME_NLP1
        # The object containing solveME methods--composite that uses a ME model object 
        me_nlp = ME_NLP1(me, growth_key='mu')
        # Use bisection for now (until the NLP formulation is worked out)
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu)
        me.solution.f = me.solution.x_dict['biomass_dilution']
        return me
        
def show_escher_map(me, solution=None):
    import escher
    view = escher.Builder("iJO1366.Central metabolism")
    view.reaction_data = me.get_metabolic_flux(solution=solution)
    return view




# In[ ]:
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run COBRAME knockouts')
    parser.add_argument('kofile', type=argparse.FileType('r'))
    parser.add_argument('outdir')
    args = parser.parse_args()
    for rxnId in open(args.kofile):
        with open('/home/meuser/me_models/iJL1678b.pickle', 'rb') as f:
            me = pickle.load(f)
            rxn = me.reactions.get_by_id( rxnId )
            rxn.lower_bound = 0
            rxn.upper_bound = 0
            solve_me_model(me, 1., min_mu = .1, precision=1e-2, using_soplex=False)
            print("Biomass dilution for {} KO: {} ".format(rxnId,  me.solution.x_dict['biomass_dilution']))
            metabolic_fluxes = pd.DataFrame({'MetabolicFlux': me.get_metabolic_flux()})
            expression_fluxes = pd.DataFrame({'TranscriptionFlux': me.get_transcription_flux(),
                                        'TranslationFlux': me.get_translation_flux()})
            metabolic_fluxes.to_csv(os.path.join(args.outdir, 'MetabolicFlux.{}.csv'.format(rxn.id)))
            expression_fluxes.to_csv(os.path.join(args.outdir('ExpressionFlux.{}.csv'.format(rxn.id))))
            



