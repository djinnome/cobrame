

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


def get_genes( me ):
    loci = set()
    for reaction in me.reactions:
        if isinstance(reaction, cobrame.core.reaction.TranscriptionReaction):
            for rna_id in reaction.transcription_data.RNA_products:
                locus_id = rna_id.replace("RNA_", "", 1)
                loci.add( locus_id )
    return loci

# In[ ]:
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run COBRAME knockouts')
    parser.add_argument('outdir', help='Output directory')
    parser.add_argument('start',type=int, help='Start Gene (0-based)')
    parser.add_argument('stop', type=int, help='End Gene (does not include this gene)')
    args = parser.parse_args()
    start, stop  = args.start, args.stop
    if start > stop:
        start, stop = stop, start
    os.makedirs(args.outdir, exist_ok=True)
    with open('/home/meuser/me_models/iJL1678b.pickle', 'rb') as f:
            me = pickle.load(f)
    genes = get_genes( me )
    for gene in genes[start:stop]:
        
        with open('/home/meuser/me_models/iJL1678b.pickle', 'rb') as f:
            me = pickle.load(f)
            me.remove_genes_from_model( [gene_id] )
            solve_me_model(me, 1., min_mu = .1, precision=1e-2, using_soplex=False)
            print("Biomass dilution for {} KO: {} ".format(gene_id,  me.solution.x_dict['biomass_dilution']))
            metabolic_fluxes = pd.DataFrame({'MetabolicFlux': me.get_metabolic_flux()})
            expression_fluxes = pd.DataFrame({'TranscriptionFlux': me.get_transcription_flux(),
                                        'TranslationFlux': me.get_translation_flux()})
            metabolic_fluxes.to_csv(os.path.join(args.outdir, 'MetabolicFlux.{}.csv'.format(gene_id)))
            expression_fluxes.to_csv(os.path.join(args.outdir('ExpressionFlux.{}.csv'.format(gene_id))))
            



