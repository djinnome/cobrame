

import argparse
import pickle
import cobrame
from cobrame.io.json import load_reduced_json_me_model, load_json_me_model
import os, re
from os.path import abspath, dirname

from cobrame.core.reaction import MetabolicReaction, TranscriptionReaction, TranslationReaction
import pandas as pd

import json
import ecolime
from cobrame.core.reaction import TranscriptionReaction
from qminospy.me1 import ME_NLP1

def compute_gene_essentiality_at_growth_rate(me, gr, out_location, start, stop):
    """Thanks to Colton Lloyd for this function: https://github.com/SBRG/cobrame/issues/29#issuecomment-417487803"""
    me_nlp = ME_NLP1(me, growth_key='mu')
    expressions = me_nlp.compile_expressions()
    me_nlp.compiled_expressions = expressions

    hs = None

    all_genes = sorted( me.metabolites.query( re.compile("^RNA_b[0-9]" )),
                        key=lambda gene: gene.id)
    results = {}
    for gene_RNA in all_genes[start:stop]:

        default_bounds = {}
        for r in gene_RNA.reactions:
            if not r.id.startswith("DM") and not \
                    isinstance(r, TranscriptionReaction):
                default_bounds[r] = (r.lower_bound, r.upper_bound)
                r.knock_out()
        x, status, hs = me_nlp.solvelp(gr, basis=hs)

        if status == 'optimal':
            results[gene_RNA.id] = 'NONESSENTIAL'
            print("Biomass dilution for {} KO: {} ".format(gene_RNA.id,  me.solution.x_dict['biomass_dilution']))
            metabolic_fluxes = pd.DataFrame( {'MetabolicFlux': me.get_metabolic_flux() } )
            expression_fluxes = pd.DataFrame( {'TranscriptionFlux': me.get_transcription_flux(),
                                               'TranslationFlux': me.get_translation_flux() })
            metabolic_fluxes.to_csv(  os.path.join( out_location,
                                                    'MetabolicFlux.{}.csv'.format( gene_RNA.id )))
            expression_fluxes.to_csv( os.path.join( out_location,
                                                    'ExpressionFlux.{}.csv'.format(gene_RNA.id)))
        else:
            results[gene_RNA.id] = 'ESSENTIAL'

        print("%s\t%s" % (gene_RNA.id.split("_")[1], str(status)))

        with open("%s/iJL1678b_essentiality_%.2f_gr_%d_%d.json" % (out_location, gr, start, stop),
                  "w") as outfile:
            json.dump(results, outfile, indent=True)

        # Reset bounds
        for r in default_bounds:
            r.lower_bound = default_bounds[r][0]
            r.upper_bound = default_bounds[r][1]

            

    





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
    return sorted(loci)

def knock_out_reactions_from_model( me, rxn_list ):
    rxns = []
    for rxn_id in rxn_list:
        rxn = me.reactions.get_by_id( rxn_id )
        rxn.bounds = (0, 0)
        rxns.add( rxn )
    return sorted(rxns, key = lambda rxn: rxn.id)
def knock_out_genes_from_model( me, gene_list ): 
    rxns = set()
    for gene in gene_list:
            # Find all complexes the gene product is part of and knock out the associated reactions
            protein = me.metabolites.get_by_id('protein_'+gene)
            for cplx in protein.complexes:
                print('Complex (%s) knocked out in model' % cplx.id)
                for rxn in cplx.metabolic_reactions:
                    rxn.bounds = (0,0)
                    rxns.add(rxn)
    return sorted(rxns, key=lambda rxn: rxn.id)
def get_metabolic_rxns( me ):
    return [rxn for rxn in me.reactions if isinstance(rxn, MetabolicReaction)]
# In[ ]:
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run COBRAME knockouts')
    parser.add_argument('outdir', help='Output directory')
    #parser.add_argument('ko', choices=['metabolic_rxns','all_rxns', 'transcription_rxns', 'translation_rxns','gene'])
    parser.add_argument('start',type=int, help='Start Gene (0-based)')
    parser.add_argument('stop', type=int, help='End Gene (does not include this gene)')
    parser.add_argument('growth_rate', type=float, help='Growth rate to consider essential')
    args = parser.parse_args()
    start, stop  = args.start, args.stop
    if start > stop:
        start, stop = stop, start
    os.makedirs(args.outdir, exist_ok=True)
    with open('/home/meuser/me_models/iJL1678b.pickle', 'rb') as f:
         me = pickle.load(f)
    compute_gene_essentiality_at_growth_rate( me, args.growth_rate, args.outdir, start, stop )

            



