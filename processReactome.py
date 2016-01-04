#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:  program.py <customer>

  program.py  <customer>
  program.py (-h | --help)
  program.py --version

Options:
  -h --help     Show help.
  --version     Show version.
"""

# TODO
#
# Ty found these errors
# Group-20 and Group-30 are identical, including rxnId.  Another pair is Group-1383 and Group-1417.  In total there are 80 pairs with identical rxnIds, and three sets of three identical rxnIds.  I spot checked a few, all the info (SET, bel statements, etc) were the same within a pair that have the same rxnId.
#   WSH Should be fixed by creating a set from ReactionList- Looks like the human, mouse and rat share rxnIds.  I’ll fix that before the next drop.

# There are cases of identical “=>” relationships that are not due to identical rxnIds.  I look into one of them (Group-65 and Group-66).  In this case, the evidence lines for these groups are the same except for GYS1-a is in one and GYS1-b is in the other.  A google search tells me that a is the unphosphorylated form, and b is the phosphorylated form.  Everything other than the evidence lines for these groups are the same, which seems like a reasonable interpretation.  I looked up a few more of these cases and they all were cases where there were separate reactions for the phosphorylated and unphosphorylated forms of an enzyme catalyzing the same reaction.
#   WSH - Not sure how to handle this.  Need a more sophisticated parse of Reactome at the very least to capture these subtleties.

# There is some weird nesting of the location with complexes.  See Group-10.  Essentially, there are complexes with loc() of this form:
# complex(complex(p(x,loc(Y)), loc(y)), loc(y)), with a single hasComponent statement showing that it is composed of complex(p(x,loc(Y)), loc(y)).  I’m still working out exactly how we will be using the hasComponent statements – right now I don’t foresee needing this to be fixed, but I may run into scenarios where this structure causes problems.
#    WSH - The highlighted complex presentation is a bug – I’ll fix that too.  I de-dup the dimers, but didn’t sort out the locations at the same time.


import os
import copy
import time
import re
from jinja2 import Environment, FileSystemLoader
import click

from toBel import toBel, dedup, dedupList, escapeBelString, setBelVersion
from reactome_webservice import getEntityData, getReactions

# Overwrite logs on each run -> filemode = 'w'
import log_setup
log = log_setup.getLogger(name='root')

# log.debug('This message should go to the log file')
# log.info('So should this')
# log.warning('And this, too')

# Jinja template
PATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(os.path.join(PATH, 'templates')),
    trim_blocks=False)

template_filename = 'belscript.jinja2'

speciesList = ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']


def convertSpeciesNameToTaxId(name):
    species = {
        'Homo sapiens': 9606,
        'Mus musculus': 10090,
        'Rattus norvegicus': 10116,
    }
    return species[name]


def render_template(template_filename, evidences, pathways=None):

    context = {}
    # Todo  add the following to a configuration file and automate the date
    if pathways:
        setname = '{} Pathway Reactome Reactions'.format(' and '.join(pathways))
    else:
        setname = 'All Reactome Reactions'
    context['BEL_DOCUMENT_NAME'] = setname
    context['BEL_DOCUMENT_DESCRIPTION'] = '{} converted to BEL 1.0'.format(setname)
    context['AUTHORS'] = 'Selventa; Nimisha Schneider; Natalie Catlett; William Hayes'
    context['DATE'] = time.strftime("%Y-%m-%d")
    context['VERSION'] = '0.1'
    context['CONTACT_EMAIL'] = 'whayes@selventa.com'
    context['evidences'] = evidences

    return TEMPLATE_ENVIRONMENT.get_template(template_filename).render(context)


def buildStatements(catalysts, inputs, outputs):

    input_bel = ', '.join([bel for input in inputs for bel in input])
    output_bel = ', '.join([bel for output in outputs for bel in output])

    rxn = 'rxn(reactants({}), products({}))'.format(input_bel, output_bel)
    statements = []
    if catalysts:
        for catalyst in catalysts:
            if catalyst:
                for bel in catalyst:
                    statement = 'cat({}) => {}'.format(bel, rxn)
                    statements.append(statement)
            else:
                log.error('No BEL in catalyst: {}'.format(catalysts))

    if not statements:
        statements.append(rxn)

    for results in (catalysts, inputs, outputs):
        for result in results:
            if result:
                for bel in result:
                    for statement in result[bel]:
                        statements.append(statement)
            else:
                log.error('No BEL in result: {}'.format(results))

    return statements


def buildBelEvidences(reactionList, belversion, pathways=None):
    ''' Load reactions and build BEL Evidences'''

    rxnUrlTpl = 'http://www.reactome.org/PathwayBrowser/#'
    evidences = []
    bad_namespaces_evidences = []  # collect bad evidences

    # indexCnt = 0
    for rxnId, rxnName in reactionList:

        print('rxnId: ', rxnId)

        # indexCnt += 1
        # if indexCnt > 2:
        #     break

        # Process Annotation information
        rxnData = getEntityData(rxnId)
        if 'stableIdentifier' not in rxnData:
            continue

        stableId = rxnData['stableIdentifier']['displayName']
        rxnUrl = '{}{}'.format(rxnUrlTpl, stableId)
        rxnName = escapeBelString(rxnData['displayName'])
        rxnType = rxnData['schemaClass']

        # Todo  collect all compartments and annotate
        compartment = 'Unknown'
        if 'compartment' in rxnData:
            compartment = rxnData['compartment'][0]['displayName']

        rxnAuthor = rxnDate = None
        if 'created' in rxnData:
            try:
                matches = re.search(r'(.*?),\s+(\d{4,4}-\d{2,2}-\d{2,2})', rxnData['created']['displayName'])
                if matches:
                    rxnAuthor = matches.group(1)
                    rxnDate = matches.group(2)
            except:
                log.info('Rxn - cannot find created date and author in object: {}'.format(rxnId))

        if rxnDate and rxnAuthor:
            citation = '{{"Online Resource", "{}", "{}", "{}", "{}"}}'.format(rxnName, rxnUrl, rxnDate, rxnAuthor)
        else:
            citation = '{{"Online Resource", "{}", "{}"}}'.format(rxnName, rxnUrl)

        evidence = {
            'name': rxnName,
            'rxnId': rxnId,  # TODO remove after debugging
            'rxnType': rxnType,
            'compartment': compartment,
            'species': rxnData['speciesName'],
            'species_tax_id': convertSpeciesNameToTaxId(rxnData['speciesName']),
            'summary_text': rxnName,
            'citation': citation,
        }

        # Process BEL Statement
        catalysts, inputs, outputs = [], [], []

        if 'catalystActivity' in rxnData:

            for catalyst in rxnData['catalystActivity']:
                catalysts.append(toBel(catalyst['dbId']))
                # print('Catalyst: {}'.format(catalyst['dbId']))

        if 'input' in rxnData:
            for input in dedup(rxnData['input']):
                inputs.append(toBel(input['dbId']))
        if 'output' in rxnData:
            for output in dedup(rxnData['output']):
                outputs.append(toBel(output['dbId']))

        print('Catalysts ', catalysts)
        print('Inputs ', inputs)
        print('Outputs ', outputs)
        print('\n')

        statements = buildStatements(catalysts, inputs, outputs)
        evidence['statements'] = dedupList(statements)

        bad_namespace_flag = False
        for statement in statements:
            if 'ENSEMBL' in statement or 'EMBL' in statement:
                bad_namespace_flag = True

        # with open('tmp_evidences.json', 'a') as f:
        #     json.dump(evidence, f, indent=4)
        #     f.write('\n\n')

        if bad_namespace_flag:
            bad_namespaces_evidences.append(copy.deepcopy(evidence))
        else:
            evidences.append(copy.deepcopy(evidence))

    belscript = render_template(template_filename, evidences, pathways=pathways)

    fn = 'reactome.bels'
    if belversion == '2':
        fn += '2'

    with open(fn, 'w') as f:
        f.write(belscript)

    import json
    with open('bad_evidences.json', 'w') as f:
        json.dump(bad_namespaces_evidences, f, indent=4)


@click.command()
@click.option('--belversion', '-b', default='1', type=click.Choice(['1', '2']), help="Use Bel 1 by default or select Bel 2")
@click.option('--species', '-s', multiple=True, type=click.Choice(['all', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus']))
@click.option('--pathways', '-p', default=None, multiple=True, help="Restrict to specific Reactome Pathway(s) - e.g. Metabolism - can use multiple -p Metabolism -p Pathway2 ...")
def main(belversion, species, pathways):
    """Process Reactome into BEL

    Example:  ./processReactome.py -b 1 -s "Homo sapiens" -s "Mus musculus" -p Metabolism
    Example:  ./processReactome.py -b 2 -s "Homo sapiens" -p Metabolism -p "Transmembrane transport of small molecules"
    Result: reactome.bels
    """
    setBelVersion(belversion)

    # quit and show help if no arguments are set
    if not species:
        click.help_option()
        print('Here')

    if 'all' in species:
        species = speciesList

    # Collect reactions
    reactionList = []
    for specie in species:
        reactionList.extend(getReactions(specie, pathways=pathways))

    # import json
    # with open('reactionlist.json', 'w') as f:
    #     json.dump(reactionList, f, indent=4)
    # quit()

    # human, mouse and rat share rxnIds
    reactionList = set(reactionList)

    buildBelEvidences(reactionList, belversion, pathways=pathways)

    # buildBelEvidences([('109514', 'Test')])
    # # buildBelEvidences([('450092', 'Test')])
    # quit()

    # # Get Metabolism reactions for human only
    # reactionList = getReactions('homo sapiens', pathway='Metabolism')
    # quit()

    # # Get all Reactome reactions for model organisms
    # reactionList = []
    # for species in ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']:
    #     reactionList.extend(getReactions(species))


if __name__ == '__main__':

    main()

