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

import os
import requests
from lxml import etree
import json

from bioservices import Reactome

import logging
log = logging.getLogger('root')

s = Reactome()

wsUrl = 'http://reactome.org/ReactomeRESTfulAPI/RESTfulWS'  # base Url for Reactome Webservice


def getEntityData(dbId, downloadDir=None):
    ''' Get reactome entities by database id'''

    if not downloadDir:
        downloadDir = './downloadedEntities'

    fn = '{}/{}.json'.format(downloadDir, dbId)

    if os.path.isfile(fn):
        with open(fn, 'r') as f:
            entity = json.load(f)

    else:
        try:
            r = requests.get("{}/queryById/DatabaseObject/{}".format(wsUrl, dbId))
        except:
            log.info('Reactome GET failed: {}'.format(dbId))

        entity = r.json()
        with open(fn, 'w') as f:
            json.dump(entity, f, indent=4)

    return entity


def getReactions(species, pathways=None, xmlFn=None):
    ''' Collect all reactions for specified species

    Inputs:
        xmlfn      optional xml filename for pathway hierarchy
        pathway    filter by an optional pathway

    Returns:
        reactionList   list of tuples containing (dbId, displayName) of each reaction

    Caching results of query locally in the xmlfn or (/tmp/reactome_pathway_hierarchy.xml).
    '''
    reaction_types = ['Reaction', 'BlackBoxEvent', 'Polymerisation', 'Depolymerisation', 'FailedReaction']

    if not xmlFn:
        xmlFn = '/tmp/reactome_pathway_hierarchy.xml'
    reactionList = []
    if os.path.isfile('pathway_hierarchy.xml'):
        with open('downloads/pathway_hierarchy.xml', 'r') as f:
            hierarchyXml = f.read()
    else:
        hierarchyXml = s.pathway_hierarchy(species)

        with open('downloads/hierarchy.xml', 'w') as f:
            f.write(hierarchyXml)

    root = etree.fromstring(hierarchyXml)

    docs = []

    for reaction_type in reaction_types:
        if pathways:
            for pathway in pathways:
                docs.append(root.xpath("/Pathways/Pathway[contains(@displayName, '{}')]//{}".format(pathway, reaction_type)))
        else:
            docs.append(root.xpath("//{}".format(reaction_type)))

    for doc in docs:
        for rxn in doc:
            dbId = rxn.get('dbId')
            displayName = rxn.get('displayName')
            # print('Reactions: ', displayName)
            reactionList.append((dbId, displayName))

    log.debug('getReactions  Species: {}  Cnt: {}'.format(species, len(reactionList)))
    return reactionList
    # print('DumpVar:\n', json.dumps(reactionList, indent=4))
    # print('Length ', len(reactionList))


def getSets():
    ''' Get collections of data from Reactome
        pathway participants and reference molecules and proteins
    '''
    queries = {
        'pathwayComplexes': 'pathwayComplexes/1430728',
        'pathwayParticipants': 'pathwayParticipants/1430728',
    }
    for key in queries:
        fn = 'downloads/{}.json'.format(key)

        if not os.path.isfile(fn):
            r = requests.get("{}/{}".format(wsUrl, queries[key]))
            data = r.json()
            with open(fn, 'w') as f:
                json.dump(data, f, indent=4)

    queries = {
        'refMolecules': 'getReferenceMolecules',
        'refProteins': 'getUniProtRefSeqs',
    }
    for key in queries:
        fn = 'downloads/{}.txt'.format(key)

        if not os.path.isfile(fn):
            r = requests.get("{}/{}".format(wsUrl, queries[key]))
            data = r.text
            with open(fn, 'w') as f:
                f.write(data)


def indexEntities():
    import glob
    files = glob.glob('./downloadedEntities/*.json')

    with open('./entityIndex.txt', mode='w') as o:
        for fn in files:
            print('FN: ', fn)
            with open(fn, 'r') as f:
                entity = json.load(f)
            if 'dbId' in entity:
                dbId = entity['dbId']
                schemaClass = entity['schemaClass']
                o.write('{}\t{}\n'.format(dbId, schemaClass))


def main():
    indexEntities()

if __name__ == '__main__':
    main()

