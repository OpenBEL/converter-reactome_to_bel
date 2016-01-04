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

import re

from reactome_webservice import getEntityData

import logging
log = logging.getLogger('root')

belVersion = 1

namespaces = ["AFFX", "CHEBIID", "CHEBI", "DOID", "DO", "EGID",
    "GOBPID", "GOBP", "GOCCID", "GOCC", "HGNC", "MESHPP", "MESHCS",
    "MESHC", "MESHCID", "MESHD", "MESHPPID", "MESHCSID", "MESHDID",
    "MGI", "RGD", "SCHEM", "SDIS", "SFAM", "SCOMP", "SPID", "SP"]


####################################################
# Common utilities
####################################################

def setBelVersion(version):
    """Set BEL Version to 1 or 2"""

    global belVersion
    belVersion = version


def dedupList(seq):
    """De-duplicate list"""

    seqset = set(seq)
    return list(seqset)


def dedup(objects):
    """ De-duplicate list of objects based on dbId key"""

    new = []
    seen = {}

    for object in objects:
        dbId = object['dbId']
        if dbId in seen:
            continue
        seen[dbId] = 1
        new.append(object)

    return new


def escapeBelString(belstring):
    """Escape double quotes in BEL string"""
    return belstring.replace('"', '\\"')


def processReferenceEntity(entity):

    compartment = getCompartment(entity)

    displayName = entity['referenceEntity']['displayName']
    name = entity['name'][0]

    # ChEBI ID
    matches = re.search(r'(.*?)\s+\[ChEBI:(\d+)\]', displayName, flags=0)
    if matches:
        chebiId = matches.group(2)
        if belVersion == '2' and compartment:
            bel = 'a(CHEBIID:{}, loc(REACTCOMP:"{}"))'.format(chebiId, compartment)
        else:
            bel = 'a(CHEBIID:{})'.format(chebiId)
        return bel

    # UniProt ID
    matches = re.search(r'(\w+):(\w+)[\W\s]*', displayName, flags=0)
    if matches:
        namespace = matches.group(1)
        accessionId = matches.group(2)
        if belVersion == '2' and compartment:
            bel = 'p({}:{}, loc(REACTCOMP:"{}"))'.format(namespace, accessionId, compartment)
        else:
            bel = 'p({}:{})'.format(namespace, accessionId)
        return bel

    # Default
    if belVersion == '2' and compartment:
        bel = 'p({}, loc(REACTCOMP:"{}"))'.format(name, compartment)
    else:
        bel = 'p({})'.format(name)
    return bel

    log.info('Cannot process referenceEntity: {}  displayName: {}'.format(entity['dbId'], displayName))


def convertCHEBI(chebiTerm):
    ''' Convert ChEBI ID to CHEBIID '''
    chebi = chebiTerm.replace('ChEBI', 'CHEBIID')
    return chebi


def formatBelEntity(belEntity):
    matches = re.search('(\w+):(.*)\s*', belEntity, flags=0)
    if matches:
        ns = matches.group(1)
        entityId = matches.group(2)
        if ns not in namespaces:
            log.info('Unknown namespace: {}'.format(belEntity))

    matches = re.search('\W', entityId, flags=0)
    if matches:
        entityId = '"{}"'.format(entityId)

    return '{}:{}'.format(ns, entityId)


def getCompartment(entity):
    """Extract compartment (location) name from entity"""

    # TODO - only getting the first compartment
    if 'compartment' in entity and 'displayName' in entity['compartment'][0]:
        return entity['compartment'][0]['displayName']
    return None


####################################################
# Convert to BEL Annotations
####################################################
def toBelCompartment(entity):

    if entity['displayName'] == 'smooth endoplasmic reticulum':
        annotation = 'SET CellStructure = "Endoplasmic Reticulum, Smooth"'
        return annotation

    log.info("Cannot process Compartment: {}".format(entity['dbId']))


def toBelEntityCompartment(entity):

    if 'displayName' in entity:
        annotation = 'SET EntityCompartment = "{}"'.format(entity['displayName'])
        return annotation

    log.info('Cannot process EntityCompartment: {}'.format(entity['dbId']))


####################################################
# Convert to BEL terms
####################################################
def toBelOtherEntity(entity):

    if 'name' in entity:
        bel = 'a("{}")'.format(escapeBelString(entity['name'][0]))
        log.debug('toBelOtherEntity: {}  dbId: {}'.format(bel, entity['dbId']))
        return {bel: []}

    log.info('Cannot process OtherEntity: {}'.format(entity['dbId']))


def toBelPolymer(entity):

    compartment = getCompartment(entity)

    if 'crossReference' in entity:
        for cross in entity['crossReference']:
            if 'ChEBI' in cross['displayName']:
                chebi = convertCHEBI(cross['displayName'])
                if belVersion == '2' and compartment:
                    bel = 'a({}, loc(REACTCOMP:"{}"))'.format(chebi, compartment)
                else:
                    bel = 'a({})'.format(chebi)
                log.debug('toBelPolymer crossReference: {}  dbId: {}'.format(bel, entity['dbId']))
                return {bel: []}

    if 'name' in entity:
        bel = 'a("{}")'.format(escapeBelString(entity['name'][0]))
        log.debug('toBelPolymer name: {}  dbId: {}'.format(bel, entity['dbId']))
        return {bel: []}

    log.info('Cannot process Polymer: {}'.format(entity['dbId']))


def toBelSimpleEntity(entity):

    if 'referenceEntity' in entity:
        bel = processReferenceEntity(entity)
        log.debug('toBelSimpleEntity: {}  dbId: {}'.format(bel, entity['dbId']))
        return {bel: []}

    log.info('Cannot process SimpleEntity: {}'.format(entity['dbId']))


def toBelSets(entity):

    compartment = getCompartment(entity)
    # log.debug('toBelSets Compartment: {}'.format(compartment))

    if 'name' in entity:
        if belVersion == '2' and compartment:
            bel = 'a({}, loc(REACTCOMP:"{}"))'.format(escapeBelString(entity['name'][0]), compartment)
        else:
            bel = 'a("{}")'.format(escapeBelString(entity['name'][0]))
        log.debug('toBelSets: {}  dbId: {}'.format(bel, entity['dbId']))

        results = {bel: []}
        members = []
        if 'hasMember' in entity:
            members = entity['hasMember']
        elif 'hasCandidate' in entity:
            members = entity['hasCandidate']
        else:
            log.error('problem with Set: no members')

        for member in dedup(members):
            # Recursively process complex components which can have complexes inside them
            result = toBel(member['dbId'])
            for key in result:
                results[bel].append('{} hasComponent {}'.format(bel, key))
                for statement in result[key]:
                    results[bel].append(statement)

        return results

    log.info('problem with Set: '.format(entity['dbId']))


def toBelGenomeEncodedEntity(entity):

    compartment = getCompartment(entity)

    if 'name' in entity:
        if belVersion == '2' and compartment:
            bel = 'a({}, loc(REACTCOMP:"{}"))'.format(escapeBelString(entity['name'][0]), compartment)
        else:
            bel = 'a("{}")'.format(escapeBelString(entity['name'][0]))

        log.debug('toBelGenomeEncodedEntity: {}  dbId: {}'.format(bel, entity['dbId']))
        return {bel: []}

    log.info('problem with Set: '.format(entity['dbId']))


####################################################
# Recursive conversion to BEL terms
####################################################
def toBelComplexNamed(entity):
    """Version that uses Reactome name of complex for complex entity name"""

    compartment = getCompartment(entity)

    if belVersion == '2' and compartment:
        bel = 'a({}, loc(REACTCOMP:"{}"))'.format(escapeBelString(entity['name'][0]), compartment)
    else:
        bel = 'complex("{}")'.format(escapeBelString(entity['name'][0]))

    results = {bel: []}
    if 'hasComponent' in entity:
        components = dedup(entity['hasComponent'])

        for component in components:
            result = toBel(component['dbId'])
            for key in result:
                results[bel].append('{} hasComponent {}'.format(bel, key))
                for statement in result[key]:
                    results[bel].append(statement)

        log.debug('toBelComplex: {}  dbId: {}'.format(bel, entity['dbId']))
        return results

    log.info('Cannot process Complex: {}'.format(entity['dbId']))


####################################################
# Recursive conversion to BEL terms
####################################################
def toBelComplexComponents(entity):
    """Version that uses complex components for complex entity name"""

    compartment = getCompartment(entity)

    results = []
    belComponents = []
    childStatements = []
    if 'hasComponent' in entity:
        components = dedup(entity['hasComponent'])

        for component in components:
            results.append(toBel(component['dbId']))

        for result in results:
            for key in result:
                belComponents.append(key)

        belComponents = dedupList(belComponents)
        belList = ', '.join(belComponents)

        if belVersion == '2' and compartment:
            bel = 'complex({}, loc(REACTCOMP:"{}"))'.format(belList, compartment)
        else:
            bel = 'complex({})'.format(belList)

        for result in results:
            for key in result:
                childStatements.append('{} hasComponent {}'.format(bel, key))
                for statement in result[key]:
                    childStatements.append(statement)

        log.debug('toBelComplex: {}  dbId: {}'.format(bel, entity['dbId']))

        return {bel: childStatements}

    log.info('Cannot process Complex: {}'.format(entity['dbId']))


def toBelEntityWithAccessionedSequence(entity):

    if 'referenceEntity' in entity:
        bel = processReferenceEntity(entity)

        if bel:
            bel = bel.replace('UniProt', 'SPID')

        log.debug('toBelEntityWithAccessionedSequence: {}  dbId: {}'.format(bel, entity['dbId']))
        return {bel: []}

    elif 'physicalEntity' in entity:
        return toBel(entity['physicalEntity']['dbId'])

    log.info('Cannot process AccessionedSequence: {}'.format(entity['dbId']))


def toBelCatalystActivity(entity):

    if 'physicalEntity' in entity:
        return(toBel(entity['physicalEntity']['dbId']))

    log.info('Cannot process CatalystActivity: {}'.format(entity['dbId']))


####################################################
# Master conversion to BEL terms
####################################################
def toBel(dbId):
    ''' Convert to BEL formats '''

    entity = getEntityData(dbId)

    print('toBel dbId: ', dbId)

    type = entity['schemaClass']

    if (type == 'Compartment'):
        bel = toBelCompartment(entity)
    elif type == 'EntityCompartment':
        bel = toBelEntityCompartment(entity)
    elif type == 'OtherEntity':
        bel = toBelOtherEntity(entity)
    elif type == 'Polymer':
        bel = toBelPolymer(entity)
    elif type == 'Complex':
        bel = toBelComplexComponents(entity)  # alternative version - toBelComplexNamed()
    elif type == 'SimpleEntity':
        bel = toBelSimpleEntity(entity)
    elif type == 'EntityWithAccessionedSequence':
        bel = toBelEntityWithAccessionedSequence(entity)
    elif 'Set' in type:
        bel = toBelSets(entity)
    elif 'GenomeEncodedEntity' in type:
        bel = toBelGenomeEncodedEntity(entity)
    elif 'CatalystActivity' in type:
        bel = toBelCatalystActivity(entity)

    return bel


def main():
    print(escapeBelString('This "is a quoted" string'))

    # print(toBel('450089'))


if __name__ == '__main__':
    main()
