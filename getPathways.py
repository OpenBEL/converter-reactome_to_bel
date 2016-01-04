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
import json
from bioservices import Reactome
from jinja2 import Environment, FileSystemLoader

s = Reactome()

PATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(os.path.join(PATH, 'templates')),
    trim_blocks=False)


def render_template(template_filename, context):
    return TEMPLATE_ENVIRONMENT.get_template(template_filename).render(context)


def getPathways():
    # pathways = s.get_list_pathways()
    # with open('pathways.json', 'w') as f:
    #     json.dump(pathways, f, indent=4)

    # reactions = s.get_all_reactions()
    # with open('reactions.txt', 'w') as f:
    #     json.dump(reactions, f, indent=4)

    hierarchy = s.pathway_hierarchy('homo sapiens')
    with open('hierarchy.xml', 'w') as f:
        f.write(hierarchy)

    # participants = s.pathway_participants(109581)
    # with open('participants.json', 'w') as f:
    #     json.dump(participants, f, indent=4)

    # reaction = 'R-MTU-870392'
    # reactants = s.bioservices_get_reactants_from_reaction_identifier(reaction)
    # with open('reactants.json', 'w') as f:
    #     json.dump(reactants, f, indent=4)




def main():
    getPathways()


if __name__ == '__main__':
    main()

