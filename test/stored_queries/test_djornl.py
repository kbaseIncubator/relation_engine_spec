"""
Tests for the Dan Jacobson ORNL Arabidopsis stored queries.
"""
import json
import time
import unittest
import requests
import os

from test.helpers import get_config, modified_environ
from test.stored_queries.helpers import create_test_docs
from importers.djornl.parser import DJORNL_Parser

_CONF = get_config()
_NOW = int(time.time() * 1000)
_TEST_DIR = '/app/test'
_VERBOSE = 0


def print_db_update(response, collection):
    if not _VERBOSE:
        return
    print(f"Saved docs to collection {collection}!")
    print(response.text)
    print('=' * 80)


class Test_DJORNL_Stored_Queries(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # import the results file
        results_file = os.path.join(_TEST_DIR, 'djornl', 'results.json')
        with open(results_file) as fh:
            cls.json_data = json.load(fh)

        cls.no_results = {'nodes': [], 'edges': []}

        # load the DB
        root_path = os.path.join(_TEST_DIR, 'djornl', 'test_data')
        with modified_environ(RES_ROOT_DATA_PATH=root_path):
            parser = DJORNL_Parser()
            config = parser.config()

            edge_data = parser.load_edges()
            r = create_test_docs(config['_NODE_NAME'], edge_data['nodes'])
            print_db_update(r, config['_NODE_NAME'])
            r = create_test_docs(config['_EDGE_NAME'], edge_data['edges'])
            print_db_update(r, config['_EDGE_NAME'])

            node_metadata = parser.load_node_metadata()
            r = create_test_docs(config['_NODE_NAME'], node_metadata['nodes'], True)
            print_db_update(r, config['_NODE_NAME'])

            cluster_data = parser.load_cluster_data()
            r = create_test_docs(config['_NODE_NAME'], cluster_data['nodes'], True)
            print_db_update(r, config['_NODE_NAME'])


    def submit_query(self, query_name, query_data={}):
        """submit a database query"""

        q_data_str = json.dumps(query_data)
        if _VERBOSE:
            print('query data string: ' + q_data_str)
        response = requests.post(
            _CONF['re_api_url'] + '/api/v1/query_results',
            params={'stored_query': query_name},
            data=q_data_str
        ).json()

        return response


    def check_expected_results(self, description, response, expected):

        if _VERBOSE:
            print("Running test " + description)

        results = response['results'][0]
        self.assertEqual(
            set([n["_key"] for n in results['nodes']]),
            set(expected['nodes'])
        )

        self.assertEqual(
            set([e["_key"] for e in results['edges']]),
            set(expected['edges'])
        )


    def test_fetch_all(self):

        response = self.submit_query('djornl_fetch_all')
        self.check_expected_results(
            "djornl_fetch_all",
            response,
            self.json_data['fetch_all']
        )

        # ensure that all the cluster data is returned OK
        node_data = response['results'][0]['nodes']
        nodes_with_clusters = [json.dumps({
            '_key':     n['_key'],
            'clusters': n['clusters']
        }) for n in node_data if 'clusters' in n]
        self.assertEqual(
            set(nodes_with_clusters),
            set([json.dumps(this) for this in self.json_data['load_cluster_data']['nodes']])
        )

    # indexing schema in results.json
    # self.json_data[query_name][param_name][param_value]["distance"][distance_param]
    # e.g. for fetch_clusters data:
    # "fetch_clusters": {
    #   "cluster_ids": {
    #     "markov_i2:6__markov_i4:3": {
    #       "distance": {
    #         1: {
    #           "nodes": [ node IDs ],
    #           "edges": [ edge data ],
    #         }
    #       }
    #     }
    #   }
    # }
    # if param_value is an array, join the array entities with "__"
    # results are in the form {"nodes": [...], "edges": [...]}
    # nodes are represented as a list of node[_key]
    # edges are objects with keys _to, _from, edge_type and score

    def test_fetch_phenotypes(self):

        for (fetch_args, key_data) in self.json_data['fetch_phenotypes']['keys'].items():
            for (distance, distance_data) in key_data['distance'].items():
                resp = self.submit_query('djornl_fetch_phenotypes', {
                    "keys": fetch_args.split('__'),
                    "distance": int(distance),
                })
                self.check_expected_results(
                    "fetch phenotypes with args " + fetch_args + " and distance " + distance,
                    resp,
                    distance_data
                )


    def test_fetch_genes(self):

        for (fetch_args, key_data) in self.json_data['fetch_genes']['keys'].items():
            for (distance, distance_data) in key_data['distance'].items():
                resp = self.submit_query('djornl_fetch_genes', {
                    "keys": fetch_args.split('__'),
                    "distance": int(distance),
                })
                self.check_expected_results(
                    "fetch genes with args " + fetch_args + " and distance " + distance,
                    resp,
                    distance_data
                )


    def test_fetch_clusters(self):

        for (fetch_args, cluster_data) in self.json_data['fetch_clusters']['cluster_ids'].items():
            for (distance, distance_data) in cluster_data['distance'].items():
                resp = self.submit_query('djornl_fetch_clusters', {
                    "cluster_ids": fetch_args.split('__'),
                    "distance": int(distance),
                })
                self.check_expected_results(
                    "fetch clusters with args " + fetch_args + " and distance " + distance,
                    resp,
                    distance_data
                )


    def test_search_nodes(self):

        for (search_text, search_data) in self.json_data['search_nodes']['search_text'].items():
            for (distance, distance_data) in search_data['distance'].items():
                resp = self.submit_query('djornl_search_nodes', {
                    "search_text": search_text,
                    "distance": int(distance),
                })
                self.check_expected_results(
                    "search nodes with args " + search_text + " and distance " + distance,
                    resp,
                    distance_data
                )
