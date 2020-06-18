"""
Loads the Dan Jacobson/ORNL group's gene and phenotype network data into
arangodb.

Running this requires a set of source files provided by the ORNL group.
"""
import json
import requests
import os
import csv

import importers.utils.config as config
CONF = config.load_from_env()

# Path config
_BASE_PATH = os.environ['RES_base_path']
_VERT_PATH = os.path.join(_BASE_PATH, 'aranet2-aragwas-MERGED-AMW-v2_091319_nodeTable.csv')
_CLUSTER_PATHS = [
    os.path.join(_BASE_PATH, 'cluster_data/out.aranetv2_subnet_AT-CX_top10percent_anno_AF_082919.abc.I2_named.tsv'),
    os.path.join(_BASE_PATH, 'cluster_data/out.aranetv2_subnet_AT-CX_top10percent_anno_AF_082919.abc.I4_named.tsv'),
    os.path.join(_BASE_PATH, 'cluster_data/out.aranetv2_subnet_AT-CX_top10percent_anno_AF_082919.abc.I6_named.tsv'),
]
_PHENO_ASSN_PATH = os.path.join(_BASE_PATH, 'aragwas_subnet_phenoassociations_AMW_083019.tsv')
_DOMAIN_CO_OCCUR_PATH = os.path.join(_BASE_PATH, 'aranetv2_subnet_AT-DC_anno_AF_082919.tsv')
_GENE_COEXPR_PATH = os.path.join(_BASE_PATH, 'aranetv2_subnet_AT-CX_top10percent_anno_AF_082919.tsv')
_PPI_HITHRU_PATH = os.path.join(_BASE_PATH, 'aranetv2_subnet_AT-HT_anno_AF_082919.tsv')
_PPI_LIT_PATH = os.path.join(_BASE_PATH, 'aranetv2_subnet_AT-LC_anno_AF_082919.tsv')

# Collection name config
_PHENO_VERT_NAME = 'djornl_phenotype'
_GENE_VERT_NAME = 'djornl_gene'
_EDGE_NAME = 'djornl_edge'


def load_edges(path, score_type):
    # Headers and sample row:
    # node1	node2	edge	edge_descrip	layer_descrip
    # AT1G01370	AT1G57820	4.40001558779779	AraNetv2_log-likelihood-score	AraNetv2-LC_lit-curated-ppi
    with open(path) as fd:
        gene_verts = []
        edges = []
        csv_reader = csv.reader(fd, delimiter='\t')
        next(csv_reader, None)  # skip headers
        for row in csv_reader:
            cols = [c.strip() for c in row]
            gene_verts.append({'_key': cols[0]})
            gene_verts.append({'_key': cols[1]})
            edges.append({
                '_from': f'{_GENE_VERT_NAME}/{cols[0]}',
                '_to': f'{_GENE_VERT_NAME}/{cols[1]}',
                'score': float(cols[2]),
                'score_type': score_type,
            })
    save_docs(_GENE_VERT_NAME, gene_verts)
    save_docs(_EDGE_NAME, edges)


def load_pheno_assns():
    # Headers and sample row:
    # node1	node2	edge	edge_descrip	layer_descrip
    # Na23	AT4G10310	41.300822742442726	AraGWAS-Association_score	AraGWAS-Phenotype_Associations
    with open(_PHENO_ASSN_PATH) as fd:
        pheno_verts = []
        gene_verts = []
        edge_verts = []
        csv_reader = csv.reader(fd, delimiter='\t')
        next(csv_reader, None)  # skip headers
        for row in csv_reader:
            cols = [c.strip() for c in row]
            edge_doc = {
                '_from': f'{_GENE_VERT_NAME}/{cols[1]}',
                '_to': f'{_PHENO_VERT_NAME}/{cols[0]}',
                'score': float(cols[2]),
                'score_type': 'pheno_assn'
            }
            edge_verts.append(edge_doc)
            pheno_verts.append({'_key': cols[0]})
            gene_verts.append({'_key': cols[1]})
    save_docs(_EDGE_NAME, edge_verts)
    save_docs(_PHENO_VERT_NAME, pheno_verts)
    save_docs(_GENE_VERT_NAME, gene_verts)


def load_vert_metadata():
    with open(_VERT_PATH) as fd:
        genes = []
        phenos = []
        csv_reader = csv.reader(fd, delimiter=',')
        next(csv_reader, None)  # skip headers
        for row in csv_reader:
            cols = [c.strip() for c in row]
            go_terms = [c.strip() for c in cols[10].split(',')]
            node_type = cols[1]
            doc = {
                '_key': cols[0],
                'node_type': node_type,
                'transcript': cols[2],
                'gene_symbol': cols[3],
                'gene_full_name': cols[4],
                'gene_model_type': cols[5],
                'tair_computational_desc': cols[6],
                'tair_curator_summary': cols[7],
                'tair_short_desc': cols[8],
                'go_descr': cols[9],
                'go_terms': go_terms,
                'mapman_bin': cols[11],
                'mapman_name': cols[12],
                'mapman_desc': cols[13],
                'pheno_aragwas_id': cols[14],
                'pheno_desc1': cols[15],
                'pheno_desc2': cols[16],
                'pheno_desc3': cols[17],
                'pheno_ref': cols[18],
                'user_notes': cols[19],
            }
            if node_type == 'gene':
                genes.append(doc)
            elif node_type == 'pheno':
                phenos.append(doc)
            else:
                raise RuntimeError(f"invalid node type {node_type}")
    save_docs(_PHENO_VERT_NAME, phenos)
    save_docs(_GENE_VERT_NAME, genes)


def load_cluster_data():
    """Annotate genes with cluster ID fields."""
    docs = []
    for path in _CLUSTER_PATHS:
        with open(path) as fd:
            csv_reader = csv.reader(fd, delimiter='\t')
            for row in csv_reader:
                cluster_label = row[0]
                gene_keys = row[1:]
                docs += [
                    {'_key': key, 'cluster': cluster_label}
                    for key in gene_keys
                ]
    save_docs(_GENE_VERT_NAME, docs)


def main():
    load_vert_metadata()
    load_pheno_assns()
    edge_paths = [
        (_GENE_COEXPR_PATH, 'gene_coexpr'),
        (_DOMAIN_CO_OCCUR_PATH, 'domain_co_occur'),
        (_PPI_HITHRU_PATH, 'ppi_hithru'),
        (_PPI_LIT_PATH, 'ppi_liter'),
    ]
    for (path, score_type) in edge_paths:
        load_edges(path, score_type)
    load_cluster_data()


def save_docs(coll_name, docs, on_dupe='update'):
    resp = requests.put(
        CONF.api_url + '/api/v1/documents',
        params={'collection': coll_name, 'on_duplicate': on_dupe},
        headers={'Authorization': CONF.auth_token},
        data='\n'.join(json.dumps(d) for d in docs)
    )
    if not resp.ok:
        raise RuntimeError(resp.text)
    else:
        print(f"Saved docs to collection {coll_name}!")
        print(resp.text)
        print("=" * 80)


if __name__ == '__main__':
    main()