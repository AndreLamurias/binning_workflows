# use edges list of assembly graph and GAF file mapping contigs to that graph to generate a graph of contigs
import sys
import re
import argparse
import networkx as nx
import gfapy
import itertools


def get_new_edges_dict(repeat_edges):
    new_edges = {}  # start: [list of starting nodes], end: [list of ending nodes]
    for e in repeat_edges:
        if repeat_edges[e][2] not in new_edges:  # new node at end
            new_edges[repeat_edges[e][2]] = {"start": [], "end": []}
        if repeat_edges[e][3] not in new_edges:  # new node at start
            new_edges[repeat_edges[e][3]] = {"start": [], "end": []}
        new_edges[repeat_edges[e][2]]["end"].append(e)  # outgoing edge from new edge
        new_edges[repeat_edges[e][3]]["start"].append(e)  # incoming edge to new edge
    return new_edges


def generate_graph_from_edges(repeat_edges):
    """
    Generate graph based on repeat graph where the edges are nodes and vice versa
    repeat_edges: edge_id -> (len, cov, start, end)
    """
    # add all edges to graph as nodes
    G = nx.Graph()
    for e in repeat_edges:
        G.add_node(e)
    # if to edges e1 e2 have e1.end = e2.start, create edge
    # import pdb; pdb.set_trace()
    new_edges = get_new_edges_dict(repeat_edges)
    for new_edge in new_edges:
        for n1 in new_edges[new_edge]["start"]:
            for n2 in new_edges[new_edge]["end"]:
                G.add_edge(n1, n2)

    return G


def parse_flye_dot(dot_fpath, min_edge_len):
    """
    Read repeat graph from GV file
    """
    G = nx.Graph()
    dict_edges = dict()  # edge ID -> edge properties
    # contig_starts = {} #node-> edge
    # contig_ends = {} #node <- edge
    pattern = '"?(?P<start>\d+)"? -> "?(?P<end>\d+)"? \[(?P<info>.+)]'
    label_pattern = "id (?P<edge_id>\-*.+) (?P<edge_len>[0-9\.]+)k (?P<coverage>\d+)"
    with open(dot_fpath) as f:
        for line in f:
            if "label =" in line:
                # "7" -> "29" [label = "id 1\l53k 59x", color = "black"] ;
                line = line.replace("\\l", " ")
                match = re.search(pattern, line)
                if not match or len(match.groups()) < 3:
                    continue
                start, end, info = match.group("start"), match.group("end"), match.group("info")
                params_dict = dict(param.split(" = ") for param in info.split(", ") if "=" in param)
                # label = params_dict.get('label')
                color = params_dict.get("color").strip().replace('"', "")
                line = line.replace(" ,", ",")
                match = re.search(label_pattern, info)
                if match and match.group("edge_id"):
                    # edge_id = get_edge_agv_id(match.group('edge_id'))
                    edge_id = "edge_{}".format(abs(int(match.group("edge_id"))))
                    cov = max(1, int(match.group("coverage")))
                    # import pdb; pdb.set_trace()
                    edge_len = max(1, int(float(match.group("edge_len")) * 1000))
                    if edge_len < min_edge_len:
                        continue
                    edge = [edge_len, cov, abs(int(start)), abs(int(end))]
                    # edge.color = color
                    # if edge.color != "black":
                    #    edge.repetitive = True
                    # edge.start, edge.end = int(start), int(end)
                    # if 'dir = both' in line:
                    #    edge.two_way = True
                    dict_edges[edge_id] = edge
                    # G.add_node(edge_id)
                    # contig_starts[start].append(edge_id)

    # dict_edges = calculate_multiplicities(dict_edges)
    return dict_edges


def read_gfa_file(input_dir):
    graph = gfapy.Gfa.from_file(input_dir)
    contigs = {}
    singletons = []
    G = nx.Graph()

    for contig in graph.segments:
        # if contig.from == "tig00000026":
        #    import pdb; pdb.set_trace()
        if len(contig.edges) == 0:
            singletons.append(contig.name)
            # continue
        # contigs[contig.name] = [c.from_name for c in contig.edges if c.from_name != contig.name] + [c.to_name for c in contig.edges if c.to_name != contig.name]
        # contigs[contig.name] = [c.from_name for c in contig.edges] + [c.to_name for c in contig.edges]
        contigs[contig.name] = {}
        G.add_node(contig.name)
        for edge in contig.edges:
            # this contig could be on either end
            import pdb

            pdb.set_trace()
            if edge.from_name == contig.name:
                other_edge = edge.to_name
            else:
                other_edge = edge.from_name
            G.add_edge(contig.name, other_edge)
            overlap = sum([a.length for a in edge.alignment if a.code == "M"])
            # total = sum([a.length for a in edge.alignment])
            if other_edge not in contigs[contig]:
                contigs[contig][other_edge] = []
            contigs[contig][other_edge].append(overlap)
        if len(contigs[contig.name]) == 0:
            del contigs[contig.name]
            singletons.append(contig.name)

    # import pdb; pdb.set_trace()
    print("{} linked contigs, {} singletons".format(len(contigs), len(singletons)))
    print("percent linked contigs: {}".format(len(contigs) / (len(contigs) + len(singletons))))
    print("avg num of links", sum([len(edges) for edges in contigs.values()]) / len(contigs))
    print("graph density", nx.density(G))
    return G


def simplify_graph(G, second_order=True):
    # iterate nodes, if node does not start with contig, combine all edges to contigs and delete node
    # import pdb; pdb.set_trace()
    all_nodes = list(G.nodes())
    for node in all_nodes:
        if "edge_" in node:
            linked_contigs = set([n for n in G.neighbors(node) if isinstance(n, str) and not n.startswith("edge_")])
            for pair in itertools.combinations(linked_contigs, 2):
                G.add_edge(pair[0], pair[1])
            # explore second order edges: if to edge nodes are linked, contigs linked to those edges should be linked
            if second_order:
                linked_edges = [n for n in G.neighbors(node) if isinstance(n, str) and n.startswith("edge_")]
                second_order_contigs = set()
                for e in linked_edges:
                    second_order_contigs.update(
                        [n for n in G.neighbors(e) if isinstance(n, str) and not n.startswith("edge_")]
                    )
                for pair in itertools.product(linked_contigs, second_order_contigs):
                    G.add_edge(pair[0], pair[1])
                # G.remove_node(node) keep node for now, to calculate second order
    for node in all_nodes:
        if "edge_" in node:
            G.remove_node(node)
    return G


def main():
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument("--edges", type=argparse.FileType("rt"), help="Use edges list as graph input")
    parser.add_argument("--gfa", help="Use GFA as graph input")
    parser.add_argument("--gaf", help="GAF file with contig alignments to edges of graph")
    parser.add_argument("--gv", help="GV graph as input")
    parser.add_argument("--minoverlap", help="minoverlap", default=0, type=int)
    parser.add_argument("--output", help="Output name")
    parser.add_argument("--simplify", help="Simplify graph to keep only edges between contigs", action="store_true")
    parser.add_argument("--second", help="Create indirect contig edges", action="store_true")
    args = parser.parse_args()
    # read graph
    if args.gfa:
        print("reading GFA")
        graph = read_gfa_file(args.gfa)

    elif args.edges:
        print("reading edges file")
        graph = nx.read_edgelist(args.edges)

    elif args.gv:
        edges = parse_flye_dot(args.gv, 0)
        print("Loading GV:", len(edges))
        graph = generate_graph_from_edges(edges)

    # read alignments
    # 1 	string 	Query sequence name
    # 2 	int 	Query sequence length
    # 3 	int 	Query start (0-based; closed)
    # 4 	int 	Query end (0-based; open)
    # 5 	char 	Strand relative to the path: "+" or "-"
    # 6 	string 	Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
    # 7 	int 	Path length
    # 8 	int 	Start position on the path (0-based)
    # 9 	int 	End position on the path (0-based)
    # 10 	int 	Number of residue matches
    # 11 	int 	Alignment block length
    # 12 	int 	Mapping quality (0-255; 255 for missing)
    print("Adding edges between contigs and assembly edges")
    overlaps = []
    with open(args.gaf, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            contig_name = values[0]
            path = values[5]
            match_length = int(values[9])
            overlaps.append(match_length)
            if match_length < args.minoverlap:
                continue
            # assume we use everything, add edges to graph
            edges = re.split(r"<|>", path)
            edges = [e for e in edges if len(e) > 0]
            for e in edges:
                # if e not in graph:
                #    import pdb; pdb.set_trace()
                graph.add_edge(contig_name, e)

    # simplify graph to keep only contigs
    if args.simplify:
        print("simplifying graph")
        graph = simplify_graph(graph, args.second)
    # TODO: break long disconnected contigs
    # write to file
    print("writing to file")
    nx.write_edgelist(graph, args.output)
    print(args.output)
    print("graph nodes", len(graph.nodes))
    print("graph edges", len(graph.edges))
    print("graph density:", nx.density(graph))

    # TODO plot edges by minverlap


if __name__ == "__main__":
    main()
