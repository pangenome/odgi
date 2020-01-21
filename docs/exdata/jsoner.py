import sys
import json
sys.path.append("../lib")
import odgi

# invoke with $python3 jsoner.py [file name without suffix]
file_name = sys.argv[1]

json_file = open(file_name + ".json")
dat = json.load(json_file)
gr = odgi.graph()
handles = dict()

for node in dat['node']:
    id = node['id']
    seq = node['sequence']
    handles[id] = gr.create_handle(seq)

for edge in dat['edge']:
    fr = edge['from']
    to = edge['to']
    gr.create_edge(handles[fr], handles[to])

for path in dat['path']:
    name = path['name']
    path['mapping'].sort(key = lambda x: x['rank'])
    path_handle = gr.create_path_handle(name)
    for mapping in path['mapping']:
        id = mapping['position']['node_id']
        gr.append_step(path_handle, handles[id])
gr.serialize(file_name + ".odgi")

# sanity check to make sure saved file is same after reading in

# gr2 = odgi.graph()
# gr2.load(file_name + ".odgi")

#sequence = []
#gr.for_each_path_handle(lambda x: gr.for_each_step_in_path(x, lambda y: sequence.append(gr.get_sequence(gr.get_handle_of_step(y)))))
#sequence2 = []
#gr2.for_each_path_handle(lambda x: gr2.for_each_step_in_path(x, lambda y: sequence2.append(gr2.get_sequence(gr2.get_handle_of_step(y)))))
#print("".join(sequence) == "".join(sequence2) )
