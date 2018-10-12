import networkx as nx
import multiprocessing as mp
import MDAnalysis as mda

def parallel_read(file):
	h = nx.read_gml(file)
	return h


	
if __name__ == "__main__":	
	start = 0
	end = 1250
	#getting a dict of "nodes" (AA:chain:number)
	u = mda.Universe("md2.pdb")
	prot = u.select_atoms("protein")
	i = 0
	dict = {}
	for atom in prot:
		if atom.name == "CA":
			dict[str(i)] = atom.resname+":"+atom.segid+":"+str(atom.resid)
			i += 1
	i = start
	graphFiles = []
	while i <= end:
		graphFiles.append("MD-TASK-master/apoMD2_2-4_"+str(i)+"_graph.gml")
		i += 1
	#reading graphs in a parallel fashion
	pool=mp.Pool(processes=int(8)) #for multiprocessing			
	graphs=pool.map(parallel_read,(graphFiles))
	
	#generating consensus graphs
	G = nx.MultiGraph()
	for i in range(len((graphs))):
		for edges in graphs[i].edges():
			n1 = dict[edges[0]]
			n2 = dict[edges[1]]
			if G.has_edge(n1,n2) == True:
				G[n1][n2][0]["percentage"] = str(int(G[n1][n2][0]["percentage"])+1)
			else:
				G.add_edge(n1,n2, percentage="1")
	toDel = []			
	for edge in G.edges():
		n1 = edge[0]
		n2 = edge[1]
		crrntPercentage = int(G[n1][n2][0]["percentage"])
		crrntPercentage = (crrntPercentage/len(graphs))*100
		G[n1][n2][0]["percentage"] = str(crrntPercentage)
		if crrntPercentage < 75:
			toDel.append([n1,n2])
	for dele in toDel:
		G.remove_edge(dele[0],dele[1])

	nx.write_gml(G, "consensus_mdtask_"+str(start)+"_"+str(end)+".gml")
