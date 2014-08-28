public class ReduceGraph {

	@SuppressWarnings("unchecked")
	static public void main( String arg[] ) throws IllegalArgumentException, SecurityException, IllegalAccessException, NoSuchMethodException, java.io.IOException {

		java.util.Random r = new java.util.Random(Integer.parseInt(arg[2]));
		final it.unimi.dsi.webgraph.ImmutableGraph graph = it.unimi.dsi.webgraph.ImmutableGraph.load(arg[0]);
		int n = Integer.parseInt(arg[1]);
		int nodes = graph.numNodes();
		long arcs = graph.numArcs();
		int startNode = r.nextInt(nodes) + 1;
		int i = 0;
		int e = 0;
		int[] adj = new int[n * n];
		java.util.ArrayList<Integer> al = new java.util.ArrayList<Integer>(n);
		al.add(startNode);

		while (i < n) {
			int[] nbrs = graph.successorArray(al.get(i));
			for (int j = 0; j < nbrs.length; j++) {
				if (!al.contains(nbrs[j])) {
					if (al.size() < n) al.add(nbrs[j]);
					else continue;
				}
				int k = al.indexOf(nbrs[j]);
				if (adj[i * n + k] == 0) adj[i * n + k] = adj[k * n + i] = ++e;
			}
			i++;
		}

		String adjs = java.util.Arrays.toString(adj).replace('[', '{').replace(']', '}');
		System.out.println("#define N " + n + "\n#define E " + e + "\nstatic const value graph[" + n * n + "] = " + adjs + ";");
	}
}
