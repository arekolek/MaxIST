MaxIST
======

To run:

```bash
make GDB=off
./eval -z 50 -n 100
```

* -z controls the number of test cases
* -t is the graph model type (gnp or rgg with a path or mst, e.g. gnp+mst)
* -n is the size of the graph (or comma-separated list of sizes)
* -p is the parameter for the model (or comma-separated list of parameters)
* -c is the list of construction algorithms to use (bfs, dfs, fifo, ilst, rdfs, rdfs50, random)
* -i is the list of improvement algorithms to use (none, prieto, lost-light, lost, lost-ex)


License
=======

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
