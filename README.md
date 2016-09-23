SR-CFSS: a Branch and Bound algorithm for Social Ridesharing
===================
SR-CFSS is a Branch and Bound algorithm for Social Ridesharing. SR-CFSS has been presented by Filippo Bistaffa, Alessandro Farinelli, and Sarvapali D. Ramchurn in “[Sharing Rides with Friends: a Coalition Formation Algorithm for Ridesharing](http://www.aaai.org/ocs/index.php/AAAI/AAAI15/paper/download/9622/9303)”, Proceedings of the 2015 AAAI Conference on Artificial Intelligence (AAAI), pages 608–614, 2015, AAAI.

Requirements
----------
SR-CFSS requires `g++` to compile, and does not require any external library to execute. In order to employ Twitter as network topology, `java` must be installed on the system, and the [Twitter GitHub repository](https://github.com/filippobistaffa/twitter) must be `git clone`'d inside SR-CFSS's root directory.

Execution
----------
SR-CFSS must be executed by means of the [`sr.sh`](https://github.com/filippobistaffa/SR-CFSS/blob/master/sr.sh) script, i.e.,
```
./sr.sh -t <scalefree|twitter> -n <#agents> -s <seed> [-m <barabasi_m>] [-d <drivers_%>] [-p <output>]

-t	Network topology (either scalefree or twitter)
-n	Number of agents
-s	Seed
-d	Drivers' percentage (optional, default d = 20)
-m	Parameter m of the Barabasi-Albert model (optional, default m = 2)
-p	Outputs a solution file formatted for PK
```

Search Tree Render
----------
SR-CFSS can generate a [DOT](http://www.graphviz.org/content/dot-language) file that represents the search tree explored during the execution. In order to generate such file, `#define` the `TREEDOT` constant as the path of the output DOT file inside [`params.h`](https://github.com/filippobistaffa/SR-CFSS/blob/master/params.h), e.g.,
```
#define TREEDOT "tree.dot"
```
To render such file as a PNG image, `dot` (part of the `graphviz` suite) must be installed in the system. Then, execute
```
unflatten -f -l3 TREEDOT | dot -Tpng -o PNGFILE
```
where `TREEDOT` is the generated DOT file and `PNGFILE` is the desired output PNG file.
![Example search tree](http://i.imgur.com/mZNSg62.png)
The ID of each coalition structure (i.e., the labels of the nodes) are printed during the execution of SR-CFSS. The optimal solution is highlighted in green.

Acknowledgements
----------
SR-CFSS employs the [GeoLife dataset by Microsoft Research](http://research.microsoft.com/en-us/projects/geolife) presented by Yu Zheng, Quannan Li, Yukun Chen, Xing Xie, and Wei-Ying Ma in “[Understanding mobility based on GPS data](https://www.microsoft.com/en-us/research/publication/understanding-mobility-based-on-gps-data)”, Proceedings of the 10th ACM conference on Ubiquitous Computing (Ubicomp), pages 312−321, 2008, ACM press.
