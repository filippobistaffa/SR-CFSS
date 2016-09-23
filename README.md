SR-CFSS: a Branch and Bound algorithm for Social Ridesharing
===================
SR-CFSS is a Branch and Bound algorithm for Social Ridesharing. SR-CFSS has been presented by Filippo Bistaffa, Alessandro Farinelli, and Sarvapali D. Ramchurn in “[Sharing Rides with Friends: a Coalition Formation Algorithm for Ridesharing](http://www.aaai.org/ocs/index.php/AAAI/AAAI15/paper/download/9622/9303)”, Proceedings of the 2015 AAAI Conference on Artificial Intelligence (AAAI), pages 608–614, 2015, AAAI.

Requirements
----------
SR-CFSS requires `g++` to compile, and does not require any external library to execute. In order to employ Twitter as network topology, `java` must be installed on the system, and the [Twitter GitHub repository](https://github.com/filippobistaffa/twitter) must be `git clone`'d inside SR-CFSS's root directory.

Execution
----------
SR-CFSS must be executed by means of the [`sr.sh`](https://github.com/filippobistaffa/SR-CFSS/blob/bound/sr.sh) script, i.e.,
```
./sr.sh -t <scalefree|twitter> -n <#agents> -s <seed> [-m <barabasi_m>] [-d <drivers_%>]

-t	Network topology (either scalefree or twitter)
-n	Number of agents
-s	Seed
-d	Drivers' percentage (optional, default d = 20)
-m	Parameter m of the Barabasi-Albert model (optional, default m = 2)
```

Acknowledgements
----------
SR-CFSS employs the [GeoLife dataset by Microsoft Research](http://research.microsoft.com/en-us/projects/geolife) presented by Yu Zheng, Quannan Li, Yukun Chen, Xing Xie, and Wei-Ying Ma in “[Understanding mobility based on GPS data](https://www.microsoft.com/en-us/research/publication/understanding-mobility-based-on-gps-data)”, Proceedings of the 10th ACM conference on Ubiquitous Computing (Ubicomp), pages 312−321, 2008, ACM press.
