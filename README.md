SR-CFSS: a Branch and Bound algorithm for Social Ridesharing
===================
SR-CFSS is a Branch and Bound algorithm for Social Ridesharing. SR-CFSS has been presented by Filippo Bistaffa, Alessandro Farinelli, and Sarvapali D. Ramchurn in “[Sharing Rides with Friends: a Coalition Formation Algorithm for Ridesharing](http://www.aaai.org/ocs/index.php/AAAI/AAAI15/paper/download/9622/9303)”, Proceedings of the 2015 AAAI Conference on Artificial Intelligence (AAAI), pages 608–614, 2015, AAAI.

Download
----------
CFSS requires `g++` to compile, and does not require any external library to execute. In order to employ Twitter as network topology, `java` must be installed on the system, and the [Twitter GitHub repository](https://github.com/filippobistaffa/twitter) must be `git clone`'d inside SR-CFSS's root directory.

Execution
----------
SR-CFSS must be executed by means of the [`sr.sh`](https://github.com/filippobistaffa/SR-CFSS/blob/master/sr.sh) script, i.e.,
```
./sr.sh -t <scalefree|twitter> -n <#agents> -s <seed>

-t	Network topology (either scalefree or twitter)
-n	Number of agents
-s	Seed
-m	Parameter m of the Barabasi-Albert model (optional, default m = 1)
```
