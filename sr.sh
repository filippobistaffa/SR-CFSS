#!/usr/bin/env bash

red='\033[0;31m'			# Red
nc='\033[0m'				# No color
re='^[0-9]+$'				# Regular expression to detect natural numbers
d=10					# Default d = 10
m=1					# Default m = 1
basename="twitter/net/twitter-2010"	# Twitter files must be in "twitter/net" subdir and have twitter-2010.* filenames
wg="twitter/wg"				# Place WebGraph libraries in "twitter/wg" subdir

usage() { echo -e "Usage: $0 -t <scalefree|twitter> -n <#agents> -s <seed> [-m <barabasi_m>] [-d <drivers_%>]\n-t\tNetwork topology (either scalefree or twitter)\n-n\tNumber of agents\n-s\tSeed\n-d\tDrivers' percentage (optional, default d = 10)\n-m\tParameter m of the Barabasi-Albert model (optional, default m = 1)" 1>&2; exit 1; }

while getopts ":t:n:s:d:m:" o; do
	case "${o}" in
	t)
		t=${OPTARG}
		if ! [[ $t == "scalefree" || $t == "twitter" ]] ; then
			echo -e "${red}Network topology must be either scalefree or twitter!${nc}\n"
			usage
		fi
		;;
	n)
		n=${OPTARG}
		if ! [[ $n =~ $re ]] ; then
			echo -e "${red}Number of agents must be a number!${nc}\n"
			usage
		fi
		;;
	s)
		s=${OPTARG}
		if ! [[ $s =~ $re ]] ; then
			echo -e "${red}Seed must be a number!${nc}\n"
			usage
		fi
		;;
	d)
		d=${OPTARG}
		if ! [[ $d =~ $re ]] ; then
			echo -e "${red}Drivers' percentage must be a number!${nc}\n"
			usage
		fi
		;;
	m)
		m=${OPTARG}
		if ! [[ $m =~ $re ]] ; then
			echo -e "${red}Parameter m must be a number!${nc}\n"
			usage
		fi
		;;
	\?)
		echo -e "${red}-$OPTARG is not a valid option!${nc}\n"
		usage
		;;
	esac
done
shift $((OPTIND-1))

if [ -z "${t}" ] || [ -z "${n}" ] || [ -z "${s}" ]; then
	echo -e "${red}Missing one or more required options!${nc}\n"
	usage
fi

echo -e "[\033[01;33m CC \033[0m] SR-CFSS"

if [[ $t == "scalefree" ]] ; then
	g++ -DN=$n -DK=$m -DSEED=$s -DDRIVERPERC=$d -Wall -march=native -O0 -funroll-loops -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 *.c *.cpp -lm -o sr
else
	tmp=`mktemp`
	java -Xmx4000m -cp .:$wg/* ReduceGraph $basename $n $s > $tmp
	g++ -DTWITTER -Wall -march=native -Ofast -funroll-loops -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 -include types.h -include $tmp *.c *.cpp -lm -o sr
	rm $tmp
fi

./sr
#echo "" >> original.log
