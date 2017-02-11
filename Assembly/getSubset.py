#!/usr/bin/env python

#	Script: getSubset.py
#	Author: Mauro Truglio
#	Massey University, New Zealand
#	Email contact <mauro.truglio@gmail.com>
#	November 2015
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   A copy of the GNU General Public License is avilable at <http://www.gnu.org/licenses/>.


#   Usage: getSubset.py [forward reads] [reverse reads] [subset size]
#   E.g.:	getSubset.py forward.fastq reverse.fastq 10000

import random
import sys

def write_random_records(fqa, fqb, N=100000):
    """ get N random headers from a fastq file without reading the
    whole thing into memory"""
    records = sum(1 for _ in open(fqa)) / 4
    rand_records = sorted([random.randint(0, records - 1) for _ in xrange(N)])

    fha, fhb = open(fqa),  open(fqb)
    suba, subb = open(fqa + ".subset", "w"), open(fqb + ".subset", "w")
    rec_no = - 1
    for rr in rand_records:

        while rec_no < rr:
            rec_no += 1       
            for i in range(4): fha.readline()
            for i in range(4): fhb.readline()
        for i in range(4):
            suba.write(fha.readline())
            subb.write(fhb.readline())
        rec_no += 1 # (thanks @anderwo)

    print >>sys.stderr, "wrote to %s, %s" % (suba.name, subb.name)

if __name__ == "__main__":
    N = 100 if len(sys.argv) < 4 else int(sys.argv[3])
    write_random_records(sys.argv[1], sys.argv[2], N)

