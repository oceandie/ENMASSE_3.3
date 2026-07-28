#!/usr/bin/env python
"""
Version: 0.2

Create a bind list for MPI rank placement.
See `create_bind_list.py -h`
"""

import argparse
import math
import sys

def msg(name=None):
    return '''
  command line: create_bind_list.py <ppn> <cores> [threads])
  mpiexec     : --cpu-bind=list:$(bind_list.py <ppn> <cores> [threads])

The command will create a bind list that tries to spread the processes and
threads approximately evenly across the available cores and chiplets
(i.e. the algorithm is chiplet aware).

The command has two required positional arguments: (ppn) process per
nodes and (cores) maximum number of cores per node. The third
positional argument is the number of threads and it is not required
(default is 1).

The number of cores used should be 128 (Milan node) or 192 (Genoa
node). If a different value is needed the optional argument
--any-cores should be used.

Examples:
  ```
  > create_bind_list.py 8 128
  0:16:32:48:64:80:96:112

  > create_bind_list.py 6 128 2
  0,1:16,17:32,33:64,65:80,81:96,97
  ```
          '''

def createBindList(bindList: str, ppn: int, cores: int, threads:int, startRank: int) -> str:

    thisProc=0.0
    rank=0
    ranksAllocated=0
    # Calculate step between ranks
    procStep = 1.0 * cores / ppn

    while ranksAllocated < ppn :
        # Add a semicolon if the bind list is not empty
        if len(bindList) != 0:
            bindList+=":"

        # Add current rank and threads
        bindList = bindList + str(startRank+rank)
        for t in range(1,threads):
            bindList = bindList + "," + str(startRank+rank+t)

        # Calculate next rank
        thisProc = thisProc + procStep
        rank = math.floor(thisProc)

        ranksAllocated = ranksAllocated + 1

    return bindList

def createChipletsBindList(ppn: int, cores: int, threads:int, debug: bool) -> str:

    # Always allocate first process to core 0
    # TODO:   Add the possibility to not allocate the first process to core 0

    bindList=""

    # Calculate step between ranks
    procStep = 1.0 * cores / ppn

    # If the proc step is a full integer or the number of cores is not
    # divisble by 8 or the number of cores is less than 8, then the
    # simple algorithm works fine
    if procStep.is_integer() or cores % 8 != 0 or cores<=8:

        bindList = createBindList(bindList, ppn, cores, threads, 0)

    else:
        # The following algorithm is chiplet aware

        # Number of chiplets available for the number of cores
        chipletsAvail = math.ceil(cores/8)

        # Assume that there are going to be two sockets on a node. So
        # halve some of the variables.
        halfChipletsAvail = math.ceil(chipletsAvail/2)
        halfPPN           = math.ceil(ppn/2)

        # Compute a list of process per chiplet needed, one per
        # chiplet. The process are going to be divided as equally as
        # possible between two sockets.

        chiplet=0
        ppc: List[int] = []
        ppcBase = math.floor(halfPPN/halfChipletsAvail)
        ppcRem  = halfPPN - (ppcBase * halfChipletsAvail)
        remStep = math.floor(halfChipletsAvail / ppcRem) if ppcRem else 0

        for s in range(2):

            while chiplet < halfChipletsAvail :

                ppcToAdd = ppcBase

                # Devide the remainder process equally if possible.
                if ppcRem > 0 and chiplet % remStep == 0 :
                    ppcToAdd += 1
                    ppcRem   -= 1

                ppc.append(ppcToAdd)
                chiplet += 1

            # Recompute variable for next cycle
            chiplet           = 0
            halfChipletsAvail = chipletsAvail - halfChipletsAvail

            if(halfChipletsAvail==0):
                break

            halfPPN = ppn - halfPPN
            ppcBase = math.floor(halfPPN/halfChipletsAvail)
            ppcRem  = halfPPN - (ppcBase * halfChipletsAvail)
            remStep = math.floor(halfChipletsAvail / ppcRem) if ppcRem else 0

        if debug :
            print("[INFO] chiplets available  : %d" % chipletsAvail)
            print("[INFO] process per chiplet :", ppc)

        if sum(ppc) != ppn :
            print(
                "[FAIL] Internal error, sum of ppc is not the same of ppn: ", sum(ppc)
            )
            sys.exit(1)


        chipsAllocated=0
        startRank = 0

        while chipsAllocated < chipletsAvail :

            if ppc[chipsAllocated] > 0 :
                bindList = createBindList(bindList, ppc[chipsAllocated], 8, threads, startRank)

            # Calculate next chiplet
            chipsAllocated = chipsAllocated + 1
            startRank = chipsAllocated*8
            if debug :
                bindList+=" "

    return bindList

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="",
                                     usage=msg())

    # Positional required
    parser.add_argument("ppn",     help="process per node", type=int)
    parser.add_argument("cores",   help="cores per node (Milan:128 Genoa:192)", type=int)

    # Positional optional
    parser.add_argument("threads", help="threads used", type=int, nargs='?',default=1)

    # Optional
    parser.add_argument("-d",   "--debug",   help="show more information",
                        required=False, action='store_true')
    parser.add_argument("--any-cores",       help="allow the number of cores to be different than 128 or 192",
                        required=False, action='store_true',default=False)

    args = parser.parse_args()

    if args.debug :
        print("[INFO] ppn                 : %d" % args.ppn)
        print("[INFO] cores               : %d" % args.cores)
        print("[INFO] threads             : %d" % args.threads)

    if  args.ppn > args.cores :
        print(
            "[FAIL] Process per node cannot be more than available cores"
        )
        sys.exit(1)

    if not args.any_cores and args.cores != 128 and args.cores != 192 :
        print(
            "[FAIL] Number of cores (%d) is not a Milan or Genoa node" % args.cores
        )
        sys.exit(1)

    if args.ppn * args.threads > args.cores:
        print(
            "[FAIL] Number of process per nodes * threads (%d) is "
            "higher than the maximum cores requested (%d)" % (args.ppn*args.threads,args.cores)
        )
        sys.exit(1)


    # Generate the bind list on stdout.
    print ( createChipletsBindList (args.ppn, args.cores, args.threads, args.debug) )


