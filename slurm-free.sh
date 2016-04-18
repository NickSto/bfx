#!/usr/bin/env bash
# original author: Nick Stoler
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

PARTITIONS=${PARTITIONS:=general}

Usage="Usage: \$ $(basename $0) [options]
Find cluster nodes that are available, given your job's needs and your wishes
about how many nodes/cpus you want to leave free. Useful in a loop like:
  node=\$($(basename $0) -f 4 -c 32)
  while [[ \$node ]]; do
    srun -w \$node your_command_here --threads 32 input.fa &
    sleep 60
    node=\$($(basename $0) -f 4 -c 32)
  done
Note: By default this is restricted to the \"$PARTITIONS\" partition(s).
Override by setting \$PARTITIONS to something else.
Options:
-f \$nodes: Leave \$nodes free.
-F \$cpus:  Leave \$cpus free.
-c \$cpus:  Job requires \$cpus. If not given, this will assume you need an
           entire node free.
-l \$cpus:  Return \"localhost\" if \$cpus will be left free there.
           Requires -c to be provided too.
-o \$path:  If this file exists, override the command line options with any that
           are present in the file. Format: one option and argument per line:
-f 2
-c 10"

function main {
  if [[ $# -lt 1 ]] || [[ $1 == '-h' ]]; then
    fail "$Usage"
  fi

  # Read in the arguments.
  free_nodes_leave=
  free_cpus_leave=
  cpus_needed=
  local_cpus_leave=
  override_file=
  while getopts ":f:F:c:l:o:h" opt; do
  case "$opt" in
      f) free_nodes_leave=$OPTARG;;
      F) free_cpus_leave=$OPTARG;;
      c) cpus_needed=$OPTARG;;
      l) local_cpus_leave=$OPTARG;;
      o) override_file=$OPTARG;;
      h) echo "$USAGE"
         exit;;
    esac
  done

  # Check arguments.
  if [[ $local_cpus_leave ]] && ! [[ $cpus_needed ]]; then
    fail "Error: If providing -l, must also give -c."
  fi
  if [[ $free_cpus_leave ]] && ! [[ $cpus_needed ]]; then
    fail "Error: If providing -F, must also give -c."
  fi

  # Read in arguments from the override file.
  if [[ $override_file ]] && [[ -f $override_file ]]; then
    while read opt arg; do
      case "$opt" in
        -f) free_nodes_leave=$arg;;
        -F) free_cpus_leave=$arg;;
        -c) cpus_needed=$arg;;
        -l) local_cpus_leave=$arg;;
      esac
    done < $override_file
  fi

  partition_arg=
  if [[ $PARTITIONS ]]; then
    partition_arg="-p $PARTITIONS"
  fi

  # Examine the node usage.
  used_cpus_total=0
  free_cpus_total=0
  free_nodes=0
  candidates=
  while read node used_cpus free_cpus x total_cpus; do
    used_cpus_total=$((used_cpus_total+used_cpus))
    free_cpus_total=$((free_cpus_total+free_cpus))
    if [[ $cpus_needed ]]; then
      if [[ $free_cpus -lt $cpus_needed ]]; then
        candidates="$candidates $node"
      fi
    elif [[ $free_cpus == $total_cpus ]]; then
      candidates="$candidates $node"
    fi
    if [[ $free_cpus == $total_cpus ]]; then
      free_nodes=$((free_nodes+1))
    fi
  done < <(sinfo -h $partition_arg -t idle,alloc -o '%n %C' | tr ' /' '\t\t')

  # Check if there will be enough nodes or cpus after launching the job.
  free_node=
  if [[ $free_nodes_leave ]]; then
    # Subtract 1 for the job that will be launched.
    if [[ $((free_nodes-1)) -ge $free_nodes_leave ]]; then
      free_node=true
    fi
  elif [[ $free_cpus_leave ]]; then
    if [[ $((free_cpus_total-cpus_needed)) -ge $free_cpus_leave ]]; then
      free_node=true
    fi
  else
    fail "Error: Must provide -f or -F."
  fi

  # If there is a free node, then print one from the list and exit.
  if [[ $free_node ]]; then
    echo "$candidates" | tr ' ' '\n' | tail -n 1
    return
  fi

  # Otherwise, see if there is room on the localhost.
  if [[ $local_cpus_leave ]] && [[ $cpus_needed ]]; then
    total_local_cpus=$(grep -c 'core id' /proc/cpuinfo)
    used_local_cpus=$(ps aux | awk '{tot+=$3} END {print int(tot/100)}')
    free_local_cpus=$((total_local_cpus-used_local_cpus))
    if [[ $free_local_cpus -ge $local_cpus_leave ]]; then
      echo localhost
    fi
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
