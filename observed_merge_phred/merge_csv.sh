#!/bin/bash

shopt -s nullglob
for d in output/evaluation/*; do
  awk 'NR==1; FNR==1{next} 1' "$d"/*.csv > "$d"_merged.csv
done
awk 'NR==1; FNR==1{next} 1' output/evaluation/*merged.csv > output/evaluation/all_merged.csv