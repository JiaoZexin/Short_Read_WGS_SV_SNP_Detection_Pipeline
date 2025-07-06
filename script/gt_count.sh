#!/bin/bash
# 用法: gt_count.sh input_file out_missing out_00 out_01 out_11

input="$1"
out_missing="$2"
out_00="$3"
out_01="$4"
out_11="$5"

while read line; do echo "$line" | grep -o '\./\.' | wc -l; done < "$input" > "$out_missing"
while read line; do echo "$line" | grep -o '0/0'   | wc -l; done < "$input" > "$out_00"
while read line; do echo "$line" | grep -o '0/1'   | wc -l; done < "$input" > "$out_01"
while read line; do echo "$line" | grep -o '1/1'   | wc -l; done < "$input" > "$out_11"
