#!/bin/bash
srr=$1
cut -f 1 ${srr}.sars.paf|fgrep -f - <(cut -f 1 ${srr}.genome.paf) -w |fgrep -f - <(cat ${srr}.sars.paf ${srr}.transcriptome.paf) -w
