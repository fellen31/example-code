#!/bin/bash

# This file needs to be created from config/config_template
# One copy on corbicula, one on UPPMAX.
source config/config.ini

raw_krill_dir=$PROJ_DIR/data/raw/krill
outgroup_dir=$PROJ_DIR/data/raw/outgroups
longest_dir=$PROJ_DIR/data/interim/longest_isoforms
pfam_dir=$PROJ_DIR/data/interim/pfam
sprot_dir=$PROJ_DIR/data/interim/sprot
transdecoder_dir=$PROJ_DIR/data/interim/transdecoder
mnor_dir=$PROJ_DIR/data/interim/mnor_model
gene_list_dir=$PROJ_DIR/data/interim/gene_list
proteinortho_dir=$PROJ_DIR/data/interim/proteinortho
orthofinder_dir=$PROJ_DIR/data/interim/orthofinder
ogs_dir=$PROJ_DIR/data/interim/orthogroups
alignments_dir=$PROJ_DIR/data/interim/alignments
busco_dir=$PROJ_DIR/data/results/busco
trees_dir=$PROJ_DIR/data/interim/trees
pruners_dir=$PROJ_DIR/data/interim/pruned
trim_dir=$PROJ_DIR/data/interim/trimmed
paltonal_dir=$PROJ_DIR/data/interim/pal2nal
dnds_dir=$PROJ_DIR/data/interim/dnds
logs_dir=$PROJ_DIR/logs
counts_dir=$PROJ_DIR/data/results/counts

mkdir -p $counts_dir

init_busco () {
    echo "placeholder"
}

run_busco () {

        # Needs to be initialized once to download files...
        busco_db=$PROJ_DIR/data/external/busco_downloads

        step=$1

        mkdir -p $busco_dir/$step && cd $busco_dir/$step

        parallel "busco \
            -i {} \
            -o {/} \
            --out_path $busco_dir/$step \
            -m transcriptome \
            --download_path $busco_db \
            --offline \
            -l arthropoda_odb10 \
            -c 1 \
            " ::: $PROJ_DIR/data/interim/$step/*.$2
}

# Selecting isoform based on expression not ideal (src/scripts/dead-ends/)

extract_longest_isoforms () {
    
    mkdir -p $longest_dir

    # get_longest_isoform_seq_per_trinity_gene.pl in.fasta > out.fasta
    parallel --progress \
        "$CONDA_PREFIX/opt/trinity-2.11.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
        {} > $longest_dir/{/.}.longest_isoforms.fasta \
        " ::: $raw_krill_dir/*.fasta
    
    cd $longest_dir
    
    grep -r -c ">" *.fasta |\
    awk -v FS=":" '{if (NR > 0) print($1,$2);}' > \
    $PROJ_DIR/data/results/counts/longest_isoforms.txt
}

# Get open reading frames (?)

get_orfs () {
    
    mkdir -p $transdecoder_dir
    cd $transdecoder_dir

    # TransDecoder.LongOrfs -t in.fasta -O out_dir
    parallel --progress \
        "TransDecoder.LongOrfs \
        -t {} \
        -O {/.} \
        " ::: $longest_dir/*.fasta

    grep -r -c ">" **/longest_orfs.pep |\
    awk -v FS=":" '{if (NR > 0) print($1,$2);}' > \
    $PROJ_DIR/data/results/counts/get_orfs.txt

}

search_sprot () {
    
    # Copied from corbicula
    uniprot_db=$PROJ_DIR/data/external/protein_dbs/uniprot_sprot.fasta

    mkdir -p $sprot_dir

    parallel --progress \
        "blastp \
            -query $transdecoder_dir/{/.}/longest_orfs.pep \
            -db $uniprot_db \
            -max_target_seqs 1 \
            -outfmt 6 \
            -evalue 1e-5 \
            -num_threads 2 \
            > $sprot_dir/{/.}.longest_orfs.pep.blastp.outfmt6 \
        " ::: $longest_dir/*.fasta

}

search_pfam () {
  
    # Copied from corbicula
    pfam_db=$PROJ_DIR/data/external/protein_dbs/Pfam-A.hmm
    
    mkdir -p $pfam_dir

    parallel --progress \
        "hmmscan \
            --cpu 1 \
            --domtblout $pfam_dir/{/.}.longest_orfs.pep.pfam.domtblout \
            $pfam_db \
            $transdecoder_dir/{/.}/longest_orfs.pep \
        " ::: $longest_dir/*.fasta

}

# Predict transcripts based on 

predict_transcripts () {

    parallel --progress \
        "cd $transdecoder_dir/{/.}; \
        TransDecoder.Predict \
            -t {} \
            --output_dir $transdecoder_dir/{/.} \
            --single_best_only \
            --retain_pfam_hits $pfam_dir/{/.}.longest_orfs.pep.pfam.domtblout \
            --retain_blastp_hits $sprot_dir/{/.}.longest_orfs.pep.blastp.outfmt6 \
        " ::: $longest_dir/*.fasta
    
}

mnor_model () {

    mnor_db=$PROJ_DIR/data/external/longest_orfs.pep
    
    mkdir -p $mnor_dir

    # Run transdecoder predicted transcripts against mnor gene model
    parallel --progress --jobs 1 \
        "$CONDA_PREFIX/bin/diamond blastp \
            --more-sensitive \
            --db $mnor_db \
            --query {} \
            --outfmt 6 \
            --out $mnor_dir/{/.}.mnor_model.tsv \
    " ::: $PROJ_DIR/data/interim/transdecoder/**/*.transdecoder.pep

}


mnor_not_best_hits () {

    # report_transcripts_w_or_wo_best_hits_mod.pl diamond.tsv best.tsv notbest.tsv
    parallel --progress \
        "$PROJ_DIR/src/scripts/report_transcripts_w_or_wo_best_hits_mod.pl \
        {} \
        {.}.best.tsv \
        {.}.not_best.tsv \
        " ::: $mnor_dir/*.mnor_model.tsv

}

# Might be removing too much? 

filter_mnor () {

    echo "Filtering."

    # Make file containing only the names (column 1) of the 
    # transcripts that were regarded as not_best against mnor

    parallel --progress \
        "awk '{print \$1}' {} > {.}.txt \
        " ::: $mnor_dir/*.not_best.tsv

    
    for ext in pep cds; do 

        # Remove everything after transcript name in header.
        # Only `>Tras_c1000891_g0_i1.p1` is kept. 
        
        parallel --progress \
             "awk '/^>/ {print \$1;next} {print \$0}' {} > {.}.clean.$ext \
             " ::: $transdecoder_dir/**/*.transdecoder.$ext

        # Exclude not_best_hits against mnor from Transdecoder.Predict results.
        
        # Transdecoder path is a bit messy but can't be helped easily.
        # {= s{^.*/|\..*+$}{}g;  =} == basename

        # faSomeRecords in.fa listFile out.fa
        parallel --progress \
            "faSomeRecords \
            -exclude \
            $transdecoder_dir/**/{= s{^.*/|\..*+$}{}g;  =}.trinity.longest_isoforms.fasta.transdecoder.$ext\
            {} \
            {.}.removed.$ext \
            " ::: $mnor_dir/*.txt
    done
}

gene_list () {

    echo "Making gene list."

    # Rename longest isoform transcripts 

    # make_gene_list_mod.pl sample in.fasta out.fasta out.tsv 
    parallel --progress \
        "src/scripts/make_gene_list_mod.pl \
            {= s{^.*/|\..*+$}{}g; =} \
            {} \
            {.}.renamed.fasta \
            {.}.renamed.tsv \
        " ::: $longest_dir/*.fasta
    
}

rename_mnor () {

    echo "Renaming Mnor not best hits against mnor gene model filtered Trandecoder.Predict results"
   
    for ext in pep cds; do

        # rename_fasta_using_gene_list_mod.pl genelistFile in.fa out.fa
        parallel --progress \
            "src/scripts/rename_fasta_using_gene_list_mod.pl \
            $longest_dir/{= s{^.*/|\..*+$}{}g; =}.trinity.longest_isoforms.renamed.tsv \
            {} \
            {.}.$ext.renamed.fasta \
            " ::: $mnor_dir/*.removed.$ext
  
    done
    
}

# Here, we want to also rename those transcripts that did not go through mnor
# So the Transdecoder.Predict outputs

rename_transdecoder () {

    echo "Renaming Transdecoder.Predict Output"

    for ext in pep cds; do

        # rename_fasta_using_gene_list_mod.pl genelistFile in.fa out.fa
        parallel --progress \
            "src/scripts/rename_fasta_using_gene_list_mod.pl \
            $longest_dir/{= s{^.*/|\..*+$}{}g; =}.trinity.longest_isoforms.renamed.tsv \
            {} \
            {.}.$ext.renamed.fasta \
            " ::: $transdecoder_dir/**/*.transdecoder.$ext
  
    done

}

# Prepare to run proteinortho both for mnor and transdecoder

proteinortho_mnor () {

    mkdir -p $proteinortho_dir/mnor/included

    ln -s $mnor_dir/*.renamed.fasta $proteinortho_dir/mnor/included/
    # /vaults/nvme0_2tb/shared/data/check_seqs/$outgroup.$seq.fasta.fixed.fasta
    # Renamed to $outgroup.$seq.fixed.fasta
    ln -s $outgroup_dir/*.fasta $proteinortho_dir/mnor/included/

    cd $proteinortho_dir/mnor
    
    proteinortho6.pl \
        -temp=$SNIC_TMP \
        -project=mnor.nbrm \
        -cpus=16 \
        -p=diamond \
        -selfblast \
        -singles \
        $proteinortho_dir/mnor/included/*.pep.*.fasta
}

# Make data libraries from proteinortho output file. 

data_lib () {

    echo "Making data lib."

    # Rename to .pep.tsv since we ned a cds version as well.
    mv $proteinortho_dir/mnor/mnor.nbrm.proteinortho.tsv \
       $proteinortho_dir/mnor/mnor.nbrm.proteinortho.pep.tsv
    
    # Make cds copy for the script to work
    sed 's/pep/cds/g' \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.pep.tsv > \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.cds.tsv
   
    # make_proteinortho_data_library_no_tarball_mod.pl 
  
    # proteinortho.tsv lib.tsv groupcount(o for outgroups) input.fa
    for ext in pep cds; do

        $PROJ_DIR/src/scripts/make_proteinortho_data_library_no_tarball_mod.pl \
            $proteinortho_dir/mnor/mnor.nbrm.proteinortho.$ext.tsv \
            $proteinortho_dir/mnor/mnor.nbrm.proteinortho.data_lib.$ext.tsv \
            'o' \
            $proteinortho_dir/mnor/included/*.$ext.*.fasta

    done

}

# Singe IQtree with bootstrapping needs min 4 seq so set min to 4.

split_header_and_get_min_4_ogs () {

    # seq_pruning_prep.R in.tsv  all.tsv
    Rscript $PROJ_DIR/src/scripts/seq_pruning_prep.R \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.data_lib.pep.tsv \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.mc.min_4.pep.tsv \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.sc.min_4.pep.tsv \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.all.min_4.pep.tsv 

}

# Split orthogroup data libraries into single files for each OG.

single_file_ogs () {

    for og in mc sc; do

        echo "Making $og-og files."
        
        mkdir -p $ogs_dir/proteinortho/mnor/$og/ && cd $ogs_dir/proteinortho/mnor/$og
        awk '(NR > 1) {print ">"$2"\n"$3 > $1".fasta"}' $proteinortho_dir/mnor/mnor.nbrm.proteinortho.$og.min_4.pep.tsv
    done
}

# Align SC-OGs, using --auto (slow) this time instead of fast

align () {

    for og in mc sc; do

        echo "Aligning $og-OGs."
    
        mkdir -p $alignments_dir/proteinortho/mnor/$og/
        mkdir -p $PROJ_DIR/logs/alignments/proteinortho/mnor/$og
        
        cd $ogs_dir/proteinortho/mnor/$og/
        find ~+ -type f -name "*.fasta" > ../$og.txt
    
        parallel --progress --eta \
            "mafft \
                --anysymbol \
                --thread 1 \
                --auto \
                {} > \
                $alignments_dir/proteinortho/mnor/$og/{/.}.aligned.fasta \
            2>$PROJ_DIR/logs/alignments/proteinortho/mnor/$og/{/.}.aligned.fasta
            " :::: $ogs_dir/proteinortho/mnor/$og.txt
    done
}

# Only need to trim the MC-OGs really. 

trim () {
     # ClipKIT vs trimal?
     # ClipKIT for now since they used it in orthosnap paper..
     for og in mc sc; do

         echo "Trimming $og-ogs."

         mkdir -p $trim_dir/proteinortho/mnor/clipkit/$og
         cd $alignments_dir/proteinortho/mnor/$og
         find ~+ -type f -name "*.fasta" > ../$og.txt
        
        # smart-gap is default parameter
        parallel --jobs 40 --eta \
            "clipkit {} \
            -o $trim_dir/proteinortho/mnor/clipkit/$og/{/.}.clipkit.fasta \
            > /dev/null \
            2> /dev/null \
            " :::: $alignments_dir/proteinortho/mnor/$og.txt
    done 
}

# Use IQTree to make trees for orthoSNAP

iqtree_run () {

    for og in mc sc; do

        mkdir -p $trees_dir/proteinortho/mnor/iqtree/$og
        
        cd $trim_dir/proteinortho/mnor/clipkit/$og
        
        # Quick fix for min 4 - not needed but leave for now  
        grep -r -c ">" |\
            awk -v FS=":" '{if ($2>3) print $1}' |\
            # 
            xargs -I {} cp {} $trees_dir/proteinortho/mnor/iqtree/$og/
        
        cd $trees_dir/proteinortho/mnor/iqtree/$og
        find ~+ -name "*.fasta" > ../$og.min_4.txt
        
        parallel --jobs 40 --eta \
            "iqtree -s {} -B 1000 -T 1 \
            " :::: $trees_dir/proteinortho/mnor/iqtree/$og.min_4.txt
    done     
}

# Run orthosnao with min 10 occupancy since we will use min 10 krill later. 
# Otherwise we would miss ~20 OGs. If a smaller number of min species is used later this needs to be changed.

orthosnap_run () {

    for tree in iqtree; do
    
        mkdir -p $pruners_dir/proteinortho/mnor/orthosnap/$tree/mc_occupancy_10
        cd $pruners_dir/proteinortho/mnor/orthosnap/$tree/mc_occupancy_10
        
        # Should be unaligned and untrimmed!
        parallel --bar "ln -s {} {/}" :::: $ogs_dir/proteinortho/mnor/mc.txt
        
        find ~+ -name "*.fasta" > ../mc.txt
        # Change back to .tree for fasttree (fasta.contree for iqtree)
        parallel --progress --eta \
            "orthosnap \
            -t $trees_dir/proteinortho/mnor/$tree/mc/{/.}.aligned.clipkit.fasta.contree \
            -f {} \
            --occupancy 10 \
            " :::: $pruners_dir/proteinortho/mnor/orthosnap/$tree/mc.txt
    done
}

# Make a single file of SNAP-OGs for R.

orthosnap_data_library () {
   
    for tree in iqtree; do
        
        cd $pruners_dir/proteinortho/mnor/orthosnap/$tree/mc_occupancy_10
 
        ls *.fa | xargs -I {} awk '/^>/ {if(N>0) printf("\n");gsub(">","",$0); printf("%s\t%s\t", FILENAME,$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' {} > ../snap_ogs_occupancy_10.tsv
    done
}

align_snap_ogs () {

    echo "Aligning SNAP-OGs."

    cd $pruners_dir/proteinortho/mnor/orthosnap/iqtree/mc_occupancy_10
    find ~+ -name "*.fa" > ../snap_ogs_occupancy_10.txt
    
    mkdir -p $alignments_dir/proteinortho/mnor/snap_ogs_occupancy_10
    mkdir -p $PROJ_DIR/logs/alignments/proteinortho/mnor/snap_ogs_occupancy_10
   
    parallel --progress --eta \
        "mafft \
            --anysymbol \
            --thread 1 \
            --auto \
            {} > \
            $alignments_dir/proteinortho/mnor/snap_ogs_occupancy_10/{/}.aligned \
        2>$PROJ_DIR/logs/alignments/proteinortho/mnor/snap_ogs_occupancy_10/{/}.aligned
        " :::: $pruners_dir/proteinortho/mnor/orthosnap/iqtree/snap_ogs_occupancy_10.txt
}

# So we can get single file nuc OGs.

split_header_and_get_min_4_ogs_cds () {

    # seq_pruning_prep.R in.tsv  all.tsv
    Rscript $PROJ_DIR/src/scripts/seq_pruning_prep.R \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.data_lib.cds.tsv \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.mc.min_4.cds.tsv \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.sc.min_4.cds.tsv \
        $proteinortho_dir/mnor/mnor.nbrm.proteinortho.all.min_4.cds.tsv 

}

# Make single file nuc OGs. 

single_file_ogs_cds () {

    #for og in mc sc; do
    for og in mc; do

        echo "Making $og-og files."
        
        mkdir -p $ogs_dir/proteinortho/mnor/cds/$og/ && cd $ogs_dir/proteinortho/mnor/cds/$og
        awk '(NR > 1) {print ">"$2"\n"$3 > $1".fasta"}' $proteinortho_dir/mnor/mnor.nbrm.proteinortho.$og.min_4.cds.tsv
    done
}

# SNAP-OGs do have less species than the original files
# Need to remove those sequences before running pal2nal

subset_snap_cds () {

    mkdir -p $ogs_dir/proteinortho/mnor/cds/snap_ogs_occupancy_10
    
    # Make list of taxa, for each OG for faSomeRecords. 
    parallel \
        "grep \">\" {} | awk -v FS=\">\" '{print \$2}' > {}.list \
        " :::: $pruners_dir/proteinortho/mnor/orthosnap/iqtree/snap_ogs_occupancy_10.txt
   
    # pipe (|) should work
    # Subset with faSomeRecords.  
    
    parallel \
        "faSomeRecords \
        $ogs_dir/proteinortho/mnor/cds/mc/{= s{^.*/|\..*+$}{}g;  =}.fasta \
        {} \
        $ogs_dir/proteinortho/mnor/cds/snap_ogs_occupancy_10/{/.}.fasta \
        " ::: $pruners_dir/proteinortho/mnor/orthosnap/iqtree/mc_occupancy_10/*.list
}

pal2nal_all () {
   
    # Combined pal2nal functions, not tested. 

    for dataset in sc snap_ogs_occupancy_10; do 

        # This should not be one function...
        
        # Make pal2nal directory
        mkdir -p $paltonal_dir/proteinortho/mnor/$dataset
        cd $paltonal_dir/proteinortho/mnor/$dataset
    
        # Lookup cds seqs now with right number of taxa
        cd $ogs_dir/proteinortho/mnor/cds/$dataset
        find ~+ -type f -name "*.fasta" > ../$dataset.txt
    
        # Make pal2nal out, and log directory
        cd $paltonal_dir/proteinortho/mnor/$dataset
        mkdir -p $paltonal_dir/proteinortho/mnor/$dataset/out
        mkdir -p $logs_dir/pal2nal/proteinortho/mnor/$dataset
    
        # Run pal2nal 
        parallel --eta \
            "pal2nal.pl \
            $alignments_dir/proteinortho/mnor/$dataset/{/}.aligned \
            {} \
            -output fasta \
            > out/{/.}.aligned.pal2nal.fasta \
            2> $logs_dir/pal2nal/proteinortho/mnor/$dataset/{/.}.log \
            " :::: $ogs_dir/proteinortho/mnor/cds/$dataset.txt
    
        # Remove everything after pipe, leaving ">ecry"
        cd $paltonal_dir/proteinortho/mnor/$dataset/out
        mkdir -p $paltonal_dir/proteinortho/mnor/$dataset/out_clean
        parallel "sed 's/|.*//g' {} > ../out_clean/{}.clean" ::: *.fasta
    
        # Make min 10 krill dir
        cd $paltonal_dir/proteinortho/mnor/$dataset/out_clean
        mkdir -p $paltonal_dir/proteinortho/mnor/$dataset/out_clean_krill_10
       
        # Grep for krill, if match is ten or more, copy into SNAP-dir
        grep -r -c "^>[e|n|m|t]" | awk -v FS=":" '{if ($2>=10) print $1}' | xargs -I {} cp {} ../out_clean_krill_10
        
        # Combine into combined-dir
        mkdir -p $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10
        cd $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10
        cp ../../$dataset/out_clean_krill_10/* .
    
    done

}

# Tried different trees to get iqtree to work, removed N, and renamed sequences.. 
# Not sure what was the problem, worked with N. 
# iqtree -p out_clean_krill_10_renamed -B 1000 -T 20 --prefix p.krill.renamed.with.N

iqtree_combined () {

    # Make species tree, nuc, non trimmed.
    cd $paltonal_dir/proteinortho/mnor/combined/

    iqtree -p out_clean_krill_10 -B 1000 -T 20
    
    # Tree used downstream is named 'tree'.
    # This tree is the same as p.krill.renamed.with.N.contree
    # cp output contree -> tree

}

# Then subset the orthogroups/alignments 

subset_pal2nal_seqs () {

    mkdir -p $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/

    cd $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10
    
    echo "Remove outgroups from pal2nal_seqs"

    \ls *.clean |\
        xargs -P 8 -I {} sh -c \
        'subsetfasta "$1" exclude ../outgroups.txt \
        > ../out_clean_krill_10_no_outgroups/"$1".no_outgroups' -- {}
}

# Remove stop codons for hyphy (and PAML but not sure if needed)

remove_stop_codon () {	
    
	# Doesn't remove gaps..
    # 11312 no_outgroups
    # 11232 nostop - guessing there are a few without stop codons? 

    # 8, 10 removes seqsuences with stop codons, 8, 3 replaces stop with caps 
	parallel \
        "(echo 8; \
        echo 3; \
        echo 1; \
        echo {}; \
        echo 1; \
        echo {}.nostop) | hyphy \
        " ::: $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/*.no_outgroups

}

fix_stop_codon () {

    find $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/ -name "*.no_outgroups" -printf '%f\n' \
        > $paltonal_dir/proteinortho/mnor/combined/no_outgroups.lst
    find $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups -name "*.nostop" -printf '%f\n' | sed 's/\.[^.]*$//' \
        > $paltonal_dir/proteinortho/mnor/combined/nostop.lst
    
    grep -v -f $paltonal_dir/proteinortho/mnor/combined/nostop.lst \
               $paltonal_dir/proteinortho/mnor/combined/no_outgroups.lst |\
               xargs -I {} cp $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/{} \
               $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/{}.nostop


}

# Make a list of sequences in each orthogroup

list_seqs_in_nostop () {

    echo "Make list of names in each alignment.."

    cd $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups
    \ls *.nostop |\
    xargs -P 8 -I {} sh -c 'grep ">" "$1" | awk -v FS=">" '\''{print $2}'\'' > "$1".list' -- {}

}

# Make a pruned species tree matching each orthogroup

make_pruned_species_trees () {

    wd=$paltonal_dir/proteinortho/mnor/combined
    
    species_all=("ecry" "edin" "efri" "elam" "elon" "emuc" "epac" "esim" "esup" "etri" "eval" "mnor" "nmeg" "tine" "tlon" "tmac" "trac" "oham" "ohaz" "opmo" "opva" "otca" "ocqu" "odma")

    echo "Make keepfile"

    printf "%s\n" "${species_all[@]}" > $wd/keep_file.txt

    echo "Make list to prune from species tree"
    
    # Find names in list that is also in the species tree in reverse
    cd $wd/out_clean_krill_10_no_outgroups
    \ls *.list |\
        xargs -P 8 -I {} sh -c \
        'grep -v -f "$1" '$wd'/keep_file.txt > "$1".prune' -- {}
    
    # Set species tree

    species_tree=$wd/tree

    echo "Make subtrees"

    cd $wd/out_clean_krill_10_no_outgroups

    # I should learn python..
    # If prune file is not empty, prune tree, otherwise, keep species tree...
    # Should redirect log
    \ls *.prune |\
    xargs -P 8 -I {} sh -c \
    'if [ -s "$1" ]; then cat "$1" | xargs ~/phyutility/phyutility.jar -pr -in '$species_tree' -names > "$1".tree; else cp '$species_tree' "$1".tree; fi '  -- {}
            

}

mark_foreground_branches () {
    
    warm=("edin" "elam" "emuc" "esim" "nmeg")

    printf "%s\n" "${warm[@]}" > $paltonal_dir/proteinortho/mnor/combined/warm.txt
    
    # Not sure what internal-nodes labeling is appropriate 
    cd $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/
    \ls *.tree |\
    xargs -P 8 -I {} sh -c 'hyphy '$paltonal_dir'/proteinortho/mnor/combined/label-tree.bf --tree "$1" --list '$paltonal_dir'/proteinortho/mnor/combined/warm.txt --output "$1".fg' -- {}
}

hyphy_prep () {
       
    cd $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/

    find ~+ -type f -name "*.nostop" > ../hyphy_nostop.txt
    find ~+ -type f -name "*.fg" > ../hyphy_fg.txt

    mkdir -p $dnds_dir/proteinortho/mnor/combined/hyphy
    cd $dnds_dir/proteinortho/mnor/combined/hyphy

    time cat $paltonal_dir/proteinortho/mnor/combined/hyphy_nostop.txt | xargs -P 40 -I {} ln -s {} .
    time cat $paltonal_dir/proteinortho/mnor/combined/hyphy_fg.txt | xargs -P 40 -I {} ln -s {} .
}

# Run hyphy on uppmax

# Some trees will have no fg seqs -> gives error.
# Removed these beforehand manually(!);
# these are 21491 and 27925...
# Should we set  fg threshold?

hyphy_run () {
    cd $dnds_dir/proteinortho/mnor/combined/hyphy
    find ~+ -name "*.nostop" > ../hyphy_nostop.txt
    
    parallel --bar --load 40 \
        "hyphy CPU=1 meme --alignment {} --tree {}.list.prune.tree.fg > {}.MEME.stdout" :::: $dnds_dir/proteinortho/mnor/combined/hyphy_nostop.txt

    # On UPPMAX HyPhy was run like:
    # parallel --jobs 20 "hyphy busted CPU=1 --alignment {} --tree {}.list.prune.tree.fg --branches Foreground > {}.BUSTED.stdout; echo {%} {} BUSTED" ::: $@
    
    # The script submitted with:
    # ls hyphy/*.nostop | xargs -n 200 -P 500 sbatch -A snic2021-5-453 -p node -n 1 -t 24:00:00 -J hyphy-busted --output=BUSTED-%j.out ./testscript-busted.sh
}

#remove_stop_codon
#hyphy_run_run

paml_prep () {
    echo "Find files."
    cd $paltonal_dir/proteinortho/mnor/combined/out_clean_krill_10_no_outgroups/
    find ~+ -type f -name "*.nostop" > ../paml_nostop.txt
    find ~+ -type f -name "*.tree" > ../paml_tree.txt
    
    mkdir -p $dnds_dir/proteinortho/mnor/combined/paml/included
    cd $dnds_dir/proteinortho/mnor/combined/paml/included
    echo "Link files."
    time cat $paltonal_dir/proteinortho/mnor/combined/paml_nostop.txt | xargs -P 40 -I {} ln -s {} .
    time cat $paltonal_dir/proteinortho/mnor/combined/paml_tree.txt | xargs -P 40 -I {} ln -s {} .

    # Ok, so for every PAML run we need the correct warm species
    # So make a list of correct warm species for every OG

    warm=$paltonal_dir/proteinortho/mnor/combined/warm.txt
    paml=$dnds_dir/proteinortho/mnor/combined/paml
    
    cd $paml/included
    mkdir -p $paml/output_3/
    mkdir -p $paml/log_bsa_test_2/
    mkdir -p $paml/img/
    # For every alignment/OG remove > and put into list
    \ls *.nostop | xargs -P 20 -I {} sh -c 'grep ">" "$1" | sed '\''s/>//g'\'' > "$1".list' -- {}
    \ls *.list | xargs -P 20 -I {} sh -c 'grep -f '$warm' "$1" > "$1".warm' -- {}
    cd $paml
    
    conda activate ete3
    # Run ete
    time parallel --jobs 40 --bar "cat {}.list.warm | xargs | sed 's/ /=/g' | xargs ete3 evol --alg {} -t {}.list.prune.tree --clear_all --codeml_binary /usr/bin/codeml --models bsA bsA1 --tests bsA,bsA1 -o $paml/output_3/{/} --cpu 1 --mark > $paml/log_bsa_test_2/{/}" ::: $paml/included/*.nostop
}

trim_dnds_to_check_quality () {
    
    mkdir -p $trim_dir/proteinortho/mnor/clipkit/dnds_combined
    
    cd $dnds_dir/proteinortho/mnor/combined/hyphy/

    \find ~+ -name "*.nostop" | xargs -P 20 -I {} ln -s {} $trim_dir/proteinortho/mnor/clipkit/dnds_combined

    cd $trim_dir/proteinortho/mnor/clipkit/dnds_combined
    \ls *.nostop | xargs -P 20 -I {} sh -c 'clipkit "$1" > "$1".stdout' -- {}

    grep -r "sites trimmed:" | awk -v FS=".aligned|trimmed:" -v OFS="\t" '{print $1,$3 >> "../trimmed"}'
    grep -r "kept" | awk -v FS=".aligned|kept:" -v OFS="\t" '{print $1,$3 >> "../kept"}'
    grep -r "Percentage" | awk -v FS=".aligned|trimmed:" -v OFS="\t" 'sub("%", "", $3) {print $1,$3/100 >> "../percentage" }'
}

kaks () {
    
    KAKS_HOME=/vaults/nvme0_2tb/felix/tools/KaKs_Calculator1.2/src
    
    # Make kaks dir 
    mkdir -p $dnds_dir/proteinortho/mnor/combined/kaks/included
    
    # Link same sequences used as in PAML/hyphy
    cd $dnds_dir/proteinortho/mnor/combined/hyphy/
    echo "Link files."
    \find ~+ -name "*.nostop" | xargs -P 20 -I {} ln -s {} $dnds_dir/proteinortho/mnor/combined/kaks/included
    cd $dnds_dir/proteinortho/mnor/combined/kaks/

    samples=("ecry" "edin" "efri" "elam" "elon" "emuc" "epac" "esim" "esup" "etri" "eval" "mnor" "nmeg" "tine" "tlon" "tmac" "trac")
    
    # Kaks dir 
    KAKS=$dnds_dir/proteinortho/mnor/combined/kaks

    # Loop through samples 
    for i in ${!samples[@]}; do
        SPEC1=${samples[$i]}
        for j in ${samples[@]:$i+1:${#samples[@]}}; do
            SPEC2=$j
            VS=$SPEC1.vs.$SPEC2 
            # Make subset list
            rm $KAKS/list.txt
            echo $SPEC1 >> $KAKS/list.txt
            echo $SPEC2 >> $KAKS/list.txt

            mkdir -p $KAKS/subset/$VS
            mkdir -p $KAKS/datasets/$VS

            # Subset fastas
            cd $KAKS/included
            echo "Look for $SPEC1 vs $SPEC2" 
            \ls *.nostop | xargs -P 40 -I {} sh -c 'faSomeRecords "$1" "'$KAKS'"/list.txt "'$KAKS'"/subset/"'$SPEC1'".vs."'$SPEC2'"/"$1"' -- {}
            
            # Keep fastas where both occur
            cd $KAKS/subset/$VS
            mkdir -p $KAKS/correct/$VS
            echo "Copy where 2 seqs."
            grep -r ">" -c | awk -v FS=":" '{if ($2 > 1) print $1}' | xargs -P 40 -I {} cp {} $KAKS/correct/$SPEC1.vs.$SPEC2/

            cd $KAKS/correct/$VS
            mkdir -p $KAKS/axt/$VS
            # Make axt
            PL () {
                name=$1
                file=$5
                SPEC1=$2
                SPEC2=$3
                KAKS=$4
                VS=$SPEC1.vs.$SPEC2
                echo $name > $KAKS/axt/$VS/$file.axt

                cat $5 | sed "s/>$SPEC1//g" | sed "s/>$SPEC2//g" | awk NF >> $KAKS/axt/$VS/$file.axt
                echo "Make axt: $name $SPEC1 $SPEC2" 
            }; export -f PL
            # Run make axt
            parallel "PL {= s{\.aligned.*}{}g; =} $SPEC1 $SPEC2 $KAKS {}" ::: *.nostop
            
            wait
            # Remove all empty files
            #rm -r $KAKS/output/$SPEC1.vs.$SPEC2
            mkdir -p $KAKS/output/$VS
            mkdir -p $KAKS/log/$VS
            cd $KAKS/axt/$VS
            # Run KaKs Calculator
            # Needs local version of parallel installed/newer version than on corbicula. Bug in load.
            echo "Run KaKs Calc."
            parallel --eta "$KAKS_HOME/KaKs_Calculator -m GY -i {} -o $KAKS/output/$VS/{}.output > $KAKS/log/$VS/{}.log" ::: *.axt
            #sleep 5
            wait

            # Make single file

            cd $KAKS/output/$VS

            #rm $KAKS/$SPEC1.vs.$SPEC2.txt
            touch $KAKS/$VS.txt
            \ls * | xargs -P 40 -I {} awk '(NR == 2) {print $0}' {} >> $KAKS/$VS.txt
   time for file in *; do
                cut $file | awk -v OG="${file%%.*}" '(NR == 2) {$0}' >> $KAKS/$VS.txt
                echo "Single file: "$file" "$SPEC1" "$SPEC2
            done

            sleep 5
            wait

        done
        wait
    done
            # Get alignments for the two species, make .axt files
}

entap_run () {

parallel --jobs 1 "/vaults/nvme0_2tb/felix/tools/EnTAP-v0.10.8-beta/EnTAP \
    --runN \
    -d /vaults/nvme0_2tb/felix/analysis/orthology/longest_isoforms/proteinortho/entap/out_cds/bin/ALL_NON_REFSEQ_TAGGED
.dmnd \
    -d /vaults/nvme0_2tb/felix/analysis/orthology/longest_isoforms/proteinortho/entap/out_cds/bin/invertebrate.dmnd \
    -d /vaults/nvme0_2tb/felix/analysis/orthology/longest_isoforms/proteinortho/entap/out_cds/bin/lobster.dmnd \
    -d /vaults/nvme0_2tb/felix/analysis/orthology/longest_isoforms/proteinortho/entap/out_cds/bin/uniprot_sprot.dmnd \
    --out-dir /vaults/nvme2_2tb/felix/projects/orthology/data/interim/entap/no_esup/{/}         -t 10         --ini /va
ults/nvme0_2tb/felix/tools/EnTAP-v0.10.8-beta/entap_config.ini -i {}" ::: /vaults/nvme2_2tb/felix/projects/orthology/da
ta/interim/longest_isoforms/*.longest_isoforms.renamed.fasta

    # Then combine results into one file and load into R
    # Possibly as such:
    # ls *0.tsv | xargs -I {} awk -v OFS="\t" '{sub("\\.trin.*","",FILENAME); if (NR>1) print FILENAME,$0}' {} > lib.tsv
}

#extract_longest_isoforms
#get_orfs
#search_sprot
#search_pfam
#predict_transcripts
#mnor_model
#mnor_not_best_hits
#filter_mnor
#gene_list
#rename_mnor
#rename_transdecoder
#proteinortho_mnor

#data_lib
#split_header_and_get_min_4_ogs
#single_file_ogs
#align
#trim
#iqtree_run
#orthosnap_run 
#orthosnap_data_library
#align_snap_ogs
#split_header_and_get_min_4_ogs_cds
#single_file_ogs_cds
#subset_snap_cds

#pal2nal_all
#iqtree_combined

#remove_outgroups_from_species_tree
#subset_pal2nal_seqs 
#remove_stop_codon
#fix_stop_codon
#list_seqs_in_nostop
#make_pruned_species_trees
#mark_foreground_branches

#hyphy_prep
#hyphy_run
#paml_prep_and_run

#entap_run
#kaks_run

#trim_dnds_to_check_quality

#run_busco longest_isoforms renamed.fasta
#run_busco mnor_model removed.pep

