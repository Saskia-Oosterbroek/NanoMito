# Assembling mitochondrial genomes from native sequencing Nanopore data.
# Author: Saskia Oosterbroek

# Required dependencies:
# NanoFilt https://github.com/wdecoster/nanofilt
# Minimap2 https://github.com/lh3/minimap2
# Samtools https://github.com/samtools/samtools
# Flye https://github.com/fenderglass/Flye

# Optional: 
# guppy_basecaller https://community.nanoporetech.com/downloads
# Canu https://github.com/marbl/canu

# Flye was developed for long read data and uses a minimum read overlap of 1000 bases. When assembling mitochondrial genomes, especially from degraded samples, shorter reads might still be valuable. 
# The 1000 bases minimum can be circumvented by adjusting the "main.py" file in flye's installation folder
# Find the following piece of code: 
#
# parser.add_argument("-m", "--min-overlap", dest="min_overlap", metavar="int",
#                        type=lambda v: check_int_range(v, 1000, 10000),
#                        default=None, help="minimum overlap between reads [auto]")
#
# Change (v, 1000, 10000) to (v, 100, 10000)
# And thr minimum read length overlap is set to 100 bases.

# Copy sequencing data from Mk1C or MinION to your computer
rsync -avh --ignore-existing minit@MC-123456:/data/project_name/sample_name/ /home/user/Nanopore_data/project_name/sample_name/RAW_data

# Run basecalling with Guppy in super accuracy(SUP). If needed, adjust the configuration file (.cfg) used for the right basecaller and set the --barcode flag to the right expansion kit (if used)
guppy_basecaller -c dna_r10.4.1_e8.2_400bps_sup.cfg --input_path /home/user/Nanopore_data/project_name/sample_name/RAW_data --save_path /home/user/Nanopore_data/project_name/sample_name/SUP_basecalled -r -x cuda:0 --do_read_splitting --max_read_split_depth 10 --trim_adapters --barcode_kits EXP-NBD104 
# --resume can be added if more data for the same sample is generated later

# Make a back-up of the RAW data 

# Create 1 fastq file per barcode. Run in the folder containing all separate barcode folders.
mkdir data_analysis
for folder in barcode* ; do
mkdir data_analysis/"$folder"
    cat "$folder"/*.fastq > data_analysis/"$folder"/"$folder".fastq
done

#Enter data_analysis folder
cd data_analysis
# Filter data for a length that suits mitochondrial genomes.
for folder in barcode* ; do
    NanoFilt -q 8 -l 250 --maxlength 30000 "$folder"/"$folder".fastq > "$folder"/filtered_"$folder".fastq
    # Check number of reads per barcode
    echo "$folder"
    line_numbers=$(wc -l < "$folder"/filtered_"$folder".fastq)
    echo "$((line_numbers/4)) reads in filtered fastq"
done

# Optional: if no related reference mitogenome is available flye could be run native which is more computationally intensive.
# This will produce MANY contigs as genomic pieces of DNA will also be assembled
# Set the number of threads appropriately for your computer system!
# Set the read-error appropriately for your flowcell type and basecalling quality.
# NOT recommended for quick assembly of mtGenome
# $ for folder in barcode* ; do
# $   flye --nano-hq "$folder"/filtered*.fastq --genome-size 15k --threads 8 --meta --iterations 5 --read-error 0.07 --min-overlap 250 --out-dir "$folder"/native_flye
# $ done

# Map the data to one or multiple reference fastas, assemble only those reads that mapped to a reference mitogenome.
# Make sure to place your reference in the newly made barcode specific folder, assuming there are different species per barcode.
for folder in barcode* ; do
    (
    cd "$folder"
    for database in *.fasta ; do
        # Make an output directory for every fasta reference file
        mkdir ${database%%.*}
        # Map the reads to the reference using a relaxed mapping setting for the draft assembly. "-ax map-ont" "Align noisy long reads of ~10% error rate to a reference genome."
        minimap2 -ax map-ont -k10 "$database" barcode*.fastq -t 20 > ${database%%.*}/align_all_data_${database%%.*}.sam ;
        samtools view -S -b ${database%%.*}/align_all_data_${database%%.*}.sam > ${database%%.*}/align_all_data_${database%%.*}.bam
        samtools view -b -F 4 ${database%%.*}/align_all_data_${database%%.*}.bam > ${database%%.*}/mapped_${database%%.*}.bam
        # Get unique mapping sequences
        samtools view -bq 1 ${database%%.*}/mapped_${database%%.*}.bam > ${database%%.*}/unique_${database%%.*}.bam
        # Make fastq with mapping reads only
        samtools fastq ${database%%.*}/unique_${database%%.*}.bam > ${database%%.*}/unique_${database%%.*}.fastq
        # Assemble reads that mapped to the reference genome(s)
        flye --nano-hq ${database%%.*}/unique_${database%%.*}.fastq --genome-size 15k --threads 20 --meta --scaffold --iterations 5 --read-error 0.05 --min-overlap 100 --out-dir ${database%%.*}/flye

        # OPTIONAL: if sample quality is low and DNA fragments are short it could be considered to try using Canu assembler
        # canu -nanopore-corrected ${database%%.*}/unique_${database%%.*}.fastq -d ${database%%.*}/canu_assembly_${database%%.*} -p Canu_${database%%.*}_ minReadLength=100 minOverlapLength=25 genomeSize=16k -assemble contigFilter="6 0 1.0 0.5 0"
  
        # Remove some large intermediate files to prevent your system running full
        rm ${database%%.*}/align*.bam
        rm ${database%%.*}/align*.sam
    done
    )
done


# Put all assemblies together in a folder
mkdir assemblies
for folder in barcode* ; do
    #Flye
    for file in "$folder"/*/*/assembly.fasta ; do
        if [ -e "$file" ] ; then
        filename=$(echo "$file" | sed 's/\//_/g' )
        cp "$file" assemblies/"$filename" ; fi
    done
    
    #Canu
    for file in "$folder"/*/canu*/*.contigs.fasta ; do
        if [ -e "$file" ] ; then
        filename=$(echo "$file" | sed 's/\//_/g' )
        cp "$file" assemblies/"$filename" ; fi
    done

    #Native
    for file in "$folder"/native_flye/assembly.fasta ; do
        if [ -e "$file" ] ; then
        filename=$(echo "$file" | sed 's/\//_/g' )
        cp "$file" assemblies/"$filename" ; fi
    done
done

# Manually check which assemblies make sense and are actually mitochondrial, for instance with Blast.
# Place your draft assemblies in a new folder with the corresponding filtered fastq file
# Repeat the assembly steps, now aligning the data against the produced draft assemblies with stricter mapping settings
for folder in barcode* ; do
(
    cd "$folder"
    for draft in *.fasta ; do
    # make an output directory for every fasta draft assembly
    mkdir ${draft%%.*}
    # Map the reads to the draft assembly using a stricter mapping setting for the assembly. "-ax asm10" and/or "-ax asm5" might be appropriate depending on the amount and quality of the data. 
    # "asm5 	Long assembly to reference mapping. Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%."
    # "asm10 	Long assembly to reference mapping. Up to 10% sequence divergence."
    minimap2 -ax asm10 -k10 "$draft" filtered*.fastq -t 20 > ${draft%%.*}/align_all_data_${draft%%.*}.sam ;
    samtools view -S -b ${draft%%.*}/align_all_data_${draft%%.*}.sam > ${draft%%.*}/align_all_data_${draft%%.*}.bam
    samtools view -b -F 4 ${draft%%.*}/align_all_data_${draft%%.*}.bam > ${draft%%.*}/mapped_${draft%%.*}.bam
    # Get unique mapping sequences
    samtools view -bq 1 ${draft%%.*}/mapped_${draft%%.*}.bam > ${draft%%.*}/unique_${draft%%.*}.bam
    # Make fastq with mapping reads only
    samtools fastq ${draft%%.*}/unique_${draft%%.*}.bam > ${draft%%.*}/unique_${draft%%.*}.fastq

    # Assemble reads that mapped to the draft genome(s)
    flye --nano-hq ${draft%%.*}/unique_${draft%%.*}.fastq --genome-size 15k --threads 20 --meta --scaffold --iterations 5 --read-error 0.05 --min-overlap 100 --out-dir ${draft%%.*}/flye

    # OPTIONAL: if sample quality is low and DNA fragments are short it could be considered to try using Canu assembler
    #$ canu -nanopore-corrected ${draft%%.*}/unique_${draft%%.*}.fastq -d ${draft%%.*}/canu_assembly_${draft%%.*} -p Canu_${draft%%.*}_ minReadLength=100 minOverlapLength=25 genomeSize=16k -assemble contigFilter="6 0 1.0 0.5 0"
  
    # remove some large intermediate files to prevent your system running full
    rm ${draft%%.*}/align*.bam
    rm ${draft%%.*}/align*.sam
    done
)
done

# Put all assemblies together in a folder
mkdir assemblies_round2
for folder in barcode* ; do
    #Flye
    for file in "$folder"/*/*/assembly.fasta ; do
        if [ -e "$file" ] ; then
        filename=$(echo "$file" | sed 's/\//_/g' )
        cp "$file" assemblies_round2/"$filename" ; fi
    done
    
    #Canu
    for file in "$folder"/*/canu*/*.contigs.fasta ; do
        if [ -e "$file" ] ; then
        filename=$(echo "$file" | sed 's/\//_/g' )
        cp "$file" assemblies_round2/"$filename" ; fi
    done
done

# Repeat aligning and assembling against new versions of the draft genomes until no improvement is seen.

# Polish with Medaka
InputFastq=unique_Ostrea_edulis.fastq # Reads that uniquely mapped to reference
DraftReference=Flye_assembly_Ostrea_edulis.fasta # Assembly from previous steps
MedakaModel=r941_min_sup_g507
MutiThread=20
# Medaka optional: -g : don't fill gaps in consensus with draft sequence.
# Make sure to use the correct medaka model for flowcell type and basecall method! For list run medaka_consensus -h

# Align data to draft sequence
minimap2 -ax map-ont "$DraftReference" "$InputFastq" -t "$MULTITHREAD" > align_"${InputFastq%%.*}".sam ;
# Polish first round with Racon (optimized for Medaka: -m 8 -x -6 -g -8 -w 500)
racon -m 8 -x -6 -g -8 -w 500 -t "$MutiThread" "$InputFastq" align_"${InputFastq%%.*}".sam "$DraftReference" > polished_"${DraftReference%%.*}".fasta ;
# polish with Medaka
medaka_consensus -g -i "$InputFastq" -d polished_"${DraftReference%%.*}".fasta -o Medaka-polished2_"{$InputFastq%%.*}".fasta -t "$MutiThread" -m "$MedakaModel"


# After successful mitochondrial assembly the quality can be checked by annotating the genome on Mitos webserver, http://mitos.bioinf.uni-leipzig.de (Bernt et al., 2013). 
# The annotation will show missing and split- or duplicate elements of the mitochondrial genome. 