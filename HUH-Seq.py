import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import numpy as np

# Function to count unique k-mers from library
def count_unique_sequences(fasta_file):
    sequence_counts = defaultdict(int)
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Bring in read as a string
        sequence = str(record.seq)
        # Find sequence upstream of the k-mer
        ref_position = sequence.find("AAAGAAGTGC")
        # Find the 'AC' dinucleotide that comes past the cleavage motif
        past_cleavage_site = sequence[ref_position+22:ref_position+24]
        # Find the cleavage motif, i.e., the k-mer
        cleavage_motif = sequence[ref_position+15:ref_position+24]
        # If the 'AC' dinucleotide is present, and if the k-mer does not contain
        # an ambiguous base call either add it to the sequence count dictionary as 
        # a key with a value of 1 or increase the value of the existing k-mer key
        # by 1
        if past_cleavage_site == "AC":
            if 'N' not in cleavage_motif:
                sequence_counts[cleavage_motif] += 1

    return dict(sequence_counts)

# Paths to the de-multiplexed FASTA files
# Exclusively bringing in the forward reads because it extends well past the k-mer
# No need for pairing the ends or reverse complementing the reverse reads
fasta_files = {
    'ref1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC1.fasta",
    'ref2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC2.fasta",
    'ref3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC3.fasta",
    'PCV2_H1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC4.fasta",
    'PCV2_H2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC5.fasta",
    'PCV2_H3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC6.fasta",
    'PCV2_L1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC7.fasta",
    'PCV2_L2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC8.fasta",
    'PCV2_L3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC9.fasta",
    'E1_H1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC10.fasta",
    'E1_H2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC11.fasta",
    'E1_H3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC12.fasta",
    'E1_L1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC13.fasta",
    'E1_L2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC14.fasta",
    'E1_L3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC15.fasta",
    'E2_H1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC16.fasta",
    'E2_H2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC17.fasta",
    'E2_H3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC18.fasta",
    'E2_L1': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC19.fasta",
    'E2_L2': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC20.fasta",
    'E2_L3': "/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/ForwardReadChunks/BC21.fasta",
}

# Initialize an empty dataframe to contain data
df_combined = pd.DataFrame()

# Process each FASTA file and populate the dataframe
for name, fasta_file in fasta_files.items():
    sequence_counts = count_unique_sequences(fasta_file)
    df_temp = pd.DataFrame(list(sequence_counts.items()), columns=['kmers', name])
    
    if df_combined.empty:
        df_combined = df_temp
    else:
        df_combined = pd.merge(df_combined, df_temp, on='kmers', how='outer')

# Sort the dataframe by 'kmers' in alphabetical order
df_combined = df_combined.sort_values(by='kmers').reset_index(drop=True)

# Assign the k-mer # based on the index
df_combined['k-mer #'] = range(1, len(df_combined) + 1)

# Reorder columns for export
cols = ['k-mer #', 'kmers'] + list(fasta_files.keys())
df_combined = df_combined[cols]

# Create the averaged DataFrame by taking the mean of the replicates
# Note the exclusion of PCV2_L1, which had poor correlation with the other two replicates
df_averaged = pd.DataFrame()
df_averaged['kmers'] = df_combined['kmers']
df_averaged['Ref_mean'] = df_combined[['ref1', 'ref2', 'ref3']].mean(axis=1)
df_averaged['PCV2_H_mean'] = df_combined[['PCV2_H1', 'PCV2_H2', 'PCV2_H3']].mean(axis=1)
df_averaged['PCV2_L_mean'] = df_combined[['PCV2_L2', 'PCV2_L3']].mean(axis=1)
df_averaged['E1_H_mean'] = df_combined[['E1_H1', 'E1_H2', 'E1_H3']].mean(axis=1)
df_averaged['E1_L_mean'] = df_combined[['E1_L1', 'E1_L2', 'E1_L3']].mean(axis=1)
df_averaged['E2_H_mean'] = df_combined[['E2_H1', 'E2_H2', 'E2_H3']].mean(axis=1)
df_averaged['E2_L_mean'] = df_combined[['E2_L1', 'E2_L2', 'E2_L3']].mean(axis=1)

# Add the k-mer # to the averaged DataFrame
df_averaged['k-mer #'] = df_combined['k-mer #']

# Reorder columns for export
df_averaged = df_averaged[['k-mer #', 'kmers', 'Ref_mean', 'PCV2_H_mean', 'PCV2_L_mean',
                           'E1_H_mean', 'E1_L_mean', 'E2_H_mean', 'E2_L_mean',
                           'RepB_mean', 'RepBm_mean']]

# Create the percent reduction DataFrame
df_percent_reduction = pd.DataFrame()
df_percent_reduction['k-mer #'] = df_averaged['k-mer #']
df_percent_reduction['kmers'] = df_averaged['kmers']
df_percent_reduction['PCV2_H_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['PCV2_H_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['PCV2_L_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['PCV2_L_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['E1_H_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['E1_H_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['E1_L_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['E1_L_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['E2_H_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['E2_H_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['E2_L_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['E2_L_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['RepB_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['RepB_mean']) / df_averaged['Ref_mean']) * 100
df_percent_reduction['RepBm_percent_reduction'] = ((df_averaged['Ref_mean'] - df_averaged['RepBm_mean']) / df_averaged['Ref_mean']) * 100

numeric_columns = df_percent_reduction.select_dtypes(include=[np.number]).columns
df_percent_reduction[numeric_columns] = df_percent_reduction[numeric_columns].clip(lower=0)

# Save all DataFrames to CSV
df_combined.to_csv("/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/Dataframes/RepKmersAll.csv", index=False)
df_averaged.to_csv("/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/Dataframes/RepKmersAvg.csv", index=False)
df_percent_reduction.to_csv("/Users/adamsmiley/Desktop/HUH-Seq Data/00_fastq/Dataframes/RepKmersReductions.csv", index=False)
