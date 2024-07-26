import os
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gzip

def compute_read_lengths_statistics(file_path):
    open_func = gzip.open if file_path.endswith(".gz") else open
    with open_func(file_path, "rt") as handle:
        read_lengths = [len(record.seq) for record in SeqIO.parse(handle, "fastq")]
    return {
        "mean": np.mean(read_lengths),
        "median": np.median(read_lengths),
        "std_dev": np.std(read_lengths),
        "all_lengths": read_lengths
    }

def plot_kde(illumina_lengths, nanopore_lengths, save_as="read_lengths_distribution.png"):
    plt.figure(figsize=(10, 6))
    sns.kdeplot(illumina_lengths, label="Illumina", fill=True)
    sns.kdeplot(nanopore_lengths, label="Nanopore", fill=True)
    plt.xlabel('Read Length')
    plt.ylabel('Density')
    plt.title('Distribution of Read Lengths')
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_as)
    plt.show()

# Add more plotting functions here if needed

def main():
    illumina_samples_dir = 'path to folder'
    nanopore_samples_dir = 'path to folder'

    illumina_lengths = []
    nanopore_lengths = []

    for file in os.listdir(illumina_samples_dir):
        if file.endswith("_1.fastq.gz"):
            sample_name = file.rsplit("_", 2)[0]
            illumina_file_1 = os.path.join(illumina_samples_dir, f'{sample_name}_1.fastq.gz')
            illumina_file_2 = os.path.join(illumina_samples_dir, f'{sample_name}_2.fastq.gz')
            nanopore_file = os.path.join(nanopore_samples_dir, f'{sample_name}.fastq')

            if not os.path.exists(illumina_file_1) or not os.path.exists(illumina_file_2):
                continue

            illumina_lengths.extend(compute_read_lengths_statistics(illumina_file_1)["all_lengths"])
            illumina_lengths.extend(compute_read_lengths_statistics(illumina_file_2)["all_lengths"])

            if os.path.exists(nanopore_file):
                nanopore_lengths.extend(compute_read_lengths_statistics(nanopore_file)["all_lengths"])

    # Create KDE plot
    plot_kde(illumina_lengths, nanopore_lengths)

    # Call more plotting functions here if needed
    print(illumina_samples_dir)
    print(nanopore_samples_dir)
    print(file)
    print(len(illumina_lengths))
    print(len(nanopore_lengths))
    print(illumina_lengths[:10])
    print(nanopore_lengths[:10])
    print(np.isnan(illumina_lengths).any())
    print(np.isnan(nanopore_lengths).any())





if __name__ == "__main__":
    main()

