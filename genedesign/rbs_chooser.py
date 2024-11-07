from Bio import SeqIO
from collections import defaultdict
import pandas as pd
from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.Translate import Translate
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from typing import Set

# Optimized function to extract gene and CDS info in a single loop
def extract_genes_info(genbank_file):
    gene_dict = defaultdict(dict)
    for record in SeqIO.parse(genbank_file, "genbank"):
        # Create a dictionary of CDS features with locus_tag for fast lookup
        cds_features = {f.qualifiers.get("locus_tag", [None])[0]: f for f in record.features if f.type == "CDS"}
        
        # Iterate through gene features and match with CDS by locus_tag
        for feature in record.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                if locus_tag and locus_tag in cds_features:
                    cds_feature = cds_features[locus_tag]
                    gene_name = feature.qualifiers.get("gene", [None])[0]
                    
                    start, end = cds_feature.location.start, cds_feature.location.end
                    strand = cds_feature.location.strand
                    
                    # Extract UTR sequence based on strand
                    utr_seq = (record.seq[max(0, start - 50):start] 
                               if strand == 1 
                               else record.seq[end:end + 50].reverse_complement())
                    
                    cds_seq = cds_feature.extract(record.seq)
                    
                    gene_dict[locus_tag] = {"gene": gene_name, "UTR": utr_seq, "CDS": cds_seq}
    return gene_dict

# Optimized function to extract the top 5% abundance
def extract_top_abundance(file_path):
    data = pd.read_csv(file_path, sep="\t", header=None, names=["string_external_id", "abundance"])
    data["abundance"] = pd.to_numeric(data["abundance"].str.replace(",", ""), errors='coerce')
    data.dropna(subset=["abundance"], inplace=True)
    
    # Use nlargest instead of sorting for efficiency
    top_5_percent = data.nlargest(int(len(data) * 0.05), "abundance")
    
    # Process and format the results
    result = [f"{row['string_external_id'].replace('511145.', '')}: {row['abundance']}" for _, row in top_5_percent.iterrows()]
    return result

# Combined function to extract top genes with RBS
def extract_top_genes_with_rbs(abundance_file, sequence_file):
    top_abundance_pairs = extract_top_abundance(abundance_file)
    genes_info = extract_genes_info(sequence_file)
    selected_genes_data = {}

    # Process top_abundance_pairs and extract relevant genes info
    for entry in top_abundance_pairs:
        locus_tag = entry.split(': ')[0].replace('511145.', '')
        if locus_tag in genes_info:
            selected_genes_data[locus_tag] = genes_info[locus_tag]

    return selected_genes_data

# Define RBSChooser class with optimizations
class RBSChooser:
    rbs_options = []

    def initiate(self):
        top_genes = extract_top_genes_with_rbs(
            "/Users/krishpatel/UCB_BioE134_GeneDesign/genedesign/511145-WHOLE_ORGANISM-integrated.txt",
            "/Users/krishpatel/UCB_BioE134_GeneDesign/genedesign/sequence.gb"
        )
        
        translator = Translate()
        translator.initiate()

        # Create RBSOptions only once for top genes
        for locus_tag, rbs_cds_info in top_genes.items():
            utr = rbs_cds_info['UTR']
            cds = rbs_cds_info['CDS']
            gene_name = rbs_cds_info['gene']
            first_six_aas = translator.run(cds[:18])
            
            rbs_option = RBSOption(utr=utr, cds=cds, gene_name=gene_name, first_six_aas=first_six_aas)
            self.rbs_options.append(rbs_option)

    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        available_rbs = [rbs for rbs in self.rbs_options if rbs not in ignores]
        
        if not available_rbs:
            raise ValueError("No available RBS options after filtering ignores.")
        
        translator = Translate()
        translator.initiate()

        best_rbs = None
        best_score = float('inf')

        for rbs_option in available_rbs:
            hairpin_count, _ = hairpin_counter(rbs_option.utr + cds)
            peptide_similarity = calculate_edit_distance(translator.run(cds[:18]), rbs_option.first_six_aas)
            score = hairpin_count + peptide_similarity

            if score < best_score:
                best_rbs = rbs_option
                best_score = score

        return best_rbs

# Instantiate and initiate RBSChooser
chooser = RBSChooser()
chooser.initiate()

# Example usage
cds = "ATGGCTAGCAAATACGATTTTACAATATAA"
ignores = set()
selected_rbs_1 = chooser.run(cds, ignores)
print(f"Selected RBS: {selected_rbs_1}")

# Ignore the selected RBS and run again
ignores.add(selected_rbs_1)
selected_rbs_2 = chooser.run(cds, ignores)
print(f"Selected RBS after ignoring {selected_rbs_1}: {selected_rbs_2}")
