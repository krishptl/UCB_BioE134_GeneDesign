from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import random
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker 
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.gc_content_checker import gc_content_checker

class TranscriptDesigner:
    def __init__(self, seed: int = 42, use_rbs_only: bool = False):
        """
        Initialize the TranscriptDesigner with codon usage data and sequence validation checkers.
        """
        self.codon_usage_file = "/Users/krishpatel/UCB_BioE134_GeneDesign/genedesign/data/codon_usage.txt"
        self.codon_usage = None  # Initially set to None
        self.forbidden_sequence_checker = None
        self.internal_promoter_checker = None
        self.codon_checker = None
        self.use_rbs_only = use_rbs_only
        random.seed(seed)

        # Instantiate RBSChooser once outside of the loop to avoid repeated instantiation
        self.rbs_chooser = None

        # Ensure codon usage is loaded during initialization
        self.initiate()

    def initiate(self):
        """
        Initializes all necessary components and checks for the TranscriptDesigner.
        """
        if self.codon_usage is None:  # Avoid reloading if it's already loaded
            self.codon_usage = self.load_codon_usage()

        if not self.codon_usage:
            raise ValueError("Codon usage data is not loaded properly.")

        # Initialize checkers
        self.forbidden_sequence_checker = ForbiddenSequenceChecker()
        self.forbidden_sequence_checker.initiate()

        self.internal_promoter_checker = PromoterChecker()
        self.internal_promoter_checker.initiate()

        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()

        # Re-initialize the RBSChooser if necessary
        self.rbs_chooser = RBSChooser()

    def load_codon_usage(self):
        """
        Load codon usage data from the file, optimizing loading to avoid redundant dictionary calls.
        """
        amino_acid_to_codons = {}
        try:
            with open(self.codon_usage_file, 'r') as file:
                for line in file:
                    codon, amino_acid = line.split()[:2]
                    amino_acid_to_codons.setdefault(amino_acid, []).append(codon)
        except FileNotFoundError:
            print(f"Error: Codon usage file not found at {self.codon_usage_file}.")
            return None

        return amino_acid_to_codons
    def run(self, peptide: str, ignores: set, iterations: int = 1000) -> Transcript:
        """
        Generate a valid transcript using a Monte Carlo approach with sliding windows.

        Parameters:
            peptide (str): Target protein sequence.
            ignores (set): RBS options to ignore.
            iterations (int): Number of attempts to optimize sequence.

        Returns:
            Transcript: The optimized transcript with RBS and CDS.
        """
        best_transcript = None
        best_score = -float('inf')  # Start with a very low initial score

        # Monte Carlo loop to iteratively find the best sequence
        for iteration in range(iterations):
            # Generate a random codon sequence for the peptide
            codons = self.generate_random_codon_sequence(peptide)
            cds = ''.join(codons) + "TAA"  # Add stop codon

            # Apply sliding window optimization
            cds_optimized = self.sliding_window_optimization(cds)

            # Validate the current sequence with all checkers
            if self.is_valid_sequence(cds_optimized):
                score = self.evaluate_cds(cds_optimized)

                # If the score is better, update the best score and best transcript
                if score > best_score:
                    best_score = score
                    print(f"New best score: {score} (Iteration {iteration+1})")  # Debugging line
                
                   # Try to find a valid RBS for the CDS sequence
                    rbs = self.rbs_chooser.run(cds_optimized, ignores)
                    if rbs is None:
                        print(f"No available RBS options after filtering ignores for CDS: {cds_optimized}")  # Debugging line
                        continue  # Skip to next iteration if no valid RBS found

                    # Create a valid transcript with the found RBS
                    best_transcript = Transcript(rbs, peptide, cds_optimized)
                else:
                    print(f"Iteration {iteration+1}/{iterations}: Score did not improve.")  # Debugging line
            else:
                print(f"Iteration {iteration+1}/{iterations}: Invalid sequence - {cds_optimized}")  # Debugging line

        # Ensure a valid transcript was found
        if best_transcript:
            return best_transcript
        else:
            raise ValueError("Failed to generate a valid transcript after multiple iterations.")
    def generate_random_codon_sequence(self, peptide: str) -> list:
        """
        Generate a random codon sequence for each amino acid in the peptide.
        """
        return [random.choice(self.codon_usage.get(aa, ['XXX'])) for aa in peptide]

    def sliding_window_optimization(self, cds: str) -> str:
        """
        Apply sliding window optimization to the CDS sequence.
        Iterate over 9 nt windows, optimizing by selecting the best window.

        Parameters:
            cds (str): The CDS sequence to optimize.

        Returns:
            str: The optimized CDS sequence after sliding window processing.
        """
        optimized_cds = list(cds)
        peptide_length = len(cds) // 3

        # Slide through the CDS in 9-nt windows
        for i in range(0, peptide_length - 2):
            # Extract the current 9-nt window (3 amino acids)
            window_start = i * 3
            window_end = window_start + 9
            current_window = cds[window_start:window_end]

            # Rank windows and select the best one (based on forbidden sequences, secondary structure, etc.)
            best_window = self.rank_window(current_window)

            # Replace the current window with the selected one
            optimized_cds[window_start:window_end] = best_window

        return ''.join(optimized_cds)

    def rank_window(self, window: str) -> str:
        """
        Rank a 9-nt window based on forbidden sequences and secondary structure.

        Parameters:
            window (str): A 9-nt window from the CDS.

        Returns:
            str: The best-ranked window after evaluation.
        """
        # Evaluate forbidden sequences
        valid, forbidden_site = self.forbidden_sequence_checker.run(window)
        if not valid:
            print(f"Window contains forbidden sequence at {forbidden_site}")
            return window  # If invalid, return original window

        # Rank the window based on other constraints like secondary structure, hairpins, etc.
        score = self.evaluate_cds(window)  # Use the same evaluation function for a quick score

        # For simplicity, we're returning the window itself as the best option
        # Add more logic if needed to choose the best possible combination
        return window

    def is_valid_sequence(self, cds: str) -> bool:
        """
        Check for forbidden sequences in the entire transcript (CDS and RBS).
        """
        # Check CDS for forbidden sequences
        valid, forbidden_site = self.forbidden_sequence_checker.run(cds)
        if not valid:
            print(f"Invalid CDS due to forbidden sequence at {forbidden_site}")
            return False
        # Check GC content
        within_range, gc_content = gc_content_checker(cds, min_gc=0.4, max_gc=0.6)
        if not within_range:
            print(f"Invalid CDS due to GC content outside the range. GC Content: {gc_content}")
            return False

        # Check RBS sequence as well (if applicable)
        if self.use_rbs_only:
            rbs_sequence = self.rbs_chooser.run(cds, set())
        
            # Ensure rbs_sequence is a Seq object before passing to the forbidden sequence checker
            if rbs_sequence:
                rbs_seq = rbs_sequence.utr  # Extract the UTR (which is the sequence to check)
                valid, forbidden_site = self.forbidden_sequence_checker.run(rbs_seq)
                if not valid:
                    print(f"RBS contains forbidden sequence at {forbidden_site}")
                    return False
    
        return True

    def evaluate_cds(self, cds: str) -> int:
        """
        Evaluate CDS with early-exit checks for severe penalties in hairpins and promoters.
        """
        score = 100

        # Check hairpins
        hairpin_result = hairpin_checker(cds)
        hairpin_count = hairpin_result[0]
        if hairpin_count > 0:
            # Exponentially increase the penalty for high hairpin counts
            score -= (5 * (2 ** (hairpin_count - 1)))  # Double penalty for each additional hairpin
            if hairpin_count > 5:  # Early exit for excessive hairpins
                return score

        # Check internal promoters
        promoter_count, _ = self.internal_promoter_checker.run(cds)
        if promoter_count > 0:
            score -= (3 * (2 ** (promoter_count - 1)))  # Exponentially increase the penalty for more promoters
            if promoter_count > 2:  # Early exit for excessive internal promoters
                return score

        # Codon adaptation index (CAI) penalties
        _, _, rare_codon_count, cai_value = self.codon_checker.run(cds)
        # Weight the CAI penalty higher for worse CAI values
        score -= (rare_codon_count * 2)  # Doubling the penalty to make it more impactful
        return score

if __name__ == "__main__":
    # Define peptide and ignores set for testing
    peptide = "MYPFIRTARARMTEYMIARMTV"
    ignores = set()  # Specify RBS sequences to ignore if any

    # Initialize the designer and run with early-exit optimizations
    designer = TranscriptDesigner(use_rbs_only=True)
    try:
        transcript = designer.run(peptide, ignores)
        print(transcript)
    except ValueError as e:
        print(f"Error: {e}")