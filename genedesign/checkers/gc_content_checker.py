# File: gc_content_checker.py

def gc_content_checker(sequence, min_gc=0.4, max_gc=0.6):
    """
    Checks if the GC content of a DNA sequence is within a specified range.
    """
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    print(f"Sequence: {sequence}, GC Content: {gc_content}, Min GC: {min_gc}, Max GC: {max_gc}")  # Debug output
    
    is_within_range = min_gc <= gc_content <= max_gc
    print(f"Result: {is_within_range}")  # Debug output

    return is_within_range, gc_content

# Example usage
if __name__ == "__main__":
    result, gc = gc_content_checker("ATGCGCGATCGGCGCTA", min_gc=0.4, max_gc=0.6)
    print("Within Range:", result, "| GC Content:", gc)