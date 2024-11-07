import pytest
from genedesign.checkers.gc_content_checker import gc_content_checker

@pytest.fixture
def gc_content_checker_fixture():
    """
    Fixture to initialize the GC content checker.
    """
    return gc_content_checker

def test_gc_content_within_range(gc_content_checker_fixture):
    """
    Test case for GC content within the specified range (0.4 to 0.6).
    """
    sequence = "ATGCGATAGCGT"  # GC content is about 50%
    result, gc_content = gc_content_checker_fixture(sequence, min_gc=0.4, max_gc=0.6)
    
    assert result == True
    assert 0.4 <= gc_content <= 0.6

def test_gc_content_below_range(gc_content_checker_fixture):
    """
    Test case for GC content below the specified range (less than 0.4).
    """
    sequence = "ATATATATATATATA"  # GC content is low
    result, gc_content = gc_content_checker_fixture(sequence, min_gc=0.4, max_gc=0.6)
    
    assert result == False
    assert gc_content < 0.4

def test_gc_content_above_range(gc_content_checker_fixture):
    """
    Test case for GC content above the specified range (greater than 0.6).
    """
    sequence = "GCGCGCGCGCGC"  # GC content is high
    result, gc_content = gc_content_checker_fixture(sequence, min_gc=0.4, max_gc=0.6)
    
    assert result == False
    assert gc_content > 0.6

def test_empty_sequence(gc_content_checker_fixture):
    """
    Test case for an empty sequence, which should return False and a GC content of 0.
    """
    sequence = ""
    result, gc_content = gc_content_checker_fixture(sequence, min_gc=0.4, max_gc=0.6)
    
    assert result == False
    assert gc_content == 0.0