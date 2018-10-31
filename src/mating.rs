//! Mate and merge overlapping paired sequences
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::mating::{mate,merge};
//! let read1 = b"tacgattcgat";
//! let read2 = b"ttcgattacgt";
//! let overlap = mate(read1, read2, 1, 1).unwrap();
//! let contig = merge(read1, read2, overlap);
//! assert_eq!(contig, b"tacgattcgattacgt");
//! ```
//!
//! Mate pair merging proceeds in two steps:
//! * mate: identifying the optimal overlap length or rejecting the case
//! * merge: combines overlapping sequences
//!
//! Mating is governed by an objective function and merging resolves
//! conflicts between sequences.
//!
//! We also need to define the scores for which mating is acceptable. The
//! minimum overlap score and minimum overlap length could be guessed
//! from the k-mer distribution. 25bp is the apparent standard.

use std::cmp;

/// Determine the index of overlap for two reads.
pub fn mate(r1: &[u8], r2: &[u8], overlap_bound: usize, min_score: i16) -> Option<usize> {
    let max_overlap = cmp::min(r1.len(), r2.len());
    let min_overlap = overlap_bound - 1;
    let mut m: i16 = 0;
    let mut overlap: usize = 0;
    for i in min_overlap..max_overlap {
        let h = score(&r2[0..i], &r1[r1.len() - i..r1.len()]);
        if h > m {
            m = h;
            overlap = i;
        }
    }
    if m >= min_score {
        return Some(overlap);
    }
    None
}

/// Mating with the Hamming rate objective function. Complexity: O(n^2).
#[allow(dead_code)]
pub fn mate_hamming_rate(
    r1: &[u8],
    r2: &[u8],
    overlap_bound: usize,
    min_score: i16,
) -> Option<usize> {
    let max_overlap = cmp::min(r1.len(), r2.len());
    let min_overlap = overlap_bound - 1;
    let mut m: i16 = 0;
    let mut overlap: usize = 0;
    for i in min_overlap..max_overlap {
        //        let h = hamming(&r2[0..i], &r1[r1.len()-i..r1.len()]) / i;
        let h = 0;
        if h > m {
            m = h;
            overlap = i;
        }
    }
    if m >= min_score {
        return Some(overlap);
    }
    None
}

/// Mating objective function that penalizes mismatches.
fn score(r1: &[u8], r2: &[u8]) -> i16 {
    let mut s: i16 = 0;
    let len = r1.len();
    for i in 0..len {
        if r1[i] == r2[i] {
            s = s + 1;
        } else {
            s = s - 1;
        }
    }
    s
}

/// Function that replaces disagreeing reads with 'N'.
fn mend_consensus(a: u8, b: u8) -> u8 {
    if a == b {
        a
    } else {
        b'N'
    }
}

/// Given two reads and the index of overlap, merge them together.
pub fn merge(r1: &[u8], r2: &[u8], overlap: usize) -> Vec<u8> {
    let r1_end = r1.len() - overlap;
    let r2_end = r2.len() - overlap;
    let len = r1_end + overlap + r2_end;

    let mut seq = vec![0; len];
    seq[0..r1_end].copy_from_slice(&r1[0..r1_end]);

    // decide what to do for the overlapping part
    for i in 0..overlap {
        seq[r1_end + i] = mend_consensus(r1[r1_end + i], r2[i]);
    }
    seq[(r1_end + overlap)..len].copy_from_slice(&r2[overlap..r2.len()]);
    seq
}

/// Mend and return the overlapping region of two reads, given an index of overlap.
pub fn truncate(r1: &[u8], r2: &[u8], overlap: usize) -> Vec<u8> {
    let r2_end = r2.len();

    let mut seq = vec![0; overlap];
    for i in 0..overlap {
        seq[(overlap - i) - 1] = mend_consensus(r1[(overlap - i) - 1], r2[(r2_end - i) - 1]);
    }
    seq
}

#[cfg(test)]
mod tests {
    use super::{mate, merge, truncate};

    #[test]
    fn test_mate_pair() {
        let r1 = b"tacgattcgat";
        let r2 = b"ttcgattacgt";
        let offset = mate(r1, r2, 3, 3).unwrap();
        assert_eq!(offset, 6);
    }

    #[test]
    fn test_merge_pair() {
        let r1 = b"tacgattcgat";
        let r2 = b"ttcgattacgt";
        assert_eq!(merge(r1, r2, 6), b"tacgattcgattacgt");
    }

    #[test]
    fn test_merge_consensus() {
        let r1 = b"actgtagtacaccatgatg";
        let r2 = b"gtttaccatgatggattga";
        let offset = mate(r1, r2, 3, 3).unwrap();
        assert_eq!(offset, 13);
        assert_eq!(merge(r1, r2, offset), b"actgtagtNNaccatgatggattga");
    }

    #[test]
    fn test_overlap_bounds() {
        let r1 = b"ctagtcagctagcatgcatgctcgcgtatacgctagc";
        let r2 = b"gtatacgctagccgctagcatgcat";
        assert_eq!(mate(r1, r2, 12, 12).unwrap(), 12);
        assert_eq!(mate(r1, r2, 13, 13), None);
    }

    #[test]
    fn test_score_bounds() {
        let r1 = b"there is no restriction to DNA strings";
        let r2 = b"ctiin ot dNA atrings asdf qwer";
        assert_eq!(mate(r1, r2, 12, 10).unwrap(), 20);
        assert_eq!(mate(r1, r2, 12, 11), None);
    }

    #[test]
    fn test_truncation() {
        let r1 = b"cgctgtcatgc";
        let r2 = b"tcatgccgctgt";
        assert_eq!(truncate(r1, r2, 6), b"cgctgt");
    }

    #[test]
    fn test_overlap() {}

    #[test]
    fn test_readthrough() {}

    #[test]
    fn test_disjoint() {}
}
