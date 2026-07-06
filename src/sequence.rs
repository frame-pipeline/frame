//! Sequence utilities for k-mer operations and genetic code translation.

use wyhash::wyhash;

/// Convert a DNA base to its 2-bit representation (A=0, C=1, G=2, T=3)
pub fn base_to_bits(base: u8) -> u64 {
    match base.to_ascii_uppercase() {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0, // A=0, N/others=A
    }
}

/// Convert a DNA base to its bitmask representation (A=1, C=2, G=4, T=8)
pub fn base_to_bit_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

/// Convert a bitmask back to a DNA base character
pub fn bit_mask_to_base(mask: u8) -> char {
    if mask & 1 != 0 {
        'A'
    } else if mask & 2 != 0 {
        'C'
    } else if mask & 4 != 0 {
        'G'
    } else {
        'T'
    }
}

/// Unpack a 64-bit packed k-mer into its string representation
pub fn unpack_kmer(mut packed: u64, k: usize) -> String {
    let mut chars = vec![0u8; k];
    for i in (0..k).rev() {
        let val = packed & 0b11;
        chars[i] = match val {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
packed >>= 2;
}
String::from_utf8(chars).unwrap()
}

/// Pack a raw byte slice into a u64 k-mer
pub fn pack_kmer(seq: &[u8], k: usize, mask: u64) -> u64 {
let mut packed: u64 = 0;
for &b in seq.iter().take(k) {
    let bits = match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    };
    packed = ((packed << 2) | bits) & mask;
}
packed
}

/// Compute reverse complement of a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
seq.iter()
    .rev()
    .map(|b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    })
    .collect()
}

/// Compute canonical k-mer (lexicographically smaller of forward/reverse complement)
pub fn canonical_kmer(kmer: &[u8]) -> Vec<u8> {
let rc = reverse_complement(kmer);
if kmer <= rc.as_slice() {
    kmer.to_vec()
} else {
    rc
}
}

/// Generate minimizers from a sequence
pub fn generate_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<u64> {
let hashes: Vec<u64> = seq
    .windows(k)
    .map(|kmer| {
        let canon = canonical_kmer(kmer);
        wyhash(&canon, 0)
    })
    .collect();

let mut minimizers: Vec<u64> = hashes
    .windows(w)
    .map(|window| *window.iter().min().unwrap())
    .collect();

minimizers.dedup();
minimizers
}

/// Translate a codon to amino acid
pub fn translate_codon(codon: &[u8]) -> char {
match codon {
    b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
    b"TGT" | b"TGC" => 'C',
    b"GAT" | b"GAC" => 'D',
    b"GAA" | b"GAG" => 'E',
    b"TTT" | b"TTC" => 'F',
    b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
    b"CAT" | b"CAC" => 'H',
    b"ATT" | b"ATC" | b"ATA" => 'I',
    b"AAA" | b"AAG" => 'K',
    b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
    b"ATG" => 'M',
    b"AAT" | b"AAC" => 'N',
    b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
    b"CAA" | b"CAG" => 'Q',
    b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
    b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
    b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
    b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
    b"TGG" => 'W',
    b"TAT" | b"TAC" => 'Y',
    b"TAA" | b"TAG" | b"TGA" => '*',
    _ => 'X',
}
}

/// Extract ORFs from a sequence with strand and frame information
pub fn get_orfs_from_sequence(seq_id: usize, dna: &str, min_aa_len: usize) -> Vec<(String, String)> {
let mut orfs = Vec::new();
let rev_dna = reverse_complement_str(dna);
let sequences = [dna, &rev_dna];

for (strand_idx, seq) in sequences.iter().enumerate() {
    let bytes = seq.as_bytes();
    let seq_len = bytes.len();
    let strand_symbol = if strand_idx == 0 { '+' } else { '-' };

    for frame_offset in 0..3 {
        let mut current_aa_seq = String::new();

        for i in (frame_offset..seq_len).step_by(3) {
            if i + 3 <= seq_len {
                let codon = &bytes[i..i + 3];
                current_aa_seq.push(translate_codon(codon));
            }
        }

        let mut current_orf = String::new();
        let mut in_orf = false;
        let mut orf_counter = 0;

        for ch in current_aa_seq.chars() {
            if ch == 'M' && !in_orf {
                in_orf = true;
            }

            if in_orf {
                current_orf.push(ch);
                if ch == '*' {
                    if current_orf.len() >= min_aa_len {
                        let header = format!(
                            ">seq_{}_strand_{}_frame_{}_orf_{}",
                            seq_id, strand_symbol, frame_offset, orf_counter
                        );
                        orfs.push((header, current_orf.clone()));
                        orf_counter += 1;
                    }
                    current_orf.clear();
                    in_orf = false;
                }
            }
        }

        if in_orf && current_orf.len() >= min_aa_len {
            let header = format!(
                ">seq_{}_strand_{}_frame_{}_orf_{}_partial",
                seq_id, strand_symbol, frame_offset, orf_counter
            );
            orfs.push((header, current_orf));
        }
    }
}
orfs
}

/// Reverse complement a DNA string
fn reverse_complement_str(dna: &str) -> String {
dna.bytes()
    .rev()
    .map(|b| match b {
        b'A' => 'T',
        b'T' => 'A',
        b'C' => 'G',
        b'G' => 'C',
        _ => 'N',
    })
    .collect()
}

#[cfg(test)]
mod tests {
use super::*;

#[test]
fn test_base_to_bits() {
    assert_eq!(base_to_bits(b'A'), 0);
    assert_eq!(base_to_bits(b'C'), 1);
    assert_eq!(base_to_bits(b'G'), 2);
    assert_eq!(base_to_bits(b'T'), 3);
}

#[test]
fn test_kmer_packing() {
    let kmer = b"ACGT";
    let mask = (1u64 << 8) - 1;
    let packed = pack_kmer(kmer, 4, mask);
    let unpacked = unpack_kmer(packed, 4);
    assert_eq!(unpacked, "ACGT");
}
}
