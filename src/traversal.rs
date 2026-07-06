//! Unitig extraction with frame-aware traversal.

use rustc_hash::{FxHashMap, FxHashSet};
use crate::sequence::{base_to_bit_mask, base_to_bits, bit_mask_to_base, unpack_kmer};
use crate::graph::{KmerCount, KmerGraph};

/// Statistics from traversal
#[derive(Debug, Clone)]
pub struct TraversalStats {
    pub total_unitigs: usize,
    pub total_junctions: usize,
    pub frame_wins: usize,
    pub count_wins: usize,
}

/// Extract unitigs with frame-aware junction resolution
pub fn extract_unitigs_frame_aware(
    graph: &KmerGraph,
    counts: &KmerCount,
    k: usize,
    min_unitig_length: usize,
    use_frame: bool
)  -> (Vec<(String, bool)>, TraversalStats)  {
    let mask_64: u64 = if k == 32 { !0 } else { (1 << (2 * k)) - 1 };

    log::info!("Calculating in-degrees...");
    let mut in_degrees: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(graph.len(), Default::default());

    for (&kmer, &out_mask) in graph.iter() {
        for &base in &[b'A', b'C', b'G', b'T'] {
            if (out_mask & base_to_bit_mask(base)) != 0 {
                let next_kmer = ((kmer << 2) | base_to_bits(base)) & mask_64;
                *in_degrees.entry(next_kmer).or_insert(0) += 1;
            }
        }
    }

    let mut visited: FxHashSet<u64> = FxHashSet::default();
    //let mut unitigs = Vec::new();
    let mut unitigs: Vec<(String, bool)> = Vec::new();
    let mut stats = TraversalStats {
        total_unitigs: 0,
        total_junctions: 0,
        frame_wins: 0,
        count_wins: 0,
    };

    log::info!("Starting frame-aware traversal...");

    for (&start_kmer, &out_mask_start) in graph.iter() {
        if visited.contains(&start_kmer) {
            continue;
        }

        let in_deg = *in_degrees.get(&start_kmer).unwrap_or(&0);
        let out_deg = out_mask_start.count_ones();

        if in_deg != 1 || out_deg != 1 {
            let mut current_path = unpack_kmer(start_kmer, k);
            let mut current_kmer = start_kmer;
            let mut current_frame = detect_frame(&current_path);
            let mut unitig_used_frame = false;
            visited.insert(current_kmer);

            loop {
                let local_mask = *graph.get(&current_kmer).unwrap_or(&0);
                let current_out_deg = local_mask.count_ones();

                if current_out_deg == 0 {
                    break;
                }

                if current_out_deg == 1 {
                    let next_base = bit_mask_to_base(local_mask);
                    let next_kmer = ((current_kmer << 2) | base_to_bits(next_base as u8)) & mask_64;

                    if visited.contains(&next_kmer) {
                        break;
                    }

                    current_path.push(next_base);
                    current_frame = (current_frame + 1) % 3;
                    visited.insert(next_kmer);
                    current_kmer = next_kmer;
                } else {
                    stats.total_junctions += 1;

                    let mut best_score = (-1, 0u32);
                    let mut best_next: Option<(u64, char)> = None;
                    let mut used_frame = false;

                    for &base_char in &['A', 'C', 'G', 'T'] {
                        let bit = base_to_bit_mask(base_char as u8);
                        if (local_mask & bit) == 0 {
                            continue;
                        }

                        let next_kmer = ((current_kmer << 2) | base_to_bits(base_char as u8)) & mask_64;

                        if visited.contains(&next_kmer) {
                            continue;
                        }
                        let count = *counts.get(&next_kmer).unwrap_or(&0);
                        let branch_phase_score = if use_frame {
                            measure_codon_phase_consistency(
                            &current_path,
                            base_char as char,
                            next_kmer,
                            graph,
                            current_frame,
                            mask_64,
                        )
                        } else {
                            0
                        };

                        let score = (branch_phase_score, count);

                        if score.0 > best_score.0
                            || (score.0 == best_score.0 && score.1 > best_score.1)
                        {
                            best_score = score;
                            best_next = Some((next_kmer, base_char as char));
                            used_frame = score.0 > 0;
                        }
                    }

                    match best_next {
                        None => break,
                        Some((next_kmer, next_base)) => {
                            if used_frame {
                                stats.frame_wins += 1;
                                unitig_used_frame = true;
                            } else {
                                stats.count_wins += 1;
                            }
                            current_path.push(next_base);
                            current_frame = (current_frame + 1) % 3;
                            visited.insert(next_kmer);
                            current_kmer = next_kmer;
                        }
                    }
                }
            }

            if current_path.len() >= min_unitig_length {
                stats.total_unitigs += 1;
                //unitigs.push(current_path);
                unitigs.push((current_path, unitig_used_frame));
            }
        }
    }

    log::info!("Traversal complete: {} unitigs extracted", stats.total_unitigs);
    (unitigs, stats)
}

/// Detect dominant reading frame from a sequence
fn detect_frame(seq: &str) -> u8 {
    const STOPS: &[&str] = &["TAA", "TAG", "TGA"];
    let bytes = seq.as_bytes();

    let mut best_frame = 0u8;
    let mut best_score = f64::INFINITY;

    for frame in 0..3 {
        let mut stop_count = 0;
        let mut codon_count = 0;
        let mut i = frame;

        while i + 3 <= bytes.len() {
            codon_count += 1;
            let codon = std::str::from_utf8(&bytes[i..i + 3]).unwrap_or("NNN");
            if STOPS.contains(&codon) {
                stop_count += 1;
            }
            i += 3;
        }

        if codon_count == 0 {
            continue;
        }

        let stop_ratio = stop_count as f64 / codon_count as f64;
        if stop_ratio < best_score {
            best_score = stop_ratio;
            best_frame = frame as u8;
        }
    }

    best_frame
}

/// Measure codon phase consistency for a candidate branch
fn measure_codon_phase_consistency(
    current_path: &str,
    next_base: char,
    next_kmer: u64,
    graph: &KmerGraph,
    current_frame: u8,
    mask_64: u64,
) -> i32 {
    let mut lookahead = String::new();
    let mut curr = next_kmer;

    for _ in 0..15 {
        match graph.get(&curr) {
            Some(&out_edges) if out_edges != 0 => {
                let next_b = bit_mask_to_base(out_edges);
                lookahead.push(next_b);
                curr = ((curr << 2) | base_to_bits(next_b as u8)) & mask_64;
            }
            _ => break,
        }
    }

    let mut candidate = String::new();
    if current_path.len() >= 2 {
        candidate.push_str(&current_path[current_path.len() - 2..]);
    }
    candidate.push(next_base);
    candidate.push_str(&lookahead);

    if candidate.len() < 9 {
        return 0;
    }

    const STOPS: &[&str] = &["TAA", "TAG", "TGA"];
    let bytes = candidate.as_bytes();
    let next_frame = (current_frame + 1) % 3;

    let mut stop_count = 0;
    let mut total_codon_count = 0;
    let mut i = next_frame as usize;

    while i + 3 <= bytes.len() {
        total_codon_count += 1;
        let codon = std::str::from_utf8(&bytes[i..i + 3]).unwrap_or("NNN");
        if STOPS.contains(&codon) {
            stop_count += 1;
        }
        i += 3;
    }

    if total_codon_count == 0 {
        return 0;
    }

    let stop_ratio = stop_count as f64 / total_codon_count as f64;

    if stop_ratio < 0.1 {
        2
    } else if stop_ratio < 0.3 {
        1
    } else {
        0
    }
}

/// Calculate ORF heuristic score for a sequence
pub fn orf_heuristic_score(seq: &str) -> f64 {
    const STOPS: &[&str] = &["TAA", "TAG", "TGA"];
    let bytes = seq.as_bytes();
    let len = bytes.len();

    let mut best_frame_score: f64 = 0.0;

    for frame in 0..3usize {
        let mut stop_count = 0usize;
        let mut codon_count = 0usize;
        let mut i = frame;

        while i + 3 <= len {
            codon_count += 1;
            let codon = &seq[i..i + 3];
            if STOPS.contains(&codon) {
                stop_count += 1;
            }
            i += 3;
        }

        if codon_count == 0 {
            continue;
        }

        let frame_score = 1.0 - (stop_count as f64 / codon_count as f64);
        if frame_score > best_frame_score {
            best_frame_score = frame_score;
        }
    }

    best_frame_score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orf_heuristic() {
        let coding = "ATGATGATGATGATGATGATGA";
        let score = orf_heuristic_score(coding);
        assert!(score > 0.5);

        let random = "NNNNNNNNNNNNNNNNNNNNNNNN";
        let score = orf_heuristic_score(random);
        assert!(score < 0.5);
    }
}