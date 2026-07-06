//! De Bruijn graph construction

use rustc_hash::FxHashMap;
use seq_io::fastq::{Reader, Record};
use crate::sequence::{base_to_bits, base_to_bit_mask};
use std::path::Path;
use rayon::prelude::*;
use dashmap::DashMap;
use std::sync::Arc;
use std::time::Instant;
use log::info;
use crate::reader::open_sequence_file;


pub type KmerCount = FxHashMap<u64, u32>;
pub type KmerGraph = FxHashMap<u64, u8>;

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct PathInfo {
    converge_kmer: u64,
    length: usize,
    avg_coverage: u64,
}

/// De Bruijn graph structure
pub struct Graph {
    pub kmers: KmerGraph,
    pub counts: KmerCount,
    pub k: usize,
    pub mask: u64,
}

impl Graph {
    /// Create a new empty graph
    pub fn new(k: usize) -> Self {
        let mask = if k == 32 { !0 } else { (1u64 << (2 * k)) - 1 };
        Graph {
            kmers: FxHashMap::default(),
            counts: FxHashMap::default(),
            k,
            mask,
        }
    }

    pub fn count_kmers<P: AsRef<Path>>(&mut self, path: P, path2: Option<P>, min_count: u32) -> std::io::Result<usize> {
        log::info!("Counting k-mers in parallel...");

        let k = self.k;
        let mask = self.mask;
        let shared_counts: DashMap<u64, u32> = DashMap::with_capacity(1 << 26);
        
        const CHUNK_SIZE: usize = 5_000_000;
        const EARLY_FILTER_THRESHOLD: u32 = 2;

        let start_time = Instant::now();
        
        let count_kmers_from_file = |file_path: &dyn AsRef<Path>| -> std::io::Result<usize> {
            let reader_io = open_sequence_file(file_path)?;
            let mut reader = Reader::new(reader_io);
            let mut sequences: Vec<Vec<u8>> = Vec::new();
            let mut file_read_count = 0usize;

            while let Some(result) = reader.next() {
                let record = result.expect("Error reading record");
                sequences.push(record.seq().to_vec());
                file_read_count += 1;

                if sequences.len() >= CHUNK_SIZE {
                    log::info!("Processing chunk: reads {}-{}", file_read_count - sequences.len(), file_read_count);

                    sequences.par_iter().for_each(|seq| {
                        if seq.len() < k {
                            return;
                        }
                        let mut current_packed: u64 = 0;
                        for i in 0..k {
                            current_packed = (current_packed << 2) | base_to_bits(seq[i]);
                        }
                        *shared_counts.entry(current_packed).or_insert(0) += 1;

                        for i in k..seq.len() {
                            current_packed = ((current_packed << 2) | base_to_bits(seq[i])) & mask;
                            *shared_counts.entry(current_packed).or_insert(0) += 1;
                        }
                    });
                    log::info!("Pre-filtering: {} k-mers in DashMap", shared_counts.len());
                    shared_counts.retain(|_, count| *count >= EARLY_FILTER_THRESHOLD);
                    log::info!("After early filter: {} k-mers retained", shared_counts.len());

                    sequences.clear();
                }
            }

            if !sequences.is_empty() {
                log::info!("Processing final chunk: {} reads", sequences.len());
                
                sequences.par_iter().for_each(|seq| {
                    if seq.len() < k {
                        return;
                    }
                    let mut current_packed: u64 = 0;
                    for i in 0..k {
                        current_packed = (current_packed << 2) | base_to_bits(seq[i]);
                    }
                    *shared_counts.entry(current_packed).or_insert(0) += 1;
        
                    for i in k..seq.len() {
                        current_packed = ((current_packed << 2) | base_to_bits(seq[i])) & mask;
                        *shared_counts.entry(current_packed).or_insert(0) += 1;
                    }
                });

            }
            Ok(file_read_count)
        };
        
        
        // Process first file
        let mut total_read_count = count_kmers_from_file(&path)?;

        // Process second file if provided (paired-end)
        if let Some(p) = path2 {
            let paired_count = count_kmers_from_file(&p)?;
            total_read_count += paired_count;
            log::info!("Paired-end mode: processed {} reads from R2", paired_count);
        }

        let count_time = start_time.elapsed();
        
        log::info!("Read {} total sequences, counting k-mers in parallel...", total_read_count);
        info!(
            "Counting completed in {:.2}s ",
            count_time.as_secs_f64(),
        );
        
        // Convert and filter
        let filter_start = std::time::Instant::now();
        log::info!("Total k-mers before filtering: {}", shared_counts.len());

        let entries: Vec<(u64, u32)> = shared_counts.into_iter().collect();
        log::info!("Collected {} entries, now filtering in parallel...", entries.len());

        // Parallel filter and fold into FxHashMap chunks, then reduce
            self.counts = entries
            .into_par_iter()
            .filter(|(_, count)| *count >= min_count)
            .fold(
                || FxHashMap::default(),
                |mut map, (k, v)| {
                    map.insert(k, v);
                    map
                }
            )
            .reduce(
                || FxHashMap::default(),
                |mut map1, map2| {
                    map1.extend(map2);
                    map1
                }
            );

        let filter_time = filter_start.elapsed();
        log::info!("K-mers after filtering (min_count={}): {} [took {:.2}s]", 
            min_count, self.counts.len(), filter_time.as_secs_f64());
            
            Ok(total_read_count)  
    }

    /// Build de Bruijn graph from FASTQ file (Pass 2)
    pub fn build_graph<P: AsRef<Path>>(&mut self, path: P, path2: Option<P>, min_count: u32) -> std::io::Result<()> {
        log::info!("Pass 2: Building de Bruijn graph");

         // Helper function to build graph from a single file
         let build_graph_from_file = |file_path: &dyn AsRef<Path>, kmers: &mut KmerGraph, counts: &KmerCount, k: usize, mask: u64, min_count: u32| -> std::io::Result<()> {
            let reader_io = open_sequence_file(file_path)?;
            let mut reader = Reader::new(reader_io);
            let mut sequence_count = 0usize;

            while let Some(result) = reader.next() {
                let record = result.expect("Error reading record");
                let seq = record.seq();
                sequence_count += 1;

            if seq.len() < k + 1 {
                continue;
            }

            let mut current_packed: u64 = 0;
                for i in 0..k {
                    current_packed = (current_packed << 2) | base_to_bits(seq[i]);
                }
                
                for i in k..seq.len() {
                    let next_base = seq[i];
                    let next_packed = ((current_packed << 2) | base_to_bits(next_base)) & mask;

                    if *counts.get(&current_packed).unwrap_or(&0) >= min_count
                        && *counts.get(&next_packed).unwrap_or(&0) >= min_count
                    {
                        let entry = kmers.entry(current_packed).or_insert(0);
                        *entry |= base_to_bit_mask(next_base);
                    }
                    current_packed = next_packed;
                }
            }
            log::info!("Processed {} sequences from file", sequence_count);
            Ok(())
        };

    // Process first file
    build_graph_from_file(&path, &mut self.kmers, &self.counts, self.k, self.mask, min_count)?;

    // Process second file if provided (paired-end)
    if let Some(p) = path2 {
        log::info!("Processing paired-end R2 file...");
        build_graph_from_file(&p, &mut self.kmers, &self.counts, self.k, self.mask, min_count)?;

    }
    log::info!("Graph size: {} nodes", self.kmers.len());
    Ok(())
    }

    /// Prune tips from the graph (short dead-end branches)
    pub fn prune_tips(&mut self) -> usize {
        log::info!("Pruning tips from graph...");

        let max_tip_length = 2 * self.k;
        let mut tips_removed = 0;
        let mut to_update = Vec::new();

        for (&kmer, &out_mask) in self.kmers.iter() {
            if out_mask.count_ones() > 1 {
                let mut new_mask = out_mask;

                for &base in &[b'A', b'C', b'G', b'T'] {
                    let base_bit = base_to_bit_mask(base);
                    if (out_mask & base_bit) != 0 && self.is_tip(kmer, base, max_tip_length) {
                        new_mask &= !base_bit;
                        tips_removed += 1;
                    }
                }

                if new_mask != out_mask {
                    to_update.push((kmer, new_mask));
                }
            }
        }

        for (kmer, new_mask) in to_update {
            self.kmers.insert(kmer, new_mask);
        }

        log::info!("Removed {} tips", tips_removed);
        tips_removed
    }

    /// Check if a branch from a k-mer is a tip (dead-end)
    fn is_tip(&self, start_kmer: u64, first_base: u8, max_len: usize) -> bool {
        use crate::sequence::{base_to_bits, bit_mask_to_base};

        let mut curr = ((start_kmer << 2) | base_to_bits(first_base)) & self.mask;
        for _ in 0..max_len {
            let mask = *self.kmers.get(&curr).unwrap_or(&0);
            let out_deg = mask.count_ones();

            if out_deg == 0 {
                return true;
            }
            if out_deg > 1 {
                return false;
            }

            let next_base = bit_mask_to_base(mask);
            curr = ((curr << 2) | base_to_bits(next_base as u8)) & self.mask;
        }
        false
    }

    /// Calculate in-degrees for all k-mers
    pub fn calculate_in_degrees(&self) -> FxHashMap<u64, u32> {
        use crate::sequence::{base_to_bit_mask, base_to_bits};

        log::debug!("Calculating in-degrees...");
        let mut in_degrees: FxHashMap<u64, u32> =
            FxHashMap::with_capacity_and_hasher(self.kmers.len(), Default::default());

        for (&kmer, &out_mask) in self.kmers.iter() {
            for &base in &[b'A', b'C', b'G', b'T'] {
                if (out_mask & base_to_bit_mask(base)) != 0 {
                    let next_kmer = ((kmer << 2) | base_to_bits(base)) & self.mask;
                    *in_degrees.entry(next_kmer).or_insert(0) += 1;
                }
            }
        }

        in_degrees
    }

/// Remove bubbles from the graph (small divergences that reconverge)
pub fn remove_bubbles(&mut self) -> usize {
    use crate::sequence::{base_to_bit_mask, base_to_bits};

    log::info!("Removing bubbles from graph...");

    let mut bubbles_removed = 0;
    let mut to_remove = Vec::new();

    // Find nodes with multiple outgoing edges
    for (&start_kmer, &out_mask) in self.kmers.iter() {
        let out_deg = out_mask.count_ones();
        if out_deg <= 1 {
            continue;
        }

        // Collect coverage of each branch
        let mut branches = Vec::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            let bit = base_to_bit_mask(base);
            if (out_mask & bit) == 0 {
                continue;
            }

            let next_kmer = ((start_kmer << 2) | base_to_bits(base)) & self.mask;
            let coverage = *self.counts.get(&next_kmer).unwrap_or(&0) as u64;
            branches.push((bit, coverage));
        }

        // If significant coverage difference, remove low-coverage branch
        if branches.len() >= 2 {
            let max_coverage = branches.iter().map(|b| b.1).max().unwrap_or(0);
            if max_coverage > 0 {
                for (bit, cov) in branches {
                    // Remove if coverage is <30% of max and count is low
                    if cov < max_coverage / 3 && cov < 10 {
                        to_remove.push((start_kmer, bit));
                        bubbles_removed += 1;
                    }
                }
            }
        }
    }

    // Remove marked edges
    for (kmer, bit) in to_remove {
        if let Some(entry) = self.kmers.get_mut(&kmer) {
            *entry &= !bit;
        }
    }

    log::info!("Removed {} bubbles", bubbles_removed);
    bubbles_removed
}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_creation() {
        let graph = Graph::new(31);
        assert_eq!(graph.k, 31);
        assert_eq!(graph.kmers.len(), 0);
        assert_eq!(graph.counts.len(), 0);
    }
}