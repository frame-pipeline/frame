use rustc_hash::FxHashMap;

#[derive(Debug, Clone)]
pub struct ReadSummary {
    pub total_reads: usize,
    pub reads_assembled: usize,
    pub reads_rescued: usize,
}

impl ReadSummary {
    /// Get rescue rate as percentage
    pub fn rescue_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            (self.reads_rescued as f64 / self.total_reads as f64) * 100.0
        }
    }
}

/// Check if a read is represented in the graph
pub fn is_read_assembled(
    seq: &[u8],
    k: usize,
    mask: u64,
    counts: &FxHashMap<u64, u32>,
    min_count: u32,
) -> bool {
    use crate::sequence::pack_kmer;

    if seq.len() < k {
        return false;
    }

    let check_positions = [
        0,
        seq.len().saturating_sub(k) / 2,
        seq.len().saturating_sub(k),
    ];

    for &start in &check_positions {
        if start + k > seq.len() {
            continue;
        }
        let kmer = pack_kmer(&seq[start..start + k], k, mask);
        if *counts.get(&kmer).unwrap_or(&0) < min_count {
            return false;
        }
    }

    true
}

pub fn format_file_size(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB"];
    let mut size = bytes as f64;
    let mut unit_idx = 0;

    while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
        size /= 1024.0;
        unit_idx += 1;
    }

    format!("{:.2} {}", size, UNITS[unit_idx])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_summary_rescue_rate() {
        let summary = ReadSummary {
            total_reads: 1000,
            reads_assembled: 800,
            reads_rescued: 200,
        };
        assert_eq!(summary.rescue_rate(), 20.0);
    }

    #[test]
    fn test_format_file_size() {
        assert_eq!(format_file_size(1024), "1.00 KB");
        assert_eq!(format_file_size(1048576), "1.00 MB");
    }
}