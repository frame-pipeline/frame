//! Gene prediction Module

use frag_gene_scan_rs::hmm::{Global, Local};
use frag_gene_scan_rs::dna::Nuc;
use std::path::PathBuf;
use log::info;

/// Gene prediction statistics
#[derive(Debug, Clone)]
pub struct PredictionStats {
    pub sequences_processed: usize,
    pub coding_sequences: usize,
}

/// Load HMM model for gene prediction
pub fn load_hmm_model(
    train_dir: PathBuf,
    model_name: PathBuf,
) -> Result<(Box<Global>, Vec<Local>), Box<dyn std::error::Error>> {
    log::info!("Loading HMM model from {:?} / {:?}", train_dir, model_name);
    let (global, locals) = frag_gene_scan_rs::hmm::get_train_from_file(train_dir, model_name)?;
    log::info!("HMM model loaded successfully");
    Ok((global, locals))
}

/// Run gene prediction on a sequence
pub fn predict_genes(
    seq: &[u8],
    header: Vec<u8>,
    global: &Box<Global>,
    locals: &Vec<Local>,
) -> frag_gene_scan_rs::gene::ReadPrediction {
        let nseq: Vec<Nuc> = seq
        .iter()
        .copied()
        .map(Nuc::from)
        .collect();

    frag_gene_scan_rs::viterbi::viterbi(global, locals, header, nseq, false)
}

/// Run gene prediction on a unitig and write results to buffers
pub fn predict_and_write_unitig(
    unitig_id: usize,
    sequence: &str,
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
) -> Result<bool, Box<dyn std::error::Error>> {
    let header = format!("unitig_{}", unitig_id).into_bytes();

    let prediction = predict_genes(sequence.as_bytes(), header, global, locals);

    if !prediction.genes.is_empty() {
        prediction.gff(gff_buffer)?;
        prediction.protein(aa_buffer, false)?;
        prediction.dna(dna_buffer, false)?;
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Run gene prediction on a raw read
pub fn predict_and_write_read(
    read_id: usize,
    sequence: &[u8],
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
) -> Result<bool, Box<dyn std::error::Error>> {
    let header = format!("read_{}", read_id).into_bytes();

    let prediction = predict_genes(sequence, header, global, locals);

    if !prediction.genes.is_empty() {
        prediction.gff(gff_buffer)?;
        prediction.protein(aa_buffer, false)?;
        prediction.dna(dna_buffer, false)?;
        Ok(true)
    } else {
        Ok(false)
    }
}

pub fn predict_and_write_read_with_header(
    read_id: usize,
    original_header: &[u8],
    sequence: &[u8],
    global: &Box<Global>,
    locals: &Vec<Local>,
    gff_buffer: &mut Vec<u8>,
    aa_buffer: &mut Vec<u8>,
    dna_buffer: &mut Vec<u8>,
    header_map_buffer: &mut Vec<u8>,
) -> Result<bool, Box<dyn std::error::Error>> {
    let header = format!("read_{}", read_id).into_bytes();

    let prediction = predict_genes(sequence, header, global, locals);

    if !prediction.genes.is_empty() {
        prediction.gff(gff_buffer)?;
        prediction.protein(aa_buffer, false)?;
        prediction.dna(dna_buffer, false)?;

        // Store read mapping (for benchmarking): original_header read_id
        let original_header_str = String::from_utf8_lossy(original_header);
        let mapping_line = format!("{}\t{}\n", original_header_str, read_id);
        header_map_buffer.extend_from_slice(mapping_line.as_bytes());

        Ok(true)
    } else {
        Ok(false)
    }
}