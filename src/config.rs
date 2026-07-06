use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Pipeline configuration parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    pub k: usize,
    pub min_count: u32,
    pub min_unitig_length: usize,
    pub hmm_train_dir: PathBuf,
    pub hmm_model: PathBuf,
    pub input_file: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_file_2: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub rescue_batch_size: usize,
}

impl PipelineConfig {
    /// Create default configuration
    pub fn default_for_file(input_file: PathBuf) -> Self {
        PipelineConfig {
            k: 31,
            min_count: 2,
            min_unitig_length: 100,
            hmm_train_dir: PathBuf::from("./lib/FragGeneScanRs/train"),
            hmm_model: PathBuf::from("illumina_5"),
            input_file,
            input_file_2: None,
            output_dir: PathBuf::from("./frame_output"),
            rescue_batch_size: 50_000,
        }
    }

    /// Create configuration from file
    pub fn from_file(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let contents = std::fs::read_to_string(path)?;
        let config: PipelineConfig = toml::from_str(&contents)?;
        log::info!("Configuration loaded from {}", path);
        Ok(config)
    }

    /// Save configuration to file
    pub fn save(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let contents = toml::to_string_pretty(self)?;
        std::fs::write(path, contents)?;
        log::info!("Configuration saved to {}", path);
        Ok(())
    }

    /// Validate configuration
    pub fn validate(&self) -> Result<(), String> {
        if self.k < 15 || self.k > 63 {
            return Err("k must be between 15 and 63".to_string());
        }

        if self.min_count < 1 {
            return Err("min_count must be at least 1".to_string());
        }

        if self.min_unitig_length < 50 {
            return Err("min_unitig_length should be at least 50".to_string());
        }

        if !self.input_file.exists() {
            return Err(format!("input_file not found: {}", self.input_file.display()));
        }
        if let Some(ref file2) = self.input_file_2 {
            if !file2.exists() {
                return Err(format!("input_file_2 not found: {}", file2.display()));
            }
        }

        if !self.hmm_train_dir.exists() {
            return Err(format!("hmm_train_dir not found: {}", self.hmm_train_dir.display()));
        }

        Ok(())
    }

    /// Get mask for k-mer packing
    pub fn get_mask(&self) -> u64 {
        if self.k == 32 { !0 } else { (1u64 << (2 * self.k)) - 1 }
    }
}