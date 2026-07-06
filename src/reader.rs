use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use flate2::read::GzDecoder;

pub enum SequenceReader {
    Uncompressed(BufReader<File>),
    Gzipped(GzDecoder<BufReader<File>>),
}

impl Read for SequenceReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            SequenceReader::Uncompressed(r) => r.read(buf),
            SequenceReader::Gzipped(r) => r.read(buf),
        }
    }
}

pub fn open_sequence_file<P: AsRef<Path>>(path: P) -> std::io::Result<SequenceReader> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)?;
    let buf_reader = BufReader::new(file);

    if path_ref
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s == "gz")
        .unwrap_or(false)
    {
        Ok(SequenceReader::Gzipped(GzDecoder::new(buf_reader)))
    } else {
        Ok(SequenceReader::Uncompressed(buf_reader))
    }
}