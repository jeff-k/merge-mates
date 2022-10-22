extern crate bio_streams;
extern crate flate2;

extern crate bio_seq;
extern crate clap;

use bio_streams::fastq::{Fastq, Record};
use std::io::{stdout, BufReader};

use flate2::read::MultiGzDecoder;
use std::fs::File;

use bio_seq::codec::dna::Dna;
use bio_seq::Seq;

use clap::Parser;

use std::path::{Path, PathBuf};

use std::cmp;

mod mating;
use mating::{mate, mend_consensus, merge, truncate};

#[derive(Parser)]
struct Cli {
    R1: PathBuf,
    R2: PathBuf,
    interleave: bool,
    csv: bool,
    stats: bool,
}

fn main() {
    let args = Cli::parse();

    // open output file handle
    /*
    let out_handle: Box<dyn Write> = match args.value_of("out") {
        Some(fp) => {
            let path = Path::new(fp);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(io::stdout()), // output to STDOUT
    };
    */
    let out_handle = stdout();

    /*
    let _unmerged_handles: Option<(Box<dyn Write>, Box<dyn Write>)> = match args.value_of("prefix")
    {
        Some(fp) => Some((
            Box::new(File::create(&format!("{}-R1.unmerged.fq", fp)).unwrap()),
            Box::new(File::create(&format!("{}-R2.unmerged.fq", fp)).unwrap()),
        )),
        None => None,
    };
    */
    let fq1: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(
            File::open(Path::new("/home/jknaggs/covid/ERR4659819_1.fastq.gz")).unwrap(),
        )
        .unwrap(),
    ));

    let fq2: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(
            File::open(Path::new("/home/jknaggs/covid/ERR4659819_2.fastq.gz")).unwrap(),
        )
        .unwrap(),
    ));

    //    args.value_of("R1").unwrap();
    let mates1: Fastq<BufReader<MultiGzDecoder<File>>, Seq<Dna>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.R1).unwrap()).unwrap(),
    ));

    let mates2: Fastq<BufReader<MultiGzDecoder<File>>, Seq<Dna>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.R2).unwrap()).unwrap(),
    ));

    let mut merged = out_handle;

    let mut lengths = vec![0; 602];
    let mut total_frags: usize = 0;
    let mut total_merged = 0;

    if args.interleave {
        for (r1, r2) in mates1.zip(mates2) {
            total_frags += 1;
            if let Some((r1, r2)) = interleave_records(&r1, &r2) {
                lengths[r1.seq.len()] += 1;
                lengths[r2.seq.len()] += 1;
                //merged.write_record(&r1).unwrap();
                //merged.write_record(&r2).unwrap();
                total_merged += 2;
            }
        }
    } else {
        for (r1, r2) in mates1.zip(mates2) {
            total_frags += 1;
            if let Some(r) = merge_records(&r1, &r2) {
                lengths[r.seq.len()] += 1;
                //merged.write_record(&r).unwrap();
                total_merged += 1;
            }
        }
    }
    if args.stats {
        eprintln!(
            "Total reads:\t{0}\nMerged reads\t{1}",
            total_frags, total_merged
        );
        eprintln!("{:?}", lengths);
    }
    if args.csv {
        println!("{},{},{:?}", total_frags, total_merged, lengths);
    }
}

fn merge_records(r1: &Record<Seq<Dna>>, r2: &Record<Seq<Dna>>) -> Option<Record<Seq<Dna>>> {
    let r2_rc = r2.seq.rc();
    let r1_rc = r1.seq.rc();

    match mate(&r1.seq, &r2_rc, 25, 20) {
        Some(overlap) => {
            let seq = merge(&r1.seq, &r2_rc, overlap, mend_consensus);
            let qual = merge(&r1.quality, &r2.quality, overlap, cmp::max);
            Some(Record::with_attrs(r1.id(), None, &seq, &qual))
        }
        None => match mate(&r1_rc, &r2.seq, 25, 20) {
            Some(overlap) => {
                let seq = truncate(&r1.seq, &r2_rc, overlap, mend_consensus);
                let qual = truncate(&r1.quality, &r2.quality, overlap, cmp::max);
                Some(Record::with_attrs(r1.id(), None, &seq, &qual))
            }
            None => None,
        },
    }
}

fn interleave_records(r1: &Record, r2: &Record) -> Option<(Record, Record)> {
    let r2_rc = r2.seq.rc();
    let r1_rc = r1.seq.rc();

    match mate(&r1.seq, &r2_rc, 25, 20) {
        // overlap
        Some(_overlap) => Some((
            Record::with_attrs(r1.id(), None, &r1.seq, &r1.quality),
            Record::with_attrs(r2.id(), None, &r2.seq, &r2.quality),
        )),

        // read-through
        None => match mate(&r1_rc, &r2.seq, 25, 20) {
            Some(overlap) => {
                let seq = truncate(&r1.seq, &r2_rc, overlap, mend_consensus);
                let qual_rc = truncate(&r2.quality, &r1.quality, overlap, cmp::max);
                let qual = truncate(&r1.quality, &r2.quality, overlap, cmp::max);
                let seq_rc = seq.rc();
                Some((
                    Record::with_attrs(r1.id(), None, &seq, &qual),
                    Record::with_attrs(r2.id(), None, &seq_rc, &qual_rc),
                ))
            }
            None => None,
        },
    }
}
