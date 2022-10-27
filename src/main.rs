extern crate bio_streams;
extern crate clap;
extern crate flate2;

use bio_streams::fastq::Fastq;
use bio_streams::Record;
use std::io::BufReader;

use flate2::read::MultiGzDecoder;
use std::fs::File;

use clap::Parser;

use std::path::PathBuf;

use std::cmp;

mod mating;
use mating::{mate, mend_consensus, merge, rc, truncate};

#[derive(Parser)]
struct Cli {
    r1: PathBuf,
    r2: PathBuf,
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
    //let out_handle = stdout();

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
    //let mut merged = out_handle;

    let mates1: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r1).unwrap()),
    ));

    let mates2: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r2).unwrap()),
    ));

    //let mut lengths = vec![0; 602];
    let mut total_frags: usize = 0;
    let mut total_merged = 0;

    for (r1, r2) in mates1.zip(mates2) {
        total_frags += 1;
        if let Some(r) = merge_records(&r1, &r2) {
            //lengths[r.seq.len()] += 1;
            println!("{}", r);
            total_merged += 1;
        }
    }

    eprintln!(
        "Total reads:\t{0}\nMerged reads:\t{1}",
        total_frags, total_merged
    );
    /*
    if args.csv {
        println!("{},{},{:?}", total_frags, total_merged, lengths);
    }
    */
}

#[inline]
fn merge_records(r1: &Record<Vec<u8>>, r2: &Record<Vec<u8>>) -> Option<Record<Vec<u8>>> {
    let r2_rc = rc(&r2.seq);
    let r1_rc = rc(&r1.seq);

    let r1_q = if let Some(r1_q) = &r1.quality {
        r1_q
    } else {
        panic!();
    };
    let r2_q = if let Some(r2_q) = &r2.quality {
        r2_q
    } else {
        panic!();
    };

    match mate(&r1.seq, &r2_rc, 25, 20) {
        Some(overlap) => {
            let seq = merge(&r1.seq, &r2_rc, overlap, mend_consensus);
            let qual = merge(r1_q, r2_q, overlap, cmp::max);
            Some(Record {
                fields: r1.fields.clone(),
                seq,
                quality: Some(qual),
            })
        }
        None => match mate(&r1_rc, &r2.seq, 25, 20) {
            Some(overlap) => {
                let seq = truncate(&r1.seq, &r2_rc, overlap, mend_consensus);
                let qual = truncate(r1_q, r2_q, overlap, cmp::max);
                Some(Record {
                    fields: r1.fields.clone(),
                    seq,
                    quality: Some(qual),
                })
            }
            None => None,
        },
    }
}
