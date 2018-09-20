extern crate bio;
extern crate clap;
extern crate flate2;

use bio::alphabets::dna;
use bio::io::fastq::{Reader, Record, Writer};
use bio::pattern_matching::mating::{mate, merge, truncate};

use clap::{App, Arg, ArgMatches};

use flate2::bufread::MultiGzDecoder;

use std::fs::File;

use std::io::{self, BufReader, Write};
use std::path::Path;

fn open_pair(
    r1_path: &str,
    r2_path: &str,
    gzip: bool,
) -> (Box<::std::io::Read>, Box<::std::io::Read>) {
    let (r1_b, r2_b) = (
        BufReader::new(File::open(r1_path).unwrap()),
        BufReader::new(File::open(r2_path).unwrap()),
    );

    if gzip || (r1_path.ends_with(".gz") && r2_path.ends_with(".gz")) {
        (
            Box::new(MultiGzDecoder::new(r1_b)),
            Box::new(MultiGzDecoder::new(r2_b)),
        )
    } else {
        (Box::new(r1_b), Box::new(r2_b))
    }
}

fn main() {
    let args = App::new("merge-reads")
        .version("0.1")
        .about("Merge paired Illumina reads")
        .author("jeff-k <jeff_k@fastmail.com>")
        .arg(
            Arg::with_name("R1")
                .help("Forward read")
                .required(true)
                .index(1),
        ).arg(
            Arg::with_name("R2")
                .help("Reverse read")
                .required(true)
                .index(2),
        ).arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .value_name("FILE")
                .help("Write to fastq file instead of STDOUT")
                .takes_value(true),
        ).arg(
            Arg::with_name("prefix")
                .short("p")
                .long("prefix")
                .value_name("PREFIX")
                .help("Write unmerged reads to paired fastq files at PREFIX-Rx.fq"),
        ).arg(
            Arg::with_name("gzip")
                .short("g")
                .long("gzip")
                .help("Input streams are gzipped"),
        ).arg(
            Arg::with_name("stats")
                .short("s")
                .long("stats")
                .help("Print merge statistics to STDERR when done"),
        ).arg(
            Arg::with_name("csv")
                .long("csv")
                .help("dump fragment length counts as comma separated values"),
        ).get_matches();

    // open output file handle
    let out_handle: Box<Write> = match args.value_of("out") {
        Some(fp) => {
            let path = Path::new(fp);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
        None => Box::new(io::stdout()), // output to STDOUT
    };

    let unmerged_handles: Option<(Box<Write>, Box<Write>)> = match args.value_of("prefix") {
        Some(fp) => Some((
            Box::new(File::create(&format!("{}-R1.unmerged.fq", fp)).unwrap()),
            Box::new(File::create(&format!("{}-R2.unmerged.fq", fp)).unwrap()),
        )),
        None => None,
    };

    //    args.value_of("R1").unwrap();
    let (r1, r2) = open_pair(
        args.value_of("R1").unwrap(),
        args.value_of("R2").unwrap(),
        args.is_present("gzip"),
    );
    let mates1 = Reader::new(r1).records();
    let mates2 = Reader::new(r2).records();

    let mut merged = Writer::new(out_handle);

    let mut lengths = vec![0; 502];
    let mut total_frags: usize = 0;
    let mut total_merged = 0;
    for (r1, r2) in mates1.zip(mates2) {
        total_frags += 1;
        match (r1, r2) {
            (Ok(x), Ok(y)) => match merge_records(&x, &y) {
                Some(r) => {
                    lengths[r.seq().len()] += 1;
                    merged.write_record(&r).unwrap();
                    total_merged += 1;
                }
                None => (),
            },
            _ => eprintln!("unsynced fastqs"),
        }
    }
    if args.is_present("stats") {
        eprintln!(
            "Total reads:\t{0}\nMerged reads\t{1}",
            total_frags, total_merged
        );
        eprintln!("{:?}", lengths);
    }
    if args.is_present("csv") {
        println!("{},{},{:?}", total_frags, total_merged, lengths);
    }
}

fn merge_records(r1: &Record, r2: &Record) -> Option<Record> {
    let r2_rc = dna::revcomp(r2.seq());
    let r1_rc = dna::revcomp(r1.seq());

    match mate(&r1.seq(), &r2_rc, 25, 20) {
        Some(overlap) => {
            let seq = merge(&r1.seq(), &r2_rc, overlap);
            let qual = merge(&r1.qual(), &r2.qual(), overlap);
            Some(Record::with_attrs(r1.id(), None, &seq, &qual))
        }
        None => match mate(&r1_rc, &r2.seq(), 25, 20) {
            Some(overlap) => {
                let seq = truncate(&r1.seq(), &r2_rc, overlap);
                let qual = truncate(&r1.qual(), &r2.qual(), overlap);
                Some(Record::with_attrs(r1.id(), None, &seq, &qual))
            }
            None => None,
        },
    }
}
