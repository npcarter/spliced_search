// parse_ensemble_dat: takes an ensemble dat file that describes a genome and
// parses it into DNA chunks that contain proteins along with vectors that
// describe their introns, exons, start points, and end points.


use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::io::{BufRead,BufReader, Lines};

use std::iter::Peekable;
use clap::Parser;

use parse_ensemble_dat::*;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    infile: String,
    proteinfile: String,
    gapfile: String,
}

// Because the dat file gives the DNA sequence for the contig at the end of its
// description, after we've parsed the contig, we need to go through and fill in
// the protein sequences for each gene.  This function handles that.
fn compute_contig_statistics(contig: &mut Contig, protein_writer: &mut BufWriter<File>, gap_writer: &mut BufWriter<File> ){
  
    let mut coding_regions: Vec<Gap> = Vec::new();
    for gene in contig.genes.iter(){
        let mut gstart:i64 = i64::MAX;
        let mut gend:i64 = 0;
        for protein in gene.proteins.iter(){
            // find the coding region of this protein
            let mut cstart = protein.mrna_start;
            let mut cend = protein.mrna_end;
            if cstart > cend // forget whether this can happen for complemented genes
            {
                println!("swapping {} and {}", cstart, cend);
                std::mem::swap(&mut cstart, &mut cend);}
            if cstart < gstart {gstart = cstart;}
            if cend > gend {gend = cend;}
        }
        coding_regions.push(Gap{start: gstart, end: gend, name: gene.name.clone()});
    }
    coding_regions.sort_by_key(|gap| gap.start);
    if coding_regions.len() > 0 {
        //println!("Coding regions:");
        for region in coding_regions.iter() {
            if region.start >=region.end {
                println!("huh! start {} is greater than end {}", region.start, region.end);
            }
            protein_writer.write_fmt(format_args!("{}\n", region.end - region.start +1)).unwrap();
        }
       // println!("Intergenic regions:");
        let mut prev_end: i64 = 0;
        for region in coding_regions.iter() {
            if region.start > prev_end {
                gap_writer.write_fmt(format_args!("{}\n", region.start - prev_end)).unwrap();
            }
    /*        else{
                if region.name != prev_name {
                    println!("found overlapping coding regions in different genes {} and {}", region.name, prev_name);
                }
            }*/

            if region.end > prev_end {
                prev_end = region.end;
            }

        }
    }
    protein_writer.flush().unwrap();
    gap_writer.flush().unwrap();
}

fn main() {
    let args = Cli::parse();

    let file = File::open(args.infile.clone()).expect(format!("Unable to open {}", args.infile).as_str());
    let proteinfile = OpenOptions::new().create(true).append(true).open(args.proteinfile.clone());
    let mut protein_writer = BufWriter::new(proteinfile.unwrap());
    let gapfile= OpenOptions::new().append(true).create(true).open(args.gapfile.clone());
    let mut gap_writer: BufWriter<File> = BufWriter::new(gapfile.unwrap());
    let reader = BufReader::new(file);
    let mut contig = Contig{dna_sequence: Vec::new(),
        genes: Vec::new(),
        length: 0,
        organism: "".to_string(),
        taxonomy: Vec::new(),
        name: "".to_string()};
    let mut lines:Peekable<Lines<BufReader<File>>> = reader.lines().peekable();
    while let Some(line_result) = lines.next() { // there's still data in the file
        let line = line_result.unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        match fields[0] {
            "ID" => {
                if contig.name != "" {
                    parse_contig(&mut contig);
                    compute_contig_statistics(&mut contig, & mut protein_writer, & mut gap_writer);

                }
                contig = Contig { dna_sequence: Vec::new(), genes: Vec::new(), length: 0, organism: "".to_string(), taxonomy: Vec::new(), name: "".to_string() };
                contig.name = fields[1].to_string();
                contig.length = fields[5].parse::<i64>().unwrap();
                //           println!("contig length is {}", contig.length);
            },
            "OS" => {
                contig.organism = fields[1..].join(" ").to_string();
            },
            "FT" => {
                if fields[1] == "gene" {
                    let temp_gene = parse_gene(&mut lines);
                    if temp_gene.proteins.len() > 0 {
                        contig.genes.push(temp_gene);
                    }
                } else {};
            }
            "SQ" => {
                contig.dna_sequence = parse_sequence(fields, &mut lines);
            },
            "OC" => {
                contig.taxonomy = parse_taxonomy(&line, &mut lines);
            },
            "//" => {
                parse_contig(&mut contig);
                compute_contig_statistics(&mut contig, & mut protein_writer, & mut gap_writer);
            },
            "XX" => {}, // spacer line
            "AC" => {}, // accession line, we don't care about this
            "SV" => {}, // sequence version line, we don't care about this
            "DT" => {}, // date line, we don't care about this
            "DE" => {}, // description line, we don't care about this
            "KW" => {}, // keywords line, we don't care about this
            "CC" => {}, // comment line, we don't care about this
            "FH" => {}, // feature header line, we don't care about this
            _ => println!("Unexpected line: {}", line),
        }
    }

    // ok, we've processed the data file, now extract what we need from it.

}
