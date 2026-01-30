use std::fs::File;
use std::io::{self, BufRead,BufReader, Lines};
use std::path::Path;
use std::iter::Peekable;
use clap::Parser;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    infile: String,
}

struct exon{
    start: i64,
    end: i64,
    complement: bool,
}

struct intron{
    start: i64,
    end: i64,
}

struct protein{
    introns: Vec<intron>,
    exons: Vec<exon>,
    dna_sequence: Vec<char>,
    protein_sequence: Vec<char>,
}
struct gene{
    start: i64, // start position of the gene in the contig
    end: i64, // end position of the gene in the contig.
    // gene runs from start to end, inclusive
    name: String,
    dna_sequence: Vec<char>,
    proteins: Vec<protein>,
}
struct contig{
    length: i64,
    name: String,
    organism: String,
    dna_sequence: Vec<char>,
    genes: Vec<gene>,
}

fn process_contig(contig: &contig){
    println!("Contig {} has length {} and organism {}", contig.name, contig.length, contig.organism);
}

fn parse_sequence(fields: Vec<&str>, lines: &mut Peekable<Lines<BufReader<File>>>) -> Vec<char>{
    let mut sequence:Vec<char> = Vec::new();
    let length:i64 = fields[2].parse::<i64>().unwrap();
    let mut residue_count:i64 = 0;
    let mut next_line = match lines.peek(){
        Some(Ok(line)) => line.to_string(),
        _ => panic!("Unexpected end of file")
    };
    while next_line != "//".to_string(){
        let line = (&mut *lines).next().unwrap().unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        let end:i64 = fields.last().unwrap().parse::<i64>().unwrap();
        for residues in fields[0..(fields.len()-1)].iter(){
            for residue in residues.chars() {
                sequence.push(residue);
                residue_count += 1;
            }
        }
        if residue_count != end{println!("Error: residue count {} does not match end position {}", residue_count, end);}
        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
    }
    println!("Sequence length is {}", sequence.len());
    sequence
}

fn parse_gene(fields: Vec<&str>, lines: &mut Peekable<Lines<BufReader<File>>>) -> gene{
    let mut the_gene = gene{start:0, end:0, dna_sequence: Vec::new(), proteins: Vec::new(), name: "".to_string()};


    let mut next_line = match lines.peek(){
        Some(Ok(line)) => line.to_string(),
        _ => panic!("Unexpected end of file")
    };
    let mut next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    while next_fields[0] != "XX" && next_fields[0] != "SQ" && next_fields[0] != "//" &&
        !(next_fields[0] == "FT" && next_fields[1] == "gene"){
        // We're still in the gene

        let line = (&mut *lines).next().unwrap().unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();

        next_line = match lines.peek(){
            Some(Ok(line)) => line.to_string(),
            _ => panic!("Unexpected end of file")
        };
        next_fields = next_line.split_whitespace().collect::<Vec<&str>>();
    }
    the_gene
}
fn main() {
    let args = Cli::parse();

    let file = File::open(args.infile.clone()).expect(format!("Unable to open {}", args.infile).as_str());
    let reader = io::BufReader::new(file);
    let mut contig = contig{dna_sequence: Vec::new(),
        genes: Vec::new(),
        length: 0,
        organism: "".to_string(),
        name: "".to_string()};
    let mut lines:Peekable<Lines<BufReader<File>>> = reader.lines().peekable();
    let mut gene_count:i64 = 0;
    while let Some(line_result) = lines.next(){ // there's still data in the file
        let line = line_result.unwrap();
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        match fields[0] {
            "ID" => {
                println! {"Processing contig {}", fields[1]};
                if contig.name != "" { process_contig(&contig); }
                contig = contig { dna_sequence: Vec::new(), genes: Vec::new(), length: 0, organism: "".to_string(), name: "".to_string() };
                contig.name = fields[1].to_string();
                contig.length = fields[5].parse::<i64>().unwrap();
                println!("contig length is {}", contig.length);
            },
            "OS" => {
                contig.organism = fields[1..].join(" ").to_string();
            },
            "FT" => {

                if fields[1] == "gene" {
                    gene_count += 1;
                    contig.genes.push(parse_gene(fields, &mut lines))
                } else {};
            }
            "SQ" => {
                contig.dna_sequence = parse_sequence(fields, &mut lines);
            },
            "//" => {
                process_contig(&contig);
            },
            "XX" => {}, // spacer line
            "AC" => {}, // accession line, we don't care about this
            "SV" => {}, // sequence version line, we don't care about this
            "DT" => {}, // date line, we don't care about this
            "DE" => {}, // description line, we don't care about this
            "KW" => {}, // keywords line, we don't care about this
            "OC" => {}, // organism classification line, we don't care about this
            "CC" => {}, // comment line, we don't care about this
            "FH" => {}, // feature header line, we don't care about this
            _ => println!("Unexpected line: {}", line),
        }
    }
    println!("Processed {} genes", gene_count);
}
